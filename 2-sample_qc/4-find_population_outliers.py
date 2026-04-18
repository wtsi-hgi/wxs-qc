# perform hail sample QC stratified by superpopulation and identify outliers
import os.path

import hail as hl
import pandas as pd
from typing import Optional

from wes_qc import hail_utils, constants,filtering
from utils.utils import parse_config, path_spark
from gnomad.sample_qc.filtering import determine_nearest_neighbors, compute_qc_metrics_residuals, compute_stratified_metrics_filter
import bokeh.plotting as bkplot
import bokeh.layouts as bklayouts
from bokeh.models import Div, Span, Range1d, Label
import numpy as np
from collections import defaultdict
import math
import json

#######################################
#for all
#######################################
def run_filtering_and_sample_qc(
    mt: hl.MatrixTable,
    control_list: list,
    min_depth: int,
    min_genotype_quality: float,
    min_vaf: float,
    **kwargs,
) -> hl.MatrixTable:
    mt = mt.filter_rows(mt.locus.in_autosome())
    """
    Filter MT entries based on DP/GQ/VAF and run sample QC
    param hl.MatrixTable mt: input MatrixTable
    param int min_depth: minimum depth threshold for filtering
    param float min_genotype_quality: minimum genotype quality threshold
    param float min_vaf: minimum variant allele fraction threshold for het calls

    ### Config fields
    step2.stratified_sample_qc.min_depth : int : minimum DP threshold
    step2.stratified_sample_qc.min_genotype_quality : float : minimum GQ threshold
    step2.stratified_sample_qc.min_vaf : float : minimum VAF threshold
    """
    # removing control samples
    mt = filtering.remove_samples(mt, control_list)
    # filter MT by depth/gq/vaf
    if min_depth > 0 or min_genotype_quality > 0 or min_vaf > 0:
        vaf = mt.AD[1] / hl.sum(mt.AD)
        print(
            f"Filtering input MT by depth: DP={min_depth}, genotype quality: GQ={min_genotype_quality}, VAF: VAF={min_vaf}"
        )
        filter_condition = (
            (mt.GT.is_het() & (vaf > min_vaf) & (mt.DP > min_depth) & (mt.GQ > min_genotype_quality))
            | (mt.GT.is_hom_ref() & (mt.DP > min_depth) & (mt.GQ > min_genotype_quality))
            | (mt.GT.is_hom_var() & (mt.DP > min_depth) & (mt.GQ > min_genotype_quality))
        )
        fraction_filtered = mt.aggregate_entries(hl.agg.fraction(~filter_condition))
        print(f"Filtering {fraction_filtered * 100:.2f}% entries out of downstream analysis.")
        filtered_mt = mt.filter_entries(filter_condition)
    print("=== Runnig sample QC ===")
    mt_with_sampleqc=hl.sample_qc(filtered_mt)
    mt_with_sampleqc = mt_with_sampleqc.annotate_cols(
        sample_qc=mt_with_sampleqc.sample_qc.annotate(
            heterozygosity_rate=mt_with_sampleqc.sample_qc.n_het / mt_with_sampleqc.sample_qc.n_called
        )
    )
    return mt_with_sampleqc

def modify_metric_dict(compute_stratified_metrics_filter_args: dict) -> dict:
    """
    Convert metric threshold values from strings to numeric values
    param dict compute_stratified_metrics_filter_args: arguments for stratified metrics filter

    ### Config fields
    step2.stratified_sample_qc.metric_threshold : dict : metric-specific thresholds
    step2.stratified_sample_qc.upper_threshold : float : default upper threshold
    step2.stratified_sample_qc.lower_threshold : float : default upper threshold

    """
    metric_dict=compute_stratified_metrics_filter_args.get("metric_threshold", None)
    if metric_dict!=None:
        compute_stratified_metrics_filter_args["metric_threshold"] = {
            k: tuple(eval(v_i, {"math": math}) if isinstance(v_i, str) else v_i for v_i in v)
            for k, v in metric_dict.items()
        }
    compute_stratified_metrics_filter_args["upper_threshold"] = math.inf if compute_stratified_metrics_filter_args.get("upper_threshold") == "math.inf" else compute_stratified_metrics_filter_args.get("upper_threshold")
    compute_stratified_metrics_filter_args["lower_threshold"] = math.inf if compute_stratified_metrics_filter_args.get("lower_threshold") == "math.inf" else compute_stratified_metrics_filter_args.get("lower_threshold")
    return compute_stratified_metrics_filter_args

def modify_metric_names(compute_stratified_metrics_filter_args: dict) -> dict:
    """
    Append '_residual' suffix to metric names in metric_threshold to run compute_stratified_metrics_filter with residuals
    param dict compute_stratified_metrics_filter_args: arguments for stratified metrics filter

    ### Config fields
    step2.stratified_sample_qc.metric_threshold : dict : metric-specific thresholds

    ### Indirect config fields
    """
    metric_dict=compute_stratified_metrics_filter_args.get("metric_threshold", None)
    if metric_dict!=None:
        for key in list(compute_stratified_metrics_filter_args['metric_threshold'].keys()):
            compute_stratified_metrics_filter_args['metric_threshold'][key + '_residual'] = compute_stratified_metrics_filter_args['metric_threshold'].pop(key)
    return compute_stratified_metrics_filter_args

def get_threshold_dict (compute_stratified_metrics_filter_args: dict, qc_metrics: list) -> dict:
    """
    Generate threshold dictionary for plotting
    param dict compute_stratified_metrics_filter_args: arguments with threshold settings
    param list qc_metrics: list of QC metrics

    ### Config fields
    step2.stratified_sample_qc.lower_threshold : float : default lower threshold
    step2.stratified_sample_qc.upper_threshold : float : default upper threshold
    step2.stratified_sample_qc.metric_threshold : dict : metric-specific overrides
    """
    threshold_dict={}
    for metric in qc_metrics:
        threshold_dict[metric]=(-compute_stratified_metrics_filter_args['lower_threshold'], compute_stratified_metrics_filter_args['upper_threshold'])
    if compute_stratified_metrics_filter_args.get('metric_threshold', None) != None:
        for metric in compute_stratified_metrics_filter_args['metric_threshold'].keys():
            threshold_dict[metric]=(-compute_stratified_metrics_filter_args['metric_threshold'][metric][0], compute_stratified_metrics_filter_args['metric_threshold'][metric][1])
    return threshold_dict

def compute_mad_score(metric, median, mad):
    """
    Compute MAD-normalized score for a metric
    param metric: metric value
    param float median: median of the metric
    param float mad: median absolute deviation
    """
    return hl.if_else(
        mad == 0,
        hl.if_else(metric == median, 0.0, hl.float64('inf')),#recheck how this will work
        (metric - median) / mad
    )

def get_mad_ht (
        metric_stat_ht: hl.Table,
        qc_ht: hl.Table,
        filter_ht: hl.Table,
        run_with_strata: bool,
        qc_type: str,
        qc_metrics: list,
) -> hl.Table:
    """
    Compute MAD-normalized QC metrics table for plotting
    param hl.Table metric_stat_ht: table with raw QC metrics or residuals
    param hl.Table qc_ht: sample QC table with annotations
    param hl.Table filter_ht: output from stratified metrics filter
    param bool run_with_strata: whether to compute MAD per strata
    param str qc_type: type of QC ('pop', 'nn', 'lr')
    param list qc_metrics: list of QC metrics
    """
    #get dicts with median and mad for each strata or per sample if nn used
    if qc_type=="pop" or qc_type=="lr":
        qc_stats_raw = filter_ht.globals.collect()[0].qc_metrics_stats
        qc_metrics_mad_info = {}
        if run_with_strata:
            for entry in qc_stats_raw:
                stratum = entry[0]
                qc_metrics_mad_info[stratum] = {}
                for metric_name, stats in qc_stats_raw[entry].items():
                    qc_metrics_mad_info[stratum][metric_name] = {
                        'median': float(stats.median),
                        'mad': float(stats.mad),
                    }
        else:
            stratum = "all"
            qc_metrics_mad_info[stratum] = {}
            for metric_name in qc_stats_raw:
                for stats in qc_stats_raw[metric_name]:
                    qc_metrics_mad_info[stratum][metric_name] = {
                    'median': float(qc_stats_raw[metric_name]["median"]),
                    'mad': float(qc_stats_raw[metric_name]["mad"]),
                    }
    elif qc_type=="nn":
        qc_metrics_mad_info = {}
        for row in filter_ht.key_by().select('s', 'qc_metrics_stats').collect():
            stratum = row.s
            qc_metrics_mad_info[stratum] = {}
            for metric in qc_metrics:
                metric_data = getattr(row.qc_metrics_stats, metric)
                qc_metrics_mad_info[stratum][metric] = {
                    'median': metric_data.median,
                    'mad': metric_data.mad
                }
    stats_dict = hl.literal(qc_metrics_mad_info)
    #run stats annotation with strata if it's used
    if run_with_strata:
        if qc_type=="pop":
            metric_stat_ht = metric_stat_ht.annotate(strata=qc_ht[metric_stat_ht.s].assigned_pop)
        else:
            metric_stat_ht = metric_stat_ht.annotate(strata=qc_ht[metric_stat_ht.s].batch)
    #rename QC metrics for lr
    if qc_type=="lr":
        qc_metrics=[metric + "_residual" for metric in qc_metrics]
    #computing MAD-normalized QC metrics per sample
    if qc_type=="pop" or (qc_type=="lr" and run_with_strata):
        mad_scores_ht = metric_stat_ht.annotate(
            **{
                f'{metric.removesuffix("_residual")}_mad': compute_mad_score(
                    metric_stat_ht[metric],
                    stats_dict[metric_stat_ht.strata][metric]['median'],  # Change 'batch' to your field
                    stats_dict[metric_stat_ht.strata][metric]['mad']
                )
                for metric in qc_metrics
            }
        )
    elif qc_type=="lr" and not run_with_strata:
        mad_scores_ht = metric_stat_ht.annotate(
            **{
                f'{metric.removesuffix("_residual")}_mad': compute_mad_score(
                    metric_stat_ht[metric],
                    stats_dict['all'][metric]['median'],
                    stats_dict['all'][metric]['mad']
                )
                for metric in qc_metrics
            }
        )
    elif qc_type=="nn":
        mad_scores_ht = metric_stat_ht.annotate(
            **{
                f'{metric}_mad': compute_mad_score(
                    metric_stat_ht[metric],
                    stats_dict[metric_stat_ht.s][metric]['median'],
                    stats_dict[metric_stat_ht.s][metric]['mad']
                )
                for metric in qc_metrics
            }
        )
    #cleaning final table
    mad_scores_ht=mad_scores_ht.drop(*qc_metrics)
    mad_scores_ht = mad_scores_ht.select(**{col.removesuffix('_mad'): mad_scores_ht[col] for col in mad_scores_ht.row_value})
    return mad_scores_ht

#######################################
#pop
#######################################
# TODO: rename to annotate_with_pop
def annotate_mt(raw_mt_file: str, pop_ht_file: str) -> hl.MatrixTable:
    """
    Annotate mt with superpopulation and sequencing runid
    :param str raw_mt_file: raw mt file
    :param str pop_ht_file: file with population annotations
    :param str runid_file: metadata file to annotate run ids TODO
    :param str annotated_mt_file: annotated mt file
    :param str pop_pandas_file: tsv from pandas df - EGA and pop
    :param dict config: A config object. No effect.

    ### Config fields
    None

    ### Indirect config fields
    step1.gatk_mt_outfile : input path : used in main
    step2.predict_populations.pop_ht_file : input path : used in main
    step2.annotate_with_pop.annotated_mt_file : output path : used in main
    """
    mt = hl.read_matrix_table(path_spark(raw_mt_file))
    pop_ht = hl.read_table(path_spark(pop_ht_file))
    mt = mt.annotate_cols(assigned_pop=pop_ht[mt.s].pop)
    return mt


def stratified_sample_qc_pop(
    mt: hl.MatrixTable,
    qc_metrics: list,
    control_list: list,
    min_depth: int,
    min_genotype_quality: float,
    min_vaf: float,
    compute_stratified_metrics_filter_args: dict,
    **kwargs,
)-> tuple[hl.MatrixTable, hl.Table, hl.Table, dict]:
    """
    Run sample QC and stratify by population
    param str annotated_mt_file: population and run id annotated MT file
    param str mt_qc_outfile: sample QC MT file
    param str ht_qccols_outfile: sample QC columns HT file
    param str qc_filter_file: output file for stratified sample QC HT
    param dict config: config object

    TODO: note about `if min_depth > 0 or min_genotype_quality > 0 or min_vaf > 0`

    ### Config fields
    step2.stratified_sample_qc.min_depth : float : TODO
    step2.stratified_sample_qc.min_genotype_quality : float : TODO
    step2.stratified_sample_qc.min_vaf : float : TODO
    step2.stratified_sample_qc.output_text_file : output path : TODO
    step2.stratified_sample_qc.output_globals_json_file : output path : TODO

    ### Indirect config fields
    step2.annotate_with_pop.annotated_mt_file : input path : used in main
    step2.stratified_sample_qc.mt_qc_outfile : output path : used in main
    step2.stratified_sample_qc.ht_qc_cols_outfile : output path : used in main
    step2.stratified_sample_qc.qc_filter_file : output path : used in main
    """
    #filtering variants and running sample QC
    mt_with_sampleqc=run_filtering_and_sample_qc(mt, control_list, min_depth, min_genotype_quality, min_vaf)
    sample_qc_ht=mt_with_sampleqc.cols()

    #modifying compute_stratified_metrics_filter_args
    compute_stratified_metrics_filter_args=modify_metric_dict(compute_stratified_metrics_filter_args)
    threshold_dict=get_threshold_dict(compute_stratified_metrics_filter_args, qc_metrics)
    if "comparison_sample_expr" in compute_stratified_metrics_filter_args.keys():
        print ("=== Warning! 'comparison_sample_expr' was provided but will be ignored in 'compute_stratified_metrics_filter'=== ")
        compute_stratified_metrics_filter_args.pop("comparison_sample_expr", None)

    print("=== Runnig stratified metrics filter ===")
    # Using gnomAD function to calculate stratified metrics
    filter_ht = compute_stratified_metrics_filter(
        sample_qc_ht,
        qc_metrics={metric: sample_qc_ht.sample_qc[metric] for metric in qc_metrics},
        strata={"qc_pop": sample_qc_ht.assigned_pop},
        **compute_stratified_metrics_filter_args,
    )
    #generating table with MADs for plots
    metrics_ht = sample_qc_ht.select(**{metric: sample_qc_ht.sample_qc[metric] for metric in qc_metrics})
    mad_ht=get_mad_ht (metrics_ht, sample_qc_ht, filter_ht, True, "pop", qc_metrics)
    #annotating sample_qc_ht with qc_metrics_stats and stratified_metrics_filter results
    globals = hl.eval(filter_ht.globals.qc_metrics_stats)
    sample_qc_ht = sample_qc_ht.annotate_globals(qc_metrics_stats=globals)
    sample_qc_ht = sample_qc_ht.annotate(**filter_ht[sample_qc_ht.key]).persist()
    checkpoint = sample_qc_ht.aggregate(hl.agg.count_where(hl.len(sample_qc_ht.qc_metrics_filters) == 0))
    print(f"=== Samples passing pop filtering: {checkpoint}")

    # return the mad_ht for plotting in the next step
    return mt_with_sampleqc, sample_qc_ht, mad_ht, threshold_dict

#######################################
#nn
#######################################
def stratified_sample_qc_nn(
    raw_mt_file: str,
    pca_score_file: str,
    qc_metrics: list,
    control_list: list,
    use_batch: bool,
    min_depth: int,
    min_genotype_quality: float,
    min_vaf: float,
    compute_stratified_metrics_filter_args: dict,
    determine_nearest_neighbors_args: dict,

    **kwargs,
)-> tuple[hl.MatrixTable, hl.Table, hl.Table, hl.Table, hl.Table, bool, dict]:
    """
    Run sample QC using nearest neighbors
    param str raw_mt_file: input MT file
    param str pca_score_file: PCA scores file
    param list qc_metrics: QC metrics to analyse
    param bool use_batch: whether to stratify by batch
    param int min_depth: minimum DP threshold
    param float min_genotype_quality: minimum GQ threshold
    param float min_vaf: minimum VAF threshold
    param dict compute_stratified_metrics_filter_args: arguments for stratified filtering
    param dict determine_nearest_neighbors_args: arguments for Determing nearest neighbors

    ### Config fields
    step2.stratified_sample_qc.use_batch : bool : use batch for stratification
    step2.stratified_sample_qc.min_depth : int : minimum DP threshold
    step2.stratified_sample_qc.min_genotype_quality : float : minimum GQ threshold
    step2.stratified_sample_qc.min_vaf : float : minimum VAF threshold

    ### Indirect config fields
    step2.prune_plot_pca.union_pca_scores_file : input path : PCA scores
    """
    pca_scores=hl.read_table(path_spark(pca_score_file))
    mt=hl.read_matrix_table(path_spark(raw_mt_file))
    #filtering variants and running sample QC
    mt_with_sampleqc=run_filtering_and_sample_qc(mt, control_list, min_depth, min_genotype_quality, min_vaf)
    sample_qc_ht=mt_with_sampleqc.cols()#sample_qc_ht=mt_with_sampleqc.cols()

    #checking if batch column exists when a user wants to use it for stratification
    if use_batch and 'batch' not in sample_qc_ht.row:
        print("=== WARNING! Batch annotation is not found, batch won't be used for stratified qc ===")
        use_batch=False
    print(f"=== Determing nearest neighbors ===")
    if determine_nearest_neighbors_args is None:
        determine_nearest_neighbors_args={}
    determine_nearest_neighbors_args["strata"]={"batch": sample_qc_ht.batch} if use_batch else None
    determine_nearest_neighbors_args["add_neighbor_distances"]=True
    # Using gnomAD function to determine nearest neighbors
    nn_ht = determine_nearest_neighbors(
        sample_qc_ht,
        pca_scores[sample_qc_ht.key].scores,
        **determine_nearest_neighbors_args
    )

    #annotating sample_qc_ht with nearest neighbors information
    sample_qc_ht = sample_qc_ht.annotate(
        nearest_neighbors=nn_ht[sample_qc_ht.key].nearest_neighbors
    )
    #modifying compute_stratified_metrics_filter_args
    compute_stratified_metrics_filter_args=modify_metric_dict(compute_stratified_metrics_filter_args)
    threshold_dict=get_threshold_dict(compute_stratified_metrics_filter_args, qc_metrics)
    if "comparison_sample_expr" in compute_stratified_metrics_filter_args.keys():
        print ("=== Warning! 'comparison_sample_expr' was provided but will be ignored in 'compute_stratified_metrics_filter' because this pipeline sets that argument internally to 'sample_qc_ht.nearest_neighbors' === ")
        compute_stratified_metrics_filter_args.pop("comparison_sample_expr", None)

    print("=== Runnig stratified metrics filter ===")
    # Using gnomAD function to calculate stratified metrics
    filter_ht = compute_stratified_metrics_filter(
        sample_qc_ht,
        qc_metrics={metric: sample_qc_ht.sample_qc[metric] for metric in qc_metrics},
        comparison_sample_expr=sample_qc_ht.nearest_neighbors,
        **compute_stratified_metrics_filter_args
    )

    #generating table with MADs for plots
    metrics_ht = sample_qc_ht.select(**{metric: sample_qc_ht.sample_qc[metric] for metric in qc_metrics})
    mad_ht=get_mad_ht (metrics_ht, sample_qc_ht, filter_ht, use_batch, "nn", qc_metrics)
    #annotating sample_qc_ht with stratified_metrics_filter results and extracting qc_metrics_stats 
    qc_metrics_stats = filter_ht.key_by().select('s', 'qc_metrics_stats').collect()
    filter_ht = filter_ht.drop("qc_metrics_stats")
    sample_qc_ht = sample_qc_ht.annotate(**filter_ht[sample_qc_ht.key]).persist()
    checkpoint = sample_qc_ht.aggregate(hl.agg.count_where(hl.len(sample_qc_ht.qc_metrics_filters) == 0))
    print(f"=== Samples passing nn filtering: {checkpoint}")
    # return the mad_ht for plotting in the next step
    return mt_with_sampleqc, sample_qc_ht, qc_metrics_stats, mad_ht, nn_ht, use_batch, threshold_dict

#######################################
#lr
#######################################
def stratified_sample_qc_lr(
    raw_mt_file: str,
    pca_score_file: str, # Scored
    qc_metrics: list,
    control_list: list,
    use_batch: bool,
    min_depth: int,
    min_genotype_quality: float,
    min_vaf: float,
    compute_stratified_metrics_filter_args: dict,
    compute_qc_metrics_residuals_args: dict,
    **kwargs,
) -> tuple[hl.MatrixTable, hl.Table, hl.Table, hl.Table, bool, dict]:
    """
    Run sample QC using regression-based residuals
    param str raw_mt_file: input MT file
    param str pca_score_file: PCA scores file from the step 2
    param list qc_metrics: QC metrics to analyse
    param bool use_batch: whether to stratify by batch
    param int min_depth: minimum DP threshold
    param float min_genotype_quality: minimum GQ threshold
    param float min_vaf: minimum VAF threshold
    param dict compute_stratified_metrics_filter_args: arguments for stratified filtering
    param dict compute_qc_metrics_residuals_args: arguments for residual computing

    ### Config fields
    step2.stratified_sample_qc.use_batch : bool : use batch for stratification
    step2.stratified_sample_qc.min_depth : int : minimum DP threshold
    step2.stratified_sample_qc.min_genotype_quality : float : minimum GQ threshold
    step2.stratified_sample_qc.min_vaf : float : minimum VAF threshold

    ### Indirect config fields
    step2.prune_plot_pca.union_pca_scores_file : input path : PCA scores
    """
    pca_scores=hl.read_table(path_spark(pca_score_file))
    mt=hl.read_matrix_table(path_spark(raw_mt_file))
    #filtering variants and running sample QC
    mt_with_sampleqc=run_filtering_and_sample_qc(mt, control_list, min_depth, min_genotype_quality, min_vaf)
    sample_qc_ht=mt_with_sampleqc.cols()#sample_qc_ht=mt_with_sampleqc.cols()

    #checking if batch column exists when a user wants to use it for stratification
    if use_batch and 'batch' not in sample_qc_ht.row:
        print("=== WARNING! Batch annotation is not found, batch won't be used for stratified qc ===")
        use_batch=False
    #annotating sample_qc_ht with PCs
    sample_qc_ht = sample_qc_ht.annotate(scores=pca_scores[sample_qc_ht.key].scores)
    #modifying compute_qc_metrics_residuals_args
    print("=== Running residual calculation ===")
    if compute_qc_metrics_residuals_args is None:
        compute_qc_metrics_residuals_args={}
    compute_qc_metrics_residuals_args["strata"]={"batch": sample_qc_ht.batch} if use_batch else None
    if "regression_sample_inclusion_expr" in compute_qc_metrics_residuals_args.keys():
        print ("=== Warning! 'regression_sample_inclusion_expr' was provided but will be ignored in 'compute_qc_metrics_residuals'=== ")
        compute_qc_metrics_residuals_args.pop("regression_sample_inclusion_expr", None)

    # Using gnomAD function to calculate residuals
    sample_qc_res_ht=compute_qc_metrics_residuals(
        sample_qc_ht,
        sample_qc_ht.scores,
        qc_metrics={metric: sample_qc_ht.sample_qc[metric] for metric in qc_metrics},
        **compute_qc_metrics_residuals_args
    )

    #modifying compute_stratified_metrics_filter_args_lr
    compute_stratified_metrics_filter_args=modify_metric_dict(compute_stratified_metrics_filter_args)
    threshold_dict=get_threshold_dict(compute_stratified_metrics_filter_args, qc_metrics)
    compute_stratified_metrics_filter_args_lr=modify_metric_names(compute_stratified_metrics_filter_args)

    if "comparison_sample_expr" in compute_stratified_metrics_filter_args_lr.keys():
        print ("=== Warning! 'comparison_sample_expr' was provided but will be ignored in 'compute_stratified_metrics_filter'=== ")
        compute_stratified_metrics_filter_args_lr.pop("comparison_sample_expr", None)

    print("=== Running stratified metrics filter ===")
    # Using gnomAD function to calculate stratified metrics
    filter_ht = compute_stratified_metrics_filter(
        sample_qc_res_ht,
        qc_metrics=dict(sample_qc_res_ht.row_value),
        strata={"batch": sample_qc_ht[sample_qc_res_ht.key].batch} if use_batch else None,
        **compute_stratified_metrics_filter_args_lr,
    )

    #generating table with MADs for plots

    mad_ht=get_mad_ht (sample_qc_res_ht, sample_qc_ht, filter_ht, use_batch, 'lr', qc_metrics)
    #annotating sample_qc_ht with qc_metrics_stats and stratified_metrics_filter results
    globals = hl.eval(filter_ht.globals.qc_metrics_stats)
    sample_qc_ht = sample_qc_ht.annotate_globals(qc_metrics_stats=globals)
    sample_qc_ht = sample_qc_ht.annotate(**filter_ht[sample_qc_ht.key]).persist()
    checkpoint = sample_qc_ht.aggregate(hl.agg.count_where(hl.len(sample_qc_ht.qc_metrics_filters) == 0))
    print(f"=== Samples passing lr filtering: {checkpoint}")
    # return the mad_ht for plotting in the next step
    return mt_with_sampleqc, sample_qc_ht, mad_ht, sample_qc_res_ht, use_batch, threshold_dict

#######################################
#plots
#######################################

def plot_sampleqc_metric(
    df: pd.Series,
    metric: str,
    lower_threshold=None,
    upper_threshold=None,
    n_bins=150,
    color="blue",
    plot_width: int = constants.width,
    plot_height: int = constants.height,
    plot_title: Optional[str] = None,
    text_size=constants.plots_text_size,
):
    """
    Plots a histogram of the given metric values for a specific population, with options to add
    thresholds, annotations, captions, and configure visual aspects of the plot.

    Parameters:
        df: Input data series for the metric to be visualized.
        Should be already filtered to contain only single pop and a single metric
        metric: Metric name to be displayed on the x-axis and caption of the plot.
        pop: Name of the population to be included in the plot's caption.
        lower_threshold: A numeric value to mark the lower boundary on the plot.
           Default is None.
        upper_threshold: A numeric value to mark the upper boundary on the plot.
           Default is None.
        n_bins: Number of bins to use for the histogram.
        Default is 150.
        color: Color used for the histogram and plot elements.
        Default is 'blue'.
        plot_width: Width of the plot.
        Default is 800.
        plot_height: Height of the plot.
        Default is 600.

    Returns:
        bokeh.plotting.figure: A Bokeh plot object that visualizes the specified metric.
    """

    p = bkplot.figure(width=plot_width, height=plot_height)
    hist, edges = np.histogram(df[df.notna()], bins=n_bins)
    interval = edges[1] - edges[0]
    padding = max(1, n_bins // 10) * interval
    p.x_range.start = edges[0] - padding
    p.x_range.end = edges[-1] + padding

    p.quad(
        top=hist,
        bottom=0,
        left=edges[:-1],
        right=edges[1:],
        fill_color=color,
        line_color=color,
    )

    if lower_threshold is not None:
        n_below = len(df[df < lower_threshold])
        hline_lower = Span(location=lower_threshold, dimension="height", line_color="red", line_width=2)
        ann_below = Label(
            x=lower_threshold,
            y=hist.max(),
            text=f"{n_below}",
            text_font_size=text_size,
            background_fill_color="white",
            background_fill_alpha=0.5,
            text_align="right",
            x_offset=-5,
        )
        p.add_layout(hline_lower)
        p.add_layout(ann_below)

    if upper_threshold is not None:
        n_above = len(df[df > upper_threshold])
        hline_upper = Span(location=upper_threshold, dimension="height", line_color="red", line_width=2)

        ann_above = Label(
            x=upper_threshold,
            y=hist.max(),
            text=f"{n_above}",
            text_font_size=text_size,
            background_fill_color="white",
            background_fill_alpha=0.5,
            text_align="left",
            x_offset=5,
        )
        p.add_layout(hline_upper)
        p.add_layout(ann_above)
    p.axis.axis_label_text_font_size = text_size
    p.axis.major_label_text_font_size = text_size
    if plot_title is not None:
        p.title.text = plot_title
        p.title.text_font_size = text_size
    p.xaxis.axis_label = metric
    return p

def plot_mad_metrics(
    ht_mad: hl.Table,
    qc_metrics: list,
    plot_outdir: str,
    plot_width: int = 400,
    plot_height: int = 400,
    text_size: str = "12pt",
    n_bins: int = 100,
    default_lower_threshold: float = -4.0,
    default_upper_threshold: float = 4.0,
    use_strata: bool = True,
    metric_thresholds: dict = None,
    color_scheme: str = None,
    **kwargs,
):
    """
    Plot MAD-normalized metrics from ht_mad table.
    
    Parameters:
        ht_mad: Hail Table with MAD-normalized metrics
        plot_outdir: Output directory for plots
        use_strata: Whether to use strata stratification if strata column exists (default True)
        metric_thresholds: Dict mapping metric names to (lower_threshold, upper_threshold) tuples.
                          Use math.inf to disable a threshold.
                          Example: {'n_snp': (-3.0, 5.0), 'r_ti_tv': (-math.inf, 4.0)}
        color_scheme: Color scheme to use. If "1KG", uses predefined colors for 1000 Genomes populations.
                     Otherwise uses default color cycling.
        plot_width: Width of each plot
        plot_height: Height of each plot
        text_size: Font size for labels
        n_bins: Number of histogram bins
        default_lower_threshold: Default lower threshold if not specified in metric_thresholds (default -4.0 MADs)
        default_upper_threshold: Default upper threshold if not specified in metric_thresholds (default 4.0 MADs)
    """
    
    os.makedirs(plot_outdir, exist_ok=True)
    
    # Initialize metric_thresholds if not provided
    if metric_thresholds is None:
        metric_thresholds = {}
    
    # Check if strata column exists AND if user wants to use it
    has_strata = 'strata' in ht_mad.row
    use_strata_stratification = has_strata and use_strata
    
    # If use_strata=False and strata exists, remove it before converting to pandas
    if has_strata and not use_strata:
        print("Strata column exists but use_strata=False - ignoring strata stratification")
        ht_mad = ht_mad.drop('strata')
    
    pd_ht = ht_mad.to_pandas()
    
    # Define metrics
    # Filter to only existing metrics
    metrics = [m for m in qc_metrics if m in pd_ht.columns]
    
    if not metrics:
        raise ValueError("No valid metrics found in table")
    
    # Helper function to get thresholds for a metric
    def get_thresholds(metric):
        """Get lower and upper thresholds for a metric, handling math.inf"""
        if metric in metric_thresholds:
            lower, upper = metric_thresholds[metric]
            # Convert inf to None (no threshold)
            lower = None if math.isinf(lower) else lower
            upper = None if math.isinf(upper) else upper
            return lower, upper
        else:
            # Use defaults
            return default_lower_threshold, default_upper_threshold
    
    if use_strata_stratification:
        print(f"Creating strata-stratified plots for {len(pd_ht['strata'].unique())} strata")
        
        # Colors - use 1KG scheme if specified
        all_strata = sorted(pd_ht['strata'].unique())
        
        if color_scheme == "1KG":
            # 1000 Genomes population color mapping
            all_pops = ["EUR", "EAS", "AFR", "AMR", "SAS", "oth"]
            colours = ["#F0E442", "#D55E00", "#0072B2", "#009E73", "#E69F00", "#56B4E9"]
            pop_color_map = dict(zip(all_pops, colours))
            
            # Assign colors based on population names
            strata_colors = {}
            for stratum in all_strata:
                if stratum in pop_color_map:
                    strata_colors[stratum] = pop_color_map[stratum]
                else:
                    # Fallback for strata not in the predefined list
                    strata_colors[stratum] = "#808080"  # Gray for unknown
            
            print(f"Using 1KG color scheme for populations: {all_pops}")
        else:
            # Default color cycling
            default_colors = ["#F0E442", "#D55E00", "#0072B2", "#009E73", "#E69F00", "#56B4E9"]
            strata_colors = dict(zip(
                all_strata,
                (default_colors * (len(all_strata) // len(default_colors) + 1))[:len(all_strata)]
            ))
        
        # Create plots for grid (don't save these individually)
        grid_plots = []
        plots_by_metric = defaultdict(list)
        
        for stratum in all_strata:
            subdf = pd_ht[pd_ht['strata'] == stratum]
            
            for metric in metrics:
                data = subdf[metric]
                color = strata_colors.get(stratum, "#0072B2")
                
                # Get thresholds for this metric
                lower_thresh, upper_thresh = get_thresholds(metric)
                
                # Create plot for grid
                p_grid = plot_sampleqc_metric(
                    data,
                    metric,
                    lower_threshold=lower_thresh,
                    upper_threshold=upper_thresh,
                    n_bins=n_bins,
                    color=color,
                    plot_width=plot_width,
                    plot_height=plot_height,
                    text_size=text_size,
                )
                
                grid_plots.append(p_grid)
                plots_by_metric[metric].append(p_grid)
                
                # Create SEPARATE plot for individual save
                p_individual = plot_sampleqc_metric(
                    data,
                    metric,
                    lower_threshold=lower_thresh,
                    upper_threshold=upper_thresh,
                    n_bins=n_bins,
                    color=color,
                    plot_width=plot_width,
                    plot_height=plot_height,
                    plot_title=f"{metric} (MAD-normalized) - {stratum}",
                    text_size=text_size,
                )
                
                # Save individual plot
                plot_name = f"MAD_hist_{metric}_{stratum}.html"
                bkplot.output_file(os.path.join(plot_outdir, plot_name))
                bkplot.save(p_individual)
            
            # Add stratum label
            grid_plots.append(
                Div(
                    text=f"<b>{stratum}</b>",
                    width=10,
                    height=plot_height,
                    styles={"writing-mode": "vertical-rl", "text-align": "center"},
                )
            )
        
        # Sync x-ranges across strata for each metric
        for k, p in plots_by_metric.items():
            x_lower = min([pl.x_range.start for pl in p])
            x_upper = max([pl.x_range.end for pl in p])
            for pl in p:
                pl.x_range = Range1d(start=x_lower, end=x_upper)
        
        # Create grid with column titles
        col_titles = [
            Div(text=f"<b>{metric}</b>", width=plot_width, height=15, styles={"text-align": "center"}) 
            for metric in metrics
        ] + [None]
        
        grid = bklayouts.gridplot(col_titles + grid_plots, ncols=len(metrics) + 1)
        plot_outfile = os.path.join(plot_outdir, "mad_metrics_by_strata.html")
        
    else:
        # Plot all samples together (no strata stratification)
        print(f"Creating global plots for all {len(pd_ht)} samples")
        
        grid_plots = []
        
        for metric in metrics:
            data = pd_ht[metric]
            
            # Get thresholds for this metric
            lower_thresh, upper_thresh = get_thresholds(metric)
            
            # Create plot for grid
            p_grid = plot_sampleqc_metric(
                data,
                metric,
                lower_threshold=lower_thresh,
                upper_threshold=upper_thresh,
                n_bins=n_bins,
                color="#0072B2",
                plot_width=plot_width,
                plot_height=plot_height,
                text_size=text_size,
            )
            
            grid_plots.append(p_grid)
            
            # Create SEPARATE plot for individual save
            p_individual = plot_sampleqc_metric(
                data,
                metric,
                lower_threshold=lower_thresh,
                upper_threshold=upper_thresh,
                n_bins=n_bins,
                color="#0072B2",
                plot_width=plot_width,
                plot_height=plot_height,
                plot_title=f"{metric} (MAD-normalized)",
                text_size=text_size,
            )
            
            # Save individual plot
            plot_name = f"MAD_hist_{metric}_all.html"
            bkplot.output_file(os.path.join(plot_outdir, plot_name))
            bkplot.save(p_individual)
        
        grid = bklayouts.gridplot(grid_plots, ncols=3)
        plot_outfile = os.path.join(plot_outdir, "mad_metrics_all.html")
    
    # Save combined grid
    bkplot.output_file(plot_outfile)
    bkplot.save(grid)
    
    print(f"Plots saved to {plot_outdir}")
    print(f"Combined plot: {plot_outfile}")

def main():
    # = STEP SETUP = #
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #
    runmode=config["step2"]["stratified_sample_qc"]["sample_qc_method"]
    control_list=config["general"]["metadata"]["control_samples"]
    if control_list is None:
        control_list=[]

    # = STEP DEPENDENCIES = #
    raw_mt_file = config["step1"]["validate_gtcheck"]["mt_gtcheck_validated"]
    pop_ht_file = config["step2"]["predict_pops"]["pop_ht_outfile"]
    pc_scores_file = config["step2"]["prune_plot_pca"]["union_pca_scores_file"]

    # = STEP OUTPUTS = #
    annotated_mt_file = config["step2"]["annotate_with_pop"]["annotated_mt_file"]
    ht_qc_cols_outfile = config["step2"]["stratified_sample_qc"]["ht_qc_cols_outfile"]
    qc_filter_file = config["step2"]["stratified_sample_qc"]["qc_filter_file"]
    mad_file = config["step2"]["stratified_sample_qc"]["mad_file"]
    output_text_file = config["step2"]["stratified_sample_qc"]["output_text_file"]
    output_stratified_metrics_json = config["step2"]["stratified_sample_qc"]["output_stratified_metrics_json_file"]
    output_nn_file = config["step2"]["stratified_sample_qc"]["output_nn_file"]
    output_residuals_file = config["step2"]["stratified_sample_qc"]["output_residuals_file"]
    output_lms_json = config["step2"]["stratified_sample_qc"]["output_lms_json_file"]

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)
    qc_metrics = [
        "heterozygosity_rate",
        "n_snp",
        "r_ti_tv",
        "n_transition",
        "n_transversion",
        "r_insertion_deletion",
        "n_insertion",
        "n_deletion",
        "r_het_hom_var",
    ]
    if runmode=="pop":
        # annotate mt with runid and pop
        mt = annotate_mt(raw_mt_file, pop_ht_file)
        mt.write(path_spark(annotated_mt_file), overwrite=True)
        # run sample QC and stratify by population
        mt_with_sampleqc, sample_qc_ht, mad_ht, threshold_dict = stratified_sample_qc_pop(mt, qc_metrics, control_list, **config["step2"]["stratified_sample_qc"])
        mt_with_sampleqc.cols().write(path_spark(ht_qc_cols_outfile), overwrite=True)
        sample_qc_ht.write(path_spark(qc_filter_file), overwrite=True)
        sample_qc_ht.export(path_spark(output_text_file), delimiter="\t")
        sample_qc_ht.globals.export(path_spark(output_stratified_metrics_json))
        mad_ht.write(path_spark(mad_file), overwrite=True)
        # plot metrics
        print("=== Plotting metrics ===")
        #all_pops = ["EUR", "EAS", "AFR", "AMR", "SAS", "oth"]
        #colors = ["#F0E442", "#D55E00", "#0072B2", "#009E73", "#E69F00", "#56B4E9"]
        plot_mad_metrics(mad_ht, qc_metrics, **config["step2"]["plot_sample_qc_metrics"], use_strata=True, metric_thresholds=threshold_dict, color_scheme="1KG")

    elif runmode=="nn":
        # run sample QC and stratify by batch if specified
        mt_with_sampleqc, sample_qc_ht, qc_metrics_stats, mad_ht, nn_ht, use_batch, threshold_dict=stratified_sample_qc_nn(raw_mt_file, pc_scores_file, qc_metrics, control_list, **config["step2"]["stratified_sample_qc"])
        mt_with_sampleqc.cols().write(path_spark(ht_qc_cols_outfile), overwrite=True)
        sample_qc_ht.write(path_spark(qc_filter_file), overwrite=True)
        sample_qc_ht.export(path_spark(output_text_file), delimiter="\t")
        nn_ht.export(path_spark(output_nn_file), delimiter="\t")
        metrics_stats_output = {
            "qc_metrics_stats": [
                {
                    "key": [row.s],
                    "value": row.qc_metrics_stats
                }
                for row in qc_metrics_stats
            ]
        }
        with open(output_stratified_metrics_json, 'w') as f:
            json.dump(metrics_stats_output, f, default=str)
        mad_ht.write(path_spark(mad_file), overwrite=True)
        # plot metrics
        print("=== Plotting metrics ===")
        plot_mad_metrics(mad_ht, qc_metrics, color_scheme=None, **config["step2"]["plot_sample_qc_metrics"], use_strata=use_batch, metric_thresholds=threshold_dict)

    elif runmode=="lr":
        # run sample QC and stratify by batch if specified
        mt_with_sampleqc, qc_metrics_ht, mad_ht, sample_qc_res_ht, use_batch, threshold_dict=stratified_sample_qc_lr(raw_mt_file, pc_scores_file, qc_metrics, control_list, **config["step2"]["stratified_sample_qc"])
        mt_with_sampleqc.cols().write(path_spark(ht_qc_cols_outfile), overwrite=True)
        qc_metrics_ht.write(path_spark(qc_filter_file), overwrite=True)
        qc_metrics_ht.export(path_spark(output_text_file), delimiter="\t")
        qc_metrics_ht.globals.export(path_spark(output_stratified_metrics_json))
        sample_qc_res_ht.export(path_spark(output_residuals_file), delimiter="\t")
        sample_qc_res_ht.globals.export(path_spark(output_lms_json))
        mad_ht.write(path_spark(mad_file), overwrite=True)
        # plot metrics
        print("=== Plotting metrics ===")
        plot_mad_metrics(mad_ht, qc_metrics, **config["step2"]["plot_sample_qc_metrics"], use_strata=use_batch, metric_thresholds=threshold_dict)
    else:
        raise ValueError(f"Unknown sample QC method: {runmode}")

if __name__ == "__main__":
    main()
