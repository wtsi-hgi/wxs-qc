# perform hail sample QC stratified by superpopulation and identify outliers
import os.path

import hail as hl
import pandas as pd
from typing import Optional

from wes_qc import hail_utils, constants
from utils.utils import parse_config, path_spark
from gnomad.sample_qc.filtering import determine_nearest_neighbors, compute_qc_metrics_residuals, compute_stratified_metrics_filter
import bokeh.plotting as bkplot
import bokeh.layouts as bklayouts
from bokeh.models import Div, Span, Range1d, Label
import numpy as np
from collections import defaultdict
import math
from sklearn.linear_model import LinearRegression
import json

#######################################
#for all
#######################################
def run_filtering_and_sample_qc(
    mt: hl.MatrixTable,
    min_depth: int,
    min_genotype_quality: float,
    min_vaf: float,
    **kwargs,
):
    mt = mt.filter_rows(mt.locus.in_autosome())

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

def modify_metric_dict(compute_stratified_metrics_filter_args: dict):
    metric_dict=compute_stratified_metrics_filter_args["metric_threshold"]
    if metric_dict!=None:
        compute_stratified_metrics_filter_args["metric_threshold"] = {
            k: tuple(eval(v_i, {"math": math}) if isinstance(v_i, str) else v_i for v_i in v)
            for k, v in metric_dict.items()
        }
    return compute_stratified_metrics_filter_args

def modify_metric_names(compute_stratified_metrics_filter_args: dict):
    metric_dict=compute_stratified_metrics_filter_args["metric_threshold"]
    if metric_dict!=None:
        for key in list(compute_stratified_metrics_filter_args['metric_threshold'].keys()):
            compute_stratified_metrics_filter_args['metric_threshold'][key + '_residual'] = compute_stratified_metrics_filter_args['metric_threshold'].pop(key)
    return compute_stratified_metrics_filter_args

def get_threshold_dict (compute_stratified_metrics_filter_args: dict, qc_metrics: list):
    threshold_dict={}
    for metric in qc_metrics:
        threshold_dict[metric]=(-compute_stratified_metrics_filter_args['lower_threshold'], compute_stratified_metrics_filter_args['upper_threshold'])
    if compute_stratified_metrics_filter_args['metric_threshold'] != None:
        for metric in compute_stratified_metrics_filter_args['metric_threshold'].keys():
            threshold_dict[metric]=(-compute_stratified_metrics_filter_args['metric_threshold'][metric][0], compute_stratified_metrics_filter_args['metric_threshold'][metric][1])
    return threshold_dict

def compute_mad_score(metric, median, mad):
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
):
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

    if run_with_strata:
        if qc_type=="pop":
            metric_stat_ht = metric_stat_ht.annotate(strata=qc_ht[metric_stat_ht.s].assigned_pop)
        else:
            metric_stat_ht = metric_stat_ht.annotate(strata=qc_ht[metric_stat_ht.s].batch)

    if qc_type=="lr":
        qc_metrics=[metric + "_residual" for metric in qc_metrics]

    stats_dict = hl.literal(qc_metrics_mad_info)

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
    mad_scores_ht=mad_scores_ht.drop(*qc_metrics)
    mad_scores_ht = mad_scores_ht.select(**{col.removesuffix('_mad'): mad_scores_ht[col] for col in mad_scores_ht.row_value})
    return mad_scores_ht

#######################################
#pop
#######################################
# TODO: rename to annotate_with_pop
def annotate_mt(raw_mt_file: str, pop_ht_file: str):
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
    min_depth: int,
    min_genotype_quality: float,
    min_vaf: float,
    compute_stratified_metrics_filter_args: dict,
    **kwargs,
):
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
    mt_with_sampleqc=run_filtering_and_sample_qc(mt, min_depth, min_genotype_quality, min_vaf)
    sample_qc_ht=mt_with_sampleqc.cols()

    #modifying compute_stratified_metrics_filter_args
    compute_stratified_metrics_filter_args=modify_metric_dict(compute_stratified_metrics_filter_args)
    threshold_dict=get_threshold_dict(compute_stratified_metrics_filter_args, qc_metrics)
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
    n_neighbors: int,
    use_batch: bool,
    min_depth: int,
    min_genotype_quality: float,
    min_vaf: float,
    compute_stratified_metrics_filter_args: dict,
    **kwargs,
):
    pca_scores=hl.read_table(path_spark(pca_score_file))
    mt=hl.read_matrix_table(path_spark(raw_mt_file))
    #filtering variants and running sample QC
    mt_with_sampleqc=run_filtering_and_sample_qc(mt, min_depth, min_genotype_quality, min_vaf)
    sample_qc_ht=mt_with_sampleqc.cols()#sample_qc_ht=mt_with_sampleqc.cols()

    #checking if batch column exists when a user wants to use it for stratification
    if use_batch and 'batch' not in sample_qc_ht.row:
        print("=== WARNING! Batch annotation is not found, batch won't be used for stratified qc ===")
        use_batch=False
    print(f"=== Determing {n_neighbors} nearest neighbors ===")
    # Using gnomAD function to determine nearest neighbors
    nn_ht = determine_nearest_neighbors(
        sample_qc_ht,
        pca_scores[sample_qc_ht.key].scores,
        n_neighbors=n_neighbors,
        strata={"batch": sample_qc_ht.batch} if use_batch else None
    )

    #annotating sample_qc_ht with nearest neighbors information
    sample_qc_ht = sample_qc_ht.annotate(
        nearest_neighbors=nn_ht[sample_qc_ht.key].nearest_neighbors
    )
    #modifying compute_stratified_metrics_filter_args
    compute_stratified_metrics_filter_args=modify_metric_dict(compute_stratified_metrics_filter_args)
    threshold_dict=get_threshold_dict(compute_stratified_metrics_filter_args, qc_metrics)

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
    pca_score_file: str,
    qc_metrics: list,
    use_batch: bool,
    min_depth: int,
    min_genotype_quality: float,
    min_vaf: float,
    compute_stratified_metrics_filter_args: dict,
    use_pc_square: bool,
    **kwargs,
):
    pca_scores=hl.read_table(path_spark(pca_score_file))
    mt=hl.read_matrix_table(path_spark(raw_mt_file))
    #filtering variants and running sample QC
    mt_with_sampleqc=run_filtering_and_sample_qc(mt, min_depth, min_genotype_quality, min_vaf)
    sample_qc_ht=mt_with_sampleqc.cols()#sample_qc_ht=mt_with_sampleqc.cols()

    #checking if batch column exists when a user wants to use it for stratification
    if use_batch and 'batch' not in sample_qc_ht.row:
        print("=== WARNING! Batch annotation is not found, batch won't be used for stratified qc ===")
        use_batch=False
    #annotating sample_qc_ht with PCs
    sample_qc_ht = sample_qc_ht.annotate(scores=pca_scores[sample_qc_ht.key].scores)
    print("=== Runnig residual calculation ===")
    # Using gnomAD function to calculate residuals
    sample_qc_res_ht=compute_qc_metrics_residuals(
        sample_qc_ht,
        sample_qc_ht.scores,
        qc_metrics={metric: sample_qc_ht.sample_qc[metric] for metric in qc_metrics},
        strata={"batch": sample_qc_ht.batch} if use_batch else None,
        use_pc_square=use_pc_square)

    #modifying compute_stratified_metrics_filter_args
    compute_stratified_metrics_filter_args=modify_metric_dict(compute_stratified_metrics_filter_args)
    threshold_dict=get_threshold_dict(compute_stratified_metrics_filter_args, qc_metrics)
    compute_stratified_metrics_filter_args_lr=modify_metric_names(compute_stratified_metrics_filter_args)

    print("=== Runnig stratified metrics filter ===")
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
#lm
#######################################
def safe_eval(s):
    """Convert string like '-math.inf' to actual number; return numbers as-is."""
    if isinstance(s, str):
        s=eval(s, {"math": math})
    return s

def run_linear_model(
        qc_ht: hl.Table,
        pca_ht: hl.Table,
        use_batch: bool,
        compute_stratified_metrics_filter_args: dict,
        **kwargs,
    ):
    """
    Returns:
      residuals_df: raw residuals (wide format)
      madnorm_df: MAD-normalized residuals (wide format)
    
    Note: 'batch' column is optional. If present, batch effects will be regressed out.
    """
    qc_df = qc_ht.to_pandas()
    pc_df = pca_ht.to_pandas()
    # -----------------------------
    # Expand PCA scores dynamically
    # -----------------------------
    pc_expanded = pd.DataFrame(
        pc_df["scores"].tolist(),
        index=pc_df.index
    )
    pc_cols = [f"PC{i}" for i in range(1, pc_expanded.shape[1] + 1)]
    pc_expanded.columns = pc_cols
    pc_expanded["s"] = pc_df["s"]
    
    # -----------------------------
    # Merge QC + PCA
    # -----------------------------
    qc_df = qc_df.merge(pc_expanded, on="s", how="left")
    
    # -----------------------------
    # Check if batch column exists
    # -----------------------------
    has_batch = "batch" in qc_df.columns
    
    # -----------------------------
    # Metrics
    # -----------------------------
    metrics = [
        "heterozygosity_rate",
        "n_snp",
        "r_ti_tv",
        "n_transition",
        "n_transversion",
        "r_insertion_deletion",
        "n_insertion",
        "n_deletion",
        "r_het_hom_var"
    ]
    metrics_full = [f"sample_qc.{m}" for m in metrics]
    
    # -----------------------------
    # Output containers (wide)
    # -----------------------------
    if use_batch and not has_batch:
        print("no batch information, batch won't be used for qc")

    if has_batch and use_batch:
        #residuals_df = qc_df[["s", "batch"]].copy()
        madnorm_df = qc_df[["s", "batch"]].copy()
    else:
        #residuals_df = qc_df[["s"]].copy()
        madnorm_df = qc_df[["s"]].copy()
    
    # -----------------------------
    # Loop over metrics
    # -----------------------------
    for m, m_full in zip(metrics, metrics_full):
        valid = qc_df[m_full].notna()
        data = qc_df.loc[valid].copy()
        
        if data.empty:
            continue
        
        # Design matrix - start with PCs
        X_pc = data[pc_cols].fillna(0)
        
        # Add batch dummies if batch column exists
        if has_batch:
            batch_dummies = pd.get_dummies(
                data["batch"],
                prefix="batch",
                drop_first=True
            )
            X = pd.concat([X_pc, batch_dummies], axis=1)
        else:
            X = X_pc
        
        y = data[m_full]
        
        # Fit model
        model = LinearRegression()
        model.fit(X, y)
        residuals = y - model.predict(X)
        
        # MAD normalization
        med = np.median(residuals)
        mad = np.median(np.abs(residuals - med))
        if mad == 0:
            madnorm = residuals - med
        else:
            madnorm = (residuals - med) / (1.4826 * mad)
        
        # Insert back into wide tables
        #residuals_df.loc[valid, m] = residuals
        madnorm_df.loc[valid, m] = madnorm
        madnorm_df.loc[data.index, m] = madnorm.values

    g_lower_threshold=compute_stratified_metrics_filter_args["lower_threshold"]
    g_upper_threshold=compute_stratified_metrics_filter_args["upper_threshold"]
    metric_threshold=compute_stratified_metrics_filter_args["metric_threshold"]
    #print(type(metric_threshold))
    #specific_thresholds=eval(metric_threshold, {"math": math})
    s_lower_thresholds = {k: v[0] for k, v in metric_threshold.items()}
    s_upper_thresholds = {k: v[1] for k, v in metric_threshold.items()}
    thresholds = {
        m: {
            "lower": safe_eval(s_lower_thresholds.get(m, g_lower_threshold)),
            "upper": safe_eval(s_upper_thresholds.get(m, g_upper_threshold))
        }
        for m in metrics
    }
    for metric, t in thresholds.items():
        madnorm_df[f"fail_{metric}"] = (madnorm_df[metric] < t["lower"]) | (madnorm_df[metric] > t["upper"])
    fail_cols = [f"fail_{m}" for m in metrics]
    madnorm_df["qc_metrics_filters"] = madnorm_df.apply(
        lambda row: {col.replace("fail_", "") for col in fail_cols if row[col]},
        axis=1
    )
    fail_cols.append("qc_metrics_filters")
    mad_ht = hl.Table.from_pandas(madnorm_df, key='s')
    mad_ht_subset = mad_ht.select(*fail_cols)
    qc_ht = qc_ht.annotate(**mad_ht_subset[qc_ht.key])
    mad_ht = mad_ht.drop(*fail_cols)
    return mad_ht, qc_ht

def stratified_sample_qc_lm(
    raw_mt_file: str,
    pca_file: str,
    use_batch: bool,
    min_depth: int,
    min_genotype_quality: float,
    min_vaf: float,
    compute_stratified_metrics_filter_args: dict,
    **kwargs,
):
    pca_ht=hl.read_table(path_spark(pca_file))
    mt=hl.read_matrix_table(path_spark(raw_mt_file))
    mt_with_sampleqc=run_filtering_and_sample_qc(mt, min_depth, min_genotype_quality, min_vaf)
    qc_metrics_ht=mt_with_sampleqc.cols()#sample_qc_ht=mt_with_sampleqc.cols()
    ########
    mad_ht, qc_metrics_ht=run_linear_model(qc_metrics_ht, pca_ht, use_batch, compute_stratified_metrics_filter_args)
    #mad_ht=modify_mad_ht(mad_ht, compute_stratified_metrics_filter_args)#made this function
    return mt_with_sampleqc, qc_metrics_ht, mad_ht

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

import math

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
        mt_with_sampleqc, sample_qc_ht, mad_ht, threshold_dict = stratified_sample_qc_pop(mt, qc_metrics, **config["step2"]["stratified_sample_qc"])
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
        mt_with_sampleqc, sample_qc_ht, qc_metrics_stats, mad_ht, nn_ht, use_batch, threshold_dict=stratified_sample_qc_nn(raw_mt_file, pc_scores_file, qc_metrics, **config["step2"]["stratified_sample_qc"])
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
        mt_with_sampleqc, qc_metrics_ht, mad_ht, sample_qc_res_ht, use_batch, threshold_dict=stratified_sample_qc_lr(raw_mt_file, pc_scores_file, qc_metrics, **config["step2"]["stratified_sample_qc"])
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

    elif runmode=="lm":
        mt_with_sampleqc, qc_metrics_ht, mad_ht=stratified_sample_qc_lm(raw_mt_file, pc_scores_file, **config["step2"]["stratified_sample_qc"])
        mt_with_sampleqc.cols().write(path_spark(ht_qc_cols_outfile), overwrite=True)
        qc_metrics_ht.write(path_spark(qc_filter_file), overwrite=True)
        qc_metrics_ht.export(path_spark(output_text_file), delimiter="\t")
        # plot population metrics
        print("=== Plotting population metrics ===")
        use_batch=config["step2"]["stratified_sample_qc"]["use_batch"]
        plot_mad_metrics(mad_ht, use_batch, **config["step2"]["plot_sample_qc_metrics"])


if __name__ == "__main__":
    main()
