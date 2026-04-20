# identify and prune related samples prior to PCA
# use mt with hard filters and sex annotation from 2-sample_qc/1-hard_filters_sex_annotation.py
import hail as hl
import os
from utils.utils import parse_config, path_local, path_spark
import bokeh.plotting as bkplt
import bokeh.layouts as bklayouts
from wes_qc import hail_utils, hail_patches, constants, filtering
from gnomad.sample_qc.ancestry import pc_project

def prune_mt(mt: hl.MatrixTable, ld_prune_args, **kwargs) -> hl.MatrixTable:
    """
    Splits multiallelic sites and runs ld pruning
    Filter to autosomes before LD pruning to decrease sample size - autosomes only wanted for later steps
    :param MatrixTable mt: input MT containing variants to be pruned
    :param dict config:
    :return: Pruned MatrixTable
    :rtype: hl.MatrixTable

    `hl.ld_prune` returns a maximal subset of variants that are nearly uncorrelated within each window.
    Requires the dataset to contain only diploid genotype calls and no multiallelic variants.
    See also: https://hail.is/docs/0.2/methods/genetics.html#hail.methods.ld_prune
    """

    print("=== Filtering to autosomes")
    mt = mt.filter_rows(mt.locus.in_autosome())
    print("=== Splitting multiallelic sites")
    mt = hail_patches.split_multi_hts(
        mt, recalculate_gq=False
    )  # this shouldn't do anything as only biallelic sites are used
    print("=== Performing LD pruning")
    pruned_ht = hl.ld_prune(mt.GT, **ld_prune_args)
    pruned_mt = mt.filter_rows(hl.is_defined(pruned_ht[mt.row_key]))
    pruned_mt = pruned_mt.select_entries(GT=hl.unphased_diploid_gt_index_call(pruned_mt.GT.n_alt_alleles()))
    return pruned_mt

def run_filtering(mt: hl.MatrixTable, long_range_ld_file, call_rate_threshold, af_threshold, hwe_threshold) -> hl.MatrixTable:
    print("=== Filtering")
    filtered_mt = filtering.filter_matrix_for_ldprune(
        mt, path_spark(long_range_ld_file), call_rate_threshold, af_threshold, hwe_threshold
    )
    return filtered_mt

def run_king(mt: hl.MatrixTable, king_args: dict, prune_args: dict) -> (hl.MatrixTable, hl.MatrixTable):
    print("=== LD pruning before KING")
    pruned_mt= prune_mt(mt, prune_args["ld_prune_args"])
    pruned_mt.write(path_spark(prune_args["pruned_mt_file"]), overwrite=True)
    print("=== Running KING")
    king_mt = hl.king(pruned_mt.GT)
    king_ht=king_mt.entries()
    king_ht.write(path_spark(king_args["king_output_file"]), overwrite=True)
    related_pairs_ht = king_ht.filter((king_ht.phi > king_args["kinship_threshold"]) & (king_ht.s_1!=king_ht.s))
    print("=== Identifying related samples")
    samples_to_remove=hl.maximal_independent_set(related_pairs_ht.s_1, related_pairs_ht.s, keep=False)
    unrelated_mt = mt.filter_cols(hl.is_defined(samples_to_remove[mt.col_key]), keep=False)
    related_mt = mt.filter_cols(hl.is_defined(samples_to_remove[mt.col_key]))
    print("=== LD pruning unrelated samples")
    pruned_unrelated_mt = prune_mt(unrelated_mt, prune_args["ld_prune_args"])
    pruned_unrelated_mt.write(path_spark(king_args["unrelated_king_file"]), overwrite=True)
    related_mt = related_mt.semi_join_rows(pruned_unrelated_mt.rows())
    related_mt.write(path_spark(king_args["related_king_file"]), overwrite=True)
    return related_mt, pruned_unrelated_mt, pruned_mt

def run_pc_project(mt_ref, mt_study, pca_components):
    print("=== Running PCA")
    pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(mt_ref.GT, k=pca_components, compute_loadings=True)
    pca_af_ht = mt_ref.annotate_rows(pca_af=hl.agg.mean(mt_ref.GT.n_alt_alleles()) / 2).rows()
    pca_loadings = pca_loadings.annotate(pca_af=pca_af_ht[pca_loadings.key].pca_af)
    print("=== Running PC projection")
    projection_pca_scores = pc_project(mt_study, pca_loadings, loading_location="loadings", af_location="pca_af")
    union_pca_scores = pca_scores.union(projection_pca_scores)

    return union_pca_scores, pca_scores, pca_loadings

def prune_pc_relate(
    mt: hl.MatrixTable, prune_args: dict, king_args: dict, pca_components: int, pc_relate_args: dict, relatedness_column: str, relatedness_threshold: float, **kwargs
) -> (hl.Table, hl.Table, hl.Table, hl.Table, hl.Table):
    print("=== Running KING")
    related_mt, unrelated_mt, pruned_mt= run_king(mt, king_args, prune_args)
    print("=== Running PCA")
    union_pca_scores, pca_scores, pca_loadings = run_pc_project(unrelated_mt, related_mt, pca_components)
    print("=== Calculating relatedness")
    relatedness_ht = hl.pc_relate(pruned_mt.GT, scores_expr=union_pca_scores[pruned_mt.col_key].scores, **pc_relate_args)
    # prune individuals to be left with unrelated - creates a table containing one column - samples to remove
    pairs = relatedness_ht.filter(relatedness_ht[relatedness_column] > relatedness_threshold)
    related_samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, keep=False)
    return related_samples_to_remove, union_pca_scores, pca_scores, pca_loadings, relatedness_ht


def plot_relatedness(
    relatedness_ht: hl.Table, relatedness_plotfile: str, text_size=constants.plots_text_size, **kwargs
) -> None:
    """Plot relatedness scores."""
    p1 = hl.plot.histogram(relatedness_ht.kin, title="Kinship distribution")
    p1.axis.axis_label_text_font_size = text_size
    p1.axis.major_label_text_font_size = text_size
    p1.title.text_font_size = text_size
    p2 = hl.plot.scatter(relatedness_ht.kin, relatedness_ht.ibd2, xlabel="Kinship", ylabel="IBD2")
    p2.axis.axis_label_text_font_size = text_size
    p2.axis.major_label_text_font_size = text_size
    layout = bklayouts.gridplot([p1, p2], ncols=2)
    bkplt.output_file(relatedness_plotfile)
    bkplt.save(layout)


# TODO: How is this step different from 2-sample_qc/3-population_pca_prediction.py/run_pca()? Only PCA plotting?
def run_population_pca(
    filtered_mt: hl.MatrixTable,
    samples_to_remove: hl.Table,
    prune_args: dict,
#    plink_outfile,
    pca_components,
    plot_outfile,
    **kwargs,
) -> (hl.MatrixTable, hl.Table, hl.Table):
    """
    Runs PCA and creates a matrix table of non-related individuals with PCA scores
    Remove related samples from PC relate from pruned MT and run PCA
    """
    print("=== Running population PCA")
    print("=== Spliting samples")
    unrelated_mt = filtered_mt.filter_cols(hl.is_defined(samples_to_remove[filtered_mt.col_key]), keep=False)
    related_mt = filtered_mt.filter_cols(hl.is_defined(samples_to_remove[filtered_mt.col_key]))
    variants, samples = unrelated_mt.count()
    print(f"=== {samples} samples for PCA")
    variants, samples = related_mt.count()
    print(f"=== {samples} samples for PC projection")
    #Why do we need plink files?
    #plink_mt = unrelated_mt.annotate_cols(uid=unrelated_mt.s).key_cols_by("uid")
    #hl.export_plink(dataset=plink_mt, output=path_spark(plink_outfile), fam_id=plink_mt.uid, ind_id=plink_mt.uid)
    print("=== Running LD pruning before PCA")
    unrelated_mt = prune_mt(unrelated_mt, prune_args["ld_prune_args"])
    unrelated_mt.write(path_spark(prune_args["pruned_unrelated_pcrelate_file"]), overwrite=True)
    related_mt = related_mt.semi_join_rows(unrelated_mt.rows())

    union_pca_scores, pca_scores, pca_loadings=run_pc_project(unrelated_mt, related_mt, pca_components)
    
    pca_mt = filtered_mt.annotate_cols(scores=union_pca_scores[filtered_mt.col_key].scores)
    print("=== Plotting PC1 vs PC2")
    os.makedirs(os.path.dirname(plot_outfile), exist_ok=True)
    p = hl.plot.scatter(pca_mt.scores[0], pca_mt.scores[1], title="PCA", xlabel="PC1", ylabel="PC2")
    print(f"=== Saving Relatedness PCA plot to {plot_outfile}")
    bkplt.output_file(plot_outfile)
    bkplt.save(p)
    return pca_mt, union_pca_scores, pca_scores, pca_loadings

def main():
    # = STEP SETUP = #
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #
    control_list=config["general"]["metadata"]["control_samples"]

    # = STEP DEPENDENCIES = #
    mt_infile = config["step2"]["impute_sex"]["sex_mt_outfile"]

    # = STEP OUTPUTS = #
    mtoutfile = config["step2"]["prune"]["pruned_mt_file"]

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)

    # ensure plotting directory exists
    pltdir = path_local(config["general"]["plots_dir"])
    if not os.path.exists(pltdir):
        os.makedirs(pltdir)
    # load input mt
    mt = hl.read_matrix_table(path_spark(mt_infile))
    #removing control samples
    mt= filtering.remove_samples(mt, control_list)
    #filter matrix to have good variants
    filtered_mt = run_filtering(mt, **config["step2"]["filter_before_pruning"])
    filtered_mt.write(path_spark(config["step2"]["filtered_mt_outfile"]), overwrite=True)
    # run pcrelate
    related_samples_to_remove_ht, scores, unrelated_scores, loadings, relatedness_ht = prune_pc_relate(
        filtered_mt, config["step2"]["prune"], config["step2"]["king"], **config["step2"]["prune_pc_relate"]
    )
    scores_file1 = config["step2"]["prune_pc_relate"]["scores_file"]
    scores_file2 = config["step2"]["prune_pc_relate"]["unrelated_samples_scores_file"]
    loadings_file = config["step2"]["prune_pc_relate"]["pca_loadings_file_pc_relate"]
    scores.write(path_spark(scores_file1), overwrite=True)  # output
    unrelated_scores.write(path_spark(scores_file2), overwrite=True)  # output
    loadings.write(path_spark(loadings_file), overwrite=True)  # output
    related_samples_to_remove_ht.write(
        path_spark(config["step2"]["prune_pc_relate"]["samples_to_remove_file"]), overwrite=True
    )
    related_samples_to_remove_ht.export(path_spark(config["step2"]["prune_pc_relate"]["samples_to_remove_tsv"]))
    relatedness_ht.export(path_spark(config["step2"]["prune_pc_relate"]["relatedness_outfile"]))

    # plot relatedness
    plot_relatedness(relatedness_ht, **config["step2"]["prune_pc_relate"])
    # run PCA
    pca_mt, union_pca_scores, pca_scores, pca_loadings = run_population_pca(
        filtered_mt, related_samples_to_remove_ht, config["step2"]["prune"], **config["step2"]["prune_plot_pca"]
    )
    pca_mt.write(path_spark(config["step2"]["prune_plot_pca"]["pca_mt_file"]), overwrite=True)
    union_pca_scores.write(path_spark(config["step2"]["prune_plot_pca"]["union_pca_scores_file"]), overwrite=True)
    pca_scores.write(path_spark(config["step2"]["prune_plot_pca"]["pca_scores_file"]), overwrite=True)
    pca_loadings.write(path_spark(config["step2"]["prune_plot_pca"]["pca_loadings_file"]), overwrite=True)  # output


if __name__ == "__main__":
    main()
