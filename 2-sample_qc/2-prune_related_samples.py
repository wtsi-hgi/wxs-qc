# identify and prune related samples prior to PCA
# use mt with hard filters and sex annotation from 2-sample_qc/1-hard_filters_sex_annotation.py
import hail as hl
import os
from utils.utils import parse_config, path_local, path_spark
import bokeh.plotting as bkplt
import bokeh.layouts as bklayouts
from wes_qc.pca_utils import prune_mt, run_pc_project
from wes_qc.compute_relatedness import prune_pc_relate
from wes_qc.filtering import filter_matrix_for_ldprune
from wes_qc import hail_utils, constants, filtering

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
    pruned_unrelated_pcrelate_file: str,
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
    unrelated_mt=unrelated_mt.checkpoint(path_spark(pruned_unrelated_pcrelate_file), overwrite=True)
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
    mtoutfile = config["step2"]["prune_params"]["pruned_mt_file"]

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
    filtered_mt = filter_matrix_for_ldprune(mt, config["step2"]["long_range_ld_file"], **config["step2"]["filter_params"])
    filtered_mt.write(path_spark(config["step2"]["filtered_mt_outfile"]), overwrite=True)

    # run pcrelate
    related_samples_to_remove_ht, relatedness_ht = prune_pc_relate(
        filtered_mt, config["step2"]["prune_params"], config["step2"]["king_params"], config["step2"]["pc_relate_params"], "study"
    )
    
    related_samples_to_remove_ht.write(
        path_spark(config["step2"]["relatedness_output"]["samples_to_remove_file"]), overwrite=True
    )
    related_samples_to_remove_ht.export(path_spark(config["step2"]["relatedness_output"]["samples_to_remove_tsv"]))
    relatedness_ht.export(path_spark(config["step2"]["relatedness_output"]["relatedness_outfile"]))

    # plot relatedness
    plot_relatedness(relatedness_ht, **config["step2"]["relatedness_output"]["relatedness_plotfile"])
    # run PCA
    pca_mt, union_pca_scores, pca_scores, pca_loadings = run_population_pca(
        filtered_mt, related_samples_to_remove_ht, config["step2"]["prune_params"], **config["step2"]["prune_plot_pca"]
    )
    pca_mt.write(path_spark(config["step2"]["prune_plot_pca"]["pca_mt_file"]), overwrite=True)
    union_pca_scores.write(path_spark(config["step2"]["prune_plot_pca"]["union_pca_scores_file"]), overwrite=True)
    pca_scores.write(path_spark(config["step2"]["prune_plot_pca"]["pca_scores_file"]), overwrite=True)
    pca_loadings.write(path_spark(config["step2"]["prune_plot_pca"]["pca_loadings_file"]), overwrite=True)  # output


if __name__ == "__main__":
    main()
