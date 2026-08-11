# identify and prune related samples prior to PCA
# use mt with hard filters and sex annotation from stage2_sample_qc/sqc1_sex_annotation.py
import argparse
import os

import bokeh.plotting as bkplt
import bokeh.layouts as bklayouts
import hail as hl

from wxs_qc.config import get_config
from wxs_qc.hail_utils import path_local, path_spark
from wxs_qc.pca_utils import prune_mt, run_pc_project
from wxs_qc.compute_relatedness import prune_pc_relate
from wxs_qc.filtering import filter_matrix_for_ldprune
from wxs_qc import hail_utils, constants, filtering
from wxs_qc.constants import PEDIGREE_THRESHOLDS


def plot_relatedness(
    relatedness_ht: hl.Table, relatedness_plotfile: str, text_size=constants.plots_text_size, **kwargs
) -> None:
    """
    Plot relatedness scores.

    Input contract:
        Recommended a materialized Table.
        Plotting triggers Hail actions.
    """
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


# TODO: How is this step different from stage2_sample_qc/sqc3_pca_population_prediction.py/run_pca()? Only PCA plotting?
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

    Input contract:
        Requires materialized input MatrixTable.
        The function branches `filtered_mt`, counts columns, performs LD pruning,
        runs PCA, and plots PCA scores.
        The `samples_to_remove` Table may be lazy.
    """
    print("=== Running population PCA")
    print("=== Spliting samples")
    unrelated_mt = filtered_mt.filter_cols(hl.is_defined(samples_to_remove[filtered_mt.col_key]), keep=False)
    related_mt = filtered_mt.filter_cols(hl.is_defined(samples_to_remove[filtered_mt.col_key]))
    samples = unrelated_mt.count_cols()
    print(f"=== {samples} samples for PCA")
    samples = related_mt.count_cols()
    print(f"=== {samples} samples for PC projection")

    # Why do we need plink files?
    # plink_mt = unrelated_mt.annotate_cols(uid=unrelated_mt.s).key_cols_by("uid")
    # hl.export_plink(dataset=plink_mt, output=path_spark(plink_outfile), fam_id=plink_mt.uid, ind_id=plink_mt.uid)
    print("=== Running LD pruning before PCA")
    unrelated_mt = prune_mt(unrelated_mt, prune_args["ld_prune_args"])
    unrelated_mt = unrelated_mt.checkpoint(path_spark(pruned_unrelated_pcrelate_file), overwrite=True)
    related_mt = related_mt.semi_join_rows(unrelated_mt.rows())

    union_pca_scores, pca_scores, pca_loadings, _ = run_pc_project(unrelated_mt, related_mt, pca_components)

    pca_mt = filtered_mt.annotate_cols(scores=union_pca_scores[filtered_mt.col_key].scores)
    print("=== Plotting PC1 vs PC2")
    os.makedirs(os.path.dirname(plot_outfile), exist_ok=True)
    p = hl.plot.scatter(pca_mt.scores[0], pca_mt.scores[1], title="PCA", xlabel="PC1", ylabel="PC2")
    print(f"=== Saving Relatedness PCA plot to {plot_outfile}")
    bkplt.output_file(plot_outfile)
    bkplt.save(p)
    return pca_mt, union_pca_scores, pca_scores, pca_loadings


def order_samples(ht: hl.Table) -> hl.Table:
    """
    Alphabelicatty reorders sample names in a Hail Table based on specified conditions.

    The function updates the `i` and `j` fields in the input Hail Table (`ht`)
    by reordering their values based on their relative size. Specifically,
    it assigns the smaller value to `i` and the larger value to `j` for each row.

    Parameters
    ----------
    ht : hl.Table
        A Hail Table containing fields `i` and `j` that define the sample indices
        to be reordered.

    Returns
    -------
    hl.Table
        A new Hail Table where the `i` and `j` fields are reordered such that
        `i` always contains the smaller value and `j` always contains the larger value.
    """
    return ht.transmute(
        i=hl.if_else(ht.i < ht.j, ht.i, ht.j),
        j=hl.if_else(ht.i > ht.j, ht.i, ht.j),
    )


def filter_trios(ped_ht, relatedness_ht):
    """
    Filters parent-child trios based on relatedness criteria.

    This function processes a pedigree table and a relatedness table to identify
    parent-child trios based on kinship and identity-by-descent (IBD) thresholds.
    It annotates the pedigree table with relatedness statistics and filters it
    to retain only the rows that satisfy the criteria for being true
    parent-child pairs.

    Parameters:
    ped_ht : Table
        The input pedigree table containing sample and parental information.
    relatedness_ht : Table
        The relatedness matrix providing kinship and IBD values for pairs of
        individuals.

    Returns:
    Table
        A filtered table containing annotated pedigree data for individuals
        that meet the parent-child relatedness criteria.
    """

    ped_ht = ped_ht.annotate(
        kin=relatedness_ht[ped_ht.i, ped_ht.j].kin,
        ibd0=relatedness_ht[ped_ht.i, ped_ht.j].ibd0,
        ibd1=relatedness_ht[ped_ht.i, ped_ht.j].ibd1,
        ibd2=relatedness_ht[ped_ht.i, ped_ht.j].ibd2,
    )

    pc_expr = (
        (ped_ht.kin < PEDIGREE_THRESHOLDS.kin_max)
        & (ped_ht.kin > PEDIGREE_THRESHOLDS.kin_min)
        & (ped_ht.ibd2 < PEDIGREE_THRESHOLDS.ibd2_max)
        & (ped_ht.ibd0 < PEDIGREE_THRESHOLDS.ibd0_max)
        & (ped_ht.ibd1 > PEDIGREE_THRESHOLDS.ibd1_min)
    )
    ped_ht = ped_ht.annotate(
        true_pc=pc_expr,
    )
    ped_ht = ped_ht.filter(ped_ht.true_pc)
    return ped_ht


def validate_trios(relatedness_ht, ped_ht):
    """
    Filters the pedigree table and keeps only trios supported by relatedness data.

    The function splits paternal and maternal data,
    keeps only pairs that are supported by relatedness data,
    that only those trios with common family members are retained.

    Parameters:
    relatedness_ht : hail.Table
        A Hail Table containing relatedness information.
    ped_ht : hail.Table
        A Hail Table containing pedigree information.

    Returns:
    hail.Table
        Pedigree with filtered family trios based on relatedness data.
    """
    ht_pat = ped_ht.select(
        ped_ht.fam_id,
        ped_ht.s,
        ped_ht.pat_id,
    )
    ht_mat = ped_ht.select(
        ped_ht.fam_id,
        ped_ht.s,
        ped_ht.mat_id,
    )
    # The relatedness table codex sample names as `i` and `j`
    # To match it with samples in pedigree we ensure that i and j goes in the alphabetic order
    relatedness_ht = order_samples(relatedness_ht)
    relatedness_ht = relatedness_ht.key_by("i", "j")

    # Extracting 'maternal' and 'paternal' parts of pedigree
    ht_mat = ht_mat.rename({"s": "i", "mat_id": "j"})
    ht_pat = ht_pat.rename({"s": "i", "pat_id": "j"})
    # Applying sorting to match relatedness_ht
    ht_mat = order_samples(ht_mat)
    ht_pat = order_samples(ht_pat)

    ht_mat = filter_trios(ht_mat, relatedness_ht)
    ht_pat = filter_trios(ht_pat, relatedness_ht)
    ht_pat = ht_pat.key_by("fam_id")
    ht_mat = ht_mat.key_by("fam_id")
    ped_ht = ped_ht.key_by("fam_id")
    common_cp = ht_mat.semi_join(ht_pat.key_by("fam_id"))
    ped_ht = ped_ht.semi_join(common_cp.key_by("fam_id"))
    return ped_ht


def get_options() -> argparse.Namespace:
    """
    Get options from the command line
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--filter-mt", help="Filter matrix table for LD pruning", action="store_true")
    parser.add_argument("--pc-relate", help="Run PC-Relate and identify related samples", action="store_true")
    parser.add_argument("--plot-pca", help="Plot relatedness and run population PCA", action="store_true")
    parser.add_argument(
        "--validate-trios",
        help="Filter provided pedigree to keep only trios supported by relatedness data",
        action="store_true",
    )
    parser.add_argument("--all", help="Run all steps", action="store_true")
    args = parser.parse_args()
    return args


def main():
    # = STEP SETUP = #
    config = get_config()
    args = get_options()
    if args.all:
        args.filter_mt = True
        args.pc_relate = True
        args.plot_pca = True
        args.validate_trios = True

    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #
    control_list = config["general"]["metadata"]["control_samples"]
    n_partitions = config["general"]["n_partitions"]

    # = STEP DEPENDENCIES = #
    mt_infile = config["stage2"]["impute_sex"]["sex_mt_outfile"]

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)

    # ensure plotting directory exists
    pltdir = path_local(config["general"]["plots_dir"])
    if not os.path.exists(pltdir):
        os.makedirs(pltdir)

    filtered_mt = None
    related_samples_to_remove_ht = None
    relatedness_ht = None

    if args.filter_mt:
        # load input mt
        mt = hl.read_matrix_table(path_spark(mt_infile))
        print(f"=== Variant count before filtering: {mt.count_rows()}")
        # removing control samples
        mt = filtering.remove_samples(mt, control_list)
        # filter matrix to have good variants
        filtered_mt = filter_matrix_for_ldprune(
            mt, config["stage2"]["long_range_ld_file"], **config["stage2"]["filter_params"]
        )
        # We potentially removed a lot of variants, so we need to coalesce partitions
        # No shuffle because subsequent LD pruning is sensitive to the number of partitions
        filtered_mt = filtered_mt.repartition(n_partitions, shuffle=False)
        filtered_mt = filtered_mt.checkpoint(path_spark(config["stage2"]["filtered_mt_outfile"]), overwrite=True)
        print(f"=== Variant count after filtering: {filtered_mt.count_rows()}")

    # -----------------------------------------------------------------------------
    if args.pc_relate:
        if filtered_mt is None:
            filtered_mt = hl.read_matrix_table(path_spark(config["stage2"]["filtered_mt_outfile"]))
        # run pcrelate
        related_samples_to_remove_ht, relatedness_ht = prune_pc_relate(
            filtered_mt,
            config["stage2"]["prune_params"],
            config["stage2"]["king_params"],
            config["stage2"]["pc_relate_params"],
        )

        related_samples_to_remove_ht = related_samples_to_remove_ht.checkpoint(
            path_spark(config["stage2"]["relatedness_output"]["samples_to_remove_file"]), overwrite=True
        )
        related_samples_to_remove_ht.export(path_spark(config["stage2"]["relatedness_output"]["samples_to_remove_tsv"]))
        relatedness_ht = relatedness_ht.checkpoint(
            path_spark(config["stage2"]["relatedness_output"]["relatedness_ht"]), overwrite=True
        )
        relatedness_ht.export(path_spark(config["stage2"]["relatedness_output"]["relatedness_outfile"]))
    # -----------------------------------------------------------------------------

    if args.plot_pca:
        if filtered_mt is None:
            filtered_mt = hl.read_matrix_table(path_spark(config["stage2"]["filtered_mt_outfile"]))
        if related_samples_to_remove_ht is None:
            related_samples_to_remove_ht = hl.read_table(
                path_spark(config["stage2"]["relatedness_output"]["samples_to_remove_file"])
            )
        if relatedness_ht is None:
            relatedness_ht = hl.import_table(
                path_spark(config["stage2"]["relatedness_output"]["relatedness_outfile"]), impute=True, force=True
            )
        # plot relatedness
        plot_relatedness(relatedness_ht, config["stage2"]["relatedness_output"]["relatedness_plotfile"])

        # run PCA
        pca_mt, union_pca_scores, pca_scores, pca_loadings = run_population_pca(
            filtered_mt,
            related_samples_to_remove_ht,
            config["stage2"]["prune_params"],
            **config["stage2"]["prune_plot_pca"],
        )
        pca_mt.write(path_spark(config["stage2"]["prune_plot_pca"]["pca_mt_file"]), overwrite=True)
        union_pca_scores.write(path_spark(config["stage2"]["prune_plot_pca"]["union_pca_scores_file"]), overwrite=True)
        pca_scores.write(path_spark(config["stage2"]["prune_plot_pca"]["pca_scores_file"]), overwrite=True)
        pca_loadings.write(
            path_spark(config["stage2"]["prune_plot_pca"]["pca_loadings_file"]), overwrite=True
        )  # output

    if args.validate_trios and config["stage2"]["relatedness_output"]["revised_pedigree"] is not None:
        print("=== Validating and correcting trios === ")
        if relatedness_ht is None:
            relatedness_ht = hl.read_table(path_spark(config["stage2"]["relatedness_output"]["relatedness_ht"]))
        print(f"=== Relatedness contains {relatedness_ht.count()} records ===")

        ped_ht = hl.import_table(path_spark(config["general"]["metadata"]["pedigree"]), no_header=True, delimiter="\t")
        print(f"=== Pedigree contains {ped_ht.count()} records ===")
        ped_ht = ped_ht.rename(
            {"f0": "fam_id", "f1": "s", "f2": "pat_id", "f3": "mat_id", "f4": "is_female", "f5": "phenotype"}
        )
        relatedness_ht = relatedness_ht.key_by()
        relatedness_ht = relatedness_ht.transmute(i=relatedness_ht.i.s, j=relatedness_ht.j.s)
        ped_ht = validate_trios(relatedness_ht, ped_ht)
        ped_ht.export(path_spark(config["stage2"]["relatedness_output"]["revised_pedigree"]), header=False)


if __name__ == "__main__":
    main()
