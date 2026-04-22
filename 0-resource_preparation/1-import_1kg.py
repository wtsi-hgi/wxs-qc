"""
This script prepares the 1000 Genome matrixtable,
starting from the published 1000G VCF files

You need to run this script only once to prepare the 1000 Genome matrix
Then you can reuse it with other datasets

To prepare data for this script
(you can use the script in scripts/1kg_download):
* Download 1000G VCFs
* Remove SV and leave only SNP and indels
"""

import argparse
from typing import Any

import hail as hl  # type: ignore

from wes_qc.hail_utils import path_spark
from wes_qc.config import get_config
from wes_qc import filtering, hail_utils, hail_patches
from gnomad.sample_qc.ancestry import pc_project

def create_1kg_mt(vcf_indir: str, kg_pop_file: str, **kwargs: dict) -> hl.MatrixTable:
    """
    Create matrixtable of 1kg data
    :param str vcf_indir: the directory with 1KG VCF files
    :param str kg_pop_file: Assigns superpopulations
    """
    print(f"Loading VCFs from {vcf_indir}")
    objects = hl.utils.hadoop_ls(vcf_indir)
    vcfs = [vcf["path"] for vcf in objects if (vcf["path"].startswith("file") and vcf["path"].endswith("vcf.gz"))]
    # create MT
    kg_unprocessed_mt = hl.import_vcf(vcfs, array_elements_required=False, force_bgz=True)
    # Annotating known populations
    kg_pop_file = path_spark(kg_pop_file)
    cohorts_pop = hl.import_table(kg_pop_file, delimiter="\t").key_by("Sample name")
    kg_unprocessed_mt = kg_unprocessed_mt.annotate_cols(
        known_pop=cohorts_pop[kg_unprocessed_mt.s]["Superpopulation code"]
    )
    # Renaming samples to avoid name clashes with cohorts samples
    kg_unprocessed_mt = kg_unprocessed_mt.key_cols_by()
    kg_unprocessed_mt = kg_unprocessed_mt.transmute_cols(s=hl.str("1kg-for-pop-pca_") + kg_unprocessed_mt.s)
    kg_unprocessed_mt = kg_unprocessed_mt.key_cols_by(kg_unprocessed_mt.s)

    return kg_unprocessed_mt


def kg_filter(
    kg_unprocessed_mt: hl.MatrixTable,
    filter_params: dict,
    long_range_ld_file: str,
    **kwargs,
) -> hl.MatrixTable:
    """
    Filter and prune the 1kg data
    :param kg_unprocessed_mt: The KG MT to filter and prune
    :param long_range_ld_file: The long range LD file
    :param call_rate_threshold: The call rate threshold
    :param af_threshold: The allele frequency threshold
    :param hwe_threshold: The HWE threshold
    :param r2_threshold: Squared correlation threshold: https://hail.is/docs/0.2/methods/genetics.html#hail.methods.ld_prune
    :return: The filtered and pruned 1kg MT
    """
    long_range_ld_file = path_spark(long_range_ld_file)
    # Filtering for good variations to make LD prune
    kg_mt_filtered = filtering.filter_matrix_for_ldprune(
        kg_unprocessed_mt,
        long_range_ld_file,
        **filter_params
    )
    # LD pruning - removing variation regions that are related to each other
    return kg_mt_filtered

def kg_remove_related_samples(kg_mt: hl.MatrixTable, related_samples_to_remove: hl.Table) -> hl.MatrixTable:
    variants, samples = kg_mt.count()
    print(f"=== Loaded form initial table: {samples} samples, {variants} variants.")
    print("=== Removing related samples")
    kg_mt_remove_related = kg_mt.filter_cols(hl.is_defined(related_samples_to_remove[kg_mt.col_key]), keep=False)
    variants, samples = kg_mt_remove_related.count()
    print(f"=== Remains after removing related samples: {samples} samples, {variants} variants.")
    return kg_mt_remove_related


def get_options() -> Any:
    """
    Get options from the command line
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--kg-to-mt", help="Convert 1kg data to matrixtable", action="store_true")
    parser.add_argument(
        "--kg-filter-and-prune",
        help="Prune related variants from 1kg matrix",
        action="store_true",
    )
    parser.add_argument("--kg-pc-relate", help="Run PC relate for 1KG ", action="store_true")
    parser.add_argument("--kg-remove-related-samples", help="Run PC relate for 1KG", action="store_true")
    parser.add_argument("--all", help="Run All steps", action="store_true")
    args = parser.parse_args()
    return args

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

def run_king(mt: hl.MatrixTable, king_args: dict, prune_args: dict) -> (hl.MatrixTable, hl.MatrixTable):
    print("=== LD pruning before KING")
    pruned_mt= prune_mt(mt, prune_args["ld_prune_args"])
    pruned_mt.write(path_spark(prune_args["pruned_kg_file"]), overwrite=True)
    print("=== Running KING")
    king_mt = hl.king(pruned_mt.GT)
    king_ht=king_mt.entries()
    king_ht.write(path_spark(king_args["king_kg_file"]), overwrite=True)
    related_pairs_ht = king_ht.filter((king_ht.phi > king_args["kinship_threshold"]) & (king_ht.s_1!=king_ht.s))
    print("=== Identifying related samples")
    samples_to_remove=hl.maximal_independent_set(related_pairs_ht.s_1, related_pairs_ht.s, keep=False)
    unrelated_mt = mt.filter_cols(hl.is_defined(samples_to_remove[mt.col_key]), keep=False)
    related_mt = mt.filter_cols(hl.is_defined(samples_to_remove[mt.col_key]))
    print("=== LD pruning unrelated samples")
    pruned_unrelated_mt = prune_mt(unrelated_mt, prune_args["ld_prune_args"])
    pruned_unrelated_mt.write(path_spark(king_args["unrelated_kg_king_file"]), overwrite=True)
    related_mt = related_mt.semi_join_rows(pruned_unrelated_mt.rows())
    related_mt.write(path_spark(king_args["related_kg_king_file"]), overwrite=True)
    return related_mt, pruned_unrelated_mt

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
    mt: hl.MatrixTable, prune_params: dict, king_params: dict, pc_relate_params: dict,  **kwargs
) -> hl.Table:
    print("=== Running KING")
    related_mt, unrelated_mt= run_king(mt, king_params, prune_params)
    print("=== Running PCA")
    union_pca_scores, pca_scores, pca_loadings = run_pc_project(unrelated_mt, related_mt, pc_relate_params["pca_components"])
    union_pca_scores.write(path_spark(pc_relate_params["scores_file"]), overwrite=True)
    pca_scores.write(path_spark(pc_relate_params["unrelated_samples_scores_file"]), overwrite=True)
    pca_loadings.write(path_spark(pc_relate_params["pca_loadings_file_pc_relate"]), overwrite=True)
    print("=== Calculating relatedness")
    pruned_mt = related_mt.union_cols(unrelated_mt)#check If i should use this or just prune the whole filtered mt
    relatedness_ht = hl.pc_relate(pruned_mt.GT, scores_expr=union_pca_scores[pruned_mt.col_key].scores, **pc_relate_params["pc_relate_args"])
    relatedness_ht.write(path_spark(pc_relate_params["relatedness_ht_file"]), overwrite=True)
    # prune individuals to be left with unrelated - creates a table containing one column - samples to remove
    pairs = relatedness_ht.filter(relatedness_ht[pc_relate_params["relatedness_column"]] > pc_relate_params["relatedness_threshold"])
    related_samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, keep=False)
    return related_samples_to_remove

def main() -> None:
    # = STEP SETUP = #
    config = get_config()
    args = get_options()
    if args.all:
        args.kg_to_mt = True
        args.kg_filter_and_prune = True
        args.kg_pc_relate = True
        args.kg_remove_related_samples = True

    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #
    conf = config["step0"]["create_1kg_mt"]

    # = STEP DEPENDENCIES = #
    vcf_indir = path_spark(conf["indir"])

    # = STEP OUTPUTS = #
    kg_unprocessed_mt_file = path_spark(conf["kg_unprocessed"])
    filtered_kg_file= path_spark(conf["filtered_kg_file"])
    samples_to_remove_file = path_spark(conf["samples_to_remove_file"])
    kg_mt_file = path_spark(conf["kg_out_mt"])

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)

    if args.kg_to_mt:
        kg_mt = create_1kg_mt(vcf_indir, **conf)
        print(f"Saving as hail mt to {kg_unprocessed_mt_file}")
        kg_mt.write(kg_unprocessed_mt_file, overwrite=True)

    if args.kg_filter_and_prune:
        kg_unprocessed_mt = hl.read_matrix_table(kg_unprocessed_mt_file)
        filtered_kg_mt = kg_filter(kg_unprocessed_mt, **conf)
        filtered_kg_mt.write(filtered_kg_file, overwrite=True)

    if args.kg_pc_relate:
        filtered_kg_mt = hl.read_matrix_table(filtered_kg_file)
        related_samples_to_remove = prune_pc_relate(
            filtered_kg_mt, **conf
        )

        related_samples_to_remove.write(samples_to_remove_file, overwrite=True)

    if args.kg_remove_related_samples:
        kg_mt = hl.read_matrix_table(kg_unprocessed_mt_file)
        related_samples_to_remove = hl.read_table(samples_to_remove_file)
        kg_mt_remove_related = kg_remove_related_samples(kg_mt, related_samples_to_remove)
        kg_mt_remove_related.write(kg_mt_file, overwrite=True)


if __name__ == "__main__":
    main()
