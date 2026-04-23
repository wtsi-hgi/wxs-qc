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
from wes_qc.compute_relatedness import prune_pc_relate
from wes_qc import filtering, hail_utils

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
        filtered_kg_mt = filtering.filter_matrix_for_ldprune(kg_unprocessed_mt, conf["long_range_ld_file"], **conf["filter_params"])
        filtered_kg_mt.write(filtered_kg_file, overwrite=True)

    if args.kg_pc_relate:
        filtered_kg_mt = hl.read_matrix_table(filtered_kg_file)
        related_samples_to_remove, _ = prune_pc_relate(
            filtered_kg_mt, conf["prune_params"], conf["king_params"], conf["pc_relate_params"], "ref"
        )

        related_samples_to_remove.write(samples_to_remove_file, overwrite=True)

    if args.kg_remove_related_samples:
        kg_mt = hl.read_matrix_table(kg_unprocessed_mt_file)
        related_samples_to_remove = hl.read_table(samples_to_remove_file)
        kg_mt_remove_related = kg_remove_related_samples(kg_mt, related_samples_to_remove)
        kg_mt_remove_related.write(kg_mt_file, overwrite=True)


if __name__ == "__main__":
    main()
