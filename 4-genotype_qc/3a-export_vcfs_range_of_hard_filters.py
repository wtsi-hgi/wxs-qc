# export to VCF after annotating with a range of hard filters
from ast import pattern
import os.path
import re
from typing import Union, Optional

import hail as hl
from utils.utils import parse_config, path_spark
from wes_qc import hail_utils, vcf_utils

# TODO move to utils constants
fail_string = "FAIL"
outside_bait_string = "OUTSIDE_BAIT"
pass_medium_string = "PASS_MEDIUM"
pass_stringent_string = "PASS_STRINGENT"


def export_vcfs(
    mtfile: str, filtered_vcf_dir: str, hard_filters: dict[str, dict[str, dict[str : Union[int, float]]]], model_id: str, csq_file: Optional[str] = None, header_file: Optional[str] = None
):
    """
    Export VCFs annotated with a range of hard filters
    :param str mtfile: matrixtable file
    :param str filtered_vcf_dir: output directory for VCFs
    :param dict hard_filters: details of all sets of filters
    :param str model_id: random forest run hash used
    """
    mt = hl.read_matrix_table(path_spark(mtfile))

    # filter to remove rows where all variants fail the most relaxed filters
    mt = mt.filter_rows(mt.info.fraction_pass_relaxed_filters > 0)

    # drop unwanted fields
    mt = mt.annotate_rows(info=mt.info.drop("AC", "AN", "AF"))
    for filter_level in ("relaxed", "medium", "stringent"):
        for metric in ("AN", "AC", "AC_Hom", "AC_Het"):
            filter_metric = f"{filter_level}_{metric}"
            mt = mt.annotate_rows(info=mt.info.annotate(**{filter_metric: mt[filter_metric]}))
            mt = mt.drop(mt[filter_metric])
    #remove variants if their relaxed_AC equals 0 after filtering
    mt=mt.filter_rows(mt.info.relaxed_AC == [0], keep=False)
    mt = mt.drop(
        mt.a_index,
        mt.was_split,
        mt.variant_qc,
        mt.stringent_pass_count,
        mt.medium_pass_count,
        mt.relaxed_pass_count,
        mt.adj,
        mt.sum_AD,
    )
    if "assigned_pop" in mt.col:
        mt = mt.drop(mt.assigned_pop)

    # info for header
    stringent_filters = (
        "SNPs: RF bin<="
        + str(hard_filters["snp"]["stringent"]["bin"])
        + " & DP>="
        + str(hard_filters["snp"]["stringent"]["dp"])
        + " & GQ>="
        + str(hard_filters["snp"]["stringent"]["gq"])
        + " & HetAB>="
        + str(hard_filters["snp"]["stringent"]["ab"])
        + " & Call_Rate>="
        + str(hard_filters["snp"]["stringent"]["call_rate"])
        + ", Indels: RF bin<="
        + str(hard_filters["indel"]["stringent"]["bin"])
        + " & DP>="
        + str(hard_filters["indel"]["stringent"]["dp"])
        + " & GQ>="
        + str(hard_filters["indel"]["stringent"]["gq"])
        + " & HetAB>="
        + str(hard_filters["indel"]["stringent"]["ab"])
        + " & Call_Rate>="
        + str(hard_filters["indel"]["stringent"]["call_rate"])
    )

    medium_filters = (
        "SNPs: RF bin<="
        + str(hard_filters["snp"]["medium"]["bin"])
        + " & DP>="
        + str(hard_filters["snp"]["medium"]["dp"])
        + " & GQ>="
        + str(hard_filters["snp"]["medium"]["gq"])
        + " & HetAB>="
        + str(hard_filters["snp"]["medium"]["ab"])
        + " & Call_Rate>="
        + str(hard_filters["snp"]["medium"]["call_rate"])
        + ", Indels: RF bin<="
        + str(hard_filters["indel"]["medium"]["bin"])
        + " & DP>="
        + str(hard_filters["indel"]["medium"]["dp"])
        + " & GQ>="
        + str(hard_filters["indel"]["medium"]["gq"])
        + " & HetAB>="
        + str(hard_filters["indel"]["medium"]["ab"])
        + " & Call_Rate>="
        + str(hard_filters["indel"]["medium"]["call_rate"])
    )

    relaxed_filters = (
        "SNPs: RF bin<="
        + str(hard_filters["snp"]["relaxed"]["bin"])
        + " & DP>="
        + str(hard_filters["snp"]["relaxed"]["dp"])
        + " & GQ>="
        + str(hard_filters["snp"]["relaxed"]["gq"])
        + " & HetAB>="
        + str(hard_filters["snp"]["relaxed"]["ab"])
        + " & Call_Rate>="
        + str(hard_filters["snp"]["relaxed"]["call_rate"])
        + ", Indels: RF bin<="
        + str(hard_filters["indel"]["relaxed"]["bin"])
        + " & DP>="
        + str(hard_filters["indel"]["relaxed"]["dp"])
        + " & GQ>="
        + str(hard_filters["indel"]["relaxed"]["gq"])
        + " & HetAB>="
        + str(hard_filters["indel"]["relaxed"]["ab"])
        + " & Call_Rate>="
        + str(hard_filters["indel"]["relaxed"]["call_rate"])
    )

    metadata = {
        "info": {
            "stringent_AN": {
                "Description": "Total number of alleles in called genotypes",
                "Number": "1",
                "Type": "Integer",
            },
            "stringent_AC": {"Description": "Allele count in genotypes", "Number": "A", "Type": "Integer"},
            "stringent_AC_Hom": {
                "Description": "Allele counts in homozygous genotypes",
                "Number": "A",
                "Type": "Integer",
            },
            "stringent_AC_Het": {
                "Description": "Allele counts in heterozygous genotypes",
                "Number": "A",
                "Type": "Integer",
            },
            "medium_AN": {
                "Description": "Total number of alleles in called genotypes",
                "Number": "1",
                "Type": "Integer",
            },
            "medium_AC": {"Description": "Allele count in genotypes", "Number": "A", "Type": "Integer"},
            "medium_AC_Hom": {"Description": "Allele counts in homozygous genotypes", "Number": "A", "Type": "Integer"},
            "medium_AC_Het": {
                "Description": "Allele counts in heterozygous genotypes",
                "Number": "A",
                "Type": "Integer",
            },
            "fraction_pass_stringent_filters": {
                "Description": "Fraction of genotypes which pass stringent hard filters " + stringent_filters,
                "Number": "A",
                "Type": "Float",
            },
            "fraction_pass_medium_filters": {
                "Description": "Fraction of genotypes which pass medium hard filters " + medium_filters,
                "Number": "A",
                "Type": "Float",
            },
            "fraction_pass_relaxed_filters": {
                "Description": "Fraction of genotypes which pass relaxed hard filters " + relaxed_filters,
                "Number": "A",
                "Type": "Float",
            },
            "rf_bin": {
                "Description": "Variant QC random forest bin, model id " + model_id,
                "Number": "A",
                "Type": "Integer",
            },
            "rf_score": {
                "Description": "Variant QC random forest score, model id " + model_id,
                "Number": "A",
                "Type": "Float",
            },
        },
        "format": {
            "HetAB": {"Description": "Hetrozygous allele balance", "Number": "A", "Type": "Float"},
            "stringent_filters": {
                "Description": "Pass/fail stringent hard filters " + stringent_filters,
                "Number": "A",
                "Type": "String",
            },
            "medium_filters": {
                "Description": "Pass/fail hard medium filters " + medium_filters,
                "Number": "A",
                "Type": "String",
            },
            "relaxed_filters": {
                "Description": "Pass/fail relaxed hard filters " + relaxed_filters,
                "Number": "A",
                "Type": "String",
            },
        },
    }

    metadata =vcf_utils.modify_vcf_metadata(metadata, csq_file, header_file)

    # export per chromosome
    chroms = [*range(1, 23), "X", "Y"]
    chromosomes = ["chr" + str(chr) for chr in chroms]
    for chromosome in chromosomes:
        print("Exporting " + chromosome)
        mt_chrom = mt.annotate_globals(chromosome=chromosome)
        mt_chrom = mt_chrom.filter_rows(mt_chrom.locus.contig == mt_chrom.chromosome)
        outfile = os.path.join(path_spark(filtered_vcf_dir), f"{chromosome}_hard_filters.vcf.bgz")
        hl.export_vcf(mt_chrom, outfile, metadata=metadata)


def main():
    # = STEP SETUP = #
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #
    hard_filters = config["step4"]["apply_hard_filters"]["hard_filters"]  # set in step 4.1a
    model_id = config["general"]["rf_model_id"]

    # = STEP DEPENDENCIES = #
    mtfile = config["step4"]["export_vcfs_a"]["mtfile"]
    csq_file = config["step4"]["annotate_cq_rf"]["cqfile"]
    csq_header_file =config["step4"]["annotate_cq_rf"]["csq_header"]

    # = STEP OUTPUTS = #
    filtered_vcf_dir = config["step4"]["export_vcfs_a"]["vcf_output_dir"]

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)
    export_vcfs(mtfile, filtered_vcf_dir, hard_filters, model_id, csq_file, csq_header_file)


if __name__ == "__main__":
    main()
