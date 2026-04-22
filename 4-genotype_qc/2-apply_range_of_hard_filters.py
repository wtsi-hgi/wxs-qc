# apply a range of different hard filters (RF bin and genotype) to SNPs and indels, also add gene and consequence
# removes samples which fail identity checks
from typing import Union

import hail as hl
import os.path
from utils.utils import parse_config, path_spark
from wes_qc import hail_utils


def remove_samples(mt: hl.MatrixTable, exclude_file: str):
    """
    Remove samples in file of samples which fail identity checks
    :param hl.MatrixTable: Input MatrixTable
    :param str exclude_file: path of file with samples to exclude
    :return: hl.MatrixTable
    """
    excl_ht = hl.import_table(path_spark(exclude_file), impute=True, key="exomeID")
    mt = mt.filter_cols(hl.is_defined(excl_ht[mt.s]), keep=False)

    return mt


def annotate_cq_rf(mt: hl.MatrixTable, rf_htfile: str, cqfile: str) -> hl.MatrixTable:
    """
    Annotate with RF bin, consequence, gene name, hgnc id, and pass/fail for 3 combinations of filters for SNPs,
    3 combinations of filters for indels and missingness for each set of filters
    :param hl.MatrixTable mtfile: Input matrixtable
    :param str rf_htfile: random forest hail table file
    :param str cqfile: consequence file
    :return hl.matrixTable:
    """
    rf_ht = hl.read_table(path_spark(rf_htfile))
    # keep only vars with rank_id = rank
    rf_ht = rf_ht.filter(rf_ht.rank_id == "rank")
    # annotate mt with score and bin
    mt = mt.annotate_rows(info=mt.info.annotate(rf_score=rf_ht[mt.row_key].score))
    mt = mt.annotate_rows(info=mt.info.annotate(rf_bin=rf_ht[mt.row_key].bin))

    cq_ht = hl.import_table(
        path_spark(cqfile),
        types={
            "f0": "str",
            "f1": "int32",
            "f2": "str",
            "f3": "str",
            "f4": "str",
            "f5": "str",
            "f6": "str",
            "f7": "str",
            "f8": "str",
        },
        no_header=True,
    )
    cq_ht = cq_ht.annotate(chr=cq_ht.f0)
    cq_ht = cq_ht.annotate(pos=cq_ht.f1)
    cq_ht = cq_ht.annotate(rs=cq_ht.f2)
    cq_ht = cq_ht.annotate(ref=cq_ht.f3)
    cq_ht = cq_ht.annotate(alt=cq_ht.f4)
    cq_ht = cq_ht.annotate(CSQ=cq_ht.f5)
    cq_ht = cq_ht.annotate(consequence=cq_ht.f6)
    cq_ht = cq_ht.annotate(gene=cq_ht.f7)
    cq_ht = cq_ht.annotate(hgnc_id=cq_ht.f8)
    cq_ht = cq_ht.key_by(locus=hl.locus(cq_ht.chr, cq_ht.pos), alleles=[cq_ht.ref, cq_ht.alt])
    cq_ht = cq_ht.drop(cq_ht.f0, cq_ht.f1, cq_ht.f2, cq_ht.f3, cq_ht.f4, cq_ht.chr, cq_ht.pos, cq_ht.ref, cq_ht.alt)
    cq_ht = cq_ht.key_by(cq_ht.locus, cq_ht.alleles)
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            CSQ=cq_ht[mt.row_key].CSQ,
            consequence=cq_ht[mt.row_key].consequence,
            gene=cq_ht[mt.row_key].gene,
            hgnc_id=cq_ht[mt.row_key].hgnc_id,
        )
    )

    return mt


def annotate_ac(mt: hl.MatrixTable, filter_name: str) -> hl.MatrixTable:
    assert filter_name in ("stringent", "medium", "relaxed")
    filter_field = f"{filter_name}_filters"

    mt_filtered = mt.filter_entries(mt[filter_field] == "Pass")
    mt_filtered = hl.variant_qc(mt_filtered, name="vqc")

    vqc = mt_filtered.index_rows(mt.row_key).vqc
    annotation = {
        f"{filter_name}_AN": vqc.AN,
        f"{filter_name}_AC": vqc.AC[1:],
        f"{filter_name}_AC_Hom": 2 * vqc.homozygote_count[1:],
        f"{filter_name}_AC_Het": vqc.AC[1:] - 2 * vqc.homozygote_count[1:],
    }
    mt = mt.annotate_rows(**annotation)
    return mt


# Build filter conditions for autosomes
def build_autosomal_condition(hard_filters: dict, filter_level: str):
    """Build filtering condition for all chromosomes"""
    condition = (
        (hl.is_snp(mt.alleles[0], mt.alleles[1]))
        & (mt.info.rf_bin <= hard_filters["snp"][filter_level]["bin"])
        & (mt.DP >= hard_filters["snp"][filter_level]["dp"])
        & (mt.GQ >= hard_filters["snp"][filter_level]["gq"])
        & (
            (mt.GT.is_het() & (mt.HetAB >= hard_filters["snp"][filter_level]["ab"]))
            | (mt.GT.is_hom_ref())
            | (mt.GT.is_hom_var())
        )
    ) | (
        (hl.is_indel(mt.alleles[0], mt.alleles[1]))
        & (mt.info.rf_bin <= hard_filters["indel"][filter_level]["bin"])
        & (mt.DP >= hard_filters["indel"][filter_level]["dp"])
        & (mt.GQ >= hard_filters["indel"][filter_level]["gq"])
        & (
            (mt.GT.is_het() & (mt.HetAB >= hard_filters["indel"][filter_level]["ab"]))
            | (mt.GT.is_hom_ref())
            | (mt.GT.is_hom_var())
        )
    )
    return condition

def build_sex_chr_condition(hard_filters: dict, filter_level: str):
    """
    Build filtering condition using different thresholds for sex chromosomes
    For males on chrX non-PAR and chrY: exclude heterozygous calls
    For females on chrX: use standard filtering
    """
    condition = (
        (hl.is_snp(mt.alleles[0], mt.alleles[1]))
        & (mt.info.rf_bin <= hard_filters["snp"][filter_level]["bin"])
        & (
            ((mt.DP >= hard_filters["snp"][filter_level]["dp"]) & (mt.sex == 'female'))
            | ((mt.DP >= hard_filters["snp"][filter_level]["dp"]) & (mt.sex == 'male') & mt.locus.in_autosome_or_par())
            | ((mt.DP >= (hard_filters["snp"][filter_level]["dp"] / 2)) & (mt.sex == 'male') & (mt.locus.in_x_nonpar() | mt.locus.in_y_nonpar()))
        )
        & (mt.GQ >= hard_filters["snp"][filter_level]["gq"])
        & (
            (mt.GT.is_het() & (mt.HetAB >= hard_filters["snp"][filter_level]["ab"]))
            | (mt.GT.is_hom_ref())
            | (mt.GT.is_hom_var())
        )
    ) | (
        (hl.is_indel(mt.alleles[0], mt.alleles[1]))
        & (mt.info.rf_bin <= hard_filters["indel"][filter_level]["bin"])
        & (
            ((mt.DP >= hard_filters["indel"][filter_level]["dp"]) & (mt.sex == 'female'))
            | ((mt.DP >= hard_filters["indel"][filter_level]["dp"]) & (mt.sex == 'male') & mt.locus.in_autosome_or_par())
            | ((mt.DP >= (hard_filters["indel"][filter_level]["dp"] / 2)) & (mt.sex == 'male') & (mt.locus.in_x_nonpar() | mt.locus.in_y_nonpar()))
        )
        & (mt.GQ >= hard_filters["indel"][filter_level]["gq"])
        & (
            (mt.GT.is_het() & (mt.HetAB >= hard_filters["indel"][filter_level]["ab"]))
            | (mt.GT.is_hom_ref())
            | (mt.GT.is_hom_var())
        )
    )
    return condition

def apply_hard_filters(
    mt: hl.MatrixTable, 
    hard_filters: dict[str, dict[str, dict[str, Union[int, float]]]],
    diff_sex_chromosome_filter: bool = False
) -> hl.MatrixTable:
    """
    Apply hard filters and annotate missingness
    :param hl.MatrixTable mt: Input MatrixTable
    :param dict hard_filters: Filters to apply
    :param bool diff_sex_chromosome_filter: If True, apply different filtering for sex chromosomes
    :return hl.MatrixTable:
    """

    # Build final conditions combining autosomal and sex chromosome logic
    if diff_sex_chromosome_filter:
        stringent_condition = build_sex_chr_condition(hard_filters, "stringent")
        medium_condition = build_sex_chr_condition(hard_filters, "medium")
        relaxed_condition = build_sex_chr_condition(hard_filters, "relaxed")
    else:
        # Use standard filtering for all chromosomes
        stringent_condition = build_autosomal_condition(hard_filters, "stringent")
        medium_condition = build_autosomal_condition(hard_filters, "medium")
        relaxed_condition = build_autosomal_condition(hard_filters, "relaxed")
    
    # Adding pass/fail tags to all genotypes
    mt = mt.annotate_entries(
        stringent_filters=hl.if_else(stringent_condition, "Pass", "Fail"),
        medium_filters=hl.if_else(medium_condition, "Pass", "Fail"),
        relaxed_filters=hl.if_else(relaxed_condition, "Pass", "Fail"),
    )

    mt = apply_missingness(
        mt,
        snv_call_rate_stringent=hard_filters["snp"]["stringent"]["call_rate"],
        snv_call_rate_medium=hard_filters["snp"]["medium"]["call_rate"],
        snv_call_rate_relaxed=hard_filters["snp"]["relaxed"]["call_rate"],
        indel_call_rate_stringent=hard_filters["indel"]["stringent"]["call_rate"],
        indel_call_rate_medium=hard_filters["indel"]["medium"]["call_rate"],
        indel_call_rate_relaxed=hard_filters["indel"]["relaxed"]["call_rate"],
    )

    # annotate variants with fraction passing/failing each set of filters
    n_samples = mt.count_cols()
    mt = mt.annotate_rows(stringent_pass_count=hl.agg.count_where(mt.stringent_filters == "Pass"))
    mt = mt.annotate_rows(info=mt.info.annotate(fraction_pass_stringent_filters=mt.stringent_pass_count / n_samples))

    mt = mt.annotate_rows(medium_pass_count=hl.agg.count_where(mt.medium_filters == "Pass"))
    mt = mt.annotate_rows(info=mt.info.annotate(fraction_pass_medium_filters=mt.medium_pass_count / n_samples))

    mt = mt.annotate_rows(relaxed_pass_count=hl.agg.count_where(mt.relaxed_filters == "Pass"))
    mt = mt.annotate_rows(info=mt.info.annotate(fraction_pass_relaxed_filters=mt.relaxed_pass_count / n_samples))

    for filter_name in ("stringent", "medium", "relaxed"):
        mt = annotate_ac(mt, filter_name=filter_name)

    return mt


def apply_missingness(
    mt: hl.MatrixTable,
    snv_call_rate_stringent: float,
    snv_call_rate_medium: float,
    snv_call_rate_relaxed: float,
    indel_call_rate_stringent: float,
    indel_call_rate_medium: float,
    indel_call_rate_relaxed: float,
) -> hl.MatrixTable:
    n = mt.count_cols()  # Count of samples

    # Calculating for each variant the counts of passed/failed genopytes
    mt = mt.annotate_rows(
        stringent_pass_count=hl.agg.count_where(mt.stringent_filters == "Pass"),
        medium_pass_count=hl.agg.count_where(mt.medium_filters == "Pass"),
        relaxed_pass_count=hl.agg.count_where(mt.relaxed_filters == "Pass"),
    )

    # Annotating filters for SNVs
    is_stringent_filters_pass = (
        (hl.is_snp(mt.alleles[0], mt.alleles[1])) & (mt.stringent_pass_count > n * snv_call_rate_stringent)
    ) | ((hl.is_indel(mt.alleles[0], mt.alleles[1])) & (mt.stringent_pass_count > n * indel_call_rate_stringent))

    is_medium_filters_pass = (
        (hl.is_snp(mt.alleles[0], mt.alleles[1])) & (mt.medium_pass_count > n * snv_call_rate_medium)
    ) | ((hl.is_indel(mt.alleles[0], mt.alleles[1])) & (mt.medium_pass_count > n * indel_call_rate_medium))

    is_relaxed_filters_pass = (
        (hl.is_snp(mt.alleles[0], mt.alleles[1])) & (mt.relaxed_pass_count > n * snv_call_rate_relaxed)
    ) | ((hl.is_indel(mt.alleles[0], mt.alleles[1])) & (mt.relaxed_pass_count > n * indel_call_rate_relaxed))

    mt = mt.annotate_entries(
        stringent_filters=hl.if_else(is_stringent_filters_pass, mt.stringent_filters, "Fail"),
        medium_filters=hl.if_else(is_medium_filters_pass, mt.medium_filters, "Fail"),
        relaxed_filters=hl.if_else(is_relaxed_filters_pass, mt.relaxed_filters, "Fail"),
    )

    mt = mt.drop(mt.stringent_pass_count, mt.medium_pass_count, mt.relaxed_pass_count)

    return mt

def sex_annotation(mt: hl.MatrixTable, sex_metadata_file: str, **kwargs) -> hl.MatrixTable:
    """
    Annotates samples in the matrix-table with sex
    """
    metadata_ht = hl.import_table(path_spark(sex_metadata_file), delimiter="\t").key_by("sample_id")
    metadata_ht = metadata_ht.transmute(self_reported_sex=metadata_ht.self_reported_sex.lower())
    mt_sex_annotated = mt.annotate_cols(sex=metadata_ht[mt.s].self_reported_sex)

    # Checking for the samples without self-reported sex
    samples = mt_sex_annotated.cols()
    samples_without_reported_sex = samples.filter(~hl.is_defined(samples.sex))
    n_samples_no_reported_sex = samples_without_reported_sex.count()
    if n_samples_no_reported_sex > 0:
        print(f"=== WARNING: Detected {n_samples_no_reported_sex} samples without reported sex: ", end="")
        print(" ".join(samples_without_reported_sex.s.collect()))
    else:
        print("=== OK: All samples are annotated with reported sex")

    return mt_sex_annotated

def main():
    # = STEP SETUP = #
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #
    model_id = config["general"]["rf_model_id"]
    hard_filters = config["step4"]["apply_hard_filters"]["hard_filters"]
    
    # Get sex chromosome filtering flag from config (default to False if not present)
    diff_sex_chromosome_filter = config["step4"]["apply_hard_filters"]["diff_sex_chromosome_filter"]

    # = STEP DEPENDENCIES = #
    rf_dir = path_spark(
        config["general"]["var_qc_rf_dir"]
    )  # TODO: use adapters inside the functions to enhance robustness
    mtfile = config["step4"]["annotate_cq_rf"]["mtfile"]
    cqfile = config["step4"]["annotate_cq_rf"]["cqfile"]
    sex_annotation_file= config["step4"]["apply_hard_filters"]["sex_annotation_file"]
    # = STEP OUTPUTS = #
    mtfile_annot = config["step4"]["mtoutfile_annot"]

    # = STEP LOGIC = #
    hail_utils.init_hl(tmp_dir)
    # TODO: implement this
    # exclude_file = annotdir + "to_be_excluded_exome.txt"
    rf_htfile = os.path.join(rf_dir, model_id, "_gnomad_score_binning_tmp.ht")  # move runhash to config
    mt = hl.read_matrix_table(path_spark(mtfile))

    # remove unwanted samples
    # mt = remove_samples(mt, exclude_file)

    # annotate mt with consequence, gene, rf bin
    print("=== Annotating mt with consequence, gene, rf bin ===")
    mt_annot = annotate_cq_rf(mt, rf_htfile, cqfile)
    print("=== Applying hard filters ===")
    if diff_sex_chromosome_filter:
        print("=== Sex chromosome-specific filtering ENABLED ===")
        mt_annot=sex_annotation(mt_annot, path_spark(sex_annotation_file))
    # annotate with all combinations of filters (pass/fail) and add missingness
    mt_annot = apply_hard_filters(mt_annot, hard_filters, diff_sex_chromosome_filter=diff_sex_chromosome_filter)
    mt_annot.write(path_spark(mtfile_annot), overwrite=True)


if __name__ == "__main__":
    main()