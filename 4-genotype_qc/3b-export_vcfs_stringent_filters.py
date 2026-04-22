# export to VCF after annotating with a range of hard filters - flexible filter level
from typing import Union
import argparse
import hail as hl
import os.path
from utils.utils import parse_config, path_spark
from wes_qc import hail_utils


def export_vcfs(
    mtfile: str,
    filtered_vcf_dir: str,
    hard_filters: dict[str, dict[str, dict[str, Union[int, float]]]],
    model_id: str,
    filter_level: str = "stringent"
):
    """
    Export VCFs annotated with a range of hard filters
    :param str mtfile: matrixtable file
    :param str filtered_vcf_dir: output directory for VCFs
    :param dict hard_filters: details of all sets of filters
    :param str model_id: random forest run hash used
    :param str filter_level: filter level to apply - 'relaxed', 'medium', or 'stringent'
    """
    # Validate filter level
    valid_levels = ["relaxed", "medium", "stringent"]
    if filter_level not in valid_levels:
        raise ValueError(f"filter_level must be one of {valid_levels}, got '{filter_level}'")
    
    mt = hl.read_matrix_table(path_spark(mtfile))

    # Define field names based on filter level
    fraction_field = f"fraction_pass_{filter_level}_filters"
    filter_field = f"{filter_level}_filters"
    an_field = f"{filter_level}_AN"
    ac_field = f"{filter_level}_AC"
    ac_hom_field = f"{filter_level}_AC_Hom"
    ac_het_field = f"{filter_level}_AC_Het"

    # Filter to remove rows where all variants fail the selected filters
    mt = mt.filter_rows(mt.info[fraction_field] > 0)
    mt = mt.filter_entries(mt[filter_field] == "Pass")

    # Drop unwanted fields
    mt = mt.drop(
        mt.a_index,
        mt.was_split,
        mt.variant_qc,
        mt.stringent_pass_count,
        mt.medium_pass_count,
        mt.relaxed_pass_count,
        mt.adj,
        mt.assigned_pop,
        mt.sum_AD,
    )
    
    # Drop the fraction fields for other filter levels
    other_levels = [level for level in valid_levels if level != filter_level]
    info_drops = [f"fraction_pass_{level}_filters" for level in other_levels]
    mt = mt.annotate_rows(info=mt.info.drop(*info_drops))
    
    # Recalculating AC, AN, AC_Hom, AC_Het based on filter results
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            AN=mt[an_field],
            AC=mt[ac_field],
            AC_Hom=mt[ac_hom_field],
            AC_Het=mt[ac_het_field]
        )
    )
    
    # Remove variants with AC equals 0
    mt = mt.filter_rows(mt.info.AC == [0], keep=False)
    
    # Remove all filter columns
    mt = mt.drop(mt.relaxed_filters, mt.medium_filters, mt.stringent_filters)
    
    # Calculate AF
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            AF=mt.info.AC / mt.info.AN
        )
    )

    # Info for header - construct filter description
    filter_desc = (
        f"SNPs: RF bin<={hard_filters['snp'][filter_level]['bin']}"
        f" & DP>={hard_filters['snp'][filter_level]['dp']}"
        f" & GQ>={hard_filters['snp'][filter_level]['gq']}"
        f" & HetAB>={hard_filters['snp'][filter_level]['ab']}"
        f" & Call_Rate>={hard_filters['snp'][filter_level]['call_rate']}"
        f", Indels: RF bin<={hard_filters['indel'][filter_level]['bin']}"
        f" & DP>={hard_filters['indel'][filter_level]['dp']}"
        f" & GQ>={hard_filters['indel'][filter_level]['gq']}"
        f" & HetAB>={hard_filters['indel'][filter_level]['ab']}"
        f" & Call_Rate>={hard_filters['indel'][filter_level]['call_rate']}"
    )

    metadata = {
        "format": {
            "HetAB": {"Description": "Hetrozygous allele balance", "Number": "A", "Type": "Float"},
            f"{filter_level}_filters": {
                "Description": f"Pass/fail {filter_level} hard filters {filter_desc}",
                "Number": "A",
                "Type": "String",
            },
        },
        "info": {
            f"fraction_pass_{filter_level}_filters": {
                "Description": f"Fraction of genotypes which pass {filter_level} hard filters {filter_desc}",
                "Number": "A",
                "Type": "Float",
            },
            "rf_score": {
                "Description": f"Variant QC random forest score, model id {model_id}",
                "Number": "A",
                "Type": "Float",
            },
            "rf_bin": {
                "Description": f"Variant QC random forest bin, model id {model_id}",
                "Number": "A",
                "Type": "Integer",
            },
            "AN": {
                "Description": "Total number of alleles in called genotypes after QC",
                "Number": "1",
                "Type": "Integer",
            },
            "AC": {"Description": "Allele count in genotypes after QC", "Number": "A", "Type": "Integer"},
            "AF": {"Description": "Allele Frequency for ALT allele after QC", "Number": "A", "Type": "Integer"},
            "AC_Hom": {"Description": "Allele counts in homozygous genotypes after QC", "Number": "A", "Type": "Integer"},
            "AC_Het": {
                "Description": "Allele counts in heterozygous genotypes after QC",
                "Number": "A",
                "Type": "Integer",
            },
            "CSQ": {
                "Description": "Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomADe_AF|gnomADe_AFR_AF|gnomADe_AMR_AF|gnomADe_ASJ_AF|gnomADe_EAS_AF|gnomADe_FIN_AF|gnomADe_NFE_AF|gnomADe_OTH_AF|gnomADe_SAS_AF|gnomADg_AF|gnomADg_AFR_AF|gnomADg_AMI_AF|gnomADg_AMR_AF|gnomADg_ASJ_AF|gnomADg_EAS_AF|gnomADg_FIN_AF|gnomADg_MID_AF|gnomADg_NFE_AF|gnomADg_OTH_AF|gnomADg_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|PHENOTYPES|Conservation",
                "Number": "A",
                "Type": "String",
            },
            "consequence": {"Description": "Most severe consequence from VEP104", "Number": "A", "Type": "String"},
            "gene": {
                "Description": "Gene affected by the most severe consequence from VEP104",
                "Number": "A",
                "Type": "String",
            },
            "hgnc_id": {
                "Description": "HGNC id of the gene affected by the most severe consequence from VEP104",
                "Number": "A",
                "Type": "String",
            },
        },
    }

    # Export per chromosome
    chroms = [*range(1, 23), "X", "Y"]
    chromosomes = [f"chr{str(chr)}" for chr in chroms]
    for chromosome in chromosomes:
        print(f"Exporting {chromosome}")
        mt_chrom = mt.annotate_globals(chromosome=chromosome)
        mt_chrom = mt_chrom.filter_rows(mt_chrom.locus.contig == mt_chrom.chromosome)
        outfile = os.path.join(path_spark(filtered_vcf_dir), f"{chromosome}_{filter_level}_filters.vcf.bgz")
        hl.export_vcf(mt_chrom, outfile, metadata=metadata)


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Export VCFs with specified filter level (relaxed, medium, or stringent)"
    )
    parser.add_argument(
        "--filter-level",
        type=str,
        choices=["relaxed", "medium", "stringent"],
        default="stringent",
        help="Filter level to apply (default: stringent)"
    )
    args = parser.parse_args()

    # = STEP SETUP = #
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #
    model_id = config["general"]["rf_model_id"]
    hard_filters = config["step4"]["apply_hard_filters"]["hard_filters"]  # set during step 4.1a

    # = STEP DEPENDENCIES = #
    mtfile = config["step4"]["export_vcfs_b"]["mtfile"]

    # = STEP OUTPUTS = #
    filtered_vcf_dir = config["step4"]["export_vcfs_b"]["filtered_vcf_dir"]

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)
    export_vcfs(mtfile, filtered_vcf_dir, hard_filters, model_id, filter_level=args.filter_level)


if __name__ == "__main__":
    main()