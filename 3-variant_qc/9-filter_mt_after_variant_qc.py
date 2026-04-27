# filter split matrixtable after variant QC
import hail as hl
import argparse
import os.path
from utils.utils import parse_config, path_spark
from wes_qc import hail_utils


def get_options():
    """
    Get options from the command line
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--snv", type=int, help="SNV threshold")
    parser.add_argument("--indel", type=int, help="Indel threshold")
    args = parser.parse_args()
    if not args.snv and args.indel:
        print("--snv and --indel must be specified")
        exit(1)

    return args


def annotate_mt_with_cq_rf_score_and_bin(
    mtfile: str, rf_htfile: str, snv_threshold: int, indel_threshold: int, filtered_mtfile: str
):
    """
    Annotate matrixtable with RF score and bin then filter SNVs and indels according to threshold
    :param str mtfile: Inputmtfile
    :param str rf_htfile: random forest hail table file
    :param int snv_threshold: bin threshold for SNVs
    :param int indel_threshold: bin threshold for indels
    :param str filtered_mtfile: random forest score annotated mtfile
    """
    mt = hl.read_matrix_table(path_spark(mtfile))
    rf_ht = hl.read_table(path_spark(rf_htfile))

    # keep only vars with rank_id = rank
    rf_ht = rf_ht.filter(rf_ht.rank_id == "rank")

    # annotate mt with score and bin
    mt = mt.annotate_rows(info=mt.info.annotate(rf_score=rf_ht[mt.row_key].score))
    mt = mt.annotate_rows(info=mt.info.annotate(rf_bin=rf_ht[mt.row_key].bin))

    # filter by SNV and indel thresholds
    mt_filtered = mt.filter_rows(
        ((hl.is_snp(mt.alleles[0], mt.alleles[1])) & (mt.info.rf_bin <= snv_threshold))
        | ((hl.is_indel(mt.alleles[0], mt.alleles[1])) & (mt.info.rf_bin <= indel_threshold))
    )

    mt_filtered.write(path_spark(filtered_mtfile), overwrite=True)

    nvars = mt.count_rows()
    nvar_filtered = mt_filtered.count_rows()
    print(f"{str(nvars)} variants before filtering, {str(nvar_filtered)} variants after filtering")


def main():
    # = STEP SETUP =
    args = get_options()
    config = parse_config()
    tmp_dir = config["general"]["tmp_dir"]

    # = STEP PARAMETERS = #
    model_id = config["general"]["rf_model_id"]

    # = STEP DEPENDENCIES = #
    rf_dir = path_spark(config["general"]["var_qc_rf_dir"])
    htfile = os.path.join(rf_dir, model_id, "_gnomad_score_binning_tmp.ht")
    mtfile = config["step3"]["annotate_mt_with_cq_rf_score_and_bin"]["mtfile"]

    # = STEP OUTPUTS = #
    mtoutfile_after_varqc = config["step3"]["annotate_mt_with_cq_rf_score_and_bin"]["mtoutfile_after_varqc"]

    # = STEP LOGIC = #
    _ = hail_utils.init_hl(tmp_dir)
    annotate_mt_with_cq_rf_score_and_bin(mtfile, htfile, args.snv, args.indel, mtoutfile_after_varqc)


if __name__ == "__main__":
    main()
