import json
from pathlib import Path
from typing import Any, cast

import hail as hl
import pytest

from tests.integration_tests.integration_stub import IntegrationTestsStub
from wxs_qc.config import parse_config_file
from wxs_qc.hail_utils import path_spark
from wxs_qc.teszt import assert_saved_tables_match

# /path/to/wxs_qc must be in PYTHONPATH

PEDIGREE_FILE_PATH_TRIOS = """ '${cvars.metadir}/control_set_small.trios.fam' """

GNOMAD_TABLE_SKIP = pytest.mark.skip(
    reason="The test depends on gnomAD table, not present in downloaded test data.\n"
    "Instead we download the resulting table from the bucket and use it as a resource.\n"
)


def _load_expected_results(suite: str, step: str) -> dict[str, Any]:
    expected_results_path = Path(__file__).with_name("expected_integration_test_results.json")
    with open(expected_results_path) as f:
        expected_results = json.load(f)
    return cast(dict[str, Any], expected_results[suite][step])


def _assert_hail_table_count_matches(actual_path: str, expected: dict[str, Any]) -> None:
    actual = hl.read_table(path_spark(actual_path))
    actual_count = actual.count()
    print(f"== VALIDATE: Actual Hail table size: {actual_count} ==")
    assert actual_count == expected["count"]


def _collect_table_samples(ht: hl.Table) -> set[str]:
    return set(ht.aggregate(hl.agg.collect_as_set(ht.s)))


def _collect_table_score_lengths(ht: hl.Table) -> set[int]:
    return set(ht.aggregate(hl.agg.collect_as_set(hl.len(ht.scores))))


def _assert_pca_scores_match_expected(actual_path: str, expected: dict[str, Any]) -> None:
    actual = hl.read_table(path_spark(actual_path))

    assert actual.count() == expected["count"]
    table_samples = _collect_table_samples(actual)
    print(f"== VALIDATE: Actual PCA scores table samples: {table_samples} ==")
    assert table_samples == set(expected["sample_ids"])
    scores_lengths = _collect_table_score_lengths(actual)
    print(f"== VALIDATE: Actual PCA scores table score lengths: {scores_lengths} ==")
    assert scores_lengths == set(expected["score_lengths"])


def _assert_pca_matrix_matches_expected(actual_path: str, expected: dict[str, Any]) -> None:
    actual = hl.read_matrix_table(path_spark(actual_path))
    actual_cols = actual.cols()

    assert actual.count() == (expected["row_count"], expected["col_count"])
    assert _collect_table_samples(actual_cols) == set(expected["sample_ids"])
    assert _collect_table_score_lengths(actual_cols) == set(expected["score_lengths"])


def assert_step_2_2_outputs_match_expected(config_path: str) -> None:
    config = parse_config_file(config_path)
    expected = _load_expected_results("trios", "step_2_2_sample_qc")

    actual_relatedness_output = config["step2"]["relatedness_output"]

    validation_dir = Path(__file__).with_name("validation")
    assert_saved_tables_match(
        validation_dir,
        {
            "related_samples_to_remove.tsv": actual_relatedness_output["samples_to_remove_tsv"],
            "relatedness.tsv": actual_relatedness_output["relatedness_outfile"],
        },
    )

    hail_outputs = expected["hail_outputs"]
    _assert_hail_table_count_matches(
        actual_relatedness_output["samples_to_remove_file"], hail_outputs["samples_to_remove"]
    )

    actual_pc_relate = config["step2"]["pc_relate_params"]
    _assert_pca_scores_match_expected(actual_pc_relate["scores_file"], hail_outputs["pc_relate_scores"])
    _assert_pca_scores_match_expected(
        actual_pc_relate["unrelated_samples_scores_file"], hail_outputs["pc_relate_unrelated_scores"]
    )
    _assert_hail_table_count_matches(
        actual_pc_relate["pca_loadings_file_pc_relate"], hail_outputs["pc_relate_pca_loadings"]
    )

    actual_pca = config["step2"]["prune_plot_pca"]
    _assert_pca_scores_match_expected(actual_pca["union_pca_scores_file"], hail_outputs["population_pca_scores"])
    _assert_pca_scores_match_expected(actual_pca["pca_scores_file"], hail_outputs["population_pca_unrelated_scores"])
    _assert_hail_table_count_matches(actual_pca["pca_loadings_file"], hail_outputs["population_pca_loadings"])
    _assert_pca_matrix_matches_expected(actual_pca["pca_mt_file"], hail_outputs["population_pca_mt"])


def assert_step_4_1_outputs_match_expected(config_path: str) -> None:
    evaluation_config = parse_config_file(config_path)["step4"]["evaluation"]
    validation_dir = Path(__file__).with_name("validation")
    assert_saved_tables_match(
        validation_dir,
        {
            "hard_filter_evaluation.snv.tsv": evaluation_config["snp_tsv"],
            "hard_filter_evaluation.indel.tsv": evaluation_config["indel_tsv"],
        },
    )


@pytest.mark.usefixtures("WES_CONFIG")
class TestIntegration(IntegrationTestsStub):
    pedigree_file_path = PEDIGREE_FILE_PATH_TRIOS

    def test_trios_0_0_create_data_folder(self) -> None:
        self.stub_0_0_create_data_folder()

    def test_trios_0_1_import_data(self) -> None:
        self.stub_0_1_import_data()

    def test_trios_0_2_import_data(self) -> None:
        self.stub_0_2_import_data()

    @GNOMAD_TABLE_SKIP
    def test_trios_0_3_import_data(self) -> None:
        self.stub_0_3_import_data()

    def test_trios_1_1_import_data(self) -> None:
        self.stub_1_1_import_data()

    def test_trios_1_2_import_data(self) -> None:
        self.stub_1_2_import_data()

    def test_trios_1_3_import_data(self) -> None:
        self.stub_1_3_import_data()

    def test_trios_1_4_import_data(self) -> None:
        self.stub_1_4_import_data()

    def test_trios_2_1_sample_qc(self) -> None:
        self.stub_2_1_sample_qc()

    def test_trios_2_2_sample_qc(self, WES_CONFIG: str) -> None:
        self.stub_2_2_sample_qc()
        assert_step_2_2_outputs_match_expected(WES_CONFIG)

    def test_trios_2_3_sample_qc(self) -> None:
        self.stub_2_3_sample_qc()

    def test_trios_2_4_sample_qc(self) -> None:
        self.stub_2_4_sample_qc()

    def test_trios_2_5_sample_qc(self) -> None:
        self.stub_2_5_sample_qc()

    def test_trios_3_1_variant_qc(self) -> None:
        self.stub_3_1_variant_qc()

    def test_trios_3_2_variant_qc(self) -> None:
        self.stub_3_2_variant_qc()

    def test_trios_3_3_variant_qc(self) -> None:
        self.stub_3_3_variant_qc()

    def test_trios_3_4_variant_qc(self) -> None:
        self.stub_3_4_variant_qc()

    def test_trios_3_5_variant_qc(self) -> None:
        self.stub_3_5_variant_qc()

    def test_trios_3_6_variant_qc(self) -> None:
        self.stub_3_6_variant_qc()

    def test_trios_3_7_variant_qc(self) -> None:
        self.stub_3_7_variant_qc()

    def test_trios_3_8_variant_qc(self) -> None:
        self.stub_3_8_variant_qc()

    def test_trios_3_9_variant_qc(self) -> None:
        self.stub_3_9_variant_qc()

    def test_trios_4_1_genotype_qc(self, WES_CONFIG: str) -> None:
        self.stub_4_1_genotype_qc()
        assert_step_4_1_outputs_match_expected(WES_CONFIG)

    def test_trios_4_2_genotype_qc(self) -> None:
        self.stub_4_2_genotype_qc()

    def test_trios_4_3a_genotype_qc(self) -> None:
        self.stub_4_3a_genotype_qc()

    def test_trios_4_3b_genotype_qc(self) -> None:
        self.stub_4_3b_genotype_qc()

    def test_trios_4_4_genotype_qc(self) -> None:
        self.stub_4_4_genotype_qc()

    def test_trios_4_5_genotype_qc(self) -> None:
        self.stub_4_5_genotype_qc()
