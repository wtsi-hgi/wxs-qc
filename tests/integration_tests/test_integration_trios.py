import csv
import gzip
import json
import math
from pathlib import Path
from typing import Any, cast

import hail as hl
import pytest

from tests.integration_tests.integration_stub import IntegrationTestsStub
from wes_qc.config import parse_config_file
from wes_qc.hail_utils import path_spark

# /path/to/wes_qc must be in PYTHONPATH

PEDIGREE_FILE_PATH_TRIOS = """ '${cvars.metadir}/control_set_small.trios.fam' """

GNOMAD_TABLE_SKIP = pytest.mark.skip(
    reason="The test depends on gnomAD table, not present in downloaded test data.\n"
    "Instead we download the resulting table from the bucket and use it as a resource.\n"
)


def _read_related_samples(path: str) -> set[str]:
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        return {json.loads(row["node"])["s"] for row in reader}


def _read_relatedness(path: str) -> dict[tuple[str, str], tuple[float, float]]:
    with gzip.open(path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return {
            (json.loads(row["i"])["s"], json.loads(row["j"])["s"]): (float(row["kin"]), float(row["ibd2"]))
            for row in reader
        }


def _load_expected_results(suite: str, step: str) -> dict[str, Any]:
    expected_results_path = Path(__file__).with_name("expected_integration_test_results.json")
    with open(expected_results_path) as f:
        expected_results = json.load(f)
    return cast(dict[str, Any], expected_results[suite][step])


def _expected_relatedness(records: list[dict[str, Any]]) -> dict[tuple[str, str], tuple[float, float]]:
    return {(record["i"], record["j"]): (record["kin"], record["ibd2"]) for record in records}


def _assert_relatedness_matches(actual_path: str, expected_records: list[dict[str, Any]]) -> None:
    actual_relatedness = _read_relatedness(actual_path)
    print(f"== VALIDATE: Actual relatedness results: {actual_relatedness} ==")
    expected_relatedness = _expected_relatedness(expected_records)

    assert set(actual_relatedness) == set(expected_relatedness)
    for pair, actual_values in actual_relatedness.items():
        expected_values = expected_relatedness[pair]
        assert math.isclose(actual_values[0], expected_values[0], rel_tol=1e-4, abs_tol=1e-6)
        assert math.isclose(actual_values[1], expected_values[1], rel_tol=1e-4, abs_tol=1e-6)


def _assert_hail_table_count_matches(actual_path: str, expected: dict[str, Any]) -> None:
    actual = hl.read_table(path_spark(actual_path))
    print(f"== VALIDATE: Actual Hail table size: {actual.count()} ==")
    assert actual.count() == expected["count"]


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

    related_samples_path = config["step2"]["relatedness_output"]["samples_to_remove_tsv"]
    assert _read_related_samples(related_samples_path) == set(expected["related_samples_to_remove"])
    _assert_relatedness_matches(
        config["step2"]["relatedness_output"]["relatedness_outfile"],
        expected["relatedness"],
    )

    actual_relatedness_output = config["step2"]["relatedness_output"]
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

    def test_trios_4_1_genotype_qc(self) -> None:
        self.stub_4_1_genotype_qc()

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
