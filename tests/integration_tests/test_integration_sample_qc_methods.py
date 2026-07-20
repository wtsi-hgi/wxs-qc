from pathlib import Path

import pytest

from tests.integration_tests.integration_stub import IntegrationTestsStub
from wxs_qc.config import parse_config_file

PEDIGREE_FILE_PATH_TRIOS = """ '${cvars.metadir}/control_set_small.trios.fam' """


def assert_sample_qc_outputs(
    config_path: str,
    sample_qc_method: str,
    output_keys: list[str],
) -> None:
    config = parse_config_file(config_path)
    assert config["general"]["sample_qc_method"] == sample_qc_method
    assert config["stage2"]["stratified_sample_qc"]["sample_qc_method"] == sample_qc_method

    for output_key in output_keys:
        output_path = Path(config["stage2"]["stratified_sample_qc"][output_key])
        assert output_path.exists(), f"Expected {sample_qc_method} output does not exist: {output_path}"
        assert output_path.stat().st_size > 0, f"Expected {sample_qc_method} output is empty: {output_path}"


@pytest.mark.usefixtures("WES_CONFIG")
class TestIntegrationSampleQcMethodNn(IntegrationTestsStub):
    pedigree_file_path = PEDIGREE_FILE_PATH_TRIOS
    sample_qc_method = "nn"

    def test_sample_qc_method_nn_2_4(self, WES_CONFIG: str) -> None:
        self.stub_2_4_sample_qc()
        assert_sample_qc_outputs(
            WES_CONFIG,
            "nn",
            [
                "output_text_file",
                "output_stratified_metrics_json_file",
                "output_nn_file",
            ],
        )


@pytest.mark.usefixtures("WES_CONFIG")
class TestIntegrationSampleQcMethodPop(IntegrationTestsStub):
    pedigree_file_path = PEDIGREE_FILE_PATH_TRIOS
    sample_qc_method = "pop"

    def test_sample_qc_method_pop_2_4(self, WES_CONFIG: str) -> None:
        self.stub_2_4_sample_qc()
        assert_sample_qc_outputs(
            WES_CONFIG,
            "pop",
            [
                "output_text_file",
                "output_stratified_metrics_json_file",
            ],
        )
