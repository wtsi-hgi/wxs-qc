from pathlib import Path

import pytest

from tests.integration_tests.integration_stub import (
    INTEGRATION_TESTS_CONFIG_RENDERED_SAVEFILE,
    INTEGRATION_TESTS_CONFIG_TEMPLATE,
    TEST_FILES_LIST,
    render_config,
)
from wes_qc import teszt


@pytest.fixture(scope="session")
def integration_tests_dir() -> Path:
    return Path(__file__).parent.resolve()


@pytest.fixture(scope="session")
def integration_tests_data_dir(integration_tests_dir: Path) -> Path:
    return integration_tests_dir.parent / "data"


@pytest.fixture(scope="session")
def test_data_dir(integration_tests_dir: Path, integration_tests_data_dir: Path) -> Path:
    test_data_download_path = integration_tests_data_dir / "test_source_data"
    test_files_list = (integration_tests_dir / TEST_FILES_LIST).resolve()

    print(f"Downloading data from the bucket using files list {test_files_list}")
    teszt.download_test_data_using_files_list(str(test_files_list), str(test_data_download_path))

    return test_data_download_path


@pytest.fixture(scope="class")
def rendered_config(
    request: pytest.FixtureRequest,
    integration_tests_dir: Path,
    integration_tests_data_dir: Path,
    test_data_dir: Path,
    tmp_path_factory: pytest.TempPathFactory,
) -> Path:
    test_class = request.cls
    rendered_config_dir = tmp_path_factory.mktemp(test_class.__name__ if test_class else "integration")
    rendered_config_savefile = rendered_config_dir / INTEGRATION_TESTS_CONFIG_RENDERED_SAVEFILE
    pedigree_file_path = getattr(test_class, "pedigree_file_path", None)
    sample_qc_method = getattr(test_class, "sample_qc_method", "lr")

    render_config(
        str(integration_tests_dir / INTEGRATION_TESTS_CONFIG_TEMPLATE),
        str(integration_tests_data_dir),
        str(test_data_dir / "control_set_small"),
        str(test_data_dir / "resources"),
        str(test_data_dir / "metadata"),
        str(test_data_dir / "training_sets"),
        str(test_data_dir / "variant_qc_random_forest"),
        pedigree_file_name=pedigree_file_path,
        sample_qc_method=sample_qc_method,
        savefile=str(rendered_config_savefile),
    )

    return rendered_config_savefile


@pytest.fixture
def WES_CONFIG(monkeypatch: pytest.MonkeyPatch, rendered_config: Path) -> str:
    config_path = str(rendered_config)
    monkeypatch.setenv("WES_CONFIG", config_path)
    return config_path
