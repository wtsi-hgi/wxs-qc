import argparse
import os
from typing import Any, Optional
from unittest.mock import patch

import pytest
from stage0_resource_preparation import (
    res0_create_project_folder as qc_step_0_0,
    res1_import_1kg as qc_step_0_1,
    res2_generate_truthset_ht as qc_step_0_2,
    res3_prepare_gnomad_ht as qc_step_0_3,
)
from stage1_import_data import (
    imt_4_mutation_spectra_preqc as qc_step_1_4,
    imt1_import_vcfs as qc_step_1_1,
    imt2_import_annotations as qc_step_1_2,
    imt3_validate_gtcheck as qc_step_1_3,
)
from stage2_sample_qc import (
    sqc1_sex_annotation as qc_step_2_1,
    sqc2_identify_related_samples as qc_step_2_2,
    sqc3_pca_population_prediction as qc_step_2_3,
    sqc4_find_outliers as qc_step_2_4,
    sqc5_filter_fail_samples as qc_step_2_5,
)
from stage3_variant_qc import (
    vqc1_split_and_family_annotate as qc_step_3_1,
    vqc2_create_rf_ht as qc_step_3_2,
    vqc3_train_rf_model as qc_step_3_3,
    vqc4_apply_rf_model as qc_step_3_4,
    vqc5_annotate_ht_after_rf as qc_step_3_5,
    vqc6_rank_and_bin as qc_step_3_6,
    vqc7_plot_metrics as qc_step_3_7,
    vqc8_select_thresholds as qc_step_3_8,
    vqc9_filter_mt_manual_vqc as qc_step_3_9,
)
from stage4_genotype_qc import (
    gqc1_compare_hard_filter_combinations as qc_step_4_1,
    gqc2_apply_range_of_hard_filters as qc_step_4_2,
    gqc3a_export_vcfs_range_of_hard_filters as qc_step_4_3a,
    gqc3b_export_vcfs_single_filter as qc_step_4_3b,
    gqc4_mutation_spectra_afterqc as qc_step_4_5,
    gqc5_counts_per_sample as qc_step_4_4,
)

# /path/to/wxs_qc must be in PYTHONPATH

# List with the test files must be located in the test directory in the repo
TEST_FILES_LIST = "../test_files_list_in_bucket.txt"

# variables for test config rendering
INTEGRATION_TESTS_DIR = "{INTEGRATION_TESTS_DIR}"
INTEGRATION_TESTS_DATA_DIR = "{INTEGRATION_TESTS_DATA_DIR}"
TEST_DATA_DIR = "{TEST_DATA_DIR}"
RESOURCES_DIR = "{RESOURCES_DIR}"
METADATA_DIR = "{METADATA_DIR}"
TRAINING_SETS_DIR = "{TRAINING_SETS_DIR}"
VARIANT_QC_RANDOM_FOREST_DIR = "{VARIANT_QC_RANDOM_FOREST_DIR}"
PEDIGREE_FILE_NAME = "{PEDIGREE_FILE_NAME}"
SAMPLE_QC_METHOD = "{SAMPLE_QC_METHOD}"

# configuration file template for the integration tests
INTEGRATION_TESTS_CONFIG_TEMPLATE = "config_test_template.yaml"
INTEGRATION_TESTS_CONFIG_RENDERED_SAVEFILE = "integration_config_rendered.yaml"


# set up explicit runhash parameter for reproducible RF training
RF_RUN_TEST_HASH = "testhash"  # manually set rf run id


def render_config(
    path_to_template: str,
    integration_tests_data_dir: Optional[str],
    test_data_dir: Optional[str],
    resources_dir: Optional[str],
    metadata_dir: Optional[str],
    training_sets_dir: Optional[str],
    variant_qc_random_forest_dir: Optional[str],
    pedigree_file_name: Optional[str],
    sample_qc_method: str = "lr",
    savefile: str = "inputs_test_rendered.yaml",
) -> None:
    """
    Read the config template and fill in the paths.
    """
    integration_tests_dir = os.path.dirname(os.path.abspath(__file__))
    default_integration_tests_data_dir = os.path.join(os.path.dirname(integration_tests_dir), "data")
    integration_tests_data_dir = (
        integration_tests_data_dir if integration_tests_data_dir else default_integration_tests_data_dir
    )

    # TODO agree on test and resources folder naming convention
    test_data_dir = test_data_dir if test_data_dir else os.path.join(integration_tests_dir, "control_set")
    resources_dir = resources_dir if resources_dir else os.path.join(integration_tests_dir, "resources")
    metadata_dir = metadata_dir if metadata_dir else os.path.join(integration_tests_dir, "metadata")
    training_sets_dir = training_sets_dir if training_sets_dir else os.path.join(integration_tests_dir, "training_sets")
    variant_qc_random_forest_dir = (
        variant_qc_random_forest_dir
        if variant_qc_random_forest_dir
        else os.path.join(integration_tests_dir, "variant_qc_random_forest_dir")
    )

    with open(path_to_template, "r") as f:
        template = f.read()

    # TODO: make versatile
    template = template.replace(INTEGRATION_TESTS_DIR, integration_tests_dir)
    template = template.replace(INTEGRATION_TESTS_DATA_DIR, integration_tests_data_dir)
    template = template.replace(TEST_DATA_DIR, test_data_dir)
    template = template.replace(RESOURCES_DIR, resources_dir)
    template = template.replace(METADATA_DIR, metadata_dir)
    template = template.replace(TRAINING_SETS_DIR, training_sets_dir)
    template = template.replace(VARIANT_QC_RANDOM_FOREST_DIR, variant_qc_random_forest_dir)
    pedigree_file_name = pedigree_file_name if pedigree_file_name is not None else "null"
    template = template.replace(PEDIGREE_FILE_NAME, pedigree_file_name)
    template = template.replace(SAMPLE_QC_METHOD, sample_qc_method)

    with open(savefile, "w") as f:
        f.write(template)


# DEBUG: might be useful https://spark.apache.org/docs/latest/api/python/getting_started/testing_pyspark.html
# TODO: clean up hail logs


class IntegrationTestsStub:
    def stub_0_0_create_data_folder(self) -> None:
        try:
            qc_step_0_0.main()
        except Exception as e:
            pytest.fail(f"Step 0.0 failed with an exception: {e}")

    @patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(
            kg_to_mt=True, kg_filter_and_prune=True, kg_pc_relate=True, kg_remove_related_samples=True, all=False
        ),
    )
    def stub_0_1_import_data(self, mock_args: Any) -> None:
        try:
            qc_step_0_1.main()
        except Exception as e:
            pytest.fail(f"Step 0.1 failed with an exception: {e}")

    def stub_0_2_import_data(self) -> None:
        try:
            qc_step_0_2.main()
        except Exception as e:
            pytest.fail(f"Step 0.2 failed with an exception: {e}")

    @pytest.mark.skip(
        reason="The test depends on gnomAD table, not present in downloaded test data.\n"
        "Instead we download the resulting table from the bucket and use it as a resource.\n"
    )
    def stub_0_3_import_data(self) -> None:
        try:
            qc_step_0_3.main()
        except Exception as e:
            pytest.fail(f"Step 0.3 failed with an exception: {e}")

    def stub_1_1_import_data(self) -> None:
        try:
            qc_step_1_1.main()
        except Exception as e:
            pytest.fail(f"Step 1.1 failed with an exception: {e}")

    def stub_1_2_import_data(self) -> None:
        try:
            qc_step_1_2.main()
        except Exception as e:
            pytest.fail(f"Step 1.2 failed with an exception: {e}")

    def stub_1_3_import_data(self) -> None:
        try:
            qc_step_1_3.main()
        except Exception as e:
            pytest.fail(f"Step 1.3 failed with an exception: {e}")

    def stub_1_4_import_data(self) -> None:
        try:
            qc_step_1_4.main()
        except Exception as e:
            pytest.fail(f"Step 1_4 failed with an exception: {e}")

    ### === Sample QC === ###
    def stub_2_1_sample_qc(self) -> None:
        try:
            qc_step_2_1.main()
        except Exception as e:
            pytest.fail(f"Step 2.1 failed with an exception: {e}")

    @patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(filter_mt=True, pc_relate=True, plot_pca=True, all=False),
    )
    def stub_2_2_sample_qc(self, mock_args: Any) -> None:
        try:
            qc_step_2_2.main()
        except Exception as e:
            pytest.fail(f"Step 2.2 failed with an exception: {e}")

    @patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(
            merge_and_ldprune=True, pca=True, pca_plot=True, assign_pops=True, pca_plot_assigned=True, all=False
        ),
    )
    def stub_2_3_sample_qc(self, mock_args: Any) -> None:
        try:
            qc_step_2_3.main()
        except Exception as e:
            pytest.fail(f"Step 2.3 failed with an exception: {e}")

    def stub_2_4_sample_qc(self) -> None:
        try:
            qc_step_2_4.main()
        except Exception as e:
            pytest.fail(f"Step 2.4 failed with an exception: {e}")

    def stub_2_5_sample_qc(self) -> None:
        try:
            qc_step_2_5.main()
        except Exception as e:
            pytest.fail(f"Step 2.5 failed with an exception: {e}")

    ### === Variant QC === ###
    # mock cli arguments
    @patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(all=False, split_qc=True, trios_stats=True, inbreeding=True),
    )
    def stub_3_1_variant_qc(self, mock_args: Any) -> None:
        try:
            qc_step_3_1.main()
        except Exception as e:
            pytest.fail(f"Step 3.1 failed with an exception: {e}")

    def stub_3_2_variant_qc(self) -> None:
        try:
            qc_step_3_2.main()
        except Exception as e:
            pytest.fail(f"Step 3.2 failed with an exception: {e}")

    @patch("argparse.ArgumentParser.parse_args", return_value=argparse.Namespace(manual_model_id=RF_RUN_TEST_HASH))
    def stub_3_3_variant_qc(self, mock_args: Any) -> None:
        try:
            qc_step_3_3.main()
        except Exception as e:
            pytest.fail(f"Step 3.3 failed with an exception: {e}")

    # mock cli arguments
    def stub_3_4_variant_qc(self) -> None:
        try:
            qc_step_3_4.main()
        except Exception as e:
            pytest.fail(f"Step 3.4 failed with an exception: {e}")

    # mock cli arguments
    def stub_3_5_variant_qc(self) -> None:
        try:
            qc_step_3_5.main()
        except Exception as e:
            pytest.fail(f"Step 3.5 failed with an exception: {e}")

    # mock cli arguments
    def stub_3_6_variant_qc(self) -> None:
        try:
            qc_step_3_6.main()
        except Exception as e:
            pytest.fail(f"Step 3.6 failed with an exception: {e}")

    # mock cli arguments
    def stub_3_7_variant_qc(self) -> None:
        try:
            qc_step_3_7.main()
        except Exception as e:
            pytest.fail(f"Step 3.7 failed with an exception: {e}")

    # mock cli arguments
    @patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(snv=92, indel=68),
    )
    def stub_3_8_variant_qc(self, mock_args: Any) -> None:
        try:
            qc_step_3_8.main()
        except Exception as e:
            pytest.fail(f"Step 3.8 failed with an exception: {e}")

    # mock cli arguments
    @patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(snv=84, indel=60),
    )
    def stub_3_9_variant_qc(self, mock_args: Any) -> None:
        try:
            qc_step_3_9.main()
        except Exception as e:
            pytest.fail(f"Step 3.9 failed with an exception: {e}")

    ### === Genotype QC === ###
    @patch(
        "argparse.ArgumentParser.parse_args",
        return_value=argparse.Namespace(prepare=True, evaluate_snv=True, evaluate_indel=True, plot=True, all=False),
    )
    def stub_4_1_genotype_qc(self, mock_args: Any) -> None:
        try:
            qc_step_4_1.main()
        except Exception as e:
            pytest.fail(f"Step 4.1 failed with an exception: {e}")

    def stub_4_2_genotype_qc(self) -> None:
        try:
            qc_step_4_2.main()
        except Exception as e:
            pytest.fail(f"Step 4.2 failed with an exception: {e}")

    def stub_4_3a_genotype_qc(self) -> None:
        try:
            qc_step_4_3a.main()
        except Exception as e:
            pytest.fail(f"Step 4.3a failed with an exception: {e}")

    @patch("argparse.ArgumentParser.parse_args", return_value=argparse.Namespace(filter_level="stringent"))
    def stub_4_3b_genotype_qc(self, mock_args: Any) -> None:
        try:
            qc_step_4_3b.main()
        except Exception as e:
            pytest.fail(f"Step 4.3b failed with an exception: {e}")

    def stub_4_4_genotype_qc(self) -> None:
        try:
            qc_step_4_4.main()
        except Exception as e:
            pytest.fail(f"Step 4.4 failed with an exception: {e}")

    def stub_4_5_genotype_qc(self) -> None:
        try:
            qc_step_4_5.main()
        except Exception as e:
            pytest.fail(f"Step 4.5 failed with an exception: {e}")
