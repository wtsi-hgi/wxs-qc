# mypy: disable-error-code="attr-defined,no-untyped-call,no-untyped-def"

import pytest

from tests.integration_tests.integration_stub import IntegrationTestsStub

# TODO: this approach doesn't work. Leave as is. Find a netter way to parametrize tests: HSH-254
# INTEGRATION_TESTS_CONFIG_TEMPLATE = "config_test_template.yaml"
# INTEGRATION_TESTS_CONFIG_RENDERED_SAVEFILE = "integration_config_rendered_non_trios.yaml"

PEDIGREE_FILE_PATH_NON_TRIOS = "null"

# set up explicit runhash parameter for reproducible RF training
RF_RUN_TEST_HASH = "testhash_non-trios"  # manually set rf run id

GNOMAD_TABLE_SKIP = pytest.mark.skip(
    reason="The test depends on gnomAD table, not present in downloaded test data.\n"
    "Instead we download the resulting table from the bucket and use it as a resource.\n"
)


class IntegrationTestsNonTrios(IntegrationTestsStub):
    @classmethod
    def setUpClass(cls):
        # Store pedigree file path as class variable
        cls.pedigree_file_path = PEDIGREE_FILE_PATH_NON_TRIOS
        # Call parent's setUpClass without arguments
        super().setUpClass()

    def test_non_trios_0_0_create_data_folder(self):
        self.stub_0_0_create_data_folder()

    def test_non_trios_0_1_import_data(self):
        self.stub_0_1_import_data()

    def test_non_trios_0_2_import_data(self):
        self.stub_0_2_import_data()

    @GNOMAD_TABLE_SKIP
    def test_non_trios_0_3_import_data(self):
        self.stub_0_3_import_data()

    def test_non_trios_1_1_import_data(self):
        self.stub_1_1_import_data()

    def test_non_trios_1_2_import_data(self):
        self.stub_1_2_import_data()

    def test_non_trios_1_3_import_data(self):
        self.stub_1_3_import_data()

    def test_non_trios_1_4_import_data(self):
        self.stub_1_4_import_data()

    def test_non_trios_2_1_sample_qc(self):
        self.stub_2_1_sample_qc()

    def test_non_trios_2_2_sample_qc(self):
        self.stub_2_2_sample_qc()

    def test_non_trios_2_3_sample_qc(self):
        self.stub_2_3_sample_qc()

    def test_non_trios_2_4_sample_qc(self):
        self.stub_2_4_sample_qc()

    def test_non_trios_2_5_sample_qc(self):
        self.stub_2_5_sample_qc()

    def test_non_trios_3_1_variant_qc(self):
        self.stub_3_1_variant_qc()

    def test_non_trios_3_2_variant_qc(self):
        self.stub_3_2_variant_qc()

    def test_non_trios_3_3_variant_qc(self):
        self.stub_3_3_variant_qc()

    def test_non_trios_3_4_variant_qc(self):
        self.stub_3_4_variant_qc()

    def test_non_trios_3_5_variant_qc(self):
        self.stub_3_5_variant_qc()

    def test_non_trios_3_6_variant_qc(self):
        self.stub_3_6_variant_qc()

    def test_non_trios_3_7_variant_qc(self):
        self.stub_3_7_variant_qc()

    def test_non_trios_3_8_variant_qc(self):
        self.stub_3_8_variant_qc()

    def test_non_trios_3_9_variant_qc(self):
        self.stub_3_9_variant_qc()

    def test_non_trios_4_1_genotype_qc(self):
        self.stub_4_1_genotype_qc()

    def test_non_trios_4_2_genotype_qc(self):
        self.stub_4_2_genotype_qc()

    def test_non_trios_4_3a_genotype_qc(self):
        self.stub_4_3a_genotype_qc()

    def test_non_trios_4_3b_genotype_qc(self):
        self.stub_4_3b_genotype_qc()

    def test_non_trios_4_4_genotype_qc(self):
        self.stub_4_4_genotype_qc()

    def test_non_trios_4_5_genotype_qc(self):
        self.stub_4_5_genotype_qc()
