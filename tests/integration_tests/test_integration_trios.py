import pytest

from tests.integration_tests.integration_stub import IntegrationTestsStub

# /path/to/wes_qc must be in PYTHONPATH

PEDIGREE_FILE_PATH_TRIOS = """ '${cvars.metadir}/control_set_small.trios.fam' """

GNOMAD_TABLE_SKIP = pytest.mark.skip(
    reason="The test depends on gnomAD table, not present in downloaded test data.\n"
    "Instead we download the resulting table from the bucket and use it as a resource.\n"
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

    def test_trios_2_2_sample_qc(self) -> None:
        self.stub_2_2_sample_qc()

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
