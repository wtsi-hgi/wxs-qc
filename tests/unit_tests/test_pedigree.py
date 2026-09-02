"""
Unit tests for wxs_qc.pedigree.

Self-contained: pedigrees are built from inline .fam text, so the tests need neither
an initialized Hail backend nor the downloaded test dataset.
`hl.Pedigree.read()` is not used because it requires the Hail/Spark filesystem.
"""

import hail as hl

from wxs_qc.pedigree import collect_pedigree_samples, select_founders

# Both fixtures mirror the test dataset pedigrees
# (tests/data/test_source_data/control_set_small_v2/metadata/control_set_small_v2.trios*.fam):
# the same four trios, recorded with and without records for the parents.
FAM_PROBANDS_ONLY = """
CEU1\tNA12878\tNA12891\tNA12892\t2\t0
YRI1\tNA19240\tNA19239\tNA19238\t2\t0
AshkenazimTrio\tHG002\tHG003\tHG004\t1\t0
ChineseTrio\tHG005\tHG006\tHG007\t1\t0
"""

FAM_WITH_PARENTS = """
CEU1\tNA12892\t0\t0\t2\t0
CEU1\tNA12891\t0\t0\t1\t0
CEU1\tNA12878\tNA12891\tNA12892\t2\t0
YRI1\tNA19238\t0\t0\t2\t0
YRI1\tNA19239\t0\t0\t1\t0
YRI1\tNA19240\tNA19239\tNA19238\t2\t0
AshkenazimTrio\tHG004\t0\t0\t2\t0
AshkenazimTrio\tHG003\t0\t0\t1\t0
AshkenazimTrio\tHG002\tHG003\tHG004\t1\t0
ChineseTrio\tHG007\t0\t0\t2\t0
ChineseTrio\tHG006\t0\t0\t1\t0
ChineseTrio\tHG005\tHG006\tHG007\t1\t0
"""

EXPECTED_FOUNDERS = {"NA12891", "NA12892", "NA19238", "NA19239", "HG003", "HG004", "HG006", "HG007"}
EXPECTED_SAMPLES = EXPECTED_FOUNDERS | {"NA12878", "NA19240", "HG002", "HG005"}


def pedigree_from_fam(fam_text: str) -> hl.Pedigree:
    """
    Build a pedigree from PLINK .fam text the same way `hl.Pedigree.read()` does:
    one trio per record, with the `0` placeholder of the parent columns mapped to `None`.
    """
    trios = []
    for line in fam_text.strip().splitlines():
        fam_id, sample_id, pat_id, mat_id, sex = line.split("\t")[:5]
        trios.append(
            hl.Trio(
                s=sample_id,
                fam_id=fam_id,
                pat_id=None if pat_id == "0" else pat_id,
                mat_id=None if mat_id == "0" else mat_id,
                is_female={"1": False, "2": True}.get(sex),
            )
        )
    return hl.Pedigree(trios)


def test_pedigree_collect_samples_probands_only() -> None:
    assert collect_pedigree_samples(pedigree_from_fam(FAM_PROBANDS_ONLY)) == EXPECTED_SAMPLES


def test_pedigree_collect_samples_with_parent_records() -> None:
    assert collect_pedigree_samples(pedigree_from_fam(FAM_WITH_PARENTS)) == EXPECTED_SAMPLES


def test_pedigree_select_founders_probands_only() -> None:
    """Parents named in trios but having no record of their own are founders."""
    assert select_founders(pedigree_from_fam(FAM_PROBANDS_ONLY)) == EXPECTED_FOUNDERS


def test_pedigree_select_founders_with_parent_records() -> None:
    """
    The parentless records of the parents must not remove them from the founders:
    both pedigree layouts describe the same families and must give the same founders.
    """
    assert select_founders(pedigree_from_fam(FAM_WITH_PARENTS)) == EXPECTED_FOUNDERS


def test_pedigree_select_founders_multigenerational() -> None:
    """Only the oldest generation without parent records is a founder."""
    fam_text = (
        "fam1\tgrandfather\t0\t0\t1\t0\n"
        "fam1\tgrandmother\t0\t0\t2\t0\n"
        "fam1\tfather\tgrandfather\tgrandmother\t1\t0\n"
        "fam1\tmother\t0\t0\t2\t0\n"
        "fam1\tchild\tfather\tmother\t2\t0\n"
    )
    assert select_founders(pedigree_from_fam(fam_text)) == {"grandfather", "grandmother", "mother"}


def test_pedigree_select_founders_single_parent_recorded() -> None:
    """A sample with only one known parent is still not a founder."""
    fam_text = "fam1\tmother\t0\t0\t2\t0\nfam1\tchild\t0\tmother\t1\t0\n"
    assert select_founders(pedigree_from_fam(fam_text)) == {"mother"}


def test_pedigree_select_founders_empty_pedigree() -> None:
    assert select_founders(hl.Pedigree([])) == set()
