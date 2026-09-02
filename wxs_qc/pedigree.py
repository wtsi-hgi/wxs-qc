"""
All functions used to inspect pedigrees (PLINK .fam files loaded as `hl.Pedigree`)
"""

import hail as hl


def collect_pedigree_samples(ped: hl.Pedigree) -> set[str]:
    """
    Collect all sample IDs mentioned anywhere in the pedigree,
    either as a proband or as a parent.
    """
    samples = {getattr(trio, member) for trio in ped.trios for member in ("mat_id", "pat_id", "s")}
    samples.discard(None)
    return samples


def select_founders(ped: hl.Pedigree) -> set[str]:
    """
    Select pedigree founders, i.e. samples that have no parent recorded in the pedigree.

    Works for pedigrees that list only the probands and for pedigrees that also contain
    records for the ancestors. `hl.Pedigree.read()` creates one `hl.Trio` per .fam line,
    so an ancestor with its own record appears both as a parent of its child's trio and
    as the proband of a trio without parents (`pat_id` and `mat_id` are `None`,
    because the .fam parent columns are `0`). Therefore only the records that actually
    name a parent can disqualify a sample from being a founder.

    Samples named as a parent but having no record of their own are founders:
    the pedigree holds no information about their ancestry.
    """
    non_founders = {trio.s for trio in ped.trios if trio.pat_id is not None or trio.mat_id is not None}
    return collect_pedigree_samples(ped) - non_founders
