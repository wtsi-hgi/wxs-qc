#!/usr/bin/env python3

"""
Builds trios.fam and trios-withparents.fam for the fixed testdata.consts trio
set, in the PLINK-fam convention already used by
tests/data/test_source_data/metadata/control_set_small.trios*.fam
(family_id, individual_id, father_id, mother_id, sex[1=male,2=female], phenotype;
tab-delimited).
"""

from pathlib import Path

from consts_loader import get_consts, get_samples, get_trios

_SEX_CODE = {"male": "1", "female": "2"}


def main() -> None:
    consts = get_consts()
    trios = get_trios()
    sex_by_sample = {sample["sample_id"]: sample["sex"] for sample in get_samples()}

    metadata_dir = Path(consts["METADATA_DIR"])
    metadata_dir.mkdir(parents=True, exist_ok=True)

    trios_path = metadata_dir / f"{consts['DATASET_NAME']}.trios.fam"
    withparents_path = metadata_dir / f"{consts['DATASET_NAME']}.trios-withparents.fam"

    with (
        open(trios_path, mode="w", newline="") as trios_file,
        open(withparents_path, mode="w", newline="") as withparents_file,
    ):
        for trio in trios:
            family_id = trio["family_id"]
            child, father, mother = trio["child"], trio["father"], trio["mother"]
            child_sex_code = _SEX_CODE[sex_by_sample[child]]

            trios_file.write(f"{family_id}\t{child}\t{father}\t{mother}\t{child_sex_code}\t0\n")

            # Matches the order used in the existing control_set_small.trios-withparents.fam
            # fixture: mother, then father, then child.
            withparents_file.write(f"{family_id}\t{mother}\t0\t0\t2\t0\n")
            withparents_file.write(f"{family_id}\t{father}\t0\t0\t1\t0\n")
            withparents_file.write(f"{family_id}\t{child}\t{father}\t{mother}\t{child_sex_code}\t0\n")

    print(f"Wrote {trios_path}")
    print(f"Wrote {withparents_path}")


if __name__ == "__main__":
    main()
