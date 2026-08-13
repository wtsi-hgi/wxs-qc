#!/usr/bin/env python3

"""
Builds self_reported_sex.tsv for the fixed testdata.consts sample set.

Adapted from artifacts/public-dataset-build/4_make_sex_info (a `cut`+`sed`
one-liner) but converted to Python: filters igsr_samples.tsv down to just the
10 fixed samples using structured csv parsing instead of `sed`, since the
generic version had no need to filter by sample.
"""

import csv
from pathlib import Path

from consts_loader import get_consts, get_samples


def main() -> None:
    consts = get_consts()
    samples = get_samples()
    sample_ids = [sample["sample_id"] for sample in samples]

    with open(consts["IGSR_SAMPLES_TSV"], newline="") as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter="\t")
        sex_by_sample = {row["Sample name"]: row["Sex"] for row in reader}

    metadata_dir = Path(consts["METADATA_DIR"])
    metadata_dir.mkdir(parents=True, exist_ok=True)
    output_path = metadata_dir / f"{consts['DATASET_NAME']}.self_reported_sex.tsv"

    with open(output_path, mode="w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(["sample_id", "self_reported_sex"])
        for sample_id in sample_ids:
            writer.writerow([sample_id, sex_by_sample[sample_id]])

    print(f"Sex info written to {output_path}")


if __name__ == "__main__":
    main()
