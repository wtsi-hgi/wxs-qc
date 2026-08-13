#!/usr/bin/env python3

"""
Builds the Sarek samplesheet CSV for the fixed testdata.consts sample set.

Unlike artifacts/public-dataset-build/2_make_sarek_manifest.py, the sample
list is not passed on the command line - it is built directly from the fixed
testdata.consts SAMPLES constant (see consts_loader.get_samples()), and the
manifest is always written to testdata.consts' SAREK_SAMPLESHEET path.
"""

import csv
from pathlib import Path

from consts_loader import get_consts, get_samples


def main() -> None:
    consts = get_consts()
    samples = get_samples()
    output_filename = consts["SAREK_SAMPLESHEET"]
    Path(output_filename).parent.mkdir(parents=True, exist_ok=True)

    with open(output_filename, mode="w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["patient", "sample", "cram", "crai"])  # Header

        for sample in samples:
            sample_id = sample["sample_id"]
            full_name = f"{sample_id}.alt_bwamem_GRCh38DH.20150826.{sample['population']}.exome"
            cram = f"{consts['CRAMS_DIR']}/{full_name}.cram"
            crai = f"{cram}.crai"
            writer.writerow([sample_id, sample_id, cram, crai])

    print(f"Manifest written to {output_filename}")


if __name__ == "__main__":
    main()
