#!/usr/bin/env python3

"""
Builds the Sarek samplesheet CSV for the fixed testdata.consts sample set.

Unlike artifacts/public-dataset-build/2_make_sarek_manifest.py, the sample
list is not passed on the command line - it is built directly from the fixed
testdata.consts SAMPLES and GIAB_SAMPLES constants (see consts_loader.get_samples()/
get_giab_samples()), and the manifest is always written to testdata.consts'
SAREK_SAMPLESHEET path.

The two sample sources use different aligned-reads formats (1000 Genomes: CRAM+CRAI;
GIAB: BAM+BAI, from 1a_download_giab_trios), so the samplesheet carries both column pairs,
leaving the unused pair blank per row - nf-core/sarek 3.10.0's samplesheet schema
(assets/schema_input.json) defines cram/crai/bam/bai as independent optional columns with no
exclusivity constraint between them, so a single file can mix both row types.
"""

import csv
from pathlib import Path

from consts_loader import get_consts, get_giab_samples, get_samples


def main() -> None:
    consts = get_consts()
    samples = get_samples()
    giab_samples = get_giab_samples()
    output_filename = consts["SAREK_SAMPLESHEET"]
    Path(output_filename).parent.mkdir(parents=True, exist_ok=True)

    with open(output_filename, mode="w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["patient", "sample", "cram", "crai", "bam", "bai"])  # Header

        for sample in samples:
            sample_id = sample["sample_id"]
            full_name = f"{sample_id}.alt_bwamem_GRCh38DH.20150826.{sample['population']}.exome"
            cram = f"{consts['CRAMS_DIR']}/{full_name}.cram"
            crai = f"{cram}.crai"
            writer.writerow([sample_id, sample_id, cram, crai, "", ""])

        for sample in giab_samples:
            sample_id = sample["sample_id"]
            bam = f"{consts['GIAB_BAMS_DIR']}/{sample_id}.exome.bam"
            bai = f"{bam}.bai"
            writer.writerow([sample_id, sample_id, "", "", bam, bai])

    print(f"Manifest written to {output_filename}")


if __name__ == "__main__":
    main()
