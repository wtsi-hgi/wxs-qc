#!/usr/bin/env python3

"""
Downloads CRAM+CRAI files for the fixed testdata.consts sample set from the
public 1000 Genomes exome-alignment collection.

Sample name example:
HG00285.alt_bwamem_GRCh38DH.20150826.FIN.exome ->
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/FIN/HG00285/exome_alignment/
HG00285.alt_bwamem_GRCh38DH.20150826.FIN.exome.cram
"""

import subprocess
from pathlib import Path

from consts_loader import get_consts, get_samples

BASE_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data"


def file_to_url(fn: str) -> str:
    fn = fn.strip()
    parts = fn.split(".")
    if len(parts) < 5:
        raise ValueError(f"Unexpected filename format: {fn}")
    sample = parts[0]
    pop = parts[3]
    dtype = parts[4]
    dirname = f"{dtype}_alignment"
    return f"{BASE_URL}/{pop}/{sample}/{dirname}/{fn}.cram"


def download_file(url: str, dest: str) -> None:
    print(f"Downloading {url}")
    try:
        subprocess.run(["wget", "-c", "-P", dest, url], check=True)
    except subprocess.CalledProcessError as e:
        raise ValueError(f"Failed to download {url}. Error: {e}") from e


def main() -> None:
    consts = get_consts()
    crams_dir = consts["CRAMS_DIR"]
    Path(crams_dir).mkdir(parents=True, exist_ok=True)

    for sample in get_samples():
        full_name = f"{sample['sample_id']}.alt_bwamem_GRCh38DH.20150826.{sample['population']}.exome"
        url = file_to_url(full_name)
        download_file(url + ".crai", crams_dir)
        download_file(url, crams_dir)


if __name__ == "__main__":
    main()
