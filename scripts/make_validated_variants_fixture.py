#!/usr/bin/env python3
"""
Create a sample-level validated TP/FP fixture TSV from an annotated MatrixTable.
"""

import argparse
import os
import shutil
import tempfile
from pathlib import Path

import hail as hl


REQUIRED_ROW_FIELDS = ("TP", "FP")


def to_spark_path(path: str) -> str:
    if "://" in path:
        return path
    return "file://" + str(Path(path).resolve())


def path_local(path: str) -> str:
    if path.startswith("file://"):
        return path[7:]
    return path


def init_hail(tmp_dir: str) -> None:
    if tmp_dir.startswith("file://"):
        shutil.rmtree(path_local(tmp_dir), ignore_errors=True)
        os.makedirs(path_local(tmp_dir), exist_ok=True)
    hl.init(tmp_dir=tmp_dir, idempotent=True)
    hl.default_reference("GRCh38")


def get_options() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract sample-level validated TP/FP fixture rows from an annotated Hail MatrixTable."
    )
    parser.add_argument(
        "mt",
        help="Input annotated MatrixTable path.",
    )
    parser.add_argument(
        "tsv",
        help="Output TSV path.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=100,
        help="Maximum rows to export for each label, TP and FP. Default: 100",
    )
    parser.add_argument(
        "--tmp-dir",
        default="./hail-tmp",
        help="Hail temporary directory.",
    )
    return parser.parse_args()


def require_annotated_tp_fp(mt: hl.MatrixTable) -> None:
    missing_fields = [field for field in REQUIRED_ROW_FIELDS if field not in mt.row]
    if missing_fields:
        missing_str = ", ".join(missing_fields)
        raise ValueError(
            f"Input MatrixTable must already contain row fields {REQUIRED_ROW_FIELDS}; missing: {missing_str}. "
            "Use an annotated MatrixTable or add TP/FP annotations before running this fixture utility."
        )


def collect_label_entries(mt: hl.MatrixTable, label: str, limit: int) -> hl.Table:
    label_expr = mt[label]
    ht = mt.filter_rows(hl.or_else(label_expr, False)).entries()
    ht = ht.filter(hl.is_defined(ht.GT) & ht.GT.is_non_ref())
    ht = ht.select(
        type=label,
        sample=ht.s,
        chr=ht.locus.contig,
        pos=ht.locus.position,
        ref=ht.alleles[0],
        alt=ht.alleles[1],
    )
    ht = ht.key_by("locus", "alleles", "sample")
    ht = ht.order_by(ht.locus.contig, ht.locus.position, ht.alleles[0], ht.alleles[1], ht["sample"])
    return ht.key_by().select("type", "sample", "chr", "pos", "ref", "alt").head(limit)


def export_single_tsv(ht: hl.Table, output_tsv: str) -> None:
    output_path = Path(path_local(output_tsv))
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="validated-variants-fixture-") as temp_dir:
        temp_output = os.path.join(temp_dir, "validated_variants.tsv")
        ht.export(to_spark_path(temp_output))
        temp_output_path = Path(temp_output)

        if temp_output_path.is_file():
            shutil.move(str(temp_output_path), output_path)
            return

        shard_files = sorted(temp_output_path.glob("part-*"))
        if len(shard_files) != 1:
            raise RuntimeError(
                f"Expected exactly one Hail export shard in {temp_output_path}, found {len(shard_files)}"
            )
        shutil.move(str(shard_files[0]), output_path)


def make_validated_variants_fixture(input_mt: str, output_tsv: str, limit: int) -> None:
    if limit < 1:
        raise ValueError("--limit must be at least 1")

    mt = hl.read_matrix_table(to_spark_path(input_mt))
    require_annotated_tp_fp(mt)

    tp_ht = collect_label_entries(mt, "TP", limit)
    fp_ht = collect_label_entries(mt, "FP", limit)
    fixture_ht = tp_ht.union(fp_ht)
    export_single_tsv(fixture_ht, output_tsv)


def main() -> None:
    args = get_options()
    init_hail(args.tmp_dir)
    make_validated_variants_fixture(args.mt, args.tsv, args.limit)
    print(f"Wrote validated TP/FP fixture to {args.tsv}")


if __name__ == "__main__":
    main()
