#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Merge linear RNA and circRNA expression matrices.

Supports both the merged two-caller circRNA table (with circRNA_ID3/Source) and
single-caller CIRCexplorer2/CIRI2 CPM tables (without circRNA_ID3).
"""

import argparse
import csv
from pathlib import Path
from typing import Dict, Iterable, List

OUTPUT_FIELDS = ["GeneID", "GeneSymbol", "GeneType", "ReadCounts", "Norm_Expr"]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Merge gene and circRNA expression tables.")
    parser.add_argument("--gene_expr", required=True)
    parser.add_argument(
        "--circRNA_expr",
        help="Optional circRNA expression matrix. Omit it to normalize a gene-only matrix.",
    )
    parser.add_argument("--out_matrix", required=True)
    return parser.parse_args()


def _require_file(path: str, label: str) -> None:
    """Raise a clear error if an input file is missing or empty."""
    file_path = Path(path)
    if not file_path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")
    if file_path.stat().st_size == 0:
        raise ValueError(f"{label} is empty: {path}")


def read_gene_data(file_path: str) -> List[Dict[str, str]]:
    """Read a linear gene expression matrix."""
    _require_file(file_path, "Gene expression matrix")
    data: List[Dict[str, str]] = []
    with open(file_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"GeneID", "GeneSymbol", "GeneType", "ReadCounts", "TPM"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Gene expression matrix missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            if row.get("GeneType") == "artifact":
                continue
            data.append(
                {
                    "GeneID": row.get("GeneID", ""),
                    "GeneSymbol": row.get("GeneSymbol", ""),
                    "GeneType": row.get("GeneType", ""),
                    "ReadCounts": row.get("ReadCounts", "0"),
                    "Norm_Expr": row.get("TPM", "0"),
                }
            )
    return data


def _circ_symbol(row: Dict[str, str]) -> str:
    """Choose a stable display symbol for merged or single-caller circRNA rows."""
    for key in ("circRNA_ID3", "circRNA_ID2", "circRNA_ID1"):
        value = (row.get(key) or "").strip()
        if value:
            return value
    return "circRNA"


def read_circ_data(file_path: str) -> List[Dict[str, str]]:
    """Read circRNA CPM data from merged or single-caller tables."""
    _require_file(file_path, "circRNA expression matrix")
    data: List[Dict[str, str]] = []
    with open(file_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"circRNA_ID1", "junction_read_counts", "CPM"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"circRNA matrix missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            circ_id = (row.get("circRNA_ID1") or "").strip()
            if not circ_id:
                continue
            data.append(
                {
                    "GeneID": circ_id,
                    "GeneSymbol": _circ_symbol(row),
                    "GeneType": "circRNAs",
                    "ReadCounts": row.get("junction_read_counts", "0"),
                    "Norm_Expr": row.get("CPM", "0"),
                }
            )
    return data


def write_merged_matrix(rows: Iterable[Dict[str, str]], output_path: str) -> None:
    """Write the merged output matrix."""
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    gene_data = read_gene_data(args.gene_expr)
    circ_data = read_circ_data(args.circRNA_expr) if args.circRNA_expr else []
    write_merged_matrix(gene_data + circ_data, args.out_matrix)


if __name__ == "__main__":
    main()
