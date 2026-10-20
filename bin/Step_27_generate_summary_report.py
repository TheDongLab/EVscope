#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Generate a concise, standalone RNA-seq summary report for one EVscope sample."""

from __future__ import annotations

import argparse
import csv
import html
import math
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


NA_VALUES = {"", "na", "n/a", "nan", "none", "null"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a one-paragraph standalone EVscope RNA-seq summary report."
    )
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--qc-summary", help="Step 24 two-column QC summary TSV.")
    parser.add_argument("--expression-matrix", help="Gene expression TSV used for the top-10 list.")
    parser.add_argument(
        "--kit",
        choices=("smarter_total_rna", "smarter_smrna"),
        default="smarter_total_rna",
    )
    parser.add_argument("--strand", choices=("forward", "reverse", "unstrand"), default="unstrand")
    parser.add_argument("--paired-end", choices=("true", "false"), required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def read_qc_metrics(path: Optional[str]) -> Dict[str, str]:
    if not path or not Path(path).is_file():
        return {}
    metrics: Dict[str, str] = {}
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"Metric", "Value"}.issubset(reader.fieldnames):
            return {}
        for row in reader:
            metric = (row.get("Metric") or "").strip()
            if metric:
                metrics[metric] = (row.get("Value") or "").strip()
    return metrics


def number(value: object) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip()
    if text.lower() in NA_VALUES:
        return None
    try:
        result = float(text.replace(",", "").rstrip("%"))
    except ValueError:
        return None
    return result if math.isfinite(result) else None


def metric_number(metrics: Dict[str, str], name: str) -> Optional[float]:
    return number(metrics.get(name))


def format_count(value: Optional[float]) -> str:
    if value is None:
        return "not available"
    absolute = abs(value)
    for scale, label in ((1_000_000_000, "billion"), (1_000_000, "million")):
        if absolute >= scale:
            shown = f"{value / scale:.2f}".rstrip("0").rstrip(".")
            return f"{shown} {label}"
    if float(value).is_integer():
        return f"{int(value):,}"
    return f"{value:,.1f}"


def format_pct(value: Optional[float]) -> str:
    return "not available" if value is None else f"{value:.2f}%"


def fragment_count(read_count: Optional[float], paired_end: bool) -> Optional[float]:
    if read_count is None:
        return None
    return read_count / 2.0 if paired_end else read_count


def join_phrases(items: List[str]) -> str:
    if len(items) < 2:
        return "".join(items)
    if len(items) == 2:
        return " and ".join(items)
    return ", ".join(items[:-1]) + ", and " + items[-1]


def first_present(row: Dict[str, str], names: Iterable[str]) -> str:
    for name in names:
        value = (row.get(name) or "").strip()
        if value and value.lower() not in NA_VALUES:
            return value
    return ""


def top_genes(path: Optional[str], limit: int = 10) -> List[Tuple[str, float]]:
    if not path or not Path(path).is_file():
        return []
    genes: Dict[str, float] = {}
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene_type = first_present(row, ("GeneType", "gene_type", "Biotype")).lower()
            if gene_type in {"circrnas", "circrna", "artifact"}:
                continue
            gene = first_present(row, ("GeneSymbol", "GeneName", "gene_name", "GeneID", "gene_id"))
            read_count = number(first_present(row, ("ReadCounts", "read_count", "expected_count")))
            if not gene or read_count is None:
                continue
            genes[gene] = genes.get(gene, 0.0) + read_count
    return sorted(genes.items(), key=lambda item: (-item[1], item[0]))[:limit]


def expected_orientation(kit: str, strand: str) -> Tuple[str, str]:
    if strand in {"forward", "reverse"}:
        expected = strand
    else:
        expected = "reverse" if kit == "smarter_total_rna" else "forward"
    if kit == "smarter_total_rna" and expected == "reverse":
        label = "expected reverse/dUTP orientation"
    elif kit == "smarter_smrna" and expected == "forward":
        label = "expected forward/sense orientation"
    else:
        label = f"configured {expected} orientation"
    return expected, label


def build_summary(
    sample_name: str,
    metrics: Dict[str, str],
    genes: List[Tuple[str, float]],
    kit: str,
    strand: str,
    paired_end: bool,
) -> str:
    unit = "read pairs" if paired_end else "reads"
    complete_unit = "complete read pairs" if paired_end else "reads"
    clauses: List[str] = []

    raw = metric_number(
        metrics,
        "Total Raw Read Pairs (Fragments)" if paired_end else "Total Raw Reads (R1+R2)",
    )
    trimmed = fragment_count(metric_number(metrics, "Total Trimmed Reads (R1+R2)"), paired_end)
    if raw is not None:
        opening = f"Sample {sample_name} initially yielded {format_count(raw)} {unit}"
        if trimmed is not None:
            opening += f", with {format_count(trimmed)} {unit} remaining after trimming"
        clauses.append(opening + ".")
    elif trimmed is not None:
        clauses.append(f"After trimming, {format_count(trimmed)} {unit} remained for sample {sample_name}.")

    refined_metric = (
        "Number of Reads after UMI-deduplication (STAR Input)"
        if kit == "smarter_total_rna"
        else "Number of Reads in Refined STAR Input (No UMI Deduplication)"
    )
    refined = metric_number(metrics, refined_metric)
    unique = metric_number(metrics, "Number of Uniquely Mapped Reads (STAR)")
    multi = metric_number(metrics, "Number of Multi-mapped Reads (STAR)")
    mapped = None if unique is None and multi is None else (unique or 0.0) + (multi or 0.0)
    mapped_pct = mapped * 100.0 / refined if mapped is not None and refined and refined > 0 else None
    if refined is not None:
        if kit == "smarter_total_rna":
            alignment = (
                f"After UMI deduplication and orphan-mate removal, {format_count(refined)} "
                f"{complete_unit} remained"
            )
        else:
            alignment = f"No UMI deduplication was applied; {format_count(refined)} {unit} entered refined alignment"
        if mapped is not None:
            alignment += f", of which {format_count(mapped)} ({format_pct(mapped_pct)}) mapped to hg38"
        clauses.append(alignment + ".")
    elif mapped is not None:
        clauses.append(f"The refined STAR alignment mapped {format_count(mapped)} {unit} to hg38.")

    if mapped and mapped > 0:
        unique_pct = (unique or 0.0) * 100.0 / mapped
        multi_pct = (multi or 0.0) * 100.0 / mapped
        clauses.append(
            f"Among human-mapped {unit}, {format_pct(unique_pct)} were unique and "
            f"{format_pct(multi_pct)} were multi-mapped."
        )

    contamination_parts: List[str] = []
    for organism, count_metric, pct_metric in (
        (
            "E. coli",
            "Number of Reads Mapped to Escherichia coli (BBSplit)",
            "Percentage of Trimmed Reads Mapped to E. coli (BBSplit)",
        ),
        (
            "Mycoplasma",
            "Number of Reads Mapped to Mycoplasma (BBSplit)",
            "Percentage of Trimmed Reads Mapped to Mycoplasma (BBSplit)",
        ),
    ):
        count = fragment_count(metric_number(metrics, count_metric), paired_end)
        pct = metric_number(metrics, pct_metric)
        if count is not None:
            pct_text = f" ({format_pct(pct)} of trimmed {unit})" if pct is not None else ""
            contamination_parts.append(f"{format_count(count)} {unit}{pct_text} to {organism}")
    if contamination_parts:
        clauses.append("BBSplit assigned " + join_phrases(contamination_parts) + ".")

    region_parts: List[str] = []
    for adjective, region in (("intronic", "Intronic"), ("exonic", "Exonic"), ("intergenic", "Intergenic")):
        count = metric_number(metrics, f"Number of Reads Mapped to {region} Regions")
        pct = metric_number(metrics, f"Percentage of Reads Mapped to {region} Regions")
        if count is not None:
            pct_text = f" ({format_pct(pct)})" if pct is not None else ""
            region_parts.append(f"{format_count(count)} {unit}{pct_text} as {adjective}")
    if region_parts:
        clauses.append(f"Region assignment counted {join_phrases(region_parts)}.")

    forward = metric_number(metrics, "Percentage of Mapped Reads on Forward Strand")
    reverse = metric_number(metrics, "Percentage of Mapped Reads on Reverse Strand")
    informative = (forward or 0.0) + (reverse or 0.0)
    expected, orientation_label = expected_orientation(kit, strand)
    expected_value = forward if expected == "forward" else reverse
    if expected_value is not None and informative > 0:
        expected_pct = expected_value * 100.0 / informative
        clauses.append(
            f"{format_pct(expected_pct)} of strand-informative alignments matched the {orientation_label}."
        )

    if genes:
        gene_text = ", ".join(f"{gene} ({format_count(count)})" for gene, count in genes)
        clauses.append(f"The top {len(genes)} genes by read count were {gene_text}.")
    else:
        clauses.append("A top-gene list was not available because no compatible expression matrix was found.")

    if not clauses:
        return f"Summary metrics were not available for sample {sample_name}."
    return " ".join(clauses)


def render_html(sample_name: str, summary: str) -> str:
    safe_sample = html.escape(sample_name)
    safe_summary = html.escape(summary)
    generated = html.escape(datetime.now().astimezone().strftime("%d %B %Y, %H:%M %Z"))
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{safe_sample} — EVscope RNA-seq Summary</title>
  <style>
    * {{ box-sizing: border-box; }}
    body {{ margin: 0; font-family: Arial, Helvetica, sans-serif; color: #263746; background: #f3f8fc; }}
    main {{ width: 90%; max-width: 1100px; margin: 7vh auto; padding: 2.5rem 3rem; background: #fff;
            border: 1px solid #dce4ec; border-radius: 14px; box-shadow: 0 8px 24px rgba(0,0,0,.08); }}
    h1 {{ margin: 0 0 .35rem; color: #1a5276; font-size: 2rem; }}
    .sample {{ color: #5d6d7e; margin-bottom: 1.6rem; }}
    .summary {{ margin: 0; font-size: 1.12rem; line-height: 1.8; }}
    footer {{ margin-top: 1.8rem; padding-top: 1rem; border-top: 1px solid #e4ebf0; color: #718096; font-size: .82rem; }}
    @media (max-width: 700px) {{ main {{ width: 94%; padding: 1.5rem; }} h1 {{ font-size: 1.55rem; }} }}
  </style>
</head>
<body>
  <main>
    <h1>EVscope RNA-seq Summary</h1>
    <div class="sample">Sample: {safe_sample}</div>
    <p class="summary">{safe_summary}</p>
    <footer>Generated by EVscope on {generated}</footer>
  </main>
</body>
</html>
"""


def main() -> None:
    args = parse_args()
    metrics = read_qc_metrics(args.qc_summary)
    genes = top_genes(args.expression_matrix)
    summary = build_summary(
        sample_name=args.sample_name,
        metrics=metrics,
        genes=genes,
        kit=args.kit,
        strand=args.strand,
        paired_end=args.paired_end == "true",
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(render_html(args.sample_name, summary), encoding="utf-8")
    print(f"Concise summary report generated: {output}")


if __name__ == "__main__":
    main()
