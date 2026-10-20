#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""
Step_24_generate_QC_matrix.py
"""

import sys
import argparse
import zipfile
import gzip
import os
import pandas as pd
from statistics import mean

def safe_call(func, *args, default=None, **kwargs):
    try:
        return func(*args, **kwargs)
    except Exception:
        return default

def value_or_na(value, available=True):
    """Return NA only when the upstream source is unavailable/failed.

    A real zero is a valid QC value and must remain 0. This matters for
    contamination/rRNA/expression thresholds where 0 means "measured but absent".
    """
    if not available or value is None:
        return "NA"
    try:
        if pd.isna(value):
            return "NA"
    except Exception:
        pass
    return value


def integer_if_integral(value):
    """Return an int for whole-number numeric values; keep non-integers unchanged."""
    try:
        value_float = float(value)
    except (TypeError, ValueError):
        return value
    if value_float.is_integer():
        return int(value_float)
    return value


def sanitize_splice_per_kb_source(source, star_log=None):
    """Return a report-friendly splice/kb source label without filesystem paths."""
    if source is None:
        return "NA"
    try:
        if pd.isna(source):
            return "NA"
    except Exception:
        pass
    source = str(source).strip()
    if source in ("", "NA", "nan", "None"):
        return "NA"
    if source.startswith("STAR refined Log.final.out:"):
        return "STAR refined Log.final.out"
    if star_log and source == star_log:
        return "STAR refined Log.final.out"
    return source


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a comprehensive QC matrix from multiple pipeline output files."
    )
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit("Error: No input parameters provided.")

    parser.add_argument("--raw_fastqc_zips", nargs='+', required=False, help="List of raw FastQC zip files from Step 1.")
    parser.add_argument(
        "--kit",
        choices=("smarter_total_rna", "smarter_smrna"),
        default="smarter_total_rna",
        help="Library kit profile used for kit-aware QC labels (default: smarter_total_rna).",
    )
    parser.add_argument("--trimmed_fastqs", nargs='+', required=False, help="List of trimmed FASTQ files from Step 3.")
    parser.add_argument("--ecoli_fastqs", nargs='*', required=False, help="List of E.coli FASTQ files from BBSplit (Step 5). Empty list means Step 5 ran and found 0 reads.")
    parser.add_argument("--myco_fastqs", nargs='*', required=False, help="List of Mycoplasma FASTQ files from BBSplit (Step 5). Empty list means Step 5 ran and found 0 reads.")
    parser.add_argument("--ribo_fastqs", nargs='*', required=False, help="List of rRNA FASTQ files from RiboDetector (Step 23). Empty list means Step 23 ran and found 0 reads.")
    parser.add_argument("--ACC_motif_fraction", required=False, help="ACC motif fraction TSV file from Step 2 (SMARTer Total RNA only).")
    parser.add_argument("--kraken_report", required=False, help="Kraken2 report file from Step 19.")
    parser.add_argument("--bam2strand_file", required=False, help="bam2strandness output file from Step 7.")
    parser.add_argument("--length_strandness_file", required=False, help="Read-length-stratified strandness TSV from Step 7.")
    parser.add_argument("--picard_insert_file", required=False, help="Picard Insert Size Metrics file from Step 11.")
    parser.add_argument("--picard_rnaseq_file", required=False, help="Picard RNA-Seq Metrics file from Step 11.")
    parser.add_argument("--expression_matrix", required=False, help="Combined gene/circRNA expression matrix from Step 15/16/17.")
    parser.add_argument("--STAR_log_initial", required=False, help="STAR Log.final.out file from the first alignment (Step 4).")
    parser.add_argument("--STAR_log", required=False, help="STAR Log.final.out file from the refined alignment (Step 6).")
    parser.add_argument("--featureCounts_3UTR", required=False, help="FeatureCounts output for 3'UTR regions.")
    parser.add_argument("--featureCounts_5UTR", required=False, help="FeatureCounts output for 5'UTR regions.")
    parser.add_argument("--featureCounts_downstream_2kb", required=False, help="FeatureCounts output for downstream regions.")
    parser.add_argument("--featureCounts_exon", required=False, help="FeatureCounts output for exonic regions.")
    parser.add_argument("--featureCounts_ENCODE_blacklist", required=False, help="FeatureCounts output for ENCODE blacklist regions.")
    parser.add_argument("--featureCounts_intergenic", required=False, help="FeatureCounts output for intergenic regions.")
    parser.add_argument("--featureCounts_intron", required=False, help="FeatureCounts output for intronic regions.")
    parser.add_argument("--featureCounts_promoter_1500_500bp", required=False, help="FeatureCounts output for promoter regions.")
    parser.add_argument("--output", required=True, help="Path to the output QC matrix TSV file.")
    return parser.parse_args()

def extract_fastqc_basic_stats(fastqc_zip):
    stats = {}
    with zipfile.ZipFile(fastqc_zip, "r") as zf:
        data_file_name = [f for f in zf.namelist() if f.endswith("fastqc_data.txt")][0]
        with zf.open(data_file_name) as f:
            lines = f.read().decode("utf-8").splitlines()
    in_basic_stats_module = False
    for line in lines:
        if line.startswith(">>Basic Statistics"):
            in_basic_stats_module = True
            continue
        if in_basic_stats_module and line.startswith(">>END_MODULE"):
            break
        if in_basic_stats_module and not line.startswith("#"):
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                key, value = parts[0], parts[1]
                if key in ["Total Sequences", "%GC", "Sequence length"]:
                    stats[key] = value
    return stats

def parse_sequence_length(val_str):
    if "-" in val_str:
        start, end = val_str.split("-")
        return (float(start) + float(end)) / 2
    return float(val_str)

def extract_acc_motif_fraction(acc_file):
    df = pd.read_csv(acc_file, sep="\t")
    return df["fraction_ACC"].iloc[0]

def extract_kraken_percentages(kraken_file):
    perc_human, perc_bacteria = 0.0, 0.0
    with open(kraken_file) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 6: continue
            percentage = float(fields[0].strip())
            taxon_name = fields[5].strip()
            if taxon_name == "Homo sapiens":
                perc_human = percentage
            elif taxon_name == "Bacteria":
                perc_bacteria = percentage
    return perc_human, perc_bacteria

def count_fastq_reads(fastq_file):
    line_count = 0
    opener = gzip.open if fastq_file.endswith(".gz") else open
    with opener(fastq_file, "rt") as f:
        for _ in f:
            line_count += 1
    return line_count // 4

def compute_average_read_length(fastq_file):
    total_bases, read_count = 0, 0
    line_num = 0
    opener = gzip.open if fastq_file.endswith(".gz") else open
    with opener(fastq_file, "rt") as f:
        for line in f:
            line_num += 1
            if line_num % 4 == 2:
                total_bases += len(line.strip())
                read_count += 1
    return total_bases / read_count if read_count > 0 else 0

def parse_bam2strand(strand_file):
    df = pd.read_csv(strand_file, sep="\t")
    if df.empty:
        return 0, 0, 0, "NA", "NA"
    row = df.iloc[0]

    def get_value(column_name, fallback_index=None, default="NA"):
        if column_name in df.columns:
            return row[column_name]
        if fallback_index is not None and len(row) > fallback_index:
            return row.iloc[fallback_index]
        return default

    forward = safe_call(float, get_value("'1++,1--,2+-,2-+ (forward)'", 2, 0), default=0.0)
    reverse = safe_call(float, get_value("'1+-,1-+,2++,2-- (reverse)'", 3, 0), default=0.0)
    failed = safe_call(float, get_value("%reads failed to determine:", 4, 0), default=0.0)

    splice_raw = get_value("Splice_per_kb", 5, "NA")
    splice_per_kb = safe_call(float, splice_raw, default="NA")
    if splice_per_kb != "NA" and pd.isna(splice_per_kb):
        splice_per_kb = "NA"
    source = str(get_value("splice_per_kb_source", None, "NA"))
    if source in ("", "nan", "None", "NA") or pd.isna(source):
        source = "NA"

    return round(forward * 100, 2), round(reverse * 100, 2), round(failed * 100, 2), splice_per_kb, source

def parse_length_strandness(length_file):
    """Parse Step 7 read-length-stratified strandness TSV.

    Supports the new threshold_nt schema and falls back to the older label-only
    schema used during development.
    """
    df = pd.read_csv(length_file, sep="\t")
    required = {"n_evaluated", "pct_forward_all", "pct_reverse_all", "pct_undetermined_all"}
    if df.empty or not required.issubset(df.columns):
        return []

    rows = []
    if "threshold_nt" in df.columns:
        for _, row in df.iterrows():
            threshold = safe_call(float, row.get("threshold_nt"), default=None)
            if threshold is None:
                continue
            rows.append((threshold, row))
        rows.sort(key=lambda item: item[0])
    else:
        preferred_order = [
            "All evaluated", "Len < 50", "Len >= 50", "Len >= 100", "Len >= 150",
            "Len >= 200", "Len >= 250", "Len >= 300", "Len >= 350", "Len >= 400",
            "Len >= 450", "Len >= 500",
        ]
        order_map = {label: i for i, label in enumerate(preferred_order)}
        for _, row in df.iterrows():
            label = str(row.get("length_label", "")).strip()
            if not label or label == "nan":
                continue
            rows.append((order_map.get(label, len(order_map)), row))
        rows.sort(key=lambda item: item[0])
    return [row for _, row in rows]


def length_row_label(row):
    """Return stable metric label for a length-stratified row."""
    threshold = _row_value(row, "threshold_nt")
    if threshold != "NA":
        try:
            threshold_float = float(threshold)
            threshold_int = int(threshold_float)
            if threshold_float == threshold_int:
                return f"Minimum aligned query length >= {threshold_int} nt"
            return f"Minimum aligned query length >= {threshold_float:g} nt"
        except (TypeError, ValueError):
            pass
    label = str(_row_value(row, "length_label"))
    if label == "All evaluated":
        return "Minimum aligned query length >= 0 nt"
    return label


def _row_value(row, column, default="NA"):
    if column not in row:
        return default
    value = row[column]
    try:
        if pd.isna(value):
            return default
    except Exception:
        pass
    return value


def _row_int(row, column):
    value = _row_value(row, column)
    if value == "NA":
        return "NA"
    return integer_if_integral(value)


def _row_pct(row, column):
    value = _row_value(row, column)
    if value == "NA":
        return "NA"
    try:
        return round(float(value), 2)
    except (TypeError, ValueError):
        return "NA"

def parse_picard_metrics(file_path):
    header, data_line = None, None
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"): continue
            if header is None:
                header = line.split("\t")
            else:
                data_line = line.split("\t")
                break
    return dict(zip(header, data_line)) if header and data_line else {}

def process_expression_matrix(expr_file):
    expr_data = {}
    df = pd.read_csv(expr_file, sep="\t")
    required = {"GeneType", "ReadCounts"}
    if not required.issubset(df.columns):
        return expr_data
    expr_col = None
    for candidate in ("Norm_Expr", "TPM", "CPM"):
        if candidate in df.columns:
            expr_col = candidate
            break
    if expr_col is None:
        return expr_data

    for _, row in df.iterrows():
        gene_type = str(row.get("GeneType", "")).strip()
        if not gene_type or gene_type == "nan":
            continue
        read_count = safe_call(float, row.get("ReadCounts"), default=0.0)
        tpm_cpm = safe_call(float, row.get(expr_col), default=0.0)
        expr_data.setdefault(gene_type, []).append((read_count, tpm_cpm))
    return expr_data

def parse_star_log(star_log_file):
    metrics = {}
    with open(star_log_file) as f:
        for line in f:
            if "|" in line:
                key, value = [item.strip() for item in line.split("|")]
                metrics[key] = value.replace('%', '')
    result = {
        "Number of input reads (STAR)": metrics.get("Number of input reads", "0"),
        "Average input read length (STAR)": metrics.get("Average input read length", "0"),
        "Average mapped length (STAR)": metrics.get("Average mapped length", "0"),
        "Uniquely mapped reads number (STAR)": metrics.get("Uniquely mapped reads number", "0"),
        "Multi-mapping reads number (STAR)": int(float(metrics.get("Number of reads mapped to multiple loci", 0)) + float(metrics.get("Number of reads mapped to too many loci", 0))),
        "Number of splices from uniquely mapped reads (STAR)": metrics.get("Number of splices: Total", "0"),
        "Unmapped reads number (STAR)": int(float(metrics.get("Number of reads unmapped: too many mismatches", 0)) + float(metrics.get("Number of reads unmapped: too short", 0)) + float(metrics.get("Number of reads unmapped: other", 0))),
        "Number of chimeric reads (STAR)": metrics.get("Number of chimeric reads", "0"),
        "Mismatch rate per base (STAR)": f"{metrics.get('Mismatch rate per base, %', '0')}%"
    }
    return result

def parse_star_log_initial(star_log_file):
    metrics = {}
    with open(star_log_file) as f:
        for line in f:
            if "|" in line:
                key, value = [item.strip() for item in line.split("|")]
                metrics[key] = value.replace('%', '')
    unique_reads = int(metrics.get("Uniquely mapped reads number", 0))
    multi_reads = int(float(metrics.get("Number of reads mapped to multiple loci", 0)) + float(metrics.get("Number of reads mapped to too many loci", 0)))
    total_mapped = unique_reads + multi_reads
    return total_mapped

def compute_splice_per_kb_from_star_log(star_log_file):
    metrics = {}
    with open(star_log_file) as f:
        for line in f:
            if "|" in line:
                key, value = [item.strip() for item in line.split("|", 1)]
                metrics[key] = value.replace('%', '')
    n_junctions = safe_call(float, metrics.get("Number of splices: Total"), default=0.0)
    n_unique = safe_call(float, metrics.get("Uniquely mapped reads number"), default=0.0)
    avg_len = safe_call(float, metrics.get("Average mapped length"), default=0.0)
    if n_unique <= 0 or avg_len <= 0:
        return "NA", "NA"
    return round(n_junctions * 1000.0 / (n_unique * avg_len), 4), star_log_file

def process_featureCounts_file(file_path):
    df = pd.read_csv(file_path, sep='\t', comment='#')
    return df.iloc[:, -1].sum()

def main():
    args = parse_args()
    ordered_metrics = [("Library Kit", args.kit)]

    total_raw_reads, raw_gcs, raw_lengths = 0, [], []
    if args.raw_fastqc_zips:
        for zip_file in args.raw_fastqc_zips:
            stats = safe_call(extract_fastqc_basic_stats, zip_file, default={})
            total_raw_reads += int(stats.get("Total Sequences", 0))
            raw_gcs.append(int(stats.get("%GC", 0)))
            raw_lengths.append(safe_call(parse_sequence_length, stats.get("Sequence length", "0"), default=0))

    num_fastqc_files = len(args.raw_fastqc_zips) if args.raw_fastqc_zips else 0
    total_raw_read_pairs = 0
    if num_fastqc_files >= 2:
        total_raw_read_pairs = total_raw_reads / 2
    elif num_fastqc_files == 1:
        total_raw_read_pairs = total_raw_reads
        
    avg_raw_gc = mean(raw_gcs) if raw_gcs else 0
    avg_raw_read_length = mean(raw_lengths) if raw_lengths else 0
    acc_fraction = "NA"
    if args.kit == "smarter_total_rna":
        acc_fraction = safe_call(extract_acc_motif_fraction, args.ACC_motif_fraction, default="NA")

    total_trimmed_reads = 0
    avg_trimmed_read_length = 0
    if args.trimmed_fastqs:
        total_trimmed_reads = sum(safe_call(count_fastq_reads, f, default=0) for f in args.trimmed_fastqs)
        trimmed_lengths = [safe_call(compute_average_read_length, f, default=0) for f in args.trimmed_fastqs]
        if trimmed_lengths:
            avg_trimmed_read_length = mean(trimmed_lengths)

    percent_reads_after_trimming = round((total_trimmed_reads / total_raw_reads) * 100, 2) if total_raw_reads > 0 else 0

    initial_star_mapped_fragments = safe_call(parse_star_log_initial, args.STAR_log_initial, default=0)

    ordered_metrics.append(("Total Raw Reads (R1+R2)", int(total_raw_reads) if total_raw_reads > 0 else "NA"))
    ordered_metrics.append(("Total Raw Read Pairs (Fragments)", int(total_raw_read_pairs) if total_raw_read_pairs > 0 else "NA"))
    ordered_metrics.append(("%GC of Raw Reads", round(avg_raw_gc, 2) if avg_raw_gc > 0 else "NA"))
    ordered_metrics.append(("Average Raw Read Length", int(round(avg_raw_read_length)) if avg_raw_read_length > 0 else "NA"))
    ordered_metrics.append(("ACC motif fraction from UMI region (SMARTer Total RNA only)", acc_fraction))
    ordered_metrics.append(("Total Trimmed Reads (R1+R2)", int(total_trimmed_reads) if total_trimmed_reads > 0 else "NA"))
    ordered_metrics.append(("Average Trimmed Read Length", int(round(avg_trimmed_read_length)) if avg_trimmed_read_length > 0 else "NA"))
    ordered_metrics.append(("Percentage of Reads Remaining after Trimming", percent_reads_after_trimming if total_raw_reads > 0 else "NA"))

    perc_human, perc_bacteria = safe_call(extract_kraken_percentages, args.kraken_report, default=(0.0, 0.0))
    
    ecoli_available = args.ecoli_fastqs is not None
    ecoli_reads = sum(safe_call(count_fastq_reads, f, default=0) for f in args.ecoli_fastqs or [])
    pct_ecoli = round((ecoli_reads / total_trimmed_reads) * 100, 2) if total_trimmed_reads > 0 else 0
    
    myco_available = args.myco_fastqs is not None
    myco_reads = sum(safe_call(count_fastq_reads, f, default=0) for f in args.myco_fastqs or [])
    pct_myco = round((myco_reads / total_trimmed_reads) * 100, 2) if total_trimmed_reads > 0 else 0

    ribo_available = args.ribo_fastqs is not None
    ribo_reads = sum(safe_call(count_fastq_reads, f, default=0) for f in args.ribo_fastqs or [])
    perc_rRNA_ribodetector = round((ribo_reads / total_trimmed_reads) * 100, 2) if total_trimmed_reads > 0 else 0

    kraken_available = bool(args.kraken_report and os.path.exists(args.kraken_report))
    ordered_metrics.append(("Percentage of Reads Mapped to Human (Kraken)", value_or_na(perc_human, kraken_available)))
    ordered_metrics.append(("Percentage of Reads Mapped to Bacteria (Kraken)", value_or_na(perc_bacteria, kraken_available)))
    ordered_metrics.append(("Number of Reads Mapped to Escherichia coli (BBSplit)", value_or_na(int(ecoli_reads), ecoli_available)))
    ordered_metrics.append(("Percentage of Trimmed Reads Mapped to E. coli (BBSplit)", value_or_na(pct_ecoli, ecoli_available and total_trimmed_reads > 0)))
    ordered_metrics.append(("Number of Reads Mapped to Mycoplasma (BBSplit)", value_or_na(int(myco_reads), myco_available)))
    ordered_metrics.append(("Percentage of Trimmed Reads Mapped to Mycoplasma (BBSplit)", value_or_na(pct_myco, myco_available and total_trimmed_reads > 0)))
    ordered_metrics.append(("Percentage of Trimmed Reads Mapped to rRNAs (RiboDetector)", value_or_na(perc_rRNA_ribodetector, ribo_available and total_trimmed_reads > 0)))

    star_metrics = safe_call(parse_star_log, args.STAR_log, default={})
    star_input_reads = int(star_metrics.get("Number of input reads (STAR)", 0))
    
    perc_reads_after_dedup = "NA"
    if args.kit == "smarter_total_rna" and initial_star_mapped_fragments > 0:
        perc_reads_after_dedup = round((star_input_reads / initial_star_mapped_fragments) * 100, 2)
    
    star_unique_reads = int(star_metrics.get("Uniquely mapped reads number (STAR)", 0))
    perc_unique_vs_star_input = round((star_unique_reads / star_input_reads) * 100, 2) if star_input_reads > 0 else 0
    
    star_multi_reads = int(star_metrics.get("Multi-mapping reads number (STAR)", 0))
    perc_multi_vs_star_input = round((star_multi_reads / star_input_reads) * 100, 2) if star_input_reads > 0 else 0
    
    strand_available = bool(args.bam2strand_file and os.path.exists(args.bam2strand_file))
    (forward_strand, reverse_strand, failed_strand, splice_per_kb, splice_per_kb_source) = safe_call(
        parse_bam2strand, args.bam2strand_file, default=(0, 0, 0, "NA", "NA")
    )
    if (not strand_available or splice_per_kb == "NA") and args.STAR_log and os.path.exists(args.STAR_log):
        splice_per_kb, splice_per_kb_source = safe_call(
            compute_splice_per_kb_from_star_log, args.STAR_log, default=("NA", "NA")
        )
    elif splice_per_kb != "NA" and splice_per_kb_source == "STAR_Log.final.out" and args.STAR_log:
        splice_per_kb_source = args.STAR_log
    splice_per_kb_source = sanitize_splice_per_kb_source(splice_per_kb_source, args.STAR_log)
    splice_available = splice_per_kb != "NA"

    if args.kit == "smarter_total_rna":
        refined_input_metrics = [
            ("Number of Reads after UMI-deduplication (STAR Input)", star_input_reads if star_input_reads > 0 else "NA"),
            ("Percentage of UMI-dedup Fragments", perc_reads_after_dedup),
        ]
    else:
        refined_input_metrics = [
            ("Number of Reads in Refined STAR Input (No UMI Deduplication)", star_input_reads if star_input_reads > 0 else "NA"),
            ("Percentage of UMI-dedup Fragments (Not Applicable: No UMI Deduplication)", "NA"),
        ]

    ordered_metrics.extend([
        ("Total Fragments Mapped to Human (First STAR)", int(initial_star_mapped_fragments) if initial_star_mapped_fragments > 0 else "NA"),
        *refined_input_metrics,
        ("Average Input Read Length (STAR)", float(star_metrics.get("Average input read length (STAR)", 0))),
        ("Average Mapped Length per Unique Fragment (STAR)", float(star_metrics.get("Average mapped length (STAR)", 0))),
        ("Percentage of Mapped Reads on Forward Strand", value_or_na(forward_strand, strand_available)),
        ("Percentage of Mapped Reads on Reverse Strand", value_or_na(reverse_strand, strand_available)),
        ("Percentage of Mapped Reads with Failed Strand", value_or_na(failed_strand, strand_available)),
    ])

    length_strand_available = bool(args.length_strandness_file and os.path.exists(args.length_strandness_file))
    length_strand_rows = safe_call(parse_length_strandness, args.length_strandness_file, default=[])
    if length_strand_available and length_strand_rows:
        for row in length_strand_rows:
            label = length_row_label(row)
            ordered_metrics.extend([
                (f"Number of Alignment Records Used for Strand Detection ({label})", _row_int(row, "n_evaluated")),
                (f"Forward Strand Percentage ({label})", _row_pct(row, "pct_forward_all")),
                (f"Reverse Strand Percentage ({label})", _row_pct(row, "pct_reverse_all")),
                (f"Undetermined Percentage ({label})", _row_pct(row, "pct_undetermined_all")),
            ])

    ordered_metrics.extend([
        ("Splices per Kilobase (splice/kb, complementary gDNA-contamination QC proxy)", value_or_na(splice_per_kb, splice_available)),
        ("Number of Uniquely Mapped Reads (STAR)", star_unique_reads if star_unique_reads > 0 else "NA"),
        ("Percentage of Uniquely Mapped Reads (vs STAR Input)", perc_unique_vs_star_input if perc_unique_vs_star_input > 0 else "NA"),
        ("Number of Multi-mapped Reads (STAR)", star_multi_reads if star_multi_reads > 0 else "NA"),
        ("Percentage of Multi-mapped Reads (vs STAR Input)", perc_multi_vs_star_input if perc_multi_vs_star_input > 0 else "NA"),
        ("Number of Splices (from unique reads, STAR)", int(star_metrics.get("Number of splices from uniquely mapped reads (STAR)", 0))),
        ("Number of Unmapped Reads (STAR)", int(star_metrics.get("Unmapped reads number (STAR)", 0))),
        ("Number of Chimeric Reads (STAR)", int(star_metrics.get("Number of chimeric reads (STAR)", 0))),
        ("Mismatch Rate per Base (STAR)", star_metrics.get("Mismatch rate per base (STAR)", "NA"))
    ])
    
    fc_sources = {
        "3'UTR": args.featureCounts_3UTR,
        "5'UTR": args.featureCounts_5UTR,
        "Downstream": args.featureCounts_downstream_2kb,
        "Exonic": args.featureCounts_exon,
        "Intergenic": args.featureCounts_intergenic,
        "Intronic": args.featureCounts_intron,
        "Promoter": args.featureCounts_promoter_1500_500bp,
    }
    fc_counts = {region: safe_call(process_featureCounts_file, path, default=0) for region, path in fc_sources.items()}
    fc_available = {region: bool(path and os.path.exists(path)) for region, path in fc_sources.items()}
    total_fc_meta = sum(fc_counts.values())
    any_fc_available = any(fc_available.values())
    for region, count in fc_counts.items():
        ordered_metrics.append((f"Number of Reads Mapped to {region} Regions", value_or_na(int(count), fc_available[region])))
    for region, count in fc_counts.items():
        pct = round((count / total_fc_meta) * 100, 2) if total_fc_meta > 0 else 0
        ordered_metrics.append((f"Percentage of Reads Mapped to {region} Regions", value_or_na(pct, any_fc_available and fc_available[region])))

    picard_metrics = {}
    picard_metrics.update(safe_call(parse_picard_metrics, args.picard_rnaseq_file, default={}))
    picard_metrics.update(safe_call(parse_picard_metrics, args.picard_insert_file, default={}))
    picard_keys_order = [
        "PF_BASES", "PF_ALIGNED_BASES", "RIBOSOMAL_BASES", "CODING_BASES", "UTR_BASES",
        "INTRONIC_BASES", "INTERGENIC_BASES", "IGNORED_READS", "CORRECT_STRAND_READS",
        "INCORRECT_STRAND_READS", "PCT_RIBOSOMAL_BASES", "PCT_CODING_BASES", "PCT_UTR_BASES",
        "PCT_INTRONIC_BASES", "PCT_INTERGENIC_BASES", "PCT_MRNA_BASES",
        "PCT_USABLE_BASES", "PCT_CORRECT_STRAND_READS", "MEDIAN_CV_COVERAGE", 
        "MEDIAN_5PRIME_BIAS", "MEDIAN_3PRIME_BIAS", "MEDIAN_5PRIME_TO_3PRIME_BIAS", "MEDIAN_INSERT_SIZE"
    ]
    for key in picard_keys_order:
        ordered_metrics.append((f"PICARD_{key}", picard_metrics.get(key, "NA")))

    expr_matrix = safe_call(process_expression_matrix, args.expression_matrix, default={})
    expr_available = bool(args.expression_matrix and os.path.exists(args.expression_matrix) and expr_matrix)
    rna_types_order = [
        "protein_coding", "lncRNAs", "pseudogenes", "miRNAs", "snoRNAs", "snRNAs", "rRNAs", "tRNAs",
        "circRNAs", "ERVs", "LINEs", "SINEs", "IG_genes", "TR_genes", "misc-sncRNAs", "TEC_protein_coding",
        "scaRNAs", "vault_RNAs", "Y_RNAs", "piRNAs", "artifact"
    ]
    read_thresh_order = [1, 5, 10]
    tpm_thresh_order = [0.01, 0.1, 1.0]

    for rna in rna_types_order:
        rna_data = expr_matrix.get(rna, [])
        for thresh in read_thresh_order:
            count = sum(1 for rc, _ in rna_data if rc > thresh)
            ordered_metrics.append((f"Expressed {rna} Genes (Read Counts > {thresh})", value_or_na(count, expr_available)))
        for thresh in tpm_thresh_order:
            count = sum(1 for _, tpm in rna_data if tpm > thresh)
            ordered_metrics.append((f"Expressed {rna} Genes (TPM/CPM > {thresh})", value_or_na(count, expr_available)))

    total_expr_reads = sum(rc for rna_data in expr_matrix.values() for rc, _ in rna_data)
    for rna in rna_types_order:
        total_rc = sum(rc for rc, _ in expr_matrix.get(rna, []))
        ordered_metrics.append((f"Number of Reads Mapped to {rna}", value_or_na(integer_if_integral(total_rc), expr_available)))
    for rna in rna_types_order:
        total_rc = sum(rc for rc, _ in expr_matrix.get(rna, []))
        pct = round((total_rc / total_expr_reads) * 100, 2) if total_expr_reads > 0 else 0
        ordered_metrics.append((f"Percentage of Reads Mapped to {rna}", value_or_na(pct, expr_available and total_expr_reads > 0)))

    with open(args.output, "w") as out:
        out.write("Metric\tValue\n")
        for metric, value in ordered_metrics:
            out.write(f"{metric}\t{value}\n")
    print(f"QC matrix successfully generated: {args.output}")

if __name__ == "__main__":
    main()
