#!/usr/bin/env python3
"""
Compute QC metrics from multiple tools.
Missing or error input files are skipped.
"""

import sys
import argparse
import zipfile
import pysam
import gzip
import pandas as pd

# Safely call a function and return default on exception.
def safe_call(func, *args, default=None):
    try:
        return func(*args)
    except Exception:
        return default

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate QC matrix table from multiple QC files."
    )
    # Check if no arguments are provided; if so, print help and exit.
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit("Error: No input parameters provided. Please provide at least one parameter.")

    # FastQC input files for raw reads
    parser.add_argument("--raw_R1_fastqc_zip", required=False, help="Raw R1 FastQC zip file")
    parser.add_argument("--raw_R2_fastqc_zip", required=False, help="Raw R2 FastQC zip file")
    # Trimmed FASTQ files for trimmed reads metrics
    parser.add_argument("--trimmed_R1_fastq", required=False, help="Trimmed R1 FASTQ file")
    parser.add_argument("--trimmed_R2_fastq", required=False, help="Trimmed R2 FASTQ file")
    # ACC motif fraction file
    parser.add_argument("--ACC_motif_fraction", required=False, help="ACC motif fraction TSV file")
    # Other input files
    parser.add_argument("--kraken_report", required=False, help="Kraken report TSV file")
    parser.add_argument("--ecoli_R1_fastq", required=False, help="E.coli R1 FASTQ (BBSplit)")
    parser.add_argument("--ecoli_R2_fastq", required=False, help="E.coli R2 FASTQ (BBSplit)")
    parser.add_argument("--myco_R1_fastq", required=False, help="Mycoplasma R1 FASTQ (BBSplit)")
    parser.add_argument("--myco_R2_fastq", required=False, help="Mycoplasma R2 FASTQ (BBSplit)")
    parser.add_argument("--bam2strand_file", required=False, help="BAM2strand file (xls)")
    parser.add_argument("--picard_insert_file", required=False, help="Picard Insert Size Metrics TSV file")
    parser.add_argument("--picard_rnaseq_file", required=False, help="Picard RNA Metrics TSV file")
    parser.add_argument("--ribo_R1_fastq", required=False, help="RiboDetector R1 FASTQ (trimmed)")
    parser.add_argument("--ribo_R2_fastq", required=False, help="RiboDetector R2 FASTQ (trimmed)")
    parser.add_argument("--expression_matrix", required=False, help="Merged expression matrix TSV")
    parser.add_argument("--STAR_log", required=False, help="STAR Log file")
    # FeatureCounts input files
    parser.add_argument("--featureCounts_3UTR", required=False, help="FeatureCounts output for 3'UTR")
    parser.add_argument("--featureCounts_5UTR", required=False, help="FeatureCounts output for 5'UTR")
    parser.add_argument("--featureCounts_downstream_2kb", required=False, help="FeatureCounts output for Downstream (2 kb)")
    parser.add_argument("--featureCounts_exon", required=False, help="FeatureCounts output for Exonic Regions")
    parser.add_argument("--featureCounts_ENCODE_blacklist", required=False, help="FeatureCounts output for ENCODE Blacklist Regions")
    parser.add_argument("--featureCounts_intergenic", required=False, help="FeatureCounts output for Intergenic Regions")
    parser.add_argument("--featureCounts_intron", required=False, help="FeatureCounts output for Intronic Regions")
    parser.add_argument("--featureCounts_promoter_1500_500bp", required=False, help="FeatureCounts output for Promoter (1500-500bp)")
    # Downsampled trimmed FASTQ files for RiboDetector input reads
    parser.add_argument("--downsampled_trimmed_R1_fastq", required=False, help="Downsampled Trimmed R1 FASTQ for RiboDetector")
    parser.add_argument("--downsampled_trimmed_R2_fastq", required=False, help="Downsampled Trimmed R2 FASTQ for RiboDetector")
    parser.add_argument("--output", required=False, help="Output QC matrix TSV file")
    return parser.parse_args()

def extract_basic_stats(fastqc_zip):
    # Extract basic stats from FastQC zip file.
    stats = {}
    with zipfile.ZipFile(fastqc_zip, "r") as zf:
        data_file = [f for f in zf.namelist() if f.endswith("fastqc_data.txt")][0]
        with zf.open(data_file) as f:
            lines = f.read().decode("utf-8").splitlines()
    in_basic = False
    for line in lines:
        if line.startswith(">>Basic Statistics"):
            in_basic = True
            continue
        if in_basic and line.startswith(">>END_MODULE"):
            break
        if in_basic and not line.startswith("#"):
            parts = line.split("\t")
            if len(parts) >= 2:
                key = parts[0].strip()
                value = parts[1].strip()
                if key in ["Total Sequences", "%GC", "Sequence length"]:
                    stats[key] = value
    return stats

def parse_sequence_length(val):
    # Return average sequence length if range, else float.
    if "-" in val:
        parts = val.split("-")
        try:
            nums = [float(x) for x in parts]
            return sum(nums) / len(nums)
        except ValueError:
            return float(parts[0])
    else:
        try:
            return float(val)
        except ValueError:
            return 0

def extract_acc_motif_fraction(acc_file):
    # Extract ACC motif fraction from TSV file.
    try:
        df = pd.read_csv(acc_file, sep="\t")
        if "fraction_ACC" in df.columns:
            return df["fraction_ACC"].iloc[0]
        return "NA"
    except Exception:
        return "NA"

def extract_kraken_percentages(kraken_file):
    # Extract percentages for Homo sapiens and Bacteria.
    perc_human = 0.0
    perc_bacteria = 0.0
    with open(kraken_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                continue
            try:
                perc = float(fields[0])
            except ValueError:
                continue
            taxon = fields[-1].strip()
            if taxon == "Homo sapiens":
                perc_human = perc
            elif taxon == "Bacteria":
                perc_bacteria = perc
    return perc_human, perc_bacteria

def count_fastq_reads(fastq_file):
    # Count reads in FASTQ file (supports gzipped).
    count = 0
    opener = gzip.open if fastq_file.endswith(".gz") else open
    with opener(fastq_file, "rt") as f:
        for _ in f:
            count += 1
    return count // 4

def compute_average_read_length(fastq_file):
    # Compute average read length from FASTQ file.
    total_bases = 0
    read_count = 0
    opener = gzip.open if fastq_file.endswith(".gz") else open
    with opener(fastq_file, "rt") as f:
        line_num = 0
        for line in f:
            line_num += 1
            if line_num % 4 == 2:
                seq = line.strip()
                total_bases += len(seq)
                read_count += 1
    return total_bases / read_count if read_count > 0 else 0

def parse_bam2strand(strand_file):
    # Parse BAM2strand file and return strand percentages.
    with open(strand_file) as f:
        header = f.readline()
        line = f.readline().strip()
    fields = line.split("\t")
    try:
        forward = float(fields[2]) * 100
    except (ValueError, IndexError):
        forward = 0.0
    try:
        reverse = float(fields[3]) * 100
    except (ValueError, IndexError):
        reverse = 0.0
    try:
        failed = float(fields[4]) * 100
    except (ValueError, IndexError):
        failed = 0.0
    return round(forward, 2), round(reverse, 2), round(failed, 2)

def parse_picard_insert_size(file_path):
    # Parse Picard Insert Size Metrics and return median insert size.
    header = None
    data = None
    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            if header is None:
                header = line.split("\t")
            else:
                data = line.split("\t")
                break
    metrics = {}
    if header and data:
        try:
            idx = header.index("MEDIAN_INSERT_SIZE")
            metrics["PICARD_MEDIAN_INSERT_SIZE"] = data[idx]
        except ValueError:
            metrics["PICARD_MEDIAN_INSERT_SIZE"] = "NA"
    return metrics

def parse_picard_rnaseq_metrics(file_path):
    # Parse Picard RNA Metrics and return selected metrics.
    desired = ["PF_BASES", "PF_ALIGNED_BASES", "RIBOSOMAL_BASES", "CODING_BASES", "UTR_BASES",
               "INTRONIC_BASES", "INTERGENIC_BASES", "IGNORED_READS", "CORRECT_STRAND_READS",
               "INCORRECT_STRAND_READS", "NUM_R1_TRANSCRIPT_STRAND_READS", "NUM_R2_TRANSCRIPT_STRAND_READS",
               "NUM_UNEXPLAINED_READS", "PCT_R1_TRANSCRIPT_STRAND_READS", "PCT_R2_TRANSCRIPT_STRAND_READS",
               "PCT_RIBOSOMAL_BASES", "PCT_CODING_BASES", "PCT_UTR_BASES", "PCT_INTRONIC_BASES",
               "PCT_INTERGENIC_BASES", "PCT_MRNA_BASES", "PCT_USABLE_BASES", "PCT_CORRECT_STRAND_READS",
               "MEDIAN_CV_COVERAGE", "MEDIAN_5PRIME_BIAS", "MEDIAN_3PRIME_BIAS", "MEDIAN_5PRIME_TO_3PRIME_BIAS"]
    header = None
    data = None
    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if header is None:
                header = line.split("\t")
            else:
                data = line.split("\t")
                break
    metrics = {}
    if header and data:
        for key in desired:
            if key in header:
                idx = header.index(key)
                metrics["PICARD_" + key] = data[idx]
            else:
                metrics["PICARD_" + key] = "NA"
    return metrics

def process_expression_matrix(expr_file):
    """
    Process merged expression matrix (5 columns: GeneID, GeneSymbol, GeneType, ReadCounts, TPM/CPM).
    Returns a dict mapping standardized RNA type to list of (read_count, tpm) tuples.
    """
    expr = {}
    rna_standard = {
        "artifact": "artifact",
        "ervs": "ERVs",
        "ig_genes": "IG_genes",
        "lines": "LINEs",
        "lncrnas": "lncRNAs",
        "mirnas": "miRNAs",
        "misc-sncrnas": "misc-sncRNAs",
        "pirnas": "piRNAs",
        "protein_coding": "protein_coding",
        "pseudogenes": "pseudogenes",
        "rrnas": "rRNAs",
        "scarnas": "scaRNAs",
        "sines": "SINEs",
        "snornas": "snoRNAs",
        "snrnas": "snRNAs",
        "tec_protein_coding": "TEC_protein_coding",
        "tr_genes": "TR_genes",
        "trnas": "tRNAs",
        "vault_rnas": "vault_RNAs",
        "y_rnas": "Y_RNAs",
        "circrnas": "circRNAs"
    }
    allowed = set(rna_standard.values())
    with open(expr_file) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                continue
            raw_type = fields[2].strip()
            std_type = rna_standard.get(raw_type.lower())
            if std_type is None or std_type not in allowed:
                continue
            try:
                read_count = float(fields[3])
            except ValueError:
                read_count = 0.0
            try:
                tpm_value = float(fields[4])
            except ValueError:
                tpm_value = 0.0
            if std_type not in expr:
                expr[std_type] = []
            expr[std_type].append((read_count, tpm_value))
    return expr

def parse_star_log(star_log_file):
    """
    Parse STAR log file and extract key metrics.
    Returns a dict with various STAR metrics.
    """
    metrics = {}
    with open(star_log_file) as f:
        for line in f:
            if "|" in line:
                parts = line.split("|")
                left = parts[0].strip()
                right = parts[1].strip()
                if right.endswith("%"):
                    right = right[:-1].strip()
                try:
                    value = float(right)
                except ValueError:
                    value = right
                metrics[left] = value
    result = {}
    result["Number of input reads (STAR)"] = metrics.get("Number of input reads", 0)
    result["Average input read length (STAR)"] = metrics.get("Average input read length", 0)
    result["Uniquely mapped reads number (STAR)"] = metrics.get("Uniquely mapped reads number", 0)
    multi_mapping = metrics.get("Number of reads mapped to multiple loci", 0) + metrics.get("Number of reads mapped to too many loci", 0)
    result["Multi-mapping reads number (STAR)"] = multi_mapping
    result["Number of splices from uniquely mapped reads (STAR)"] = metrics.get("Number of splices: Total", 0)
    unmapped = (metrics.get("Number of reads unmapped: too many mismatches", 0) +
                metrics.get("Number of reads unmapped: too short", 0) +
                metrics.get("Number of reads unmapped: other", 0))
    result["Unmapped reads number (STAR)"] = unmapped
    result["Number of chimeric reads (STAR)"] = metrics.get("Number of chimeric reads", 0)
    mr = metrics.get("Mismatch rate per base, %", 0)
    result["Mismatch rate per base (STAR)"] = f"{mr}%"
    return result

def process_featureCounts_file(file_path):
    """
    Process a featureCounts file and return total read counts.
    """
    total_counts = 0
    with open(file_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if fields[0] == "Geneid":
                continue
            try:
                count = float(fields[-1])
            except ValueError:
                count = 0
            total_counts += count
    return int(total_counts)

def main():
    args = parse_args()

    # 1. Raw and Trimmed Read Metrics
    raw_stats_r1 = safe_call(extract_basic_stats, args.raw_R1_fastqc_zip, default={})
    raw_stats_r2 = safe_call(extract_basic_stats, args.raw_R2_fastqc_zip, default={})
    total_raw_reads_r1 = safe_call(lambda: float(raw_stats_r1.get("Total Sequences", 0)), default=0)
    total_raw_reads_r2 = safe_call(lambda: float(raw_stats_r2.get("Total Sequences", 0)), default=0)
    total_raw_reads = total_raw_reads_r1 + total_raw_reads_r2
    raw_gc_avg = ((safe_call(lambda: float(raw_stats_r1.get("%GC", 0)), default=0) +
                   safe_call(lambda: float(raw_stats_r2.get("%GC", 0)), default=0)) / 2)
    raw_read_length_r1 = safe_call(lambda: parse_sequence_length(raw_stats_r1.get("Sequence length", "0")), default=0)
    raw_read_length_r2 = safe_call(lambda: parse_sequence_length(raw_stats_r2.get("Sequence length", "0")), default=0)
    avg_raw_read_length = (raw_read_length_r1 + raw_read_length_r2) / 2
    acc_fraction = safe_call(extract_acc_motif_fraction, args.ACC_motif_fraction, default="NA")
    trimmed_reads_r1 = safe_call(count_fastq_reads, args.trimmed_R1_fastq, default=0)
    trimmed_reads_r2 = safe_call(count_fastq_reads, args.trimmed_R2_fastq, default=0)
    total_trimmed_reads = trimmed_reads_r1 + trimmed_reads_r2
    avg_trimmed_read_length = (safe_call(compute_average_read_length, args.trimmed_R1_fastq, default=0) +
                               safe_call(compute_average_read_length, args.trimmed_R2_fastq, default=0)) / 2

    # 2. Contamination Metrics
    perc_human, perc_bacteria = safe_call(extract_kraken_percentages, args.kraken_report, default=(0.0, 0.0))
    ecoli_reads_r1 = safe_call(count_fastq_reads, args.ecoli_R1_fastq, default=0)
    ecoli_reads_r2 = safe_call(count_fastq_reads, args.ecoli_R2_fastq, default=0)
    ecoli_total_reads = ecoli_reads_r1 + ecoli_reads_r2
    pct_ecoli = round((ecoli_total_reads / total_raw_reads) * 100, 2) if total_raw_reads > 0 else 0
    myco_reads_r1 = safe_call(count_fastq_reads, args.myco_R1_fastq, default=0)
    myco_reads_r2 = safe_call(count_fastq_reads, args.myco_R2_fastq, default=0)
    myco_total_reads = myco_reads_r1 + myco_reads_r2
    pct_myco = round((myco_total_reads / total_raw_reads) * 100, 2) if total_raw_reads > 0 else 0
    ribo_reads_r1 = safe_call(count_fastq_reads, args.ribo_R1_fastq, default=0)
    ribo_reads_r2 = safe_call(count_fastq_reads, args.ribo_R2_fastq, default=0)
    total_ribo_reads = ribo_reads_r1 + ribo_reads_r2
    down_ribo_reads_r1 = safe_call(count_fastq_reads, args.downsampled_trimmed_R1_fastq, default=0)
    down_ribo_reads_r2 = safe_call(count_fastq_reads, args.downsampled_trimmed_R2_fastq, default=0)
    total_down_ribo_reads = down_ribo_reads_r1 + down_ribo_reads_r2
    perc_rRNA_ribodetector = round((total_ribo_reads / total_down_ribo_reads) * 100, 2) if total_down_ribo_reads > 0 else 0

    # 3. STAR Metrics and UMI-deduped Strand Metrics
    star_metrics = safe_call(parse_star_log, args.STAR_log, default={})
    umi_input_reads_star = star_metrics.get("Number of input reads (STAR)", 0)
    perc_umi_reads_star = round((umi_input_reads_star / total_raw_reads) * 100, 2) if total_raw_reads > 0 else 0
    avg_umi_fragment_length_star = star_metrics.get("Average input read length (STAR)", 0)
    (umi_forward_strand,
     umi_reverse_strand,
     umi_failed_strand) = safe_call(parse_bam2strand, args.bam2strand_file, default=(0, 0, 0))
    umi_unique_reads_star = star_metrics.get("Uniquely mapped reads number (STAR)", 0)
    perc_umi_unique_reads_star = round((umi_unique_reads_star / umi_input_reads_star) * 100, 2) if umi_input_reads_star > 0 else 0
    umi_multi_reads_star = star_metrics.get("Multi-mapping reads number (STAR)", 0)
    perc_umi_multi_reads_star = round((umi_multi_reads_star / umi_input_reads_star) * 100, 2) if umi_input_reads_star > 0 else 0
    umi_splices_star = star_metrics.get("Number of splices from uniquely mapped reads (STAR)", 0)
    umi_unmapped_reads_star = star_metrics.get("Unmapped reads number (STAR)", 0)
    umi_chimeric_reads_star = star_metrics.get("Number of chimeric reads (STAR)", 0)
    umi_mismatch_rate_star = star_metrics.get("Mismatch rate per base (STAR)", "0%")

    # 3.5. PICARD Metrics
    picard_insert_metrics = safe_call(parse_picard_insert_size, args.picard_insert_file, default={})
    picard_rnaseq_metrics = safe_call(parse_picard_rnaseq_metrics, args.picard_rnaseq_file, default={})
    picard_metrics = {}
    if picard_rnaseq_metrics is not None:
        picard_metrics.update(picard_rnaseq_metrics)
    if picard_insert_metrics is not None:
        picard_metrics.update(picard_insert_metrics)

    # 4. FeatureCounts Meta Gene Metrics
    fc_3utr = safe_call(process_featureCounts_file, args.featureCounts_3UTR, default=0)
    fc_5utr = safe_call(process_featureCounts_file, args.featureCounts_5UTR, default=0)
    fc_downstream = safe_call(process_featureCounts_file, args.featureCounts_downstream_2kb, default=0)
    fc_exon = safe_call(process_featureCounts_file, args.featureCounts_exon, default=0)
    fc_intergenic = safe_call(process_featureCounts_file, args.featureCounts_intergenic, default=0)
    fc_intron = safe_call(process_featureCounts_file, args.featureCounts_intron, default=0)
    fc_promoter = safe_call(process_featureCounts_file, args.featureCounts_promoter_1500_500bp, default=0)
    fc_blacklist = safe_call(process_featureCounts_file, args.featureCounts_ENCODE_blacklist, default=0)
    total_fc_meta = fc_3utr + fc_5utr + fc_downstream + fc_exon + fc_intergenic + fc_intron + fc_promoter + fc_blacklist
    pct_fc_3utr = round((fc_3utr / total_fc_meta) * 100, 2) if total_fc_meta > 0 else 0
    pct_fc_5utr = round((fc_5utr / total_fc_meta) * 100, 2) if total_fc_meta > 0 else 0
    pct_fc_downstream = round((fc_downstream / total_fc_meta) * 100, 2) if total_fc_meta > 0 else 0
    pct_fc_exon = round((fc_exon / total_fc_meta) * 100, 2) if total_fc_meta > 0 else 0
    pct_fc_intergenic = round((fc_intergenic / total_fc_meta) * 100, 2) if total_fc_meta > 0 else 0
    pct_fc_intron = round((fc_intron / total_fc_meta) * 100, 2) if total_fc_meta > 0 else 0
    pct_fc_promoter = round((fc_promoter / total_fc_meta) * 100, 2) if total_fc_meta > 0 else 0
    pct_fc_blacklist = round((fc_blacklist / total_fc_meta) * 100, 2) if total_fc_meta > 0 else 0

    # 7. Expression Matrix Metrics - Threshold-based Counts
    rna_types_order = ["artifact", "ERVs", "IG_genes", "LINEs", "lncRNAs", "miRNAs", "misc-sncRNAs",
                       "piRNAs", "protein_coding", "pseudogenes", "rRNAs", "scaRNAs", "SINEs",
                       "snoRNAs", "snRNAs", "TEC_protein_coding", "TR_genes", "tRNAs", "vault_RNAs",
                       "Y_RNAs", "circRNAs"]
    read_thresh_order = [1, 10, 5]
    tpm_thresh_order = [0.001, 0.01, 0.1, 0.2, 0.5, 1]
    expr_matrix = safe_call(process_expression_matrix, args.expression_matrix, default={})

    # 8. Expression Matrix Metrics - FeatureCounts Totals and Percentages
    # (Counts and percentages will be calculated below)

    # Assemble final ordered metrics list
    ordered_metrics = []

    # 1. Raw and Trimmed Reads
    ordered_metrics.append(("Raw Read Counts Number", int(round(total_raw_reads)) if total_raw_reads > 0 else "NA"))
    ordered_metrics.append(("%GC of Raw Reads", round(raw_gc_avg, 2) if raw_gc_avg > 0 else "NA"))
    ordered_metrics.append(("Raw Read Length", int(round(avg_raw_read_length)) if avg_raw_read_length > 0 else "NA"))
    ordered_metrics.append(("ACC motif fraction from UMI region", acc_fraction))
    ordered_metrics.append(("Average Trimmed Read Length", int(round(avg_trimmed_read_length)) if avg_trimmed_read_length > 0 else "NA"))
    ordered_metrics.append(("Number of Trimmed Reads", int(round(total_trimmed_reads)) if total_trimmed_reads > 0 else "NA"))
    ordered_metrics.append(("Percentage of Reads after Trimming", round((total_trimmed_reads / total_raw_reads) * 100, 2) if total_raw_reads > 0 else "NA"))

    # 2. Contamination Metrics
    ordered_metrics.append(("Percentage of Clean Reads Mapped to Human (Kraken)", perc_human if perc_human > 0 else "NA"))
    ordered_metrics.append(("Percentage of Clean Reads Mapped to Non-human (Kraken)", round(100 - perc_human, 2) if perc_human > 0 else "NA"))
    ordered_metrics.append(("Percentage of Clean Reads Mapped to Bacteria (Kraken)", perc_bacteria if perc_bacteria > 0 else "NA"))
    ordered_metrics.append(("Number of Clean Reads Mapped to Escherichia coli", int(round(ecoli_total_reads)) if ecoli_total_reads > 0 else "NA"))
    ordered_metrics.append(("Percentage of Clean Reads Mapped to Escherichia coli", pct_ecoli if pct_ecoli > 0 else "NA"))
    ordered_metrics.append(("Number of Clean Reads Mapped to Mycoplasma", int(round(myco_total_reads)) if myco_total_reads > 0 else "NA"))
    ordered_metrics.append(("Percentage of Clean Reads Mapped to Mycoplasma", pct_myco if pct_myco > 0 else "NA"))
    ordered_metrics.append(("Percentage of Trimmed Reads Mapped to rRNAs (RiboDetector)", perc_rRNA_ribodetector if perc_rRNA_ribodetector > 0 else "NA"))

    # 3. STAR Metrics and UMI-deduped Strand Metrics
    ordered_metrics.append(("Number of UMI-deduped Input Reads (STAR)", umi_input_reads_star if umi_input_reads_star > 0 else "NA"))
    ordered_metrics.append(("Percentage of UMI-deduped Reads (STAR)", perc_umi_reads_star if perc_umi_reads_star > 0 else "NA"))
    ordered_metrics.append(("Average Mapped UMI-deduped Fragment Length (STAR)", avg_umi_fragment_length_star if avg_umi_fragment_length_star > 0 else "NA"))
    ordered_metrics.append(("Percentage of UMI-deduped Reads on the Forward Strand", umi_forward_strand if umi_forward_strand > 0 else "NA"))
    ordered_metrics.append(("Percentage of UMI-deduped Reads on the Reverse Strand", umi_reverse_strand if umi_reverse_strand > 0 else "NA"))
    ordered_metrics.append(("Percentage of UMI-deduped Reads with Failed Strand", umi_failed_strand if umi_failed_strand > 0 else "NA"))
    ordered_metrics.append(("Number of UMI-deduped Reads Uniquely Mapped to Human (STAR)", umi_unique_reads_star if umi_unique_reads_star > 0 else "NA"))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Uniquely Mapped to Human (STAR)", perc_umi_unique_reads_star if perc_umi_unique_reads_star > 0 else "NA"))
    ordered_metrics.append(("Number of UMI-deduped Reads Multi-mapped to Human (STAR)", umi_multi_reads_star if umi_multi_reads_star > 0 else "NA"))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Multi-mapped to Human (STAR)", perc_umi_multi_reads_star if perc_umi_multi_reads_star > 0 else "NA"))
    ordered_metrics.append(("Number of Splices from UMI-deduped Reads Uniquely Mapped to Human (STAR)", umi_splices_star if umi_splices_star > 0 else "NA"))
    ordered_metrics.append(("Number of Human Unmapped UMI-deduped Reads (STAR)", umi_unmapped_reads_star if umi_unmapped_reads_star > 0 else "NA"))
    ordered_metrics.append(("Number of Chimeric UMI-deduped Reads Mapped to Human (STAR)", umi_chimeric_reads_star if umi_chimeric_reads_star > 0 else "NA"))
    ordered_metrics.append(("Mismatch Rate per Base for UMI-deduped Reads (STAR)", umi_mismatch_rate_star if umi_mismatch_rate_star != "0%" else "NA"))

    # 4. FeatureCounts Meta Gene Metrics - Counts
    ordered_metrics.append(("Number of UMI-deduped Reads Mapped to Human 3'UTR (featureCounts)", fc_3utr if fc_3utr > 0 else "NA"))
    ordered_metrics.append(("Number of UMI-deduped Reads Mapped to Human 5'UTR (featureCounts)", fc_5utr if fc_5utr > 0 else "NA"))
    ordered_metrics.append(("Number of UMI-deduped Reads Mapped to Human Downstream (featureCounts)", fc_downstream if fc_downstream > 0 else "NA"))
    ordered_metrics.append(("Number of UMI-deduped Reads Mapped to Human Exonic Regions (featureCounts)", fc_exon if fc_exon > 0 else "NA"))
    ordered_metrics.append(("Number of UMI-deduped Reads Mapped to Human Intergenic Regions (featureCounts)", fc_intergenic if fc_intergenic > 0 else "NA"))
    ordered_metrics.append(("Number of UMI-deduped Reads Mapped to Human Intronic Regions (featureCounts)", fc_intron if fc_intron > 0 else "NA"))
    ordered_metrics.append(("Number of UMI-deduped Reads Mapped to Human Promoter Regions (featureCounts)", fc_promoter if fc_promoter > 0 else "NA"))
    ordered_metrics.append(("Number of UMI-deduped Reads Mapped to Human Blacklist Regions (featureCounts)", fc_blacklist if fc_blacklist > 0 else "NA"))

    # 5. FeatureCounts Meta Gene Metrics - Percentages
    ordered_metrics.append(("Percentage of UMI-deduped Reads Mapped to Human 3'UTR (featureCounts)", pct_fc_3utr if pct_fc_3utr > 0 else "NA"))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Mapped to Human 5'UTR (featureCounts)", pct_fc_5utr if pct_fc_5utr > 0 else "NA"))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Mapped to Human Downstream (featureCounts)", pct_fc_downstream if pct_fc_downstream > 0 else "NA"))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Mapped to Human Exonic Regions (featureCounts)", pct_fc_exon if pct_fc_exon > 0 else "NA"))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Mapped to Human Intergenic Regions (featureCounts)", pct_fc_intergenic if pct_fc_intergenic > 0 else "NA"))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Mapped to Human Intronic Regions (featureCounts)", pct_fc_intron if pct_fc_intron > 0 else "NA"))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Mapped to Human Promoter Regions (featureCounts)", pct_fc_promoter if pct_fc_promoter > 0 else "NA"))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Mapped to Human Blacklist Regions (featureCounts)", pct_fc_blacklist if pct_fc_blacklist > 0 else "NA"))

    # 6. PICARD Metrics
    picard_keys_order = [
        "PICARD_PF_BASES", "PICARD_PF_ALIGNED_BASES", "PICARD_RIBOSOMAL_BASES",
        "PICARD_CODING_BASES", "PICARD_UTR_BASES", "PICARD_INTRONIC_BASES",
        "PICARD_INTERGENIC_BASES", "PICARD_IGNORED_READS", "PICARD_CORRECT_STRAND_READS",
        "PICARD_INCORRECT_STRAND_READS", "PICARD_NUM_R1_TRANSCRIPT_STRAND_READS",
        "PICARD_NUM_R2_TRANSCRIPT_STRAND_READS", "PICARD_NUM_UNEXPLAINED_READS",
        "PICARD_PCT_R1_TRANSCRIPT_STRAND_READS", "PICARD_PCT_R2_TRANSCRIPT_STRAND_READS",
        "PICARD_PCT_RIBOSOMAL_BASES", "PICARD_PCT_CODING_BASES", "PICARD_PCT_UTR_BASES",
        "PICARD_PCT_INTRONIC_BASES", "PICARD_PCT_INTERGENIC_BASES", "PICARD_PCT_MRNA_BASES",
        "PICARD_PCT_USABLE_BASES", "PICARD_PCT_CORRECT_STRAND_READS",
        "PICARD_MEDIAN_CV_COVERAGE", "PICARD_MEDIAN_5PRIME_BIAS",
        "PICARD_MEDIAN_3PRIME_BIAS", "PICARD_MEDIAN_5PRIME_TO_3PRIME_BIAS",
        "PICARD_MEDIAN_INSERT_SIZE"
    ]
    for key in picard_keys_order:
        ordered_metrics.append((key, picard_metrics.get(key, "NA")))

    # 7. Expression Matrix Metrics - Threshold-based Counts
    for rna in rna_types_order:
        for thresh in read_thresh_order:
            value = sum(1 for (rc, tpm) in expr_matrix.get(rna, []) if rc > thresh)
            ordered_metrics.append((f"Expressed Human {rna} (UMI-deduped Read Counts > {thresh})", value if value > 0 else "NA"))
        for thresh in tpm_thresh_order:
            value = sum(1 for (rc, tpm) in expr_matrix.get(rna, []) if tpm > thresh)
            ordered_metrics.append((f"Expressed Human {rna} (UMI-deduped TPM/CPM > {thresh})", value if value > 0 else "NA"))

    # 8. Expression Matrix Metrics - FeatureCounts Totals and Percentages
    for rna in rna_types_order:
        total_fc = sum(rc for (rc, tpm) in expr_matrix.get(rna, []))
        ordered_metrics.append((f"Number of UMI-deduped Reads Mapped to {rna} (featureCounts)", total_fc if total_fc > 0 else "NA"))
    total_fc_all = sum(sum(rc for (rc, tpm) in expr_matrix.get(rna, [])) for rna in rna_types_order)
    for rna in rna_types_order:
        total_fc = sum(rc for (rc, tpm) in expr_matrix.get(rna, []))
        pct_fc = round((total_fc / total_fc_all) * 100, 2) if total_fc_all > 0 else 0
        ordered_metrics.append((f"Percentage of UMI-deduped Reads Mapped to {rna} (featureCounts)", pct_fc if pct_fc > 0 else "NA"))

    # Write ordered metrics to output TSV
    with open(args.output, "w") as out:
        out.write("Metric\tValue\n")
        for metric, value in ordered_metrics:
            if value is not None:
                out.write(f"{metric}\t{value}\n")
    print(f"{args.output} generated.")

if __name__ == "__main__":
    main()