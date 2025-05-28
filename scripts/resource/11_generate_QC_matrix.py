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

# Safely call a function and return default on exception.
def safe_call(func, *args, default=None):
    try:
        return func(*args)
    except Exception:
        return default

def parse_args():
    parser = argparse.ArgumentParser(description="Generate QC matrix table from multiple QC files.")
    # FastQC input files for raw reads
    parser.add_argument("--raw_R1_fastqc_zip", required=True, help="Raw R1 FastQC zip file")
    parser.add_argument("--raw_R2_fastqc_zip", required=True, help="Raw R2 FastQC zip file")
    # New: trimmed FASTQ files for trimmed reads metrics
    parser.add_argument("--trimmed_R1_fastq", required=True, help="Trimmed R1 FASTQ file")
    parser.add_argument("--trimmed_R2_fastq", required=True, help="Trimmed R2 FASTQ file")
    # Other input files
    parser.add_argument("--kraken_report", required=True, help="Kraken report TSV file")
    parser.add_argument("--ecoli_R1_fastq", required=True, help="E.coli R1 FASTQ (BBSplit)")
    parser.add_argument("--ecoli_R2_fastq", required=True, help="E.coli R2 FASTQ (BBSplit)")
    parser.add_argument("--myco_R1_fastq", required=True, help="Mycoplasma R1 FASTQ (BBSplit)")
    parser.add_argument("--myco_R2_fastq", required=True, help="Mycoplasma R2 FASTQ (BBSplit)")
    parser.add_argument("--umi_dedup_bam", required=True, help="UMI-deduplicated BAM file")
    parser.add_argument("--bam2strand_file", required=True, help="BAM2strand file (xls)")
    parser.add_argument("--picard_insert_file", required=True, help="Picard Insert Size Metrics TSV file")
    parser.add_argument("--picard_rnaseq_file", required=True, help="Picard RNA Metrics TSV file")
    parser.add_argument("--ribo_R1_fastq", required=True, help="RiboDetector R1 FASTQ (trimmed)")
    parser.add_argument("--ribo_R2_fastq", required=True, help="RiboDetector R2 FASTQ (trimmed)")
    parser.add_argument("--bedtools_cov_3UTR", required=True, help="bedtools 3'UTR coverage TSV")
    parser.add_argument("--bedtools_cov_5UTR", required=True, help="bedtools 5'UTR coverage TSV")
    parser.add_argument("--bedtools_cov_downstream_2kb", required=True, help="bedtools Downstream (2 kb) coverage TSV")
    parser.add_argument("--bedtools_cov_exon", required=True, help="bedtools Exon coverage TSV")
    parser.add_argument("--bedtools_cov_ENCODE_blacklist", required=True, help="bedtools ENCODE Blacklist coverage TSV")
    parser.add_argument("--bedtools_cov_intergenic", required=True, help="bedtools Intergenic coverage TSV")
    parser.add_argument("--bedtools_cov_intron", required=True, help="bedtools Intron coverage TSV")
    parser.add_argument("--bedtools_cov_promoter_1500_500bp", required=True, help="bedtools Promoter (1500-500bp) coverage TSV")
    parser.add_argument("--expression_matrix", required=True, help="Merged expression matrix TSV")
    parser.add_argument("--STAR_log", required=True, help="STAR Log file")
    # FeatureCounts input files
    parser.add_argument("--featureCounts_3UTR", required=True, help="FeatureCounts output for 3'UTR")
    parser.add_argument("--featureCounts_5UTR", required=True, help="FeatureCounts output for 5'UTR")
    parser.add_argument("--featureCounts_downstream_2kb", required=True, help="FeatureCounts output for Downstream (2 kb)")
    parser.add_argument("--featureCounts_exon", required=True, help="FeatureCounts output for Exonic Regions")
    parser.add_argument("--featureCounts_ENCODE_blacklist", required=True, help="FeatureCounts output for ENCODE Blacklist Regions")
    parser.add_argument("--featureCounts_intergenic", required=True, help="FeatureCounts output for Intergenic Regions")
    parser.add_argument("--featureCounts_intron", required=True, help="FeatureCounts output for Intronic Regions")
    parser.add_argument("--featureCounts_promoter_1500_500bp", required=True, help="FeatureCounts output for Promoter (1500-500bp)")
    # New: downsampled trimmed FASTQ files for RiboDetector input reads
    parser.add_argument("--downsampled_trimmed_R1_fastq", required=True, help="Downsampled Trimmed R1 FASTQ for RiboDetector")
    parser.add_argument("--downsampled_trimmed_R2_fastq", required=True, help="Downsampled Trimmed R2 FASTQ for RiboDetector")
    parser.add_argument("--output", required=True, help="Output QC matrix TSV file")
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

def process_bedtools_file(file_path):
    """
    Process a bedtools coverage file.
    Returns total read counts, total mapped bases, and total length.
    """
    total_counts = 0
    total_mapped = 0
    total_len = 0
    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            try:
                count = float(fields[6])
                mapped = float(fields[7])
                length = float(fields[8])
            except ValueError:
                continue
            total_counts += count
            total_mapped += mapped
            total_len += length
    return total_counts, total_mapped, total_len

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

    # ---------------------------
    # 1. Raw and Trimmed Read Metrics
    stats_R1 = safe_call(extract_basic_stats, args.raw_R1_fastqc_zip, default={})
    stats_R2 = safe_call(extract_basic_stats, args.raw_R2_fastqc_zip, default={})
    Total_Raw_Reads_R1 = safe_call(lambda: float(stats_R1.get("Total Sequences", 0)), default=0)
    Total_Raw_Reads_R2 = safe_call(lambda: float(stats_R2.get("Total Sequences", 0)), default=0)
    Total_Raw_Reads = Total_Raw_Reads_R1 + Total_Raw_Reads_R2
    GC_Total_Raw_Reads = ((safe_call(lambda: float(stats_R1.get("%GC", 0)), default=0) +
                           safe_call(lambda: float(stats_R2.get("%GC", 0)), default=0)) / 2)
    Read_Length_R1 = safe_call(lambda: parse_sequence_length(stats_R1.get("Sequence length", "0")), default=0)
    Read_Length_R2 = safe_call(lambda: parse_sequence_length(stats_R2.get("Sequence length", "0")), default=0)
    Read_Length_for_raw_reads = (Read_Length_R1 + Read_Length_R2) / 2
    Trimmed_Reads_R1 = safe_call(count_fastq_reads, args.trimmed_R1_fastq, default=0)
    Trimmed_Reads_R2 = safe_call(count_fastq_reads, args.trimmed_R2_fastq, default=0)
    Number_Of_Reads_After_TrimGalore_From_Raw_Reads = Trimmed_Reads_R1 + Trimmed_Reads_R2
    Avg_Read_Length_of_Trimmed_Reads = (safe_call(compute_average_read_length, args.trimmed_R1_fastq, default=0) +
                                        safe_call(compute_average_read_length, args.trimmed_R2_fastq, default=0)) / 2

    # ---------------------------
    # 2. Contamination Metrics
    perc_Human, perc_Bacteria = safe_call(extract_kraken_percentages, args.kraken_report, default=(0.0, 0.0))
    Ecoli_Reads_R1 = safe_call(count_fastq_reads, args.ecoli_R1_fastq, default=0)
    Ecoli_Reads_R2 = safe_call(count_fastq_reads, args.ecoli_R2_fastq, default=0)
    Number_of_Clean_Reads_Mapped_to_Escherichia_coli = Ecoli_Reads_R1 + Ecoli_Reads_R2
    pct_Ecoli = round((Number_of_Clean_Reads_Mapped_to_Escherichia_coli / Total_Raw_Reads) * 100, 2) if Total_Raw_Reads > 0 else 0
    Myco_Reads_R1 = safe_call(count_fastq_reads, args.myco_R1_fastq, default=0)
    Myco_Reads_R2 = safe_call(count_fastq_reads, args.myco_R2_fastq, default=0)
    Number_of_Clean_Reads_Mapped_to_Mycoplasma = Myco_Reads_R1 + Myco_Reads_R2
    pct_Mycoplasma = round((Number_of_Clean_Reads_Mapped_to_Mycoplasma / Total_Raw_Reads) * 100, 2) if Total_Raw_Reads > 0 else 0
    Ribo_Reads_R1 = safe_call(count_fastq_reads, args.ribo_R1_fastq, default=0)
    Ribo_Reads_R2 = safe_call(count_fastq_reads, args.ribo_R2_fastq, default=0)
    Total_Ribo_Reads = Ribo_Reads_R1 + Ribo_Reads_R2
    Down_Ribo_Reads_R1 = safe_call(count_fastq_reads, args.downsampled_trimmed_R1_fastq, default=0)
    Down_Ribo_Reads_R2 = safe_call(count_fastq_reads, args.downsampled_trimmed_R2_fastq, default=0)
    Total_Down_Ribo_Reads = Down_Ribo_Reads_R1 + Down_Ribo_Reads_R2
    Percentage_of_trimmed_Reads_Mapped_to_rRNA_RiboDetector = round((Total_Ribo_Reads / Total_Down_Ribo_Reads) * 100, 2) if Total_Down_Ribo_Reads > 0 else 0

    # ---------------------------
    # 3. STAR Metrics and UMI-deduped Strand Metrics
    star_metrics = safe_call(parse_star_log, args.STAR_log, default={})
    Number_of_UMI_deduped_Input_Reads_STAR = star_metrics.get("Number of input reads (STAR)", 0)
    Percentage_of_UMI_deduped_Reads_STAR = round((Number_of_UMI_deduped_Input_Reads_STAR / Total_Raw_Reads) * 100, 2) if Total_Raw_Reads > 0 else 0
    Average_Mapped_UMI_deduped_Fragment_Length_STAR = star_metrics.get("Average input read length (STAR)", 0)
    (UMI_deduped_Forward_Strand,
     UMI_deduped_Reverse_Strand,
     UMI_deduped_Failed_Strand) = safe_call(parse_bam2strand, args.bam2strand_file, default=(0, 0, 0))
    Number_of_UMI_deduped_Reads_Uniquely_Mapped_to_Human_STAR = star_metrics.get("Uniquely mapped reads number (STAR)", 0)
    Percentage_of_UMI_deduped_Reads_Uniquely_Mapped_to_Human_STAR = round((Number_of_UMI_deduped_Reads_Uniquely_Mapped_to_Human_STAR / Number_of_UMI_deduped_Input_Reads_STAR) * 100, 2) if Number_of_UMI_deduped_Input_Reads_STAR > 0 else 0
    Number_of_UMI_deduped_Reads_Multi_mapped_to_Human_STAR = star_metrics.get("Multi-mapping reads number (STAR)", 0)
    Percentage_of_UMI_deduped_Reads_Multi_mapped_to_Human_STAR = round((Number_of_UMI_deduped_Reads_Multi_mapped_to_Human_STAR / Number_of_UMI_deduped_Input_Reads_STAR) * 100, 2) if Number_of_UMI_deduped_Input_Reads_STAR > 0 else 0
    Number_of_Splices_from_UMI_deduped_Reads_Uniquely_Mapped_to_Human_STAR = star_metrics.get("Number of splices from uniquely mapped reads (STAR)", 0)
    Number_of_Human_Unmapped_UMI_deduped_Reads_STAR = star_metrics.get("Unmapped reads number (STAR)", 0)
    Number_of_Chimeric_UMI_deduped_Reads_Mapped_to_Human_STAR = star_metrics.get("Number of chimeric reads (STAR)", 0)
    Mismatch_Rate_per_Base_for_UMI_deduped_Reads_STAR = star_metrics.get("Mismatch rate per base (STAR)", "0%")
    
    # ---------------------------
    # 3.5. PICARD Metrics
    picard_insert = safe_call(parse_picard_insert_size, args.picard_insert_file, default={})
    picard_rnaseq = safe_call(parse_picard_rnaseq_metrics, args.picard_rnaseq_file, default={})
    picard_metrics = {}
    if picard_rnaseq is not None:
        picard_metrics.update(picard_rnaseq)
    if picard_insert is not None:
        picard_metrics.update(picard_insert)
    
    # ---------------------------
    # 4. FeatureCounts Meta Gene Metrics (Counts and Percentages)
    fc_3UTR = safe_call(process_featureCounts_file, args.featureCounts_3UTR, default=0)
    fc_5UTR = safe_call(process_featureCounts_file, args.featureCounts_5UTR, default=0)
    fc_downstream_2kb = safe_call(process_featureCounts_file, args.featureCounts_downstream_2kb, default=0)
    fc_exon = safe_call(process_featureCounts_file, args.featureCounts_exon, default=0)
    fc_intergenic = safe_call(process_featureCounts_file, args.featureCounts_intergenic, default=0)
    fc_intron = safe_call(process_featureCounts_file, args.featureCounts_intron, default=0)
    fc_promoter = safe_call(process_featureCounts_file, args.featureCounts_promoter_1500_500bp, default=0)
    fc_encode_blacklist = safe_call(process_featureCounts_file, args.featureCounts_ENCODE_blacklist, default=0)
    Total_meta_FeatureCounts = fc_3UTR + fc_5UTR + fc_downstream_2kb + fc_exon + fc_intergenic + fc_intron + fc_promoter + fc_encode_blacklist
    pct_fc_3UTR = round((fc_3UTR / Total_meta_FeatureCounts) * 100, 2) if Total_meta_FeatureCounts > 0 else 0
    pct_fc_5UTR = round((fc_5UTR / Total_meta_FeatureCounts) * 100, 2) if Total_meta_FeatureCounts > 0 else 0
    pct_fc_downstream = round((fc_downstream_2kb / Total_meta_FeatureCounts) * 100, 2) if Total_meta_FeatureCounts > 0 else 0
    pct_fc_exon = round((fc_exon / Total_meta_FeatureCounts) * 100, 2) if Total_meta_FeatureCounts > 0 else 0
    pct_fc_intergenic = round((fc_intergenic / Total_meta_FeatureCounts) * 100, 2) if Total_meta_FeatureCounts > 0 else 0
    pct_fc_intron = round((fc_intron / Total_meta_FeatureCounts) * 100, 2) if Total_meta_FeatureCounts > 0 else 0
    pct_fc_promoter = round((fc_promoter / Total_meta_FeatureCounts) * 100, 2) if Total_meta_FeatureCounts > 0 else 0
    pct_fc_blacklist = round((fc_encode_blacklist / Total_meta_FeatureCounts) * 100, 2) if Total_meta_FeatureCounts > 0 else 0

    # ---------------------------
    # 5. Bedtools Mapping Metrics
    bedtools_regions = ["3'UTR", "5'UTR", "Downstream (2 kb)", "Exonic", "Intergenic", "Intronic", "Promoter", "ENCODE Blacklist"]
    bedtools_files = {
        "3'UTR": args.bedtools_cov_3UTR,
        "5'UTR": args.bedtools_cov_5UTR,
        "Downstream (2 kb)": args.bedtools_cov_downstream_2kb,
        "Exonic": args.bedtools_cov_exon,
        "Intergenic": args.bedtools_cov_intergenic,
        "Intronic": args.bedtools_cov_intron,
        "Promoter": args.bedtools_cov_promoter_1500_500bp,
        "ENCODE Blacklist": args.bedtools_cov_ENCODE_blacklist
    }
    bedtools_counts = {}
    for region in bedtools_regions:
        bedtools_counts[region] = safe_call(process_bedtools_file, bedtools_files[region], default=(0,0,0))[0]
    Total_bedtools = sum(int(round(bedtools_counts[r])) for r in bedtools_regions)
    pct_bt_3UTR = round((int(round(bedtools_counts["3'UTR"])) / Total_bedtools) * 100, 2) if Total_bedtools > 0 else 0
    pct_bt_5UTR = round((int(round(bedtools_counts["5'UTR"])) / Total_bedtools) * 100, 2) if Total_bedtools > 0 else 0
    pct_bt_downstream = round((int(round(bedtools_counts["Downstream (2 kb)"])) / Total_bedtools) * 100, 2) if Total_bedtools > 0 else 0
    pct_bt_exonic = round((int(round(bedtools_counts["Exonic"])) / Total_bedtools) * 100, 2) if Total_bedtools > 0 else 0
    pct_bt_intergenic = round((int(round(bedtools_counts["Intergenic"])) / Total_bedtools) * 100, 2) if Total_bedtools > 0 else 0
    pct_bt_intronic = round((int(round(bedtools_counts["Intronic"])) / Total_bedtools) * 100, 2) if Total_bedtools > 0 else 0
    pct_bt_promoter = round((int(round(bedtools_counts["Promoter"])) / Total_bedtools) * 100, 2) if Total_bedtools > 0 else 0
    pct_bt_blacklist = round((int(round(bedtools_counts["ENCODE Blacklist"])) / Total_bedtools) * 100, 2) if Total_bedtools > 0 else 0

    # ---------------------------
    # 6. PICARD Metrics (固定顺序)
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

    # ---------------------------
    # 7. Expression Matrix Metrics - Threshold-based Counts
    rna_types_order = ["artifact", "ERVs", "IG_genes", "LINEs", "lncRNAs", "miRNAs", "misc-sncRNAs",
                       "piRNAs", "protein_coding", "pseudogenes", "rRNAs", "scaRNAs", "SINEs",
                       "snoRNAs", "snRNAs", "TEC_protein_coding", "TR_genes", "tRNAs", "vault_RNAs",
                       "Y_RNAs", "circRNAs"]
    read_thresh_order = [1, 10, 5]
    tpm_thresh_order = [0.001, 0.01, 0.1, 0.2, 0.5, 1]
    expr_matrix = safe_call(process_expression_matrix, args.expression_matrix, default={})
    expr_metrics_list = []
    for rna in rna_types_order:
        for thresh in read_thresh_order:
            value = sum(1 for (rc, tpm) in expr_matrix.get(rna, []) if rc > thresh)
            expr_metrics_list.append((f"Expressed Human {rna} (UMI-deduped Read Counts > {thresh})", value))
        for thresh in tpm_thresh_order:
            value = sum(1 for (rc, tpm) in expr_matrix.get(rna, []) if tpm > thresh)
            expr_metrics_list.append((f"Expressed Human {rna} (UMI-deduped TPM/CPM > {thresh})", value))

    # ---------------------------
    # 8. Expression Matrix Metrics - FeatureCounts Totals and Percentages
    expr_fc_metrics = []
    for rna in rna_types_order:
        total_fc = sum(rc for (rc, tpm) in expr_matrix.get(rna, []))
        expr_fc_metrics.append((f"Number of UMI-deduped Reads Mapped to {rna} (featureCounts)", total_fc))
    total_fc_all = sum(sum(rc for (rc, tpm) in expr_matrix.get(rna, [])) for rna in rna_types_order)
    expr_fc_pct_metrics = []
    for rna in rna_types_order:
        total_fc = sum(rc for (rc, tpm) in expr_matrix.get(rna, []))
        pct_fc = round((total_fc / total_fc_all) * 100, 2) if total_fc_all > 0 else 0
        expr_fc_pct_metrics.append((f"Percentage of UMI-deduped Reads Mapped to {rna} (featureCounts)", pct_fc))

    # ---------------------------
    # Assemble final ordered metrics list according to specified order
    ordered_metrics = []

    # 1. Raw and Trimmed Reads
    ordered_metrics.append(("Number of Total Raw Reads", int(round(Total_Raw_Reads))))
    ordered_metrics.append(("%GC Total Raw Reads", GC_Total_Raw_Reads))
    ordered_metrics.append(("Read Length for Raw Reads", int(round(Read_Length_for_raw_reads))))
    ordered_metrics.append(("Average Trimmed Read Length", int(round(Avg_Read_Length_of_Trimmed_Reads))))
    ordered_metrics.append(("Number of Trimmed Reads (TrimGalore)", int(round(Number_Of_Reads_After_TrimGalore_From_Raw_Reads))))
    ordered_metrics.append(("Percentage of Trimmed Reads (TrimGalore)", round((Number_Of_Reads_After_TrimGalore_From_Raw_Reads / Total_Raw_Reads) * 100, 2) if Total_Raw_Reads > 0 else 0))

    # 2. Contamination Metrics
    ordered_metrics.append(("Percentage of Clean Reads Mapped to Human (Kraken)", perc_Human))
    ordered_metrics.append(("Percentage of Clean Reads Mapped to Non-human (Kraken)", round(100 - perc_Human, 2)))
    ordered_metrics.append(("Percentage of Clean Reads Mapped to Bacteria (Kraken)", perc_Bacteria))
    ordered_metrics.append(("Number of Clean Reads Mapped to Escherichia coli", int(round(Number_of_Clean_Reads_Mapped_to_Escherichia_coli))))
    ordered_metrics.append(("Percentage of Clean Reads Mapped to Escherichia coli", pct_Ecoli))
    ordered_metrics.append(("Number of Clean Reads Mapped to Mycoplasma", int(round(Number_of_Clean_Reads_Mapped_to_Mycoplasma))))
    ordered_metrics.append(("Percentage of Clean Reads Mapped to Mycoplasma", pct_Mycoplasma))
    ordered_metrics.append(("Percentage of Trimmed Reads Mapped to rRNAs (RiboDetector)", Percentage_of_trimmed_Reads_Mapped_to_rRNA_RiboDetector))

    # 3. STAR Metrics and UMI-deduped Strand Metrics
    ordered_metrics.append(("Number of UMI-deduped Input Reads (STAR)", Number_of_UMI_deduped_Input_Reads_STAR))
    ordered_metrics.append(("Percentage of UMI-deduped Reads (STAR)", Percentage_of_UMI_deduped_Reads_STAR))
    ordered_metrics.append(("Average Mapped UMI-deduped Fragment Length (STAR)", Average_Mapped_UMI_deduped_Fragment_Length_STAR))
    ordered_metrics.append(("Percentage of UMI-deduped Reads on the Forward Strand", UMI_deduped_Forward_Strand))
    ordered_metrics.append(("Percentage of UMI-deduped Reads on the Reverse Strand", UMI_deduped_Reverse_Strand))
    ordered_metrics.append(("Percentage of UMI-deduped Reads with Failed Strand", UMI_deduped_Failed_Strand))
    ordered_metrics.append(("Number of UMI-deduped Reads Uniquely Mapped to Human (STAR)", Number_of_UMI_deduped_Reads_Uniquely_Mapped_to_Human_STAR))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Uniquely Mapped to Human (STAR)", Percentage_of_UMI_deduped_Reads_Uniquely_Mapped_to_Human_STAR))
    ordered_metrics.append(("Number of UMI-deduped Reads Multi-mapped to Human (STAR)", Number_of_UMI_deduped_Reads_Multi_mapped_to_Human_STAR))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Multi-mapped to Human (STAR)", Percentage_of_UMI_deduped_Reads_Multi_mapped_to_Human_STAR))
    ordered_metrics.append(("Number of Splices from UMI-deduped Reads Uniquely Mapped to Human (STAR)", Number_of_Splices_from_UMI_deduped_Reads_Uniquely_Mapped_to_Human_STAR))
    ordered_metrics.append(("Number of Human Unmapped UMI-deduped Reads (STAR)", Number_of_Human_Unmapped_UMI_deduped_Reads_STAR))
    ordered_metrics.append(("Number of Chimeric UMI-deduped Reads Mapped to Human (STAR)", Number_of_Chimeric_UMI_deduped_Reads_Mapped_to_Human_STAR))
    ordered_metrics.append(("Mismatch Rate per Base for UMI-deduped Reads (STAR)", Mismatch_Rate_per_Base_for_UMI_deduped_Reads_STAR))

    # 4. FeatureCounts Meta Gene Metrics - Counts
    ordered_metrics.append(("Number of UMI-deduped Reads Mapped to Human 3'UTR (featureCounts)", fc_3UTR))
    ordered_metrics.append(("Number of UMI-deduped Reads Mapped to Human 5'UTR (featureCounts)", fc_5UTR))
    ordered_metrics.append(("Number of UMI-deduped Reads Mapped to Human Downstream (featureCounts)", fc_downstream_2kb))
    ordered_metrics.append(("Number of UMI-deduped Reads Mapped to Human Exonic Regions (featureCounts)", fc_exon))
    ordered_metrics.append(("Number of UMI-deduped Reads Mapped to Human Intergenic Regions (featureCounts)", fc_intergenic))
    ordered_metrics.append(("Number of UMI-deduped Reads Mapped to Human Intronic Regions (featureCounts)", fc_intron))
    ordered_metrics.append(("Number of UMI-deduped Reads Mapped to Human Promoter Regions (featureCounts)", fc_promoter))
    ordered_metrics.append(("Number of UMI-deduped Reads Mapped to Human Blacklist Regions (featureCounts)", fc_encode_blacklist))

    # 5. FeatureCounts Meta Gene Metrics - Percentages
    ordered_metrics.append(("Percentage of UMI-deduped Reads Mapped to Human 3'UTR (featureCounts)", pct_fc_3UTR))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Mapped to Human 5'UTR (featureCounts)", pct_fc_5UTR))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Mapped to Human Downstream (featureCounts)", pct_fc_downstream))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Mapped to Human Exonic Regions (featureCounts)", pct_fc_exon))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Mapped to Human Intergenic Regions (featureCounts)", pct_fc_intergenic))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Mapped to Human Intronic Regions (featureCounts)", pct_fc_intron))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Mapped to Human Promoter Regions (featureCounts)", pct_fc_promoter))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Mapped to Human Blacklist Regions (featureCounts)", pct_fc_blacklist))

    # 6. Bedtools Mapping Metrics - Counts
    ordered_metrics.append(("Number of UMI-deduped Reads Mapped to Human 3'UTR (bedtools)", int(round(bedtools_counts.get("3'UTR", 0)))))
    ordered_metrics.append(("Number of UMI-deduped Reads Mapped to Human 5'UTR (bedtools)", int(round(bedtools_counts.get("5'UTR", 0)))))
    ordered_metrics.append(("Number of UMI-deduped Reads Mapped to Human Downstream (2 kb) Regions (bedtools)", int(round(bedtools_counts.get("Downstream (2 kb)", 0)))))
    ordered_metrics.append(("Number of UMI-deduped Reads Mapped to Human Exonic Regions (bedtools)", int(round(bedtools_counts.get("Exonic", 0)))))
    ordered_metrics.append(("Number of UMI-deduped Reads Mapped to Human Intergenic Regions (bedtools)", int(round(bedtools_counts.get("Intergenic", 0)))))
    ordered_metrics.append(("Number of UMI-deduped Reads Mapped to Human Intronic Regions (bedtools)", int(round(bedtools_counts.get("Intronic", 0)))))
    ordered_metrics.append(("Number of UMI-deduped Reads Mapped to Human Promoter Regions (bedtools)", int(round(bedtools_counts.get("Promoter", 0)))))
    ordered_metrics.append(("Number of UMI-deduped Reads Mapped to Human Blacklist Regions (bedtools)", int(round(bedtools_counts.get("ENCODE Blacklist", 0)))))

    # 7. Bedtools Mapping Metrics - Percentages
    ordered_metrics.append(("Percentage of UMI-deduped Reads Mapped to Human 3'UTR (bedtools)", pct_bt_3UTR))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Mapped to Human 5'UTR (bedtools)", pct_bt_5UTR))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Mapped to Human Downstream (2 kb) Regions (bedtools)", pct_bt_downstream))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Mapped to Human Exonic Regions (bedtools)", pct_bt_exonic))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Mapped to Human Intergenic Regions (bedtools)", pct_bt_intergenic))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Mapped to Human Intronic Regions (bedtools)", pct_bt_intronic))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Mapped to Human Promoter Regions (bedtools)", pct_bt_promoter))
    ordered_metrics.append(("Percentage of UMI-deduped Reads Mapped to Human Blacklist Regions (bedtools)", pct_bt_blacklist))

    # 8. PICARD Metrics (固定顺序)
    for key in picard_keys_order:
        ordered_metrics.append((key, picard_metrics.get(key, "NA")))

    # 9. Expression Matrix Metrics - Threshold-based Counts
    for rna in rna_types_order:
        for thresh in read_thresh_order:
            value = sum(1 for (rc, tpm) in expr_matrix.get(rna, []) if rc > thresh)
            ordered_metrics.append((f"Expressed Human {rna} (UMI-deduped Read Counts > {thresh})", value))
        for thresh in tpm_thresh_order:
            value = sum(1 for (rc, tpm) in expr_matrix.get(rna, []) if tpm > thresh)
            ordered_metrics.append((f"Expressed Human {rna} (UMI-deduped TPM/CPM > {thresh})", value))

    # 10. Expression Matrix Metrics - FeatureCounts Totals and Percentages
    for rna in rna_types_order:
        total_fc = sum(rc for (rc, tpm) in expr_matrix.get(rna, []))
        ordered_metrics.append((f"Number of UMI-deduped Reads Mapped to {rna} (featureCounts)", total_fc))
    total_fc_all = sum(sum(rc for (rc, tpm) in expr_matrix.get(rna, [])) for rna in rna_types_order)
    for rna in rna_types_order:
        total_fc = sum(rc for (rc, tpm) in expr_matrix.get(rna, []))
        pct_fc = round((total_fc / total_fc_all) * 100, 2) if total_fc_all > 0 else 0
        ordered_metrics.append((f"Percentage of UMI-deduped Reads Mapped to {rna} (featureCounts)", pct_fc))

    # ---------------------------
    # Write ordered metrics to output TSV
    with open(args.output, "w") as out:
        out.write("Metric\tValue\n")
        for metric, value in ordered_metrics:
            if value is not None:
                out.write(f"{metric}\t{value}\n")
    print(f"{args.output} generated.")

if __name__ == "__main__":
    main()

