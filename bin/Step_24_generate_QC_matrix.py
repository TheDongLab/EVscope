#!/usr/bin/env python3
"""
Step_24_generate_QC_matrix.py
"""

import sys
import argparse
import zipfile
import gzip
import pandas as pd
from statistics import mean

def safe_call(func, *args, default=None, **kwargs):
    try:
        return func(*args, **kwargs)
    except Exception:
        return default

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a comprehensive QC matrix from multiple pipeline output files."
    )
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit("Error: No input parameters provided.")

    parser.add_argument("--raw_fastqc_zips", nargs='+', required=False, help="List of raw FastQC zip files from Step 1.")
    parser.add_argument("--trimmed_fastqs", nargs='+', required=False, help="List of trimmed FASTQ files from Step 3.")
    parser.add_argument("--ecoli_fastqs", nargs='+', required=False, help="List of E.coli FASTQ files from BBSplit (Step 5).")
    parser.add_argument("--myco_fastqs", nargs='+', required=False, help="List of Mycoplasma FASTQ files from BBSplit (Step 5).")
    parser.add_argument("--ribo_fastqs", nargs='+', required=False, help="List of rRNA FASTQ files from RiboDetector (Step 23).")
    parser.add_argument("--ACC_motif_fraction", required=False, help="ACC motif fraction TSV file from Step 2.")
    parser.add_argument("--kraken_report", required=False, help="Kraken2 report file from Step 19.")
    parser.add_argument("--bam2strand_file", required=False, help="bam2strandness output file from Step 7.")
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
    with open(strand_file) as f:
        _ = f.readline()
        fields = f.readline().strip().split("\t")
    forward_pct = float(fields[2]) * 100
    reverse_pct = float(fields[3]) * 100
    failed_pct = float(fields[4]) * 100
    return round(forward_pct, 2), round(reverse_pct, 2), round(failed_pct, 2)

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
    with open(expr_file) as f:
        _ = f.readline()
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 5: continue
            gene_type, read_count, tpm_cpm = fields[2].strip(), fields[3], fields[4]
            if gene_type not in expr_data:
                expr_data[gene_type] = []
            expr_data[gene_type].append((safe_call(float, read_count, default=0.0), safe_call(float, tpm_cpm, default=0.0)))
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

def process_featureCounts_file(file_path):
    df = pd.read_csv(file_path, sep='\t', comment='#')
    return df.iloc[:, -1].sum()

def main():
    args = parse_args()
    ordered_metrics = []

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
    ordered_metrics.append(("ACC motif fraction from UMI region", acc_fraction))
    ordered_metrics.append(("Total Trimmed Reads (R1+R2)", int(total_trimmed_reads) if total_trimmed_reads > 0 else "NA"))
    ordered_metrics.append(("Average Trimmed Read Length", int(round(avg_trimmed_read_length)) if avg_trimmed_read_length > 0 else "NA"))
    ordered_metrics.append(("Percentage of Reads Remaining after Trimming", percent_reads_after_trimming if total_raw_reads > 0 else "NA"))

    perc_human, perc_bacteria = safe_call(extract_kraken_percentages, args.kraken_report, default=(0.0, 0.0))
    
    ecoli_reads = sum(safe_call(count_fastq_reads, f, default=0) for f in args.ecoli_fastqs or [])
    pct_ecoli = round((ecoli_reads / total_trimmed_reads) * 100, 2) if total_trimmed_reads > 0 else 0
    
    myco_reads = sum(safe_call(count_fastq_reads, f, default=0) for f in args.myco_fastqs or [])
    pct_myco = round((myco_reads / total_trimmed_reads) * 100, 2) if total_trimmed_reads > 0 else 0

    ribo_reads = sum(safe_call(count_fastq_reads, f, default=0) for f in args.ribo_fastqs or [])
    perc_rRNA_ribodetector = round((ribo_reads / total_trimmed_reads) * 100, 2) if total_trimmed_reads > 0 else 0

    ordered_metrics.append(("Percentage of Reads Mapped to Human (Kraken)", perc_human if perc_human > 0 else "NA"))
    ordered_metrics.append(("Percentage of Reads Mapped to Bacteria (Kraken)", perc_bacteria if perc_bacteria > 0 else "NA"))
    ordered_metrics.append(("Number of Reads Mapped to Escherichia coli (BBSplit)", int(ecoli_reads) if ecoli_reads > 0 else "NA"))
    ordered_metrics.append(("Percentage of Trimmed Reads Mapped to E. coli (BBSplit)", pct_ecoli if pct_ecoli > 0 else "NA"))
    ordered_metrics.append(("Number of Reads Mapped to Mycoplasma (BBSplit)", int(myco_reads) if myco_reads > 0 else "NA"))
    ordered_metrics.append(("Percentage of Trimmed Reads Mapped to Mycoplasma (BBSplit)", pct_myco if pct_myco > 0 else "NA"))
    ordered_metrics.append(("Percentage of Trimmed Reads Mapped to rRNAs (RiboDetector)", perc_rRNA_ribodetector if perc_rRNA_ribodetector > 0 else "NA"))

    star_metrics = safe_call(parse_star_log, args.STAR_log, default={})
    star_input_reads = int(star_metrics.get("Number of input reads (STAR)", 0))
    
    perc_reads_after_dedup = round((star_input_reads / initial_star_mapped_fragments) * 100, 2) if initial_star_mapped_fragments > 0 else 0
    
    star_unique_reads = int(star_metrics.get("Uniquely mapped reads number (STAR)", 0))
    perc_unique_vs_star_input = round((star_unique_reads / star_input_reads) * 100, 2) if star_input_reads > 0 else 0
    
    star_multi_reads = int(star_metrics.get("Multi-mapping reads number (STAR)", 0))
    perc_multi_vs_star_input = round((star_multi_reads / star_input_reads) * 100, 2) if star_input_reads > 0 else 0
    
    (forward_strand, reverse_strand, failed_strand) = safe_call(parse_bam2strand, args.bam2strand_file, default=(0, 0, 0))

    ordered_metrics.extend([
        ("Total Fragments Mapped to Human (First STAR)", int(initial_star_mapped_fragments) if initial_star_mapped_fragments > 0 else "NA"),
        ("Number of Reads after UMI-deduplication (STAR Input)", star_input_reads if star_input_reads > 0 else "NA"),
        ("Percentage of UMI-dedup Fragments", perc_reads_after_dedup if perc_reads_after_dedup > 0 else "NA"),
        ("Average Mapped Read Length (STAR)", float(star_metrics.get("Average input read length (STAR)", 0))),
        ("Percentage of Mapped Reads on Forward Strand", forward_strand if forward_strand > 0 else "NA"),
        ("Percentage of Mapped Reads on Reverse Strand", reverse_strand if reverse_strand > 0 else "NA"),
        ("Percentage of Mapped Reads with Failed Strand", failed_strand if failed_strand > 0 else "NA"),
        ("Number of Uniquely Mapped Reads (STAR)", star_unique_reads if star_unique_reads > 0 else "NA"),
        ("Percentage of Uniquely Mapped Reads (vs STAR Input)", perc_unique_vs_star_input if perc_unique_vs_star_input > 0 else "NA"),
        ("Number of Multi-mapped Reads (STAR)", star_multi_reads if star_multi_reads > 0 else "NA"),
        ("Percentage of Multi-mapped Reads (vs STAR Input)", perc_multi_vs_star_input if perc_multi_vs_star_input > 0 else "NA"),
        ("Number of Splices (from unique reads, STAR)", int(star_metrics.get("Number of splices from uniquely mapped reads (STAR)", 0))),
        ("Number of Unmapped Reads (STAR)", int(star_metrics.get("Unmapped reads number (STAR)", 0))),
        ("Number of Chimeric Reads (STAR)", int(star_metrics.get("Number of chimeric reads (STAR)", 0))),
        ("Mismatch Rate per Base (STAR)", star_metrics.get("Mismatch rate per base (STAR)", "NA"))
    ])
    
    fc_counts = {
        "3'UTR": safe_call(process_featureCounts_file, args.featureCounts_3UTR, default=0),
        "5'UTR": safe_call(process_featureCounts_file, args.featureCounts_5UTR, default=0),
        "Downstream": safe_call(process_featureCounts_file, args.featureCounts_downstream_2kb, default=0),
        "Exonic": safe_call(process_featureCounts_file, args.featureCounts_exon, default=0),
        "Intergenic": safe_call(process_featureCounts_file, args.featureCounts_intergenic, default=0),
        "Intronic": safe_call(process_featureCounts_file, args.featureCounts_intron, default=0),
        "Promoter": safe_call(process_featureCounts_file, args.featureCounts_promoter_1500_500bp, default=0),
        "Blacklist": safe_call(process_featureCounts_file, args.featureCounts_ENCODE_blacklist, default=0)
    }
    total_fc_meta = sum(fc_counts.values())
    for region, count in fc_counts.items():
        ordered_metrics.append((f"Number of Reads Mapped to {region} Regions", int(count) if count > 0 else "NA"))
    for region, count in fc_counts.items():
        pct = round((count / total_fc_meta) * 100, 2) if total_fc_meta > 0 else 0
        ordered_metrics.append((f"Percentage of Reads Mapped to {region} Regions", pct if pct > 0 else "NA"))

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
            ordered_metrics.append((f"Expressed {rna} Genes (Read Counts > {thresh})", count if count > 0 else "NA"))
        for thresh in tpm_thresh_order:
            count = sum(1 for _, tpm in rna_data if tpm > thresh)
            ordered_metrics.append((f"Expressed {rna} Genes (TPM/CPM > {thresh})", count if count > 0 else "NA"))

    total_expr_reads = sum(rc for rna_data in expr_matrix.values() for rc, _ in rna_data)
    for rna in rna_types_order:
        total_rc = sum(rc for rc, _ in expr_matrix.get(rna, []))
        ordered_metrics.append((f"Number of Reads Mapped to {rna}", total_rc if total_rc > 0 else "NA"))
    for rna in rna_types_order:
        total_rc = sum(rc for rc, _ in expr_matrix.get(rna, []))
        pct = round((total_rc / total_expr_reads) * 100, 2) if total_expr_reads > 0 else 0
        ordered_metrics.append((f"Percentage of Reads Mapped to {rna}", pct if pct > 0 else "NA"))

    with open(args.output, "w") as out:
        out.write("Metric\tValue\n")
        for metric, value in ordered_metrics:
            out.write(f"{metric}\t{value}\n")
    print(f"QC matrix successfully generated: {args.output}")

if __name__ == "__main__":
    main()