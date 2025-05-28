#!/usr/bin/env python3
"""
QC_matrix.py

Generates a QC matrix table from multiple QC files.
The table (two columns: Metric and Value) includes metrics from FastQC (raw and trimmed),
Kraken, BBSplit, RiboDetector, UMI-dedup BAM (and BAM2strand), Picard, and expression filtering.
For the expression matrix (merged RNA distribution), for each of 21 RNA types, the number of genes
meeting nine thresholds (3 for read counts and 6 for TPM/CPM) is computed. If an RNA type is absent,
all counts are 0.
Output file: <sample>_QC_matrix.tsv
"""

import sys
import argparse
import zipfile
import pysam
import gzip

def parse_args():
    parser = argparse.ArgumentParser(description="Generate QC matrix table from multiple QC files.")
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--raw_R1_zip", required=True, help="Raw R1 FastQC zip file")
    parser.add_argument("--raw_R2_zip", required=True, help="Raw R2 FastQC zip file")
    parser.add_argument("--trimmed_R1_zip", required=True, help="Trimmed R1 FastQC zip file")
    parser.add_argument("--trimmed_R2_zip", required=True, help="Trimmed R2 FastQC zip file")
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
    parser.add_argument("--bedtools_cov_3UTR", required=True, help="bedtools 3′UTR coverage TSV")
    parser.add_argument("--bedtools_cov_5UTR", required=True, help="bedtools 5′UTR coverage TSV")
    parser.add_argument("--bedtools_cov_downstream_2kb", required=True, help="bedtools Downstream (2 kb) coverage TSV")
    parser.add_argument("--bedtools_cov_exon", required=True, help="bedtools Exon coverage TSV")
    parser.add_argument("--bedtools_cov_GENECODE_blacklist", required=True, help="bedtools GENECODE Blacklist coverage TSV")
    parser.add_argument("--bedtools_cov_intergenic", required=True, help="bedtools Intergenic coverage TSV")
    parser.add_argument("--bedtools_cov_intron", required=True, help="bedtools Intron coverage TSV")
    parser.add_argument("--bedtools_cov_promoter_1500_500bp", required=True, help="bedtools Promoter (1500-500bp) coverage TSV")
    parser.add_argument("--expression_matrix", required=True, help="Merged expression matrix TSV (5 columns: GeneID, GeneSymbol, GeneType, ReadCounts, TPM/CPM)")
    parser.add_argument("--output", required=True, help="Output QC matrix TSV file")
    return parser.parse_args()

def extract_basic_stats(fastqc_zip):
    """Extract Total Sequences, %GC, and Sequence length from fastqc_data.txt."""
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
    """Return average if value is a range; otherwise return float."""
    if "-" in val:
        parts = val.split("-")
        try:
            nums = [float(x) for x in parts]
            return sum(nums)/len(nums)
        except ValueError:
            return float(parts[0])
    else:
        try:
            return float(val)
        except ValueError:
            return 0

def extract_kraken_percentages(kraken_file):
    """Extract percentages for Homo sapiens and Bacteria from Kraken report."""
    perc_human = 0.0
    perc_bacteria = 0.0
    with open(kraken_file) as f:
        for line in f:
            line = line.strip()
            if not line: continue
            fields = line.split("\t")
            if len(fields) < 8: continue
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
    """Count reads in FASTQ file (4 lines per read, supports gzipped)."""
    count = 0
    opener = gzip.open if fastq_file.endswith(".gz") else open
    with opener(fastq_file, "rt") as f:
        for _ in f:
            count += 1
    return count // 4

def process_umi_dedup_bam(bam_file):
    """Process UMI-deduplicated BAM; return total unique, uniquely mapped, multiply mapped, and avg fragment length."""
    unique_reads = {}
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for aln in bam.fetch(until_eof=True):
            if aln.is_secondary or aln.is_supplementary:
                continue
            qname = aln.query_name
            if qname not in unique_reads:
                unique_reads[qname] = aln
    total_unique = len(unique_reads)
    uniquely_mapped = 0
    multiply_mapped = 0
    frag_lengths = []
    for aln in unique_reads.values():
        if aln.mapping_quality >= 30:
            uniquely_mapped += 1
        else:
            multiply_mapped += 1
        tmpl = abs(aln.template_length)
        if tmpl > 0:
            frag_lengths.append(tmpl)
    avg_frag_len = int(round(sum(frag_lengths)/len(frag_lengths))) if frag_lengths else 0
    return total_unique, uniquely_mapped, multiply_mapped, avg_frag_len

def parse_bam2strand(strand_file):
    """Parse BAM2strand file; return forward, reverse, and failed percentages."""
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
    return round(forward,2), round(reverse,2), round(failed,2)

def parse_picard_insert_size(file_path):
    """Parse Picard Insert Size Metrics file; return dict with PICARD_MEDIAN_INSERT_SIZE."""
    header = None
    data = None
    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith("#"): continue
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
    """Parse Picard RNA Metrics file; return dict with keys prefixed with 'PICARD_'."""
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
            if not line: continue
            if line.startswith("#"): continue
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
    Returns: total_read_counts (sum of col7), total_mapped_bases (sum of col8), total_length (sum of col9).
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
    Process the expression matrix (5 columns: GeneID, GeneSymbol, GeneType, ReadCounts, TPM/CPM).
    Returns a dictionary mapping RNA type to a list of tuples (read_count, tpm).
    """
    expr = {}
    with open(expr_file) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            line = line.strip()
            if not line: continue
            fields = line.split("\t")
            if len(fields) < 5: continue
            rna_type = fields[2].strip()
            try:
                rc = float(fields[3])
            except ValueError:
                rc = 0.0
            try:
                tpm = float(fields[4])
            except ValueError:
                tpm = 0.0
            if rna_type not in expr:
                expr[rna_type] = []
            expr[rna_type].append((rc, tpm))
    return expr

def main():
    args = parse_args()
    
    # Raw FastQC metrics.
    stats_R1 = extract_basic_stats(args.raw_R1_zip)
    stats_R2 = extract_basic_stats(args.raw_R2_zip)
    try:
        total_R1 = float(stats_R1.get("Total Sequences", 0))
    except ValueError:
        total_R1 = 0
    try:
        total_R2 = float(stats_R2.get("Total Sequences", 0))
    except ValueError:
        total_R2 = 0
    total_raw = total_R1 + total_R2
    try:
        gc_R1 = float(stats_R1.get("%GC", 0))
    except ValueError:
        gc_R1 = 0
    try:
        gc_R2 = float(stats_R2.get("%GC", 0))
    except ValueError:
        gc_R2 = 0
    avg_gc = (gc_R1 + gc_R2) / 2
    read_len_R1 = parse_sequence_length(stats_R1.get("Sequence length", "0"))
    read_len_R2 = parse_sequence_length(stats_R2.get("Sequence length", "0"))
    avg_read_length = (read_len_R1 + read_len_R2) / 2

    # Trimmed FastQC metrics.
    stats_trim_R1 = extract_basic_stats(args.trimmed_R1_zip)
    stats_trim_R2 = extract_basic_stats(args.trimmed_R2_zip)
    try:
        trimmed_total_R1 = float(stats_trim_R1.get("Total Sequences", 0))
    except ValueError:
        trimmed_total_R1 = 0
    try:
        trimmed_total_R2 = float(stats_trim_R2.get("Total Sequences", 0))
    except ValueError:
        trimmed_total_R2 = 0
    total_clean_reads = trimmed_total_R1 + trimmed_total_R2

    # Kraken percentages.
    perc_human, perc_bacteria = extract_kraken_percentages(args.kraken_report)

    # BBSplit FASTQ counts.
    ecoli_R1_count = count_fastq_reads(args.ecoli_R1_fastq)
    ecoli_R2_count = count_fastq_reads(args.ecoli_R2_fastq)
    total_ecoli = ecoli_R1_count + ecoli_R2_count

    myco_R1_count = count_fastq_reads(args.myco_R1_fastq)
    myco_R2_count = count_fastq_reads(args.myco_R2_fastq)
    total_myco = myco_R1_count + myco_R2_count

    # RiboDetector FASTQ counts.
    ribo_R1_count = count_fastq_reads(args.ribo_R1_fastq)
    ribo_R2_count = count_fastq_reads(args.ribo_R2_fastq)
    total_ribo = ribo_R1_count + ribo_R2_count

    # UMI-dedup BAM metrics.
    total_umi, uniquely_mapped, multiply_mapped, avg_frag_len = process_umi_dedup_bam(args.umi_dedup_bam)
    umi_mapped = total_umi

    # BAM2strand percentages.
    forward_pct, reverse_pct, failed_pct = parse_bam2strand(args.bam2strand_file)

    # Picard metrics.
    picard_insert = parse_picard_insert_size(args.picard_insert_file)
    picard_rnaseq = parse_picard_rnaseq_metrics(args.picard_rnaseq_file)
    picard_metrics = {}
    picard_metrics.update(picard_rnaseq)
    picard_metrics.update(picard_insert)

    # Bedtools metrics.
    bedtools_files = {
        "3′ UTR": args.bedtools_cov_3UTR,
        "5′ UTR": args.bedtools_cov_5UTR,
        "Downstream (2 kb)": args.bedtools_cov_downstream_2kb,
        "Exonic": args.bedtools_cov_exon,
        "GENECODE Blacklist": args.bedtools_cov_GENECODE_blacklist,
        "Intergenic": args.bedtools_cov_intergenic,
        "Intronic": args.bedtools_cov_intron,
        "Promoter": args.bedtools_cov_promoter_1500_500bp
    }
    bedtools_counts = {}
    bedtools_avg_cov = {}
    for region, file_path in bedtools_files.items():
        total_counts, total_mapped, total_len = process_bedtools_file(file_path)
        bedtools_counts[region] = total_counts
        avg_cov = total_mapped / total_len if total_len > 0 else 0
        bedtools_avg_cov[region] = avg_cov

    # Expression matrix filtering.
    # Expression matrix has columns: GeneID, GeneSymbol, GeneType, ReadCounts, TPM/CPM.
    expr_data = process_expression_matrix(args.expression_matrix)
    rna_types = ["artifact", "ERVs", "IG_genes", "LINEs", "lncRNAs", "miRNAs", "misc-sncRNAs",
                 "piRNAs", "protein_coding", "pseudogenes", "rRNAs", "scaRNAs", "SINEs",
                 "snoRNAs", "snRNAs", "TEC_protein_coding", "TR_genes", "tRNAs", "vault_RNAs",
                 "Y_RNAs", "circRNAs"]
    read_thresh = [1, 10, 5]
    tpm_thresh = [0.001, 0.01, 0.1, 0.2, 0.5, 1]
    expr_metrics = []
    # For each RNA type, if not present, treat as empty list.
    for rna in rna_types:
        gene_list = expr_data.get(rna, [])
        # For read count thresholds.
        for thresh in read_thresh:
            count = sum(1 for (rc, tpm) in gene_list if rc > thresh)
            expr_metrics.append((f"Expressed Human {rna} (UMI-dedup Read Counts > {thresh})", count))
        # For TPM thresholds.
        for thresh in tpm_thresh:
            count = sum(1 for (rc, tpm) in gene_list if tpm > thresh)
            expr_metrics.append((f"Expressed Human {rna} (UMI-dedup TPM/CPM > {thresh})", count))

    # Build overall QC metrics list.
    qc_metrics = []
    qc_metrics.append(("Total Raw Reads", int(round(total_raw))))
    qc_metrics.append(("%GC Total Raw Reads", avg_gc))
    qc_metrics.append(("Read length", int(round(avg_read_length))))
    qc_metrics.append(("Total Clean Reads (Post-TrimGalore)", int(round(total_clean_reads))))
    qc_metrics.append(("Percentage of Clean Reads Mapping to Human (Kraken)", perc_human))
    qc_metrics.append(("Percentage of Clean Reads Mapping to Bacteria (Kraken)", perc_bacteria))
    qc_metrics.append(("Clean Reads Mapped to Escherichia coli", int(round(total_ecoli))))
    qc_metrics.append(("Clean Reads Mapped to Mycoplasma", int(round(total_myco))))
    qc_metrics.append(("Clean Reads Mapped to rRNA (RiboDetector)", int(round(total_ribo))))
    qc_metrics.append(("UMI-Dedup Reads Mapped to Human (BWA)", int(round(umi_mapped))))
    qc_metrics.append(("UMI-Dedup Reads Uniquely Mapped to Human (BWA)", int(round(uniquely_mapped))))
    qc_metrics.append(("UMI-Dedup Reads Multiply Mapped to Human (BWA)", int(round(multiply_mapped))))
    qc_metrics.append(("UMI-Deduplicated Alignment Average Fragment Length", int(round(avg_frag_len))))
    qc_metrics.append(("Percentage of UMI-Dedup Reads on the Forward Strand", forward_pct))
    qc_metrics.append(("Percentage of UMI-Dedup Reads on the Reverse Strand", reverse_pct))
    qc_metrics.append(("Percentage of UMI-Dedup Reads with Failed Strand", failed_pct))
    # Append bedtools counts.
    qc_metrics.append(("UMI-Dedup Reads in Human 3′ UTR (bedtools)", int(round(bedtools_counts["3′ UTR"]))))
    qc_metrics.append(("UMI-Dedup Reads in Human 5′ UTR (bedtools)", int(round(bedtools_counts["5′ UTR"]))))
    qc_metrics.append(("UMI-Dedup Reads in Human Downstream (2 kb) Regions (bedtools)", int(round(bedtools_counts["Downstream (2 kb)"])) ))
    qc_metrics.append(("UMI-Dedup Reads in Human Exonic Regions (bedtools)", int(round(bedtools_counts["Exonic"])) ))
    qc_metrics.append(("UMI-Dedup Reads in Human Intergenic Regions (bedtools)", int(round(bedtools_counts["Intergenic"])) ))
    qc_metrics.append(("UMI-Dedup Reads in Human Intronic Regions (bedtools)", int(round(bedtools_counts["Intronic"])) ))
    qc_metrics.append(("UMI-Dedup Reads in Human Promoter Regions (bedtools)", int(round(bedtools_counts["Promoter"])) ))
    qc_metrics.append(("UMI-Dedup Reads in Human Blacklist Regions (bedtools)", int(round(bedtools_counts["GENECODE Blacklist"])) ))
    # Append bedtools average coverage.
    qc_metrics.append(("UMI-Dedup Reads AVG Coverage in Human 3′ UTR (bedtools)", bedtools_avg_cov["3′ UTR"]))
    qc_metrics.append(("UMI-Dedup Reads AVG Coverage in Human 5′ UTR (bedtools)", bedtools_avg_cov["5′ UTR"]))
    qc_metrics.append(("UMI-Dedup Reads AVG Coverage in Human Downstream (2 kb) Regions (bedtools)", bedtools_avg_cov["Downstream (2 kb)"]))
    qc_metrics.append(("UMI-Dedup Reads AVG Coverage in Human Exonic Regions (bedtools)", bedtools_avg_cov["Exonic"]))
    qc_metrics.append(("UMI-Dedup Reads AVG Coverage in Human Intergenic Regions (bedtools)", bedtools_avg_cov["Intergenic"]))
    qc_metrics.append(("UMI-Dedup Reads AVG Coverage in Human Intronic Regions (bedtools)", bedtools_avg_cov["Intronic"]))
    qc_metrics.append(("UMI-Dedup Reads AVG Coverage in Human Promoter Regions (bedtools)", bedtools_avg_cov["Promoter"]))
    qc_metrics.append(("UMI-Dedup Reads AVG Coverage in Human Blacklist Regions (bedtools)", bedtools_avg_cov["GENECODE Blacklist"]))
    # Append Picard metrics.
    desired_picard = ["PF_BASES", "PF_ALIGNED_BASES", "RIBOSOMAL_BASES", "CODING_BASES", "UTR_BASES",
                      "INTRONIC_BASES", "INTERGENIC_BASES", "IGNORED_READS", "CORRECT_STRAND_READS",
                      "INCORRECT_STRAND_READS", "NUM_R1_TRANSCRIPT_STRAND_READS", "NUM_R2_TRANSCRIPT_STRAND_READS",
                      "NUM_UNEXPLAINED_READS", "PCT_R1_TRANSCRIPT_STRAND_READS", "PCT_R2_TRANSCRIPT_STRAND_READS",
                      "PCT_RIBOSOMAL_BASES", "PCT_CODING_BASES", "PCT_UTR_BASES", "PCT_INTRONIC_BASES",
                      "PCT_INTERGENIC_BASES", "PCT_MRNA_BASES", "PCT_USABLE_BASES", "PCT_CORRECT_STRAND_READS",
                      "MEDIAN_CV_COVERAGE", "MEDIAN_5PRIME_BIAS", "MEDIAN_3PRIME_BIAS", "MEDIAN_5PRIME_TO_3PRIME_BIAS"]
    for key in desired_picard:
        metric_name = "PICARD_" + key
        value = picard_metrics.get(metric_name, "NA")
        qc_metrics.append((metric_name, value))
    qc_metrics.append(("PICARD_MEDIAN_INSERT_SIZE", picard_metrics.get("PICARD_MEDIAN_INSERT_SIZE", "NA")))
    
    # Process expression matrix filtering.
    expr_data = process_expression_matrix(args.expression_matrix)
    # RNA types to consider.
    rna_types = ["artifact", "ERVs", "IG_genes", "LINEs", "lncRNAs", "miRNAs", "misc-sncRNAs",
                 "piRNAs", "protein_coding", "pseudogenes", "rRNAs", "scaRNAs", "SINEs",
                 "snoRNAs", "snRNAs", "TEC_protein_coding", "TR_genes", "tRNAs", "vault_RNAs",
                 "Y_RNAs", "circRNAs"]
    read_thresh = [1, 10, 5]
    tpm_thresh = [0.001, 0.01, 0.1, 0.2, 0.5, 1]
    # For each RNA type, count genes that meet each threshold.
    for rna in rna_types:
        gene_list = expr_data.get(rna, [])
        # If no genes of this type, counts are zero.
        for thresh in read_thresh:
            count = sum(1 for (rc, tpm) in gene_list if rc > thresh)
            qc_metrics.append((f"Expressed Human {rna} (UMI-dedup Read Counts > {thresh})", count))
        for thresh in tpm_thresh:
            count = sum(1 for (rc, tpm) in gene_list if tpm > thresh)
            qc_metrics.append((f"Expressed Human {rna} (UMI-dedup TPM/CPM > {thresh})", count))
    
    # Write QC matrix table.
    with open(args.output, "w") as out:
        out.write("Metric\tValue\n")
        for metric, value in qc_metrics:
            out.write(f"{metric}\t{value}\n")
    print(f"{args.output} generated.")

if __name__ == "__main__":
    args = parse_args()
    main()

