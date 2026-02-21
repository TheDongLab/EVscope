#!/usr/bin/env python3
"""
Usage:
    python UMIAdapterTrimR1.py \
        --input_R1_fq 20250205_egAA215d42d44_EG12799_S114_R1_adapter_trimmed.fq.gz \
        --input_R2_fq 20250205_egAA215d42d44_EG12799_S114_R2_adapter_trimmed.fq.gz \
        --output_R1_fq 20250205_egAA215d42d44_EG12799_S114_R1_adapter_UMI_trimmed.fq.gz \
        --output_R2_fq 20250205_egAA215d42d44_EG12799_S114_R2_adapter_UMI_trimmed.fq.gz \
        --output_tsv 20250205_egAA215d42d44_EG12799_S114_R1_readthrough_UMI_trimming.log \
        --min-overlap 3 \
        --min-length 10 \
        --chunk-size 100000 \
        --error-rate 0.1

This script processes strictly paired-end FASTQ files to remove UMI-derived adapter contamination
from R1 reads by extracting the UMI from the R2 header. It generates a TSV log, pie charts
summarizing the trimming statistics, and a summary file with key metrics.
"""

import sys
import os
import gzip
import argparse
import logging
import time
import numpy as np
from numba import njit
import matplotlib.pyplot as plt

# Attempt to import pyfastx for fast FASTQ parsing; if unavailable, fallback to Biopython.
try:
    import pyfastx

    def fastq_reader(filename):
        """Yield tuples (read_id, sequence, quality) from a FASTQ file using pyfastx."""
        for read in pyfastx.Fastq(filename, build_index=False):
            yield read.name, read.seq, read.qual
except ImportError:
    from Bio.SeqIO.QualityIO import FastqGeneralIterator

    def fastq_reader(filename):
        """Yield tuples (read_id, sequence, quality) from a FASTQ file using Biopython."""
        with smart_open(filename) as f:
            for rec in FastqGeneralIterator(f):
                yield rec


def smart_open(filename, mode="rt"):
    """
    Open a file in text mode, using gzip if the file extension is .gz.
    """
    if filename.endswith(".gz"):
        return gzip.open(filename, mode)
    return open(filename, mode)


def reverse_complement(seq):
    """
    Compute the reverse complement of a nucleotide sequence.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(complement.get(base, 'N') for base in seq)[::-1]

@njit(cache=True, fastmath=True)
def numba_trim(seq_arr, adapter_arr, min_overlap, error_rate):
    seq_len = len(seq_arr)
    adapter_len = len(adapter_arr)
    max_possible = seq_len if seq_len < adapter_len else adapter_len
    for L in range(max_possible, min_overlap - 1, -1):
        allowed_mismatches = int(L * error_rate)
        mismatches = 0
        for i in range(L):
            if seq_arr[seq_len - L + i] != adapter_arr[i]:
                mismatches += 1
                if mismatches > allowed_mismatches:
                    break
        if mismatches <= allowed_mismatches:
            return seq_len - L, L
    return seq_len, 0

def trim_by_adapter(seq, adapter, min_overlap, error_rate):
    seq_bytes = seq.encode('ascii')
    adapter_bytes = adapter.encode('ascii')
    seq_arr = np.frombuffer(seq_bytes, dtype=np.uint8)
    adapter_arr = np.frombuffer(adapter_bytes, dtype=np.uint8)
    new_length, trimmed_length = numba_trim(seq_arr, adapter_arr, min_overlap, error_rate)
    return seq[:new_length], trimmed_length


def modify_header(header, new_umi):
    """
    Modifies a FASTQ header to replace only the UMI portion with a new UMI (first 8 bp),
    while preserving the rest of the header exactly as in the original.
    """
    # Remove the leading '@' if present
    header_content = header.lstrip('@')
    # Separate token and any trailing description (including the space)
    if ' ' in header_content:
        token, rest = header_content.split(' ', 1)
        rest = ' ' + rest
    else:
        token = header_content
        rest = ''
    # Split token at the last underscore to separate prefix and old UMI
    if '_' in token:
        prefix, _ = token.rsplit('_', 1)
        new_token = f"{prefix}_{new_umi}"
    else:
        new_token = token
    # Reconstruct header
    return '@' + new_token + rest


def paired_chunker(r1_iter, r2_iter, chunk_size):
    while True:
        chunk = []
        for _ in range(chunk_size):
            try:
                r1 = next(r1_iter)
                r2 = next(r2_iter)
                chunk.append((r1, r2))
            except StopIteration:
                break
        if chunk:
            yield chunk
        else:
            break


def process_paired_chunk(chunk, min_overlap, error_rate, min_length):
    out1_lines = []
    out2_lines = []
    tsv_lines = []
    stats = {"count": 0, "trimmed": 0, "discarded": 0}

    for r1, r2 in chunk:
        stats["count"] += 1
        r1_id, r1_seq, r1_qual = r1
        r2_id, r2_seq, r2_qual = r2

        # Extract UMI from R2 header
        r2_header = r2_id.lstrip('@')
        token = r2_header.split()[0]
        if '_' in token:
            _, umi_full = token.split('_', 1)
        else:
            umi_full = ""

        # Discard if UMI too short
        if len(umi_full) < 14:
            stats["discarded"] += 1
            continue

        # Compute adapter and new UMI
        umi_adapter = umi_full[:14]
        computed_adapter = reverse_complement(umi_adapter)
        new_umi = umi_adapter[:8]

        # Trim R1
        trimmed_r1_seq, trimmed_length = trim_by_adapter(r1_seq, computed_adapter, min_overlap, error_rate)
        trimmed_r1_qual = r1_qual[:len(trimmed_r1_seq)]

        # Discard pairs below min_length
        if len(trimmed_r1_seq) < min_length or len(r2_seq) < min_length:
            stats["discarded"] += 1
            continue

        # Log trimming if occurred
        if trimmed_length > 0:
            modified_r1_name = modify_header(r1_id, new_umi)
            tsv_lines.append(f"{modified_r1_name.lstrip('@')}\t{computed_adapter}\t{r1_seq[-trimmed_length:]}\t{trimmed_length}\n")
            stats["trimmed"] += 1

        # Modify headers for output
        new_r1_id = modify_header(r1_id, new_umi)
        new_r2_id = modify_header(r2_id, new_umi)

        out1_lines.append(f"{new_r1_id}\n{trimmed_r1_seq}\n+\n{trimmed_r1_qual}\n")
        out2_lines.append(f"{new_r2_id}\n{r2_seq}\n+\n{r2_qual}\n")

    return out1_lines, out2_lines, stats, tsv_lines


def generate_pie_chart(stats, output_prefix):
    """
    Generates pie charts (PNG and PDF) for read statistics:
    1. No_UMI_readthrough (Untrimmed valid)
    2. UMI_readthrough (Trimmed valid)
    3. Too_short (Discarded)
    """
    total = stats['count']
    discarded = stats['discarded']      # Reads removed (too short)
    trimmed = stats['trimmed']          # Reads kept and trimmed
    no_umi = total - discarded - trimmed # Reads kept and not trimmed

    labels = ['No_UMI_readthrough', 'UMI_readthrough', 'Too_short']
    sizes = [no_umi, trimmed, discarded]
    colors = ['#66b3ff', '#ff9999', '#99ff99'] # Blue, Red, Green
    explode = (0.05, 0.05, 0.05)

    plt.figure(figsize=(8, 6))
    plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=False, startangle=140)
    plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.title(f"UMI Adapter Trimming Stats\nTotal Reads: {total}")

    # Save plots
    base_name = os.path.basename(output_prefix)
    # Remove extension if present in prefix logic, or just append
    if base_name.endswith('.log') or base_name.endswith('.tsv'):
        base_name = os.path.splitext(base_name)[0]

    # Construct output filenames
    # If the user passed a full path for TSV, we try to use that directory, otherwise local
    dir_name = os.path.dirname(output_prefix) if os.path.dirname(output_prefix) else "."

    save_base = os.path.join(dir_name, base_name)

    png_out = f"{save_base}_pie.png"
    pdf_out = f"{save_base}_pie.pdf"

    plt.savefig(png_out, bbox_inches='tight', dpi=300)
    plt.savefig(pdf_out, bbox_inches='tight')
    plt.close()

    print(f"Pie charts generated:\n  {png_out}\n  {pdf_out}")


def main():
    start_time = time.time()
    parser = argparse.ArgumentParser(
        description="Process strictly paired-end FASTQ files for UMI-based adapter trimming on R1."
    )
    parser.add_argument("--input_R1_fq", required=True,
                        help="Input FASTQ file for R1 (adapter trimmed; gzipped if applicable)")
    parser.add_argument("--input_R2_fq", required=True,
                        help="Input FASTQ file for R2 (adapter trimmed; gzipped if applicable)")
    parser.add_argument("--output_R1_fq", required=True,
                        help="Output FASTQ file for R1 (UMI trimmed)")
    parser.add_argument("--output_R2_fq", required=True,
                        help="Output FASTQ file for R2 (UMI trimmed)")
    parser.add_argument("--output_tsv", default=None,
                        help="Output TSV file for UMI adapter trimming details (default: derived from input_R1 filename)")
    parser.add_argument("--min-overlap", type=int, default=3,
                        help="Minimum overlap between the adapter and R1 sequence for trimming (default: 3)")
    parser.add_argument("--error-rate", type=float, default=0.1,
                        help="Maximum allowed mismatch rate during adapter matching (default: 0.1)")
    parser.add_argument("--min-length", type=int, default=10,
                        help="Minimum allowed sequence length after trimming (default: 10)")
    parser.add_argument("--chunk-size", type=int, default=100000,
                        help="Number of read pairs processed per chunk (default: 100000)")
    args = parser.parse_args()

    if args.output_tsv is None:
        base = os.path.basename(args.input_R1_fq)
        args.output_tsv = f"{base}_UMI_trimming.log"

    logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
    logging.info("Starting the paired-end UMI trimming process...")

    r1_iter = fastq_reader(args.input_R1_fq)
    r2_iter = fastq_reader(args.input_R2_fq)

    total_stats = {"count": 0, "trimmed": 0, "discarded": 0}
    all_tsv_lines = []

    with smart_open(args.output_R1_fq, "wt") as out1_handle, \
         smart_open(args.output_R2_fq, "wt") as out2_handle:
        for chunk in paired_chunker(r1_iter, r2_iter, args.chunk_size):
            chunk_out1, chunk_out2, stats, tsv_lines = process_paired_chunk(
                chunk, args.min_overlap, args.error_rate, args.min_length
            )
            out1_handle.write("".join(chunk_out1))
            out2_handle.write("".join(chunk_out2))
            all_tsv_lines.extend(tsv_lines)
            for key in total_stats:
                total_stats[key] += stats.get(key, 0)

    try:
        with smart_open(args.output_tsv, "wt") as tsv_handle:
            header = "read_name\tcomputed_adapter\ttrimmed_adapter\tadapter_trim_length\n"
            tsv_handle.write(header)
            tsv_handle.write("".join(all_tsv_lines))
    except Exception as e:
        logging.error(f"Failed to write TSV file: {e}")
        raise

    paired_output = total_stats["count"] - total_stats["discarded"]
    print(f"Total read pairs processed: {total_stats['count']}")
    print(f"Read pairs output (after filtering by min length): {paired_output}")
    print(f"Read pairs with trimmed R1 (adapter contamination removed): {total_stats['trimmed']}")
    print(f"Discarded read pairs (due to insufficient length or UMI issues): {total_stats['discarded']}")

    # Generate Pie Charts
    generate_pie_chart(total_stats, args.output_tsv)

    # =========================================================================
    # Generate Summary File
    # =========================================================================
    
    # 1. Define Variables based on total_stats
    Input_read_pairs_number = total_stats['count']
    Output_read_pairs_number = paired_output  # Calculated above as count - discarded
    Read_pairs_num_with_trimmed_R1_UMI_readthrough = total_stats['trimmed']
    
    # Avoid division by zero
    if Output_read_pairs_number > 0:
        R1_UMI_readthrough_ratio = Read_pairs_num_with_trimmed_R1_UMI_readthrough / Output_read_pairs_number
    else:
        R1_UMI_readthrough_ratio = 0.0

    # 2. Determine output filename based on --output_tsv logic
    # If args.output_tsv ends with .log or .tsv, replace extension with .summary
    summary_filename = args.output_tsv
    if summary_filename.endswith(".log"):
        summary_filename = summary_filename[:-4] + ".summary"
    elif summary_filename.endswith(".tsv"):
        summary_filename = summary_filename[:-4] + ".summary"
    else:
        summary_filename = summary_filename + ".summary"

    try:
        with open(summary_filename, "w") as summary_handle:
            # Header
            summary_handle.write(
                "Input_read_pairs_number\tOutput_read_pairs_number\t"
                "Read_pairs_num_with_trimmed_R1_UMI_readthrough\tR1_UMI_readthrough_ratio\n"
            )
            # Data Line
            summary_handle.write(
                f"{Input_read_pairs_number}\t{Output_read_pairs_number}\t"
                f"{Read_pairs_num_with_trimmed_R1_UMI_readthrough}\t{R1_UMI_readthrough_ratio:.6f}\n"
            )
        print(f"Summary statistics file generated: {summary_filename}")
    except Exception as e:
        logging.error(f"Failed to write summary file: {e}")

    end_time = time.time()
    runtime = end_time - start_time
    print(f"Total running time: {runtime:.2f} seconds")

if __name__ == "__main__":
    main()
