#!/usr/bin/env python3
"""
Calculate the fraction of reads in a FASTQ file with the ACC motif at specified positions.
Usage: python Step_02_calculate_ACC_motif_fraction.py --input_fastq <fastq_file> --positions <start-end> [--output_tsv <tsv_file>]
"""

import sys
import gzip
import argparse
import os
import re

def parse_positions(positions_str):
    """
    Parse the positions string (e.g., '9-11') into start and end positions (0-based).

    Args:
        positions_str (str): Position range in format 'start-end' (e.g., '9-11').

    Returns:
        tuple: (start, end) where start is 0-based and end is exclusive.

    Raises:
        ValueError: If positions format is invalid or start > end.
    """
    match = re.match(r'^(\d+)-(\d+)$', positions_str)
    if not match:
        raise ValueError(f"Invalid positions format: {positions_str}. Expected 'start-end' (e.g., '9-11').")
    start = int(match.group(1))
    end = int(match.group(2))
    if start < 1 or end < start:
        raise ValueError(f"Invalid positions: start={start}, end={end}. Start must be >= 1 and <= end.")
    return start - 1, end  # Convert to 0-based indexing

def calculate_motif_ratio(fastq_path, motif="ACC", start=8, end=11):
    """
    Calculate the number and fraction of reads with the specified motif at given positions.

    Args:
        fastq_path (str): Path to the input FASTQ file (gzipped or plain).
        motif (str): Motif to search for (default: "ACC").
        start (int): Start position (0-based, default: 8).
        end (int): End position (exclusive, default: 11).

    Returns:
        tuple: (total_reads, motif_reads, ratio) where ratio is motif_reads / total_reads.
    """
    total_reads = 0
    motif_reads = 0

    # Check if the file is gzipped and open accordingly
    if fastq_path.endswith(".gz"):
        f = gzip.open(fastq_path, "rt")  # Open gzipped file in text mode
    else:
        f = open(fastq_path, "r")

    with f:
        while True:
            header = f.readline()
            if not header:
                break  # End of file
            seq = f.readline().strip()
            plus = f.readline()
            qual = f.readline()

            total_reads += 1
            # Check if sequence is long enough and if the motif matches at the specified positions
            if len(seq) >= end:
                if seq[start:end] == motif:
                    motif_reads += 1

    ratio = motif_reads / total_reads if total_reads > 0 else 0
    return total_reads, motif_reads, ratio

def main():
    """
    Main function to parse arguments, calculate ACC motif ratio, and write output.
    """
    parser = argparse.ArgumentParser(
        description="Calculate the fraction of reads with ACC motif at specified positions in a FASTQ file."
    )
    parser.add_argument(
        "--input_fastq",
        required=True,
        help="Input FASTQ file (gzipped or plain, e.g., 1.fastq or 1.fastq.gz)"
    )
    parser.add_argument(
        "--positions",
        default="9-11",
        help="Position range for ACC motif (e.g., '9-11' for bases 9-11, default: '9-11')"
    )
    parser.add_argument(
        "--output_tsv",
        default=None,
        help="Output TSV file with columns fastq_name and fraction_ACC (optional, e.g., 1.tsv)"
    )
    args = parser.parse_args()

    # Parse positions
    try:
        start, end = parse_positions(args.positions)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    fastq_file = args.input_fastq
    total, count, ratio = calculate_motif_ratio(fastq_file, start=start, end=end)

    # Print results to stdout
    print(f"Total reads: {total}")
    print(f"Reads with bases {start+1}-{end} == 'ACC': {count}")
    print(f"Ratio: {ratio:.4f} ({ratio*100:.2f}%)")

    # Write TSV file if --output_tsv is provided
    if args.output_tsv:
        fastq_name = os.path.basename(fastq_file)
        with open(args.output_tsv, "w") as tsv_out:
            tsv_out.write("fastq_name\tfraction_ACC\n")
            tsv_out.write(f"{fastq_name}\t{ratio:.4f}\n")
        print(f"TSV output written to: {args.output_tsv}")

if __name__ == "__main__":
    main()
