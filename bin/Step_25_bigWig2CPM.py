#!/usr/bin/env python3
"""
bigWig2CPM.py: Calculate mean per-base CPM for genes from BigWig files.

Computes Counts Per Million (CPM) for each gene based on strand-specific or unstranded BigWig coverage files and a GTF or BED annotation. Outputs a TSV file with gene_id and mean CPM (MCPM).

Version: 1.0.0
"""

import os
import time
import sys
import logging
import argparse
import pyBigWig
import numpy as np
import re
import tempfile
from typing import Dict, List, Tuple, Optional
from collections import defaultdict

# Script version
__version__ = "1.0.0"

# Global logger
global_logger: logging.Logger = None

###############################################################################
# Helper Functions
###############################################################################
def setup_logging(log_path: str) -> logging.Logger:
    """Configure logging to file and console."""
    logger = logging.getLogger('bigWig2CPM')
    logger.handlers.clear()
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s', "%Y-%m-%d %H:%M:%S")
    file_handler = logging.FileHandler(log_path, mode='w')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    return logger

def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=f"bigWig2CPM v{__version__}: Calculate gene CPM from BigWig files.")
    parser.add_argument("--stranded", choices=["forward", "reverse", "no"], required=True,
                        help="Strand specificity: forward, reverse, or no (unstranded).")
    parser.add_argument("--output", required=True, help="Output TSV file for gene_id and MCPM.")
    parser.add_argument("--input_F1R2_bw", required=True, help="BigWig file for F1R2 strand coverage.")
    parser.add_argument("--input_F2R1_bw", required=True, help="BigWig file for F2R1 strand coverage.")
    parser.add_argument("--gtf", help="GTF annotation file with gene features.")
    parser.add_argument("--bed", help="BED6 annotation file with gene features.")
    parser.add_argument("--version", action="version", version=f"bigWig2CPM {__version__}")
    # Ensure at least one of --gtf or --bed is provided
    args = parser.parse_args()
    if not args.gtf and not args.bed:
        parser.error("At least one of --gtf or --bed must be provided.")
    if args.gtf and args.bed:
        parser.error("Only one of --gtf or --bed can be provided.")
    return args

def validate_file(file_path: str, logger: logging.Logger, file_type: str) -> None:
    """Validate existence and non-empty status of a file."""
    if not os.path.exists(file_path):
        logger.error(f"{file_type} file {os.path.basename(file_path)} not found")
        raise FileNotFoundError(f"{file_type} file {file_path} not found")
    if os.path.getsize(file_path) == 0:
        logger.error(f"{file_type} file {os.path.basename(file_path)} is empty")
        raise ValueError(f"{file_type} file {file_path} is empty")

def read_gene_intervals(gtf_path: Optional[str], bed_path: Optional[str], logger: logging.Logger) -> Dict[str, List[Tuple[int, int, str, str]]]:
    """Read GTF or BED file and extract gene intervals (0-based). Returns {chrom: [(start, end, strand, gene_id)]}."""
    gene_intervals = defaultdict(list)
    
    if gtf_path:
        logger.info(f"Reading GTF file {os.path.basename(gtf_path)}")
        try:
            with open(gtf_path, 'r') as gtf_file:
                for line in gtf_file:
                    if line.startswith('#'):
                        continue
                    fields = line.strip().split('\t')
                    if len(fields) < 9:
                        logger.warning(f"Invalid GTF line: {line.strip()}")
                        continue
                    if fields[2] != 'gene':
                        continue
                    chrom, start, end, strand = fields[0], fields[3], fields[4], fields[6]
                    try:
                        start = int(start) - 1  # Convert to 0-based
                        end = int(end)
                        if start < 0 or end <= start:
                            logger.warning(f"Invalid interval in GTF: {chrom}:{start+1}-{end}")
                            continue
                        if strand not in ['+', '-']:
                            logger.warning(f"Invalid strand in GTF: {strand}")
                            continue
                        # Parse attributes robustly
                        attributes = {}
                        attr_string = fields[8].strip()
                        attr_pairs = [pair.strip() for pair in attr_string.split(';') if pair.strip()]
                        for pair in attr_pairs:
                            match = re.match(r'(\w+)\s+"([^"]+)"', pair)
                            if match:
                                key, value = match.groups()
                                attributes[key] = value
                            else:
                                match = re.match(r'(\w+)\s+([^"\s]+)', pair)
                                if match:
                                    key, value = match.groups()
                                    attributes[key] = value
                                else:
                                    logger.debug(f"Skipping malformed attribute: {pair}")
                        gene_id = attributes.get('gene_id', '')
                        if not gene_id:
                            logger.warning(f"No gene_id found in GTF line: {line.strip()}")
                            continue
                        gene_intervals[chrom].append((start, end, strand, gene_id))
                    except (ValueError, IndexError) as e:
                        logger.warning(f"Error parsing GTF line: {line.strip()} - {e}")
                        continue
        except Exception as e:
            logger.error(f"Failed to read GTF file: {e}")
            raise
    
    elif bed_path:
        logger.info(f"Reading BED file {os.path.basename(bed_path)}")
        try:
            with open(bed_path, 'r') as bed_file:
                for line in bed_file:
                    if line.startswith('#'):
                        continue
                    fields = line.strip().split('\t')
                    if len(fields) < 6:
                        logger.warning(f"Invalid BED line (expected 6+ columns): {line.strip()}")
                        continue
                    chrom, start, end, gene_id, score, strand = fields[:6]
                    try:
                        start = int(start)  # BED is 0-based
                        end = int(end)
                        if start < 0 or end <= start:
                            logger.warning(f"Invalid interval in BED: {chrom}:{start}-{end}")
                            continue
                        if strand not in ['+', '-']:
                            logger.warning(f"Invalid strand in BED: {strand}")
                            continue
                        if not gene_id:
                            logger.warning(f"No gene_id found in BED line: {line.strip()}")
                            continue
                        gene_intervals[chrom].append((start, end, strand, gene_id))
                    except (ValueError, IndexError) as e:
                        logger.warning(f"Error parsing BED line: {line.strip()} - {e}")
                        continue
        except Exception as e:
            logger.error(f"Failed to read BED file: {e}")
            raise
    
    if not gene_intervals:
        logger.error("No valid gene entries found in input file")
        raise ValueError("No valid gene entries found in input file")
    logger.info(f"Extracted {sum(len(intervals) for intervals in gene_intervals.values())} genes from input")
    return gene_intervals

def compute_total_reads(bigwig_f1r2_path: str, bigwig_f2r1_path: str, logger: logging.Logger) -> float:
    """Compute total reads from F1R2 and F2R1 BigWig files."""
    logger.info("Computing total reads from BigWig files")
    total_reads = 0.0
    try:
        with pyBigWig.open(bigwig_f1r2_path) as bw_f1r2, pyBigWig.open(bigwig_f2r1_path) as bw_f2r1:
            chroms = set(bw_f1r2.chroms().keys()) | set(bw_f2r1.chroms().keys())
            for chrom in chroms:
                length = max(bw_f1r2.chroms().get(chrom, 0), bw_f2r1.chroms().get(chrom, 0))
                if length == 0:
                    continue
                f1r2_sum = bw_f1r2.stats(chrom, 0, length, type="sum", nBins=1)[0] or 0.0
                f2r1_sum = bw_f2r1.stats(chrom, 0, length, type="sum", nBins=1)[0] or 0.0
                total_reads += f1r2_sum + f2r1_sum
        if total_reads <= 0:
            logger.error("Total reads is zero or negative")
            raise ValueError("Total reads is zero or negative")
        logger.info(f"Total reads: {total_reads:.2f}")
        return total_reads
    except Exception as e:
        logger.error(f"Failed to compute total reads from BigWig files: {e}")
        raise

def calculate_gene_cpm(
    gene_intervals: Dict[str, List[Tuple[int, int, str, str]]],
    bigwig_f1r2_path: str,
    bigwig_f2r1_path: str,
    total_reads: float,
    stranded_mode: str,
    logger: logging.Logger
) -> Dict[str, float]:
    """Calculate mean per-base CPM (MCPM) for each gene."""
    logger.info("Calculating CPM values")
    gene_cpm = {}
    try:
        with pyBigWig.open(bigwig_f1r2_path) as bw_f1r2, pyBigWig.open(bigwig_f2r1_path) as bw_f2r1:
            for chrom, intervals in sorted(gene_intervals.items()):
                if chrom not in bw_f1r2.chroms() and chrom not in bw_f2r1.chroms():
                    logger.warning(f"Chromosome {chrom} not in BigWig files; skipping")
                    continue
                for start, end, strand, gene_id in intervals:
                    length = end - start
                    if length <= 0:
                        logger.warning(f"Invalid gene length for {gene_id} on {chrom}:{start}-{end}")
                        continue
                    if stranded_mode == "no":
                        f1r2_sum = bw_f1r2.stats(chrom, start, end, type="sum", nBins=1)[0] or 0.0
                        f2r1_sum = bw_f2r1.stats(chrom, start, end, type="sum", nBins=1)[0] or 0.0
                        coverage_sum = f1r2_sum + f2r1_sum
                    elif stranded_mode == "forward":
                        if strand == '+':
                            coverage_sum = bw_f1r2.stats(chrom, start, end, type="sum", nBins=1)[0] or 0.0
                        else:  # strand == '-'
                            coverage_sum = bw_f2r1.stats(chrom, start, end, type="sum", nBins=1)[0] or 0.0
                    elif stranded_mode == "reverse":
                        if strand == '+':
                            coverage_sum = bw_f2r1.stats(chrom, start, end, type="sum", nBins=1)[0] or 0.0
                        else:  # strand == '-'
                            coverage_sum = bw_f1r2.stats(chrom, start, end, type="sum", nBins=1)[0] or 0.0
                    else:
                        logger.error(f"Invalid stranded mode: {stranded_mode}")
                        raise ValueError(f"Invalid stranded mode: {stranded_mode}")
                    # Calculate CPM sum: (sum of coverage / total_reads) * 1e6
                    cpm_sum = (coverage_sum / total_reads) * 1e6
                    # Calculate MCPM: cpm_sum / length
                    mean_cpm = cpm_sum / length if length > 0 else 0.0
                    gene_cpm[gene_id] = mean_cpm
                    #logger.debug(f"Gene {gene_id} on {chrom}:{start}-{end} ({strand}): MCPM={mean_cpm:.4f}")
        return gene_cpm
    except Exception as e:
        logger.error(f"Failed to calculate gene CPM: {e}")
        raise

def write_cpm_output(gene_cpm: Dict[str, float], output_path: str, logger: logging.Logger) -> None:
    """Write gene_id and MCPM to TSV file atomically."""
    logger.info(f"Writing CPM results to {os.path.basename(output_path)}")
    try:
        with tempfile.NamedTemporaryFile(mode='w', delete=False, dir=os.path.dirname(output_path), suffix='.tmp') as temp_file:
            temp_file.write("gene_id\tMCPM\n")
            for gene_id, mcpm in sorted(gene_cpm.items()):
                temp_file.write(f"{gene_id}\t{mcpm:.6f}\n")
        os.replace(temp_file.name, output_path)
        logger.info(f"Successfully wrote {len(gene_cpm)} genes to {os.path.basename(output_path)}")
    except IOError as e:
        logger.error(f"Failed to write output TSV {os.path.basename(output_path)}: {e}")
        raise

def main():
    """Main function."""
    global global_logger
    try:
        args = parse_arguments()
        output_directory = os.path.dirname(os.path.abspath(args.output))
        os.makedirs(output_directory, exist_ok=True)
        log_path = os.path.join(output_directory, f"{os.path.splitext(os.path.basename(args.output))[0]}.log")
        global_logger = setup_logging(log_path)
        logger = global_logger
        logger.info(f"bigWig2CPM v{__version__}")
        logger.info(f"Command: {' '.join(sys.argv)}")
        start_time = time.time()

        # Validate input files
        validate_file(args.input_F1R2_bw, logger, "F1R2 BigWig")
        validate_file(args.input_F2R1_bw, logger, "F2R1 BigWig")
        if args.gtf:
            validate_file(args.gtf, logger, "GTF")
        elif args.bed:
            validate_file(args.bed, logger, "BED")

        # Read gene intervals from GTF or BED
        gene_intervals = read_gene_intervals(args.gtf, args.bed, logger)

        # Compute total reads
        total_reads = compute_total_reads(args.input_F1R2_bw, args.input_F2R1_bw, logger)

        # Calculate gene CPM
        gene_cpm = calculate_gene_cpm(
            gene_intervals, args.input_F1R2_bw, args.input_F2R1_bw, total_reads, args.stranded, logger
        )

        # Write output
        write_cpm_output(gene_cpm, args.output, logger)

        logger.info(f"Total runtime: {time.time() - start_time:.2f} seconds")
        logger.info("bigWig2CPM completed successfully")
    except Exception as e:
        if global_logger:
            global_logger.exception(f"Error in main: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()