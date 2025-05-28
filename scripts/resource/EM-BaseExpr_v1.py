#!/usr/bin/env python3
"""
Optimized EM-BaseExpr_v1.py

This script calculates per-base coverage from a BAM file using both uniquely mapped and multi-mapped reads.
It implements an Expectation-Maximization (EM) algorithm for fractional assignment of multi-mapped reads.
Modifications in this version address:
  1. Clear differentiation of fully mapped, partially mapped, and fully unmapped paired-end reads.
  2. Explicit memory release of temporary coverage arrays in EM iterations.
  3. Input file existence check.
  4. Combined bedGraph writing function to avoid duplicate code.
  5. Robust fraction_dict initialization.
  
Author: (Your Name)
Date: (Date)
"""

import os
import sys
import time
import logging
import numpy as np
import pysam
from collections import defaultdict
import argparse

###############################################################################
# 1) Logging Setup
###############################################################################
def setup_logging(log_path):
    """
    Setup logging to both console and file.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s', "%Y-%m-%d %H:%M:%S")

    # File handler with DEBUG level
    fh = logging.FileHandler(log_path, mode='w')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    # Console handler set to DEBUG so that detailed logs are shown
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    return logger

###############################################################################
# 2) ChunkedCoverage Class
###############################################################################
class ChunkedCoverage:
    CHUNK_SIZE = 5_000_000  # Each chunk covers 5 million bases

    def __init__(self, total_length):
        """
        Initialize a chunked coverage array for a reference sequence of the given total_length.
        """
        self.total_length = total_length
        self.num_chunks = (total_length + self.CHUNK_SIZE - 1) // self.CHUNK_SIZE
        self.chunks = [
            np.zeros(min(self.CHUNK_SIZE, total_length - i * self.CHUNK_SIZE), dtype=np.float32)
            for i in range(self.num_chunks)
        ]

    def __len__(self):
        return self.total_length

    def _add_value_to_range(self, start, end, value):
        """
        Add a given value to the coverage over the interval [start, end), splitting across chunks as needed.
        """
        if end <= start:
            return
        start = max(start, 0)
        end = min(end, self.total_length)
        while start < end:
            chunk_index = start // self.CHUNK_SIZE
            chunk_offset = start % self.CHUNK_SIZE
            space_in_chunk = self.CHUNK_SIZE - chunk_offset
            length_to_add = min(space_in_chunk, end - start)
            if length_to_add <= 0:
                break
            self.chunks[chunk_index][chunk_offset:chunk_offset + length_to_add] += value
            start += length_to_add

    def add_one_coverage(self, start, end):
        """
        Add a coverage value of 1 for the interval [start, end).
        """
        self._add_value_to_range(start, end, 1.0)

    def add_fractional_coverage(self, start, end, fraction):
        """
        Add a fractional coverage value for the interval [start, end).
        """
        if fraction <= 0:
            return
        self._add_value_to_range(start, end, fraction)

    def sum_range(self, start, end):
        """
        Sum the coverage values over the interval [start, end), handling chunk boundaries correctly.
        """
        if end <= start:
            return 0.0
        start = max(start, 0)
        end = min(end, self.total_length)
        first_chunk = start // self.CHUNK_SIZE
        last_chunk = (end - 1) // self.CHUNK_SIZE

        if first_chunk == last_chunk:
            arr = self.chunks[first_chunk]
            start_offset = start % self.CHUNK_SIZE
            end_offset = end % self.CHUNK_SIZE if end % self.CHUNK_SIZE != 0 else self.CHUNK_SIZE
            return float(np.sum(arr[start_offset:end_offset]))

        coverage_sum = 0.0
        # Sum partial block from the first chunk
        arr = self.chunks[first_chunk]
        coverage_sum += np.sum(arr[start % self.CHUNK_SIZE:])

        # Sum full chunks in between
        for chunk_idx in range(first_chunk + 1, last_chunk):
            coverage_sum += np.sum(self.chunks[chunk_idx])

        # Sum partial block from the last chunk
        end_offset = end % self.CHUNK_SIZE if end % self.CHUNK_SIZE != 0 else self.CHUNK_SIZE
        coverage_sum += np.sum(self.chunks[last_chunk][:end_offset])
        return float(coverage_sum)

    def copy_array(self):
        """
        Flatten and return the coverage array (used for bedGraph output).
        """
        return np.concatenate(self.chunks)

    def free_memory(self):
        """
        Release memory by clearing the chunk arrays.
        """
        self.chunks = []

###############################################################################
# 3) Gathering Alignments from BAM
###############################################################################
def gather_alignments(args, logger):
    """
    Open the BAM file, classify alignments as unique or multi-mapped,
    and record maximum coordinate per chromosome.
    Applies paired-end insert size filtering (if applicable) and uses get_blocks() for spliced reads.
    Also, explicitly counts fully mapped, partially mapped, and fully unmapped read pairs.
    """
    # Check if input BAM file exists
    if not os.path.exists(args.input_bam):
        logger.error(f"Input BAM file {args.input_bam} not found.")
        sys.exit(1)
    # Check for BAM index
    bam_index = args.input_bam + ".bai"
    if not os.path.exists(bam_index):
        logger.error(f"BAM index file {bam_index} not found. Please index the BAM file first.")
        sys.exit(1)

    try:
        bam_file = pysam.AlignmentFile(args.input_bam, "rb")
    except Exception as e:
        logger.error(f"Error opening BAM file: {e}")
        sys.exit(1)

    logger.info("Reading BAM file and gathering alignments...")
    chr2max = defaultdict(int)
    unique_read_records = []   # Each record: (chrom, blocks, is_read1, is_reverse)
    multi_map_dict = defaultdict(list)  # key: read_name, value: list of (chrom, blocks, align_score, is_read1, is_reverse)
    
    # Dictionary to track paired-end read statuses: key = query_name, value = list of booleans per end
    read_pair_status = {}
    total_alignments = 0
    total_blocks = 0  # Sum of blocks for logging

    for aln in bam_file.fetch(until_eof=True):
        total_alignments += 1
        query_name = aln.query_name

        # Apply paired-end insert size filtering (only if paired and max_insert_size > 0)
        if aln.is_paired and args.max_insert_size > 0:
            if abs(aln.template_length) > args.max_insert_size:
                continue

        # Retrieve continuous aligned blocks (exons) via get_blocks()
        try:
            blocks = aln.get_blocks()  # List of (start, end) tuples
        except Exception as e:
            logger.warning(f"Error parsing blocks for read {query_name}: {e}")
            continue
        if not blocks:
            continue
        total_blocks += len(blocks)

        chrom = bam_file.get_reference_name(aln.reference_id)
        # Update maximum coordinate for the chromosome based on blocks
        for (_, block_end) in blocks:
            if block_end > chr2max[chrom]:
                chr2max[chrom] = block_end

        # Update paired-end status
        if aln.is_paired:
            if query_name not in read_pair_status:
                read_pair_status[query_name] = [False, False]
            if not aln.is_unmapped:
                if aln.is_read1:
                    read_pair_status[query_name][0] = True
                else:
                    read_pair_status[query_name][1] = True
        else:
            read_pair_status[query_name] = [not aln.is_unmapped]

        # Skip unmapped alignments for coverage calculation
        if aln.is_unmapped:
            continue

        # Determine if alignment is unique based on aligner and MAPQ
        if is_unique_alignment(aln, args):
            unique_read_records.append((chrom, blocks, aln.is_read1, aln.is_reverse))
        else:
            if aln.has_tag("AS"):
                align_score = float(aln.get_tag("AS"))
            else:
                align_score = float(aln.mapping_quality)
            multi_map_dict[query_name].append((chrom, blocks, align_score, aln.is_read1, aln.is_reverse))
    
    bam_file.close()

    # Count paired-end mapping statistics
    fully_mapped = 0
    partially_mapped = 0
    fully_unmapped = 0
    for status in read_pair_status.values():
        if len(status) == 1:
            if status[0]:
                fully_mapped += 1
            else:
                fully_unmapped += 1
        elif len(status) == 2:
            if status[0] and status[1]:
                fully_mapped += 1
            elif not (status[0] or status[1]):
                fully_unmapped += 1
            else:
                partially_mapped += 1

    logger.info(f"Total alignments processed: {total_alignments}")
    logger.info(f"Average blocks per alignment: {total_blocks/total_alignments:.2f}")
    logger.info(f"Fully mapped fragments: {fully_mapped}")
    logger.info(f"Partially mapped fragments: {partially_mapped}")
    logger.info(f"Fully unmapped fragments: {fully_unmapped}")
    logger.info(f"Unique mapped fragment count: {len(unique_read_records)}")
    logger.info(f"Multi-mapped fragment count: {len(multi_map_dict)}")

    return (fully_unmapped, partially_mapped, fully_mapped, len(unique_read_records),
            len(multi_map_dict), multi_map_dict, chr2max, unique_read_records)

def is_unique_alignment(aln, args):
    """
    Determine if the alignment is unique.
    For STAR: MAPQ==255 indicates unique; otherwise, MAPQ>=30.
    """
    if args.bam_aligner == 'STAR':
        return (aln.mapping_quality == 255)
    else:
        return (aln.mapping_quality >= 30)

###############################################################################
# 4) Unique Coverage Calculation
###############################################################################
def add_unique_coverage(unique_cov_dict, unique_read_records, logger):
    """
    For each uniquely mapped read, add coverage = 1 for each valid aligned block.
    Blocks with invalid intervals (start >= end) are skipped with a warning.
    """
    for (chrom, blocks, is_read1, is_reverse) in unique_read_records:
        if chrom not in unique_cov_dict:
            continue
        for (block_start, block_end) in blocks:
            if block_start >= block_end:
                logger.warning(f"Invalid block interval ({block_start}, {block_end}) on {chrom}; skipping.")
                continue
            unique_cov_dict[chrom].add_one_coverage(block_start, block_end)

###############################################################################
# 5) EM Algorithm for Multi-mapped Reads
###############################################################################
def update_fractions_for_chunk(read_names, fraction_dict, multi_map_dict, coverage_all):
    """
    For each multi-mapped read in read_names, recalculate fractional assignments.
    For each candidate alignment, weight = (sum over blocks of coverage) * alignment_score.
    Returns updated fraction mapping along with cumulative differences.
    """
    new_fraction_map = {}
    diff_sum = 0.0
    base_sum = 0.0

    for read_name in read_names:
        alignments = multi_map_dict[read_name]
        old_fraction = fraction_dict[read_name]
        weight_list = []
        for (chrom, blocks, align_score, is_read1, is_reverse) in alignments:
            block_weight = 0.0
            for (st, ed) in blocks:
                block_weight += coverage_all[chrom].sum_range(st, ed)
            weight_list.append(block_weight * align_score)
        total_weight = sum(weight_list)
        if total_weight > 0:
            new_fraction = [w / total_weight for w in weight_list]
        else:
            new_fraction = [1.0 / len(weight_list)] * len(weight_list) if weight_list else []
        fraction_diff = sum(abs(nf - of) for nf, of in zip(new_fraction, old_fraction))
        fraction_base = sum(abs(nf) + abs(of) for nf, of in zip(new_fraction, old_fraction))
        diff_sum += fraction_diff
        base_sum += fraction_base
        new_fraction_map[read_name] = new_fraction

    return new_fraction_map, diff_sum, base_sum

def build_multi_coverage_unstranded(fraction_dict, multi_map_dict, coverage_all):
    """
    Build multi-mapped coverage using updated fractional assignments.
    For each candidate alignment of each multi-mapped read, add its fractional contribution over each valid block.
    """
    multi_coverage = {}
    for chrom, cov_obj in coverage_all.items():
        multi_coverage[chrom] = ChunkedCoverage(len(cov_obj))
    for read_name, frac_list in fraction_dict.items():
        alignments = multi_map_dict[read_name]
        for i, (chrom, blocks, align_score, is_read1, is_reverse) in enumerate(alignments):
            fraction_value = frac_list[i]
            if fraction_value <= 0:
                continue
            for (st, ed) in blocks:
                if st >= ed:
                    continue
                multi_coverage[chrom].add_fractional_coverage(st, ed, fraction_value)
    return multi_coverage

def run_em_unstranded(multi_map_dict, unique_cov, fraction_dict, max_iter=100, tol=1e-3, logger=None):
    """
    Run the EM algorithm for multi-mapped reads (unstranded).
    Scheme: total_coverage = unique_cov + multi_coverage_prev.
    Iteratively update fractional assignments and multi-mapped coverage.
    For performance, only multi_coverage_prev is updated and reused.
    Memory for temporary coverage arrays is explicitly released.
    """
    # Initialize multi_coverage_prev (start with zeros) and total coverage array
    multi_coverage_prev = {}
    coverage_all = {}
    for chrom, cov_obj in unique_cov.items():
        multi_coverage_prev[chrom] = ChunkedCoverage(len(cov_obj))
        coverage_all[chrom] = ChunkedCoverage(len(cov_obj))
        # Initialize total coverage as unique_cov (since multi_coverage_prev is all zeros)
        for idx in range(len(cov_obj.chunks)):
            coverage_all[chrom].chunks[idx] = cov_obj.chunks[idx].copy()

    multi_read_names = list(multi_map_dict.keys())
    if not multi_read_names:
        return fraction_dict, coverage_all

    for iteration in range(1, max_iter + 1):
        iter_start_time = time.time()
        # Update fractional assignments for multi-mapped reads using current total coverage
        chunk_size = max(1, len(multi_read_names) // 8)
        read_chunks = [multi_read_names[i:i+chunk_size] for i in range(0, len(multi_read_names), chunk_size)]
        diff_sum_total = 0.0
        base_sum_total = 0.0
        new_fraction_dict = {}
        for chunk in read_chunks:
            updated_frac, diff_sum, base_sum = update_fractions_for_chunk(chunk, fraction_dict, multi_map_dict, coverage_all)
            diff_sum_total += diff_sum
            base_sum_total += base_sum
            new_fraction_dict.update(updated_frac)
        fraction_dict = new_fraction_dict
        fraction_change = diff_sum_total / (base_sum_total + 1e-9)
        logger.info(f"[EM iteration {iteration}] Processed {len(multi_read_names)} reads, fraction change = {fraction_change:.3e}")

        # Reuse multi_coverage_prev: reset its arrays to zero
        for chrom in multi_coverage_prev:
            for arr in multi_coverage_prev[chrom].chunks:
                arr.fill(0)

        # Build new multi coverage based on updated fractions
        multi_coverage_current = build_multi_coverage_unstranded(fraction_dict, multi_map_dict, coverage_all)
        # Copy new coverage into multi_coverage_prev and free memory of multi_coverage_current
        for chrom in multi_coverage_prev:
            for idx in range(len(multi_coverage_prev[chrom].chunks)):
                multi_coverage_prev[chrom].chunks[idx] = multi_coverage_current[chrom].chunks[idx].copy()
        # Free multi_coverage_current: first free internal arrays, then clear dictionary
        for chrom in multi_coverage_current:
            multi_coverage_current[chrom].free_memory()
        multi_coverage_current.clear()

        # Update total coverage in place: unique_cov + multi_coverage_prev
        for chrom in coverage_all:
            for idx in range(len(coverage_all[chrom].chunks)):
                coverage_all[chrom].chunks[idx] = unique_cov[chrom].chunks[idx].copy() + multi_coverage_prev[chrom].chunks[idx]
        iter_end_time = time.time()
        logger.info(f"[EM iteration {iteration}] Completed in {iter_end_time - iter_start_time:.1f}s")
        if fraction_change < tol:
            logger.info(f"[EM] Converged at iteration {iteration}")
            break

    return fraction_dict, coverage_all

###############################################################################
# 6) Parted Coverage for Stranded Data
###############################################################################
def parted_coverage_unique_final(chr2max, unique_read_records):
    """
    Split unique coverage into two strands: F1R2 and F2R1.
      F1R2: (is_read1 != is_reverse)
      F2R1: (is_read1 == is_reverse)
    Process each valid block individually.
    """
    coverage_f1r2_unique = {}
    coverage_f2r1_unique = {}
    for chrom, max_pos in chr2max.items():
        coverage_f1r2_unique[chrom] = ChunkedCoverage(max_pos + 1)
        coverage_f2r1_unique[chrom] = ChunkedCoverage(max_pos + 1)
    for (chrom, blocks, is_read1, is_reverse) in unique_read_records:
        for (st, ed) in blocks:
            if st >= ed:
                continue
            if is_read1 != is_reverse:
                coverage_f1r2_unique[chrom].add_one_coverage(st, ed)
            else:
                coverage_f2r1_unique[chrom].add_one_coverage(st, ed)
    return coverage_f1r2_unique, coverage_f2r1_unique

def parted_coverage_multi_final(fraction_dict, multi_map_dict, chr2max):
    """
    Similar to unique coverage, split multi-mapped coverage into F1R2 and F2R1 based on strand.
    """
    coverage_f1r2_multi = {}
    coverage_f2r1_multi = {}
    for chrom, max_pos in chr2max.items():
        coverage_f1r2_multi[chrom] = ChunkedCoverage(max_pos + 1)
        coverage_f2r1_multi[chrom] = ChunkedCoverage(max_pos + 1)
    for read_name, frac_list in fraction_dict.items():
        alignments = multi_map_dict[read_name]
        for i, (chrom, blocks, align_score, is_read1, is_reverse) in enumerate(alignments):
            fraction_value = frac_list[i]
            if fraction_value <= 0:
                continue
            for (st, ed) in blocks:
                if st >= ed:
                    continue
                if is_read1 != is_reverse:
                    coverage_f1r2_multi[chrom].add_fractional_coverage(st, ed, fraction_value)
                else:
                    coverage_f2r1_multi[chrom].add_fractional_coverage(st, ed, fraction_value)
    return coverage_f1r2_multi, coverage_f2r1_multi

###############################################################################
# 7) Combined BedGraph Writing Function
###############################################################################
def write_bedgraph(coverage_dict, out_path, logger):
    """
    Write a bedGraph file from the given coverage dictionary.
    This function is used for both unstranded and stranded outputs.
    Lines with zero coverage are skipped.
    """
    logger.info(f"Writing bedGraph to {out_path}")
    lines = []
    for chrom in sorted(coverage_dict.keys()):
        cov_array = coverage_dict[chrom].copy_array()
        length_arr = len(cov_array)
        if length_arr == 0:
            continue
        start_idx = 0
        prev_val = cov_array[0]
        for i in range(1, length_arr):
            if cov_array[i] != prev_val:
                if prev_val != 0.0 and i > start_idx:
                    lines.append(f"{chrom}\t{start_idx}\t{i}\t{prev_val:.6f}\n")
                start_idx = i
                prev_val = cov_array[i]
        if prev_val != 0.0 and length_arr > start_idx:
            lines.append(f"{chrom}\t{start_idx}\t{length_arr}\t{prev_val:.6f}\n")
    try:
        with open(out_path, 'w') as out_file:
            out_file.writelines(lines)
    except Exception as e:
        logger.error(f"Error writing bedGraph file {out_path}: {e}")

###############################################################################
# 8) Argument Parsing and Main Function
###############################################################################
def parse_args():
    parser = argparse.ArgumentParser(description="EM-based read coverage estimation from a BAM file.")
    parser.add_argument("--input_bam", required=True, help="Input BAM file with alignments.")
    parser.add_argument("--sample_name", default=None,
                        help="Sample name prefix for output files. Default: use BAM filename.")
    parser.add_argument("--bam_aligner", choices=["STAR", "other"], required=True,
                        help="For STAR: MAPQ=255 indicates unique; otherwise, MAPQ>=30 for unique.")
    parser.add_argument("--stranded", choices=["None", "Forward", "Reverse"], required=True,
                        help="Specify if the data is unstranded or uses forward/reverse orientation.")
    parser.add_argument("--mode", choices=["uniq", "multi"], required=True,
                        help="Mode: 'uniq' for unique coverage only, 'multi' for EM-based multi coverage.")
    parser.add_argument("--max_iter", type=int, default=100,
                        help="Maximum number of EM iterations.")
    parser.add_argument("--tol", type=float, default=1e-3,
                        help="Convergence threshold for EM.")
    parser.add_argument("--max_insert_size", type=int, default=600,
                        help="Maximum allowed insert size for paired-end reads (bp).")
    return parser.parse_args()

def get_output_directory_and_prefix(args):
    """
    Determine the output directory and file prefix based on sample_name or input BAM filename.
    """
    if args.sample_name is None:
        base_name = os.path.basename(args.input_bam)
        if base_name.lower().endswith(".bam"):
            base_name = base_name[:-4]
    else:
        base_name = args.sample_name
    output_dir = base_name + "_out"
    return output_dir, base_name

def main():
    args = parse_args()
    output_dir, prefix = get_output_directory_and_prefix(args)
    os.makedirs(output_dir, exist_ok=True)
    log_path = os.path.join(output_dir, prefix + ".log")
    logger = setup_logging(log_path)
    logger.info("Command: " + " ".join(sys.argv))
    start_time = time.time()

    logger.info(f"Mode: {args.mode}, Stranded: {args.stranded}, Aligner: {args.bam_aligner}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Log file: {log_path}")

    try:
        # 1) Gather alignment information and paired-end mapping statistics
        (fully_unmapped, partially_mapped, fully_mapped, unique_frag_count,
         multi_frag_count, multi_map_dict, chr2max, unique_read_records) = gather_alignments(args, logger)

        # Check if no mapped reads exist
        if not chr2max:
            logger.error("No mapped reads found in the BAM file. Exiting.")
            sys.exit(1)

        # 2) Process Unique Coverage
        if args.mode == "uniq":
            unique_coverage = {chrom: ChunkedCoverage(length_max + 1) for chrom, length_max in chr2max.items()}
            add_unique_coverage(unique_coverage, unique_read_records, logger)
            if args.stranded == "None":
                out_bedgraph = os.path.join(output_dir, f"{prefix}.unstrand_uniq_raw.bedGraph")
                write_bedgraph(unique_coverage, out_bedgraph, logger)
                logger.info(f"Done (uniq, None) => {out_bedgraph}")
            else:
                cov_f1r2_unique, cov_f2r1_unique = parted_coverage_unique_final(chr2max, unique_read_records)
                out_f1r2 = os.path.join(output_dir, f"{prefix}_F1R2_Unique.bedGraph")
                out_f2r1 = os.path.join(output_dir, f"{prefix}_F2R1_Unique.bedGraph")
                write_bedgraph(cov_f1r2_unique, out_f1r2, logger)
                write_bedgraph(cov_f2r1_unique, out_f2r1, logger)
                logger.info(f"Done (uniq, {args.stranded}) => {out_f1r2}, {out_f2r1}")
            elapsed = time.time() - start_time
            logger.info(f"Total run time: {elapsed:.1f}s")
            return

        # 3) Process Multi-mapped Coverage with EM algorithm
        unique_coverage = {chrom: ChunkedCoverage(length_max + 1) for chrom, length_max in chr2max.items()}
        add_unique_coverage(unique_coverage, unique_read_records, logger)
        # Initialize fraction_dict: for each multi-mapped read, assign initial fraction 1/n per candidate alignment.
        fraction_dict = {}
        for read_name, align_list in multi_map_dict.items():
            if len(align_list) > 0:
                fraction_dict[read_name] = [1.0 / len(align_list)] * len(align_list)
            else:
                fraction_dict[read_name] = []
        fraction_dict, coverage_all = run_em_unstranded(multi_map_dict, unique_coverage, fraction_dict,
                                                        max_iter=args.max_iter, tol=args.tol, logger=logger)
        if args.stranded == "None":
            out_bedgraph = os.path.join(output_dir, f"{prefix}.unstrand_multi_raw.bedGraph")
            write_bedgraph(coverage_all, out_bedgraph, logger)
            logger.info(f"Done (multi, None) => {out_bedgraph}")
        else:
            cov_f1r2_unique, cov_f2r1_unique = parted_coverage_unique_final(chr2max, unique_read_records)
            cov_f1r2_multi, cov_f2r1_multi = parted_coverage_multi_final(fraction_dict, multi_map_dict, chr2max)
            out_f1r2_unique = os.path.join(output_dir, f"{prefix}_F1R2_Unique.bedGraph")
            out_f2r1_unique = os.path.join(output_dir, f"{prefix}_F2R1_Unique.bedGraph")
            out_f1r2_multi = os.path.join(output_dir, f"{prefix}_F1R2_Multi.bedGraph")
            out_f2r1_multi = os.path.join(output_dir, f"{prefix}_F2R1_Multi.bedGraph")
            write_bedgraph(cov_f1r2_unique, out_f1r2_unique, logger)
            write_bedgraph(cov_f2r1_unique, out_f2r1_unique, logger)
            write_bedgraph(cov_f1r2_multi, out_f1r2_multi, logger)
            write_bedgraph(cov_f2r1_multi, out_f2r1_multi, logger)
            # Merge unique and multi coverage for total stranded output
            cov_f1r2_total = {}
            cov_f2r1_total = {}
            for chrom in cov_f1r2_unique:
                cov_f1r2_total[chrom] = ChunkedCoverage(len(cov_f1r2_unique[chrom]))
                cov_f2r1_total[chrom] = ChunkedCoverage(len(cov_f2r1_unique[chrom]))
                for i in range(len(cov_f1r2_total[chrom].chunks)):
                    cov_f1r2_total[chrom].chunks[i] = cov_f1r2_unique[chrom].chunks[i] + cov_f1r2_multi[chrom].chunks[i]
                    cov_f2r1_total[chrom].chunks[i] = cov_f2r1_unique[chrom].chunks[i] + cov_f2r1_multi[chrom].chunks[i]
            out_f1r2_total = os.path.join(output_dir, f"{prefix}_F1R2_Total.bedGraph")
            out_f2r1_total = os.path.join(output_dir, f"{prefix}_F2R1_Total.bedGraph")
            write_bedgraph(cov_f1r2_total, out_f1r2_total, logger)
            write_bedgraph(cov_f2r1_total, out_f2r1_total, logger)
            logger.info(f"Done (multi, {args.stranded}) => Stranded bedGraph files generated.")
        elapsed = time.time() - start_time
        logger.info(f"Total run time: {elapsed:.1f}s")
    except Exception as e:
        logger.exception(f"Error in main: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

