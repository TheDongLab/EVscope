#!/usr/bin/env python3
"""
Optimized EM-BaseExpr_v2.py

This script calculates per-base coverage from a BAM file using both uniquely mapped and multi-mapped reads.
It implements an Expectation-Maximization (EM) algorithm for fractional assignment of multi-mapped reads at single-base resolution.

Key modifications in this version:
  1. Cached interval merging via an LRU cache.
  2. Mate cache is now limited using an OrderedDict with a maximum size.
  3. If the BAM file lacks a header, reference lengths are computed dynamically from mapped reads (skipping unmapped), and zero-length chromosomes are filtered.
  4. The ChunkedCoverage constructor now raises an error if total_length <= 0.
  5. In the EM loop, reusable multi-coverage objects are reset in place by zeroing existing arrays and then updating them via copy_chunks_from().
  6. The sum_range method includes stricter range checks and detailed logging.
  7. The stranded parameter is renamed to --split_by_strand with options "yes" or "no".
  
Author: (Your Name)
Date: (Date)
"""

import os
import sys
import time
import logging
import numpy as np
import pysam
from collections import defaultdict, OrderedDict
import argparse
from typing import List, Tuple, Dict, Any
from functools import lru_cache

# Module-level logger for internal debug warnings.
_module_logger = logging.getLogger(__name__)

###############################################################################
# Caching for Interval Merging
###############################################################################
@lru_cache(maxsize=10000)
def cached_merge_intervals(intervals: Tuple[Tuple[int, int], ...]) -> Tuple[Tuple[int, int], ...]:
    """
    Cached version of merge_intervals. Converts the input tuple back to a list,
    performs merging, and returns the result as a tuple.
    """
    return tuple(merge_intervals(list(intervals)))

###############################################################################
# Helper Functions
###############################################################################
def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Merge a list of intervals (each a tuple (start, end)) into non-overlapping intervals.
    
    Example:
      merge_intervals([(100,150), (250,300), (10,70), (100,120), (250,440)])
      returns [(10,70), (100,150), (250,440)]
    
    Parameters:
      intervals: List of (start, end) intervals.
      
    Returns:
      Merged non-overlapping intervals.
      
    Raises:
      ValueError if an interval with start > end is encountered.
    """
    for start, end in intervals:
        if start > end:
            raise ValueError(f"Invalid interval: start ({start}) > end ({end})")
    if not intervals:
        return []
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    merged = []
    current_start, current_end = sorted_intervals[0]
    for start, end in sorted_intervals[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = start, end
    merged.append((current_start, current_end))
    return merged

def merge_pair_blocks(read1_blocks: List[Tuple[int, int]], read2_blocks: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Merge the block intervals from two paired-end reads to form the candidate insert fragment.
    If either block list is empty, log a warning and return an empty list.
    
    Example:
      read1_blocks = [(100,150), (200,250)]
      read2_blocks = [(350,400)]
      returns [(100,150), (200,250), (350,400)]
    
    Parameters:
      read1_blocks: Blocks for the first read.
      read2_blocks: Blocks for the mate.
      
    Returns:
      Merged non-overlapping intervals.
    """
    if not read1_blocks or not read2_blocks:
        _module_logger.warning("Empty blocks detected for read pair. Skipping merge.")
        return []
    all_blocks = read1_blocks + read2_blocks
    return list(cached_merge_intervals(tuple(all_blocks)))

def validate_coverage(coverage_dict: Dict[str, "ChunkedCoverage"], logger: logging.Logger = None) -> None:
    """
    Validate that no negative coverage exists in the coverage dictionary.
    Logs the specific positions of negative coverage and raises a ValueError.
    """
    for chrom, cov_obj in coverage_dict.items():
        arr = cov_obj.copy_array()
        negative_indices = np.where(arr < 0)[0]
        if negative_indices.size > 0:
            if logger:
                logger.error(f"Negative coverage detected in {chrom} at positions: {negative_indices}")
            raise ValueError(f"Negative coverage detected in {chrom} at positions: {negative_indices}")

###############################################################################
# 1) Logging Setup
###############################################################################
def setup_logging(log_path: str) -> logging.Logger:
    """
    Set up logging to file and console at INFO level.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s', "%Y-%m-%d %H:%M:%S")
    
    fh = logging.FileHandler(log_path, mode='w')
    fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    
    return logger

###############################################################################
# 2) ChunkedCoverage Class
###############################################################################
class ChunkedCoverage:
    CHUNK_SIZE = 5_000_000  # Each chunk covers 5 million bases

    def __init__(self, total_length: int):
        """
        Initialize a chunked coverage array for a reference sequence with the given total_length.
        total_length represents the maximum coordinate (region: [0, total_length)).
        Memory is pre-allocated as a NumPy array and cumulative sums are computed.
        Raises ValueError if total_length <= 0.
        """
        if total_length <= 0:
            raise ValueError(f"Invalid total_length={total_length}. Must be > 0.")
        self.total_length = total_length
        self.num_chunks = (total_length + self.CHUNK_SIZE - 1) // self.CHUNK_SIZE
        self.chunks = [np.zeros(min(self.CHUNK_SIZE, total_length - i * self.CHUNK_SIZE), dtype=np.float32)
                       for i in range(self.num_chunks)]
        self.cumsums = [np.cumsum(chunk.copy()) for chunk in self.chunks]

    def __len__(self) -> int:
        return self.total_length

    def _add_value_to_range(self, start: int, end: int, value: float) -> None:
        """
        Add a given value to the coverage over the interval [start, end), splitting across chunks.
        After updating a chunk, its cumulative sum is recomputed.
        """
        if end <= start:
            return
        start = max(start, 0)
        end = min(end, self.total_length)
        while start < end:
            chunk_index = start // self.CHUNK_SIZE
            chunk_offset = start % self.CHUNK_SIZE
            length_to_add = min(self.CHUNK_SIZE - chunk_offset, end - start)
            self.chunks[chunk_index][chunk_offset:chunk_offset+length_to_add] += value
            self.cumsums[chunk_index] = np.cumsum(self.chunks[chunk_index])
            start += length_to_add

    def add_one_coverage(self, start: int, end: int) -> None:
        """
        Add a coverage value of 1 for the interval [start, end).
        """
        self._add_value_to_range(start, end, 1.0)

    def add_fractional_coverage(self, start: int, end: int, fraction: float) -> None:
        """
        Add a fractional coverage value for the interval [start, end).
        """
        if fraction <= 0:
            return
        self._add_value_to_range(start, end, fraction)

    def sum_range(self, start: int, end: int) -> float:
        """
        Compute the sum of coverage over the interval [start, end) using the differences of precomputed
        cumulative sums (cumsums). This optimization greatly improves performance.

        Parameters:
          start: The starting base (inclusive).
          end: The ending base (exclusive).

        Returns:
          The total coverage sum over the interval.
        """
        if self.num_chunks == 0:
            _module_logger.warning(f"Empty coverage object for range [{start}-{end}]. No chunks available.")
            return 0.0
        if start >= self.total_length or end <= 0:
            _module_logger.debug(f"Range [{start}-{end}] is outside the coverage bounds [0-{self.total_length}].")
            return 0.0
        start = max(start, 0)
        end = min(end, self.total_length)
        if end <= start:
            return 0.0
        first_chunk = start // self.CHUNK_SIZE
        last_chunk = (end - 1) // self.CHUNK_SIZE
        if first_chunk >= len(self.chunks) or last_chunk >= len(self.chunks):
            _module_logger.warning(f"Chunk indices out of bounds: chrom_length={self.total_length}, first_chunk={first_chunk}, last_chunk={last_chunk}, num_chunks={len(self.chunks)}")
            return 0.0
        total = 0.0
        if first_chunk == last_chunk:
            cs = self.cumsums[first_chunk]
            start_offset = start % self.CHUNK_SIZE
            end_offset = end % self.CHUNK_SIZE
            if end_offset == 0:
                end_offset = len(self.chunks[first_chunk])
            total = cs[end_offset - 1] - (cs[start_offset - 1] if start_offset > 0 else 0)
        else:
            cs_first = self.cumsums[first_chunk]
            start_offset = start % self.CHUNK_SIZE
            total += cs_first[-1] - (cs_first[start_offset - 1] if start_offset > 0 else 0)
            for idx in range(first_chunk + 1, last_chunk):
                total += self.cumsums[idx][-1]
            cs_last = self.cumsums[last_chunk]
            end_offset = end % self.CHUNK_SIZE
            if end_offset == 0:
                end_offset = len(self.chunks[last_chunk])
            total += cs_last[end_offset - 1]
        return float(total)

    def copy_array(self) -> np.ndarray:
        """
        Flatten and return the coverage array (used for bedGraph output).
        If no chunks exist, return an empty array.
        """
        if not self.chunks:
            return np.array([])
        return np.concatenate(self.chunks)

    def get_full_array(self) -> np.ndarray:
        """
        Return the full one-dimensional coverage array.
        """
        return self.copy_array()

    def copy_chunks(self) -> List[np.ndarray]:
        """
        Return a deep copy of the chunks list.
        """
        return [chunk.copy() for chunk in self.chunks]

    def copy_chunks_from(self, source: "ChunkedCoverage") -> None:
        """
        Copy chunks from a source ChunkedCoverage object to self, and synchronize total_length and num_chunks.
        Update cumulative sums accordingly.
        """
        self.chunks = source.copy_chunks()
        self.total_length = source.total_length  # Synchronize total_length
        self.num_chunks = source.num_chunks      # Synchronize num_chunks
        self.cumsums = [np.cumsum(chunk.copy()) for chunk in self.chunks]

    def free_memory(self) -> None:
        """
        Release memory by explicitly deleting each chunk's array and clearing the chunks and cumsums lists.
        """
        for chunk in self.chunks:
            del chunk
        self.chunks.clear()
        for cs in self.cumsums:
            del cs
        self.cumsums.clear()

###############################################################################
# 3) Gathering Alignments from BAM (Paired and Single-end)
###############################################################################
# Use an OrderedDict for mate_cache to limit its size.
MAX_MATE_CACHE_SIZE = 10000
mate_cache: Dict[str, Any] = OrderedDict()

def gather_alignments(args: argparse.Namespace, logger: logging.Logger) -> Tuple[int, int, int, int, int, Dict[Any, List[Any]], Dict[str, int], List[Any]]:
    """
    Open the input BAM file, classify alignments as unique or multi-mapped, and record reference lengths.
    Preload reference lengths from the BAM header if available; otherwise, compute the maximum coordinate per reference dynamically.
    This ensures that ChunkedCoverage is initialized with a correct length.

    For paired-end reads, process as a pair by:
      - Processing only read1 and retrieving its mate (with up to 3 retries using pysam.mate and caching).
      - If the mate is on a different chromosome, log a warning, update status as partially mapped, and process the reads separately.
      - Otherwise, merge their block intervals via merge_pair_blocks.
      - Record strand information as ("paired", read1_is_reverse, mate_is_reverse).
      - Store the candidate in multi_mapped_reads_dict with key (query_name, True).

    For single-end reads, merge the blocks once and, if the alignment is unique, store ("single", is_read1, is_reverse)
    in unique_fragment_records; otherwise, add to multi_mapped_reads_dict with the same marker.

    Returns:
      fully_unmapped: count of fragments with 0 mapped ends.
      partially_mapped: count of fragments with 1 mapped end.
      fully_mapped: count of fragments with 2 mapped ends (for paired) or 1 for single-end.
      unique_frag_count: number of unique mapped fragments.
      multi_frag_count: number of multi-mapped fragments.
      multi_mapped_reads_dict: dictionary of multi-mapped candidate alignments.
      ref_lengths: dictionary mapping reference (chromosome) to its actual length.
      unique_fragment_records: list of unique mapped records.
         For paired, each record is (chrom, pre-merged intervals, ("paired", read1_is_reverse, mate_is_reverse));
         For single-end, each record is (chrom, pre-merged intervals, ("single", is_read1, is_reverse)).
    """
    try:
        bam_file = pysam.AlignmentFile(args.input_bam, "rb")
    except (IOError, ValueError) as e:
        logger.error(f"Error opening BAM file: {e}")
        sys.exit(1)
    # If header is missing, compute ref_lengths dynamically (skip unmapped reads).
    if not bam_file.references:
        ref_lengths = defaultdict(int)
        for aln in bam_file.fetch(until_eof=True):
            if aln.is_unmapped or aln.reference_id < 0:
                continue
            chrom = getattr(aln, "reference_name", f"chr{aln.reference_id}")
            for (start, end) in aln.get_blocks():
                ref_lengths[chrom] = max(ref_lengths[chrom], end)
        bam_file.reset()
        ref_lengths = {chrom: length for chrom, length in ref_lengths.items() if length > 0}
    else:
        ref_lengths = {ref: length for ref, length in zip(bam_file.references, bam_file.lengths)}
        ref_lengths = {chrom: length for chrom, length in ref_lengths.items() if length > 0}
    if not ref_lengths:
        logger.error("No valid reference lengths found.")
        sys.exit(1)
    unique_fragment_records = []
    multi_mapped_reads_dict = defaultdict(list)
    processed_pairs = set()
    read_pair_status = {}
    total_alignments = 0
    total_blocks = 0
    global mate_cache

    def get_mate(aln: pysam.AlignedSegment, max_retries: int = 3) -> Any:
        """Retrieve mate using caching to reduce repeated I/O. Limits cache size."""
        query_name = aln.query_name
        if query_name in mate_cache:
            mate = mate_cache.pop(query_name)
            mate_cache[query_name] = mate
            return mate
        if len(mate_cache) >= MAX_MATE_CACHE_SIZE:
            mate_cache.popitem(last=False)
        mate_local = None
        for _ in range(max_retries):
            try:
                mate_local = bam_file.mate(aln)
                break
            except (ValueError, StopIteration):
                continue
        mate_cache[query_name] = mate_local
        return mate_local

    for aln in bam_file.fetch(until_eof=True):
        total_alignments += 1
        query_name = aln.query_name
        try:
            blocks = aln.get_blocks()
        except Exception as e:
            logger.error(f"Skipping read {query_name} due to error in get_blocks(): {e}")
            continue
        if not blocks:
            continue
        total_blocks += len(blocks)
        if aln.is_unmapped or aln.reference_id < 0:
            continue
        chrom = bam_file.get_reference_name(aln.reference_id) if bam_file.references else getattr(aln, "reference_name", f"chr{aln.reference_id}")
        ref_length = ref_lengths.get(chrom, 0)
        try:
            merged_blocks = list(cached_merge_intervals(tuple(blocks)))
        except ValueError as e:
            logger.error(f"Skipping read {query_name} due to invalid blocks: {e}")
            continue
        if aln.is_paired:
            if query_name in processed_pairs or not aln.is_read1:
                continue
            mate = get_mate(aln)
            if mate is None:
                logger.debug(f"Processing {query_name} as single-end due to mate retrieval failure.")
                if is_unique_alignment(aln, args):
                    unique_fragment_records.append((chrom, merged_blocks, ("single", aln.is_read1, aln.is_reverse)))
                else:
                    align_score = float(aln.get_tag("AS")) if aln.has_tag("AS") else float(aln.mapping_quality)
                    multi_mapped_reads_dict[(query_name, False)].append((chrom, merged_blocks, align_score, ("single", aln.is_read1, aln.is_reverse)))
                read_pair_status[(query_name, False)] = (1 if not aln.is_unmapped else 0, False)
                continue
            if mate.reference_id != aln.reference_id:
                logger.warning(f"Mate mapped to a different chromosome for {query_name}; treating reads separately.")
                if is_unique_alignment(aln, args):
                    unique_fragment_records.append((chrom, list(cached_merge_intervals(tuple(aln.get_blocks()))), ("single", aln.is_read1, aln.is_reverse)))
                else:
                    align_score = float(aln.get_tag("AS")) if aln.has_tag("AS") else float(aln.mapping_quality)
                    multi_mapped_reads_dict[(query_name, False)].append((chrom, list(cached_merge_intervals(tuple(aln.get_blocks()))), align_score, ("single", aln.is_read1, aln.is_reverse)))
                mate_chrom = bam_file.get_reference_name(mate.reference_id) if bam_file.references else getattr(mate, "reference_name", f"chr{mate.reference_id}")
                if (not mate.is_unmapped) and is_unique_alignment(mate, args):
                    unique_fragment_records.append((mate_chrom, list(cached_merge_intervals(tuple(mate.get_blocks()))), ("single", mate.is_read1, mate.is_reverse)))
                elif mate and not mate.is_unmapped:
                    mate_align_score = float(mate.get_tag("AS")) if mate.has_tag("AS") else float(mate.mapping_quality)
                    multi_mapped_reads_dict[(query_name, True)].append((mate_chrom, list(cached_merge_intervals(tuple(mate.get_blocks()))), mate_align_score, ("single", mate.is_read1, mate.is_reverse)))
                read_pair_status[(query_name, True)] = (1, True)
                continue
            if args.max_insert_size >= 0 and abs(aln.template_length) > args.max_insert_size:
                continue
            merged_blocks = merge_pair_blocks(aln.get_blocks(), mate.get_blocks())
            if not merged_blocks:
                continue
            total_blocks += len(merged_blocks)
            read_pair_status[(query_name, True)] = (2, True)
            strand_info = ("paired", aln.is_reverse, mate.is_reverse)
            if is_unique_alignment(aln, args) and is_unique_alignment(mate, args):
                unique_fragment_records.append((chrom, merged_blocks, strand_info))
            else:
                align_score = ((float(aln.get_tag("AS")) + float(mate.get_tag("AS"))) / 2.0
                               if aln.has_tag("AS") and mate.has_tag("AS")
                               else (float(aln.mapping_quality) + float(mate.mapping_quality)) / 2.0)
                multi_mapped_reads_dict[(query_name, True)].append((chrom, merged_blocks, align_score, strand_info))
            processed_pairs.add(query_name)
        else:
            read_pair_status[(query_name, False)] = (1 if not aln.is_unmapped else 0, False)
            if is_unique_alignment(aln, args):
                unique_fragment_records.append((chrom, merged_blocks, ("single", aln.is_read1, aln.is_reverse)))
            else:
                align_score = float(aln.get_tag("AS")) if aln.has_tag("AS") else float(aln.mapping_quality)
                multi_mapped_reads_dict[(query_name, False)].append((chrom, merged_blocks, align_score, ("single", aln.is_read1, aln.is_reverse)))
    bam_file.close()
    if total_alignments == 0:
        logger.warning("No alignments were processed. Check the BAM file for issues.")
    avg_blocks = total_blocks / total_alignments if total_alignments > 0 else 0
    fully_mapped = sum(1 for v in read_pair_status.values() if v[0] == 2)
    partially_mapped = sum(1 for v in read_pair_status.values() if v[0] == 1)
    fully_unmapped = sum(1 for v in read_pair_status.values() if v[0] == 0)
    logger.info(f"Total alignments processed: {total_alignments}")
    logger.info(f"Average blocks per alignment: {avg_blocks:.2f}")
    logger.info(f"Fully mapped fragments: {fully_mapped}")
    logger.info(f"Partially mapped fragments: {partially_mapped}")
    logger.info(f"Fully unmapped fragments: {fully_unmapped}")
    logger.info(f"Unique mapped fragment count: {len(unique_fragment_records)}")
    logger.info(f"Multi-mapped fragment count: {len(multi_mapped_reads_dict)}")
    mate_cache.clear()  # Clear mate cache to prevent memory leakage.
    return (fully_unmapped, partially_mapped, fully_mapped, len(unique_fragment_records),
            len(multi_mapped_reads_dict), multi_mapped_reads_dict, ref_lengths, unique_fragment_records)

def is_unique_alignment(aln: pysam.AlignedSegment, args: argparse.Namespace) -> bool:
    """
    Determine if an alignment is unique.
    For STAR, MAPQ==255 indicates unique; otherwise, MAPQ>=30 is considered unique.
    """
    if args.bam_aligner == 'STAR':
        return (aln.mapping_quality == 255)
    else:
        return (aln.mapping_quality >= 30)

###############################################################################
# 4) Unique Coverage Calculation
###############################################################################
def add_unique_coverage(unique_coverage_dict: Dict[str, ChunkedCoverage],
                        unique_fragment_records: List[Any],
                        logger: logging.Logger) -> None:
    """
    For each uniquely mapped fragment, add coverage = 1 for each valid interval in its pre-merged intervals.
    """
    for record in unique_fragment_records:
        chrom, merged_blocks, strand_info = record
        if chrom not in unique_coverage_dict:
            logger.warning(f"Chromosome {chrom} not found in reference lengths; skipping intervals.")
            continue
        for start, end in merged_blocks:
            if start < end:
                unique_coverage_dict[chrom].add_one_coverage(start, end)

###############################################################################
# 5) EM Algorithm for Multi-mapped Reads (Paired and Single-end)
###############################################################################
def update_fractions_for_chunk(read_keys: List[Any],
                               fraction_dict: Dict[Any, List[float]],
                               multi_mapped_reads_dict: Dict[Any, List[Any]],
                               unique_coverage_dict: Dict[str, ChunkedCoverage],
                               multi_coverage_prev: Dict[str, ChunkedCoverage]
                              ) -> Tuple[Dict[Any, List[float]], float, float]:
    """
    For each multi-mapped candidate in read_keys, recalculate fractional assignments.
    For each candidate, use the stored pre-merged intervals and compute:
         w_{r,n} = (sum over intervals [unique_coverage(p) + multi_coverage_prev(p)]) * align_score.
    Then update fraction as: f_{r,n} = w_{r,n} / sum_j w_{r,j}.

    Returns:
      new_fraction_map, diff_sum, base_sum.
    """
    new_fraction_map = {}
    diff_sum = 0.0
    base_sum = 0.0
    for key in read_keys:
        alignments = multi_mapped_reads_dict[key]
        old_fraction = fraction_dict[key]
        weight_list = []
        for (chrom, merged_blocks, align_score, strand_info) in alignments:
            interval_sum = sum(unique_coverage_dict[chrom].sum_range(start, end) +
                               multi_coverage_prev[chrom].sum_range(start, end)
                               for start, end in merged_blocks)
            interval_sum = max(interval_sum, 0)
            weight_list.append(interval_sum * align_score)
        if not weight_list:
            new_fraction = []
        else:
            total_weight = sum(weight_list)
            new_fraction = [w / total_weight if total_weight > 0 else 1.0 / len(weight_list) for w in weight_list]
        fraction_diff = sum(abs(nf - of) for nf, of in zip(new_fraction, old_fraction))
        fraction_base = sum(abs(nf) + abs(of) for nf, of in zip(new_fraction, old_fraction))
        diff_sum += fraction_diff
        base_sum += fraction_base
        new_fraction_map[key] = new_fraction
    return new_fraction_map, diff_sum, base_sum

def build_multi_coverage_unstranded(fraction_dict: Dict[Any, List[float]],
                                    multi_mapped_reads_dict: Dict[Any, List[Any]],
                                    unique_coverage_dict: Dict[str, ChunkedCoverage]
                                   ) -> Dict[str, ChunkedCoverage]:
    """
    Build the multi-mapped coverage using current fractional assignments.
    For each candidate, use the stored pre-merged intervals and add its fractional contribution.
    """
    multi_coverage = {}
    for chrom, cov_obj in unique_coverage_dict.items():
        multi_coverage[chrom] = ChunkedCoverage(cov_obj.total_length)
    for key, frac_list in fraction_dict.items():
        alignments = multi_mapped_reads_dict[key]
        for i, (chrom, merged_blocks, align_score, strand_info) in enumerate(alignments):
            if chrom not in multi_coverage:
                _module_logger.warning(f"Chromosome {chrom} not in unique_coverage_dict; skipping candidate.")
                continue
            fraction_value = frac_list[i]
            if fraction_value <= 0:
                continue
            for start, end in merged_blocks:
                if start < end:
                    multi_coverage[chrom].add_fractional_coverage(start, end, fraction_value)
    return multi_coverage

def build_multi_coverage_stranded(fraction_dict: Dict[Any, List[float]],
                                  multi_mapped_reads_dict: Dict[Any, List[Any]],
                                  ref_lengths: Dict[str, int]
                                 ) -> Tuple[Dict[str, ChunkedCoverage], Dict[str, ChunkedCoverage]]:
    """
    Build stranded multi-mapped coverage.

    For each candidate, use the stored pre-merged intervals and assign fractional coverage to F1R2 or F2R1 based on strand info.
    For single-end records, the marker ("single", is_read1, is_reverse) is used:
      - If is_read1 is True: forward → F1R2, reverse → F2R1.
      - If is_read1 is False: forward → F2R1, reverse → F1R2.
    For paired records, the marker ("paired", read1_is_reverse, mate_is_reverse) is used:
      - If (not read1_is_reverse) and mate_is_reverse then F1R2; otherwise F2R1.
      
    Returns:
      coverage_f1r2, coverage_f2r1: Two dictionaries for F1R2 and F2R1 coverage.
    """
    coverage_f1r2 = {}
    coverage_f2r1 = {}
    for chrom, length in ref_lengths.items():
        coverage_f1r2[chrom] = ChunkedCoverage(length)
        coverage_f2r1[chrom] = ChunkedCoverage(length)
    for key, frac_list in fraction_dict.items():
        alignments = multi_mapped_reads_dict[key]
        for i, (chrom, merged_blocks, align_score, strand_info) in enumerate(alignments):
            if chrom not in coverage_f1r2 or chrom not in coverage_f2r1:
                _module_logger.warning(f"Chromosome {chrom} missing in ref_lengths; skipping candidate.")
                continue
            fraction_value = frac_list[i]
            if fraction_value <= 0:
                continue
            if isinstance(strand_info, tuple):
                marker = strand_info[0]
                if marker == "paired":
                    _, read1_rev, mate_rev = strand_info
                    is_f1r2 = (not read1_rev) and mate_rev
                elif marker == "single":
                    _, is_read1, is_reverse = strand_info
                    if is_read1:
                        is_f1r2 = not is_reverse
                    else:
                        is_f1r2 = is_reverse
                else:
                    logger.warning(f"Unexpected strand_info marker: {strand_info}")
                    raise ValueError(f"Unexpected strand_info format: {strand_info}")
            else:
                logger.warning(f"Invalid strand_info format: {strand_info}")
                raise ValueError(f"Unexpected strand_info format: {strand_info}")
            for start, end in merged_blocks:
                if start < end:
                    if is_f1r2:
                        coverage_f1r2[chrom].add_fractional_coverage(start, end, fraction_value)
                    else:
                        coverage_f2r1[chrom].add_fractional_coverage(start, end, fraction_value)
    return coverage_f1r2, coverage_f2r1

def run_em_unstranded(multi_mapped_reads_dict: Dict[Any, List[Any]],
                      unique_coverage_dict: Dict[str, ChunkedCoverage],
                      fraction_dict: Dict[Any, List[float]],
                      max_iter: int = 100,
                      tol: float = 1e-3,
                      logger: logging.Logger = None
                     ) -> Tuple[Dict[Any, List[float]], Dict[str, ChunkedCoverage]]:
    """
    Run the EM algorithm for multi-mapped reads (unstranded) at single-base resolution.

    Initialization:
      - Compute the initial multi-mapping coverage from fraction_dict (with uniform 1/n).
      - Set initial total coverage: coverage_all(0) = unique_coverage_dict + multi_coverage(0).
      - Perform one warm-up iteration to update fractions.

    Iteration:
      - E-step: Update fractional assignments based on current total coverage.
      - M-step: Rebuild multi-mapping coverage using updated fractions.
      - Replace multi_coverage_prev with the new multi-coverage using in-place updates.
      - Update coverage_all = unique_coverage_dict + multi_coverage_prev using np.copyto and np.add in-place.
      - Log iteration time and fraction change (diff_ratio) every 5 iterations (and on the first iteration).
      - Reuse a reusable multi_coverage dictionary to avoid new object creation.

    Returns:
      Updated fraction_dict and total coverage (coverage_all).
    """
    multi_coverage_prev = {}
    coverage_all = {}
    for chrom, cov_obj in unique_coverage_dict.items():
        multi_coverage_prev[chrom] = ChunkedCoverage(cov_obj.total_length)
        coverage_all[chrom] = ChunkedCoverage(cov_obj.total_length)

    # Compute initial multi-mapping coverage.
    initial_multi_coverage = build_multi_coverage_unstranded(fraction_dict, multi_mapped_reads_dict, unique_coverage_dict)
    for chrom in unique_coverage_dict:
        multi_coverage_prev[chrom].chunks = initial_multi_coverage[chrom].copy_chunks()
        num_chunks = min(len(coverage_all[chrom].chunks), len(unique_coverage_dict[chrom].chunks),
                         len(multi_coverage_prev[chrom].chunks))
        for idx in range(num_chunks):
            target = coverage_all[chrom].chunks[idx]
            np.copyto(target, unique_coverage_dict[chrom].chunks[idx])
            np.add(target, multi_coverage_prev[chrom].chunks[idx], out=target)
    for chrom in initial_multi_coverage:
        initial_multi_coverage[chrom].free_memory()
    initial_multi_coverage.clear()

    # Create a reusable multi_coverage dictionary.
    reusable_multi_coverage = {chrom: ChunkedCoverage(unique_coverage_dict[chrom].total_length)
                               for chrom in unique_coverage_dict}

    # Warm-up iteration: pre-adjustment.
    multi_read_keys = list(multi_mapped_reads_dict.keys())
    if not multi_read_keys:
        logger.info("No multi-mapped reads found. Skipping EM iterations.")
        return fraction_dict, coverage_all
    else:
        chunk_size = max(1, len(multi_read_keys) // 8)
    read_chunks = [multi_read_keys[i:i+chunk_size] for i in range(0, len(multi_read_keys), chunk_size)]
    diff_sum_total = 0.0
    base_sum_total = 0.0
    new_fraction_dict = {}
    for chunk in read_chunks:
        updated_frac, diff_sum, base_sum = update_fractions_for_chunk(chunk, fraction_dict, multi_mapped_reads_dict,
                                                                        unique_coverage_dict, multi_coverage_prev)
        diff_sum_total += diff_sum
        base_sum_total += base_sum
        new_fraction_dict.update(updated_frac)
    fraction_dict = new_fraction_dict
    warmup_multi_coverage = build_multi_coverage_unstranded(fraction_dict, multi_mapped_reads_dict, unique_coverage_dict)
    for chrom in multi_coverage_prev:
        multi_coverage_prev[chrom].free_memory()
    multi_coverage_prev.clear()
    # Instead of directly updating with warmup objects, perform a deep copy.
    for chrom, warm_obj in warmup_multi_coverage.items():
        new_cov = ChunkedCoverage(warm_obj.total_length)
        new_cov.copy_chunks_from(warm_obj)
        multi_coverage_prev[chrom] = new_cov
    for chrom in warmup_multi_coverage:
        warmup_multi_coverage[chrom].free_memory()
    warmup_multi_coverage.clear()
    for chrom in coverage_all:
        num_chunks = min(len(coverage_all[chrom].chunks), len(unique_coverage_dict[chrom].chunks),
                         len(multi_coverage_prev[chrom].chunks))
        for idx in range(num_chunks):
            target = coverage_all[chrom].chunks[idx]
            np.copyto(target, unique_coverage_dict[chrom].chunks[idx])
            np.add(target, multi_coverage_prev[chrom].chunks[idx], out=target)
    logger.info("[Warm-up] Pre-adjustment iteration completed.")

    # Official EM iterations.
    for iteration in range(1, max_iter + 1):
        iter_start_time = time.time()
        chunk_size = max(1, len(multi_read_keys) // 8)
        read_chunks = [multi_read_keys[i:i+chunk_size] for i in range(0, len(multi_read_keys), chunk_size)]
        diff_sum_total = 0.0
        base_sum_total = 0.0
        new_fraction_dict = {}
        for chunk in read_chunks:
            updated_frac, diff_sum, base_sum = update_fractions_for_chunk(chunk, fraction_dict, multi_mapped_reads_dict,
                                                                            unique_coverage_dict, multi_coverage_prev)
            diff_sum_total += diff_sum
            base_sum_total += base_sum
            new_fraction_dict.update(updated_frac)
        fraction_dict = new_fraction_dict
        fraction_change = diff_sum_total / (base_sum_total + 1e-9)
        iter_end_time = time.time()
        iteration_time = iter_end_time - iter_start_time
        if iteration % 5 == 0 or iteration == 1:
            logger.info(f"[EM iteration {iteration}] Time: {iteration_time:.1f}s, fraction change: {fraction_change:.3e}")
        else:
            logger.debug(f"[EM iteration {iteration}] Time: {iteration_time:.1f}s, fraction change: {fraction_change:.3e}")

        multi_coverage_current = build_multi_coverage_unstranded(fraction_dict, multi_mapped_reads_dict, unique_coverage_dict)
        for chrom in reusable_multi_coverage:
            # Reset the reusable object: clear existing arrays (fill with zeros) and update by copying new chunks.
            for arr in reusable_multi_coverage[chrom].chunks:
                arr.fill(0)
            reusable_multi_coverage[chrom].copy_chunks_from(multi_coverage_current[chrom])
        for chrom in multi_coverage_prev:
            multi_coverage_prev[chrom].free_memory()
        multi_coverage_prev.clear()
        multi_coverage_prev.update({chrom: reusable_multi_coverage[chrom] for chrom in reusable_multi_coverage})
        for chrom in multi_coverage_current:
            multi_coverage_current[chrom].free_memory()
        multi_coverage_current.clear()
        for chrom in coverage_all:
            num_chunks = min(len(coverage_all[chrom].chunks), len(unique_coverage_dict[chrom].chunks),
                             len(multi_coverage_prev[chrom].chunks))
            for idx in range(num_chunks):
                target = coverage_all[chrom].chunks[idx]
                np.copyto(target, unique_coverage_dict[chrom].chunks[idx])
                np.add(target, multi_coverage_prev[chrom].chunks[idx], out=target)
        if fraction_change < tol:
            logger.info(f"[EM] Converged at iteration {iteration}")
            break

    for cov in reusable_multi_coverage.values():
        cov.free_memory()
    return fraction_dict, coverage_all

###############################################################################
# 6) Parted Coverage for Stranded Data
###############################################################################
def parted_coverage_unique_final(ref_lengths: Dict[str, int], unique_fragment_records: List[Any]) -> Tuple[Dict[str, ChunkedCoverage], Dict[str, ChunkedCoverage]]:
    """
    Split unique coverage into two strands: F1R2 and F2R1.
    
    For paired fragments, the stored marker is ("paired", read1_is_reverse, mate_is_reverse).
    For single-end fragments, the marker is ("single", is_read1, is_reverse) and the following logic applies:
      - If is_read1 is True: if not is_reverse, assign to F1R2; if is_reverse, assign to F2R1.
      - If is_read1 is False: if not is_reverse, assign to F2R1; if is_reverse, assign to F1R2.
    
    Returns:
      Two dictionaries for F1R2 and F2R1 coverage.
    """
    coverage_f1r2_unique = {}
    coverage_f2r1_unique = {}
    for chrom, length in ref_lengths.items():
        coverage_f1r2_unique[chrom] = ChunkedCoverage(length)
        coverage_f2r1_unique[chrom] = ChunkedCoverage(length)
    for record in unique_fragment_records:
        chrom, merged_blocks, strand_info = record
        if isinstance(strand_info, tuple) and strand_info[0] == "paired":
            _, read1_rev, mate_rev = strand_info
            if (not read1_rev) and mate_rev:
                for start, end in merged_blocks:
                    if start < end:
                        coverage_f1r2_unique[chrom].add_one_coverage(start, end)
            else:
                for start, end in merged_blocks:
                    if start < end:
                        coverage_f2r1_unique[chrom].add_one_coverage(start, end)
        elif isinstance(strand_info, tuple) and strand_info[0] == "single":
            _, is_read1, is_reverse = strand_info
            if is_read1:
                # For read1: forward -> F1R2, reverse -> F2R1.
                if not is_reverse:
                    for start, end in merged_blocks:
                        if start < end:
                            coverage_f1r2_unique[chrom].add_one_coverage(start, end)
                else:
                    for start, end in merged_blocks:
                        if start < end:
                            coverage_f2r1_unique[chrom].add_one_coverage(start, end)
            else:
                # For read2: forward -> F2R1, reverse -> F1R2.
                if not is_reverse:
                    for start, end in merged_blocks:
                        if start < end:
                            coverage_f2r1_unique[chrom].add_one_coverage(start, end)
                else:
                    for start, end in merged_blocks:
                        if start < end:
                            coverage_f1r2_unique[chrom].add_one_coverage(start, end)
        else:
            if not strand_info:
                for start, end in merged_blocks:
                    if start < end:
                        coverage_f1r2_unique[chrom].add_one_coverage(start, end)
            else:
                for start, end in merged_blocks:
                    if start < end:
                        coverage_f2r1_unique[chrom].add_one_coverage(start, end)
    return coverage_f1r2_unique, coverage_f2r1_unique

def parted_coverage_multi_final(fraction_dict: Dict[Any, List[float]],
                                multi_mapped_reads_dict: Dict[Any, List[Any]],
                                ref_lengths: Dict[str, int]
                               ) -> Tuple[Dict[str, ChunkedCoverage], Dict[str, ChunkedCoverage]]:
    """
    Similarly, split multi-mapped coverage into F1R2 and F2R1.
    
    For each candidate, use the stored pre-merged intervals and the recorded strand information.
    
    Returns:
      Two dictionaries for F1R2 and F2R1 coverage.
    """
    coverage_f1r2_multi = {}
    coverage_f2r1_multi = {}
    for chrom, length in ref_lengths.items():
        coverage_f1r2_multi[chrom] = ChunkedCoverage(length)
        coverage_f2r1_multi[chrom] = ChunkedCoverage(length)
    for key, frac_list in fraction_dict.items():
        alignments = multi_mapped_reads_dict[key]
        for i, candidate in enumerate(alignments):
            chrom, merged_blocks, align_score, strand_info = candidate
            fraction_value = frac_list[i]
            if fraction_value <= 0:
                continue
            if isinstance(strand_info, tuple) and strand_info[0] == "paired":
                _, read1_rev, mate_rev = strand_info
                is_f1r2 = (not read1_rev) and mate_rev
            elif isinstance(strand_info, tuple) and strand_info[0] == "single":
                _, is_read1, is_reverse = strand_info
                if is_read1:
                    is_f1r2 = not is_reverse
                else:
                    is_f1r2 = is_reverse
            else:
                _module_logger.warning(f"Unexpected strand_info format: {strand_info}")
                raise ValueError(f"Unexpected strand_info format: {strand_info}")
            for start, end in merged_blocks:
                if start < end:
                    if is_f1r2:
                        coverage_f1r2_multi[chrom].add_fractional_coverage(start, end, fraction_value)
                    else:
                        coverage_f2r1_multi[chrom].add_fractional_coverage(start, end, fraction_value)
    return coverage_f1r2_multi, coverage_f2r1_multi

###############################################################################
# 7) Combined BedGraph Writing Functions
###############################################################################
def generate_bedgraph_lines(coverage_dict: Dict[str, ChunkedCoverage]) -> List[str]:
    """
    Generate bedGraph lines from a coverage dictionary by processing each chunk separately.
    
    Returns:
      List of bedGraph lines with 4 decimal places.
      Each line represents a non-overlapping interval.
    """
    lines = []
    for chrom in sorted(coverage_dict.keys()):
        if not coverage_dict[chrom].chunks:
            _module_logger.warning(f"No coverage data found for chromosome {chrom}.")
            continue
        pos_offset = 0
        for chunk in coverage_dict[chrom].chunks:
            diff = np.diff(chunk, prepend=0)
            change_indices = np.where(diff != 0)[0]
            current_pos = 0
            for idx in change_indices:
                if current_pos < idx:
                    val = chunk[current_pos]
                    if val != 0:
                        lines.append(f"{chrom}\t{pos_offset + current_pos}\t{pos_offset + idx}\t{val:.4f}\n")
                current_pos = idx
            # Process any remaining interval at the end of the chunk.
            if current_pos < len(chunk):
                val = chunk[current_pos]
                if val != 0:
                    end_pos = pos_offset + len(chunk)
                    lines.append(f"{chrom}\t{pos_offset + current_pos}\t{end_pos}\t{val:.4f}\n")
            pos_offset += len(chunk)
    return lines

def write_bedgraph(coverage_dict: Dict[str, ChunkedCoverage], out_path: str, logger: logging.Logger) -> None:
    """
    Write a bedGraph file from the provided coverage dictionary.
    If writing fails, log detailed error information and propagate the exception.
    """
    logger.info(f"Writing bedGraph to {out_path}")
    lines = generate_bedgraph_lines(coverage_dict)
    try:
        with open(out_path, 'w') as out_file:
            out_file.writelines(lines)
    except IOError as e:
        logger.error(f"IOError writing bedGraph file {out_path}: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error writing bedGraph file {out_path}: {e}")
        raise

###############################################################################
# 8) Argument Parsing and Main Function
###############################################################################
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="EM-based read coverage estimation from a BAM file.")
    parser.add_argument("--input_bam", required=True,
                        help="Input BAM file with alignments. (Default: None)")
    parser.add_argument("--sample_name", default=None,
                        help="Sample name prefix for output files. (Default: use BAM filename)")
    parser.add_argument("--bam_aligner", choices=["STAR", "other"], required=True,
                        help="For STAR: MAPQ=255 indicates unique; otherwise, MAPQ>=30 for unique. (Default: None)")
    parser.add_argument("--split_by_strand", choices=["yes", "no"], required=True,
                        help="Whether to split coverage by strand. 'yes' for stranded data, 'no' for unstranded.")
    parser.add_argument("--mode", choices=["uniq", "multi"], required=True,
                        help="Mode: 'uniq' for unique coverage only, 'multi' for EM-based multi coverage. (Default: None)")
    parser.add_argument("--max_iter", type=int, default=100,
                        help="Maximum EM iterations (default: 100).")
    parser.add_argument("--tol", type=float, default=1e-3,
                        help="Tolerance for EM convergence (default: 1e-3).")
    parser.add_argument("--max_insert_size", type=int, default=600,
                        help="Maximum allowed insert size for paired-end reads (bp, must be ≥0; default: 600).")
    args = parser.parse_args()
    if args.max_insert_size < 0:
        parser.error("--max_insert_size must be ≥0.")
    return args

def get_output_directory_and_prefix(args: argparse.Namespace) -> Tuple[str, str]:
    """
    Determine the output directory and file prefix based on sample_name or input BAM filename.
    Ensures that the base name is stripped of any path separators.
    """
    if args.sample_name is None:
        base_name = os.path.basename(args.input_bam)
        if base_name.lower().endswith(".bam"):
            base_name = base_name[:-4]
    else:
        base_name = os.path.basename(args.sample_name)
    output_dir = os.path.join(os.getcwd(), base_name + "_out")
    return output_dir, base_name

def main():
    unique_coverage_dict = None
    coverage_all = None
    cov_f1r2_unique = cov_f2r1_unique = cov_f1r2_multi = cov_f2r1_multi = cov_f1r2_total = cov_f2r1_total = None
    try:
        args = parse_args()
        output_dir, prefix = get_output_directory_and_prefix(args)
        os.makedirs(output_dir, exist_ok=True)
        log_path = os.path.join(output_dir, prefix + ".log")
        logger = setup_logging(log_path)
        logger.info("Command: " + " ".join(sys.argv))
        start_time = time.time()

        logger.info(f"Mode: {args.mode}, Split by strand: {args.split_by_strand}, Aligner: {args.bam_aligner}")
        logger.info(f"Output directory: {output_dir}")
        logger.info(f"Log file: {log_path}")

        (fully_unmapped, partially_mapped, fully_mapped, unique_frag_count,
         multi_frag_count, multi_mapped_reads_dict, ref_lengths, unique_fragment_records) = gather_alignments(args, logger)

        if not ref_lengths:
            logger.error("No mapped reads found in the BAM file. Exiting.")
            sys.exit(1)
        elif sum(ref_lengths.values()) == 0:
            logger.error("Reference lengths sum to zero. Exiting.")
            sys.exit(1)

        if args.mode == "uniq":
            unique_coverage_dict = {chrom: ChunkedCoverage(length) for chrom, length in ref_lengths.items() if length > 0}
            add_unique_coverage(unique_coverage_dict, unique_fragment_records, logger)
            validate_coverage(unique_coverage_dict, logger)
            if args.split_by_strand == "no":
                out_bedgraph = os.path.join(output_dir, f"{prefix}.unstrand_uniq_raw.bedGraph")
                write_bedgraph(unique_coverage_dict, out_bedgraph, logger)
                logger.info(f"Done (uniq, unstranded) => {out_bedgraph}")
            else:  # "yes"
                cov_f1r2_unique, cov_f2r1_unique = parted_coverage_unique_final(ref_lengths, unique_fragment_records)
                out_f1r2 = os.path.join(output_dir, f"{prefix}_F1R2_Unique.bedGraph")
                out_f2r1 = os.path.join(output_dir, f"{prefix}_F2R1_Unique.bedGraph")
                write_bedgraph(cov_f1r2_unique, out_f1r2, logger)
                write_bedgraph(cov_f2r1_unique, out_f2r1, logger)
                logger.info(f"Done (uniq, stranded) => {out_f1r2}, {out_f2r1}")
            elapsed = time.time() - start_time
            logger.info(f"Total run time: {elapsed:.1f}s")
            if args.split_by_strand == "yes":
                for cov in cov_f1r2_unique.values():
                    cov.free_memory()
                for cov in cov_f2r1_unique.values():
                    cov.free_memory()
            return

        # Process Multi-mapped Coverage with the EM algorithm.
        unique_coverage_dict = {chrom: ChunkedCoverage(length) for chrom, length in ref_lengths.items() if length > 0}
        add_unique_coverage(unique_coverage_dict, unique_fragment_records, logger)
        fraction_dict = {}
        for key, align_list in multi_mapped_reads_dict.items():
            if align_list:
                fraction_dict[key] = [1.0 / len(align_list)] * len(align_list)
        fraction_dict, coverage_all = run_em_unstranded(multi_mapped_reads_dict, unique_coverage_dict, fraction_dict,
                                                        max_iter=args.max_iter, tol=args.tol, logger=logger)
        validate_coverage(coverage_all, logger)
        if args.split_by_strand == "no":
            out_bedgraph = os.path.join(output_dir, f"{prefix}.unstrand_multi_raw.bedGraph")
            write_bedgraph(coverage_all, out_bedgraph, logger)
            logger.info(f"Done (multi, unstranded) => {out_bedgraph}")
        else:
            cov_f1r2_unique, cov_f2r1_unique = parted_coverage_unique_final(ref_lengths, unique_fragment_records)
            cov_f1r2_multi, cov_f2r1_multi = parted_coverage_multi_final(fraction_dict, multi_mapped_reads_dict, ref_lengths)
            out_f1r2_unique = os.path.join(output_dir, f"{prefix}_F1R2_Unique.bedGraph")
            out_f2r1_unique = os.path.join(output_dir, f"{prefix}_F2R1_Unique.bedGraph")
            out_f1r2_multi = os.path.join(output_dir, f"{prefix}_F1R2_Multi.bedGraph")
            out_f2r1_multi = os.path.join(output_dir, f"{prefix}_F2R1_Multi.bedGraph")
            write_bedgraph(cov_f1r2_unique, out_f1r2_unique, logger)
            write_bedgraph(cov_f2r1_unique, out_f2r1_unique, logger)
            write_bedgraph(cov_f1r2_multi, out_f1r2_multi, logger)
            write_bedgraph(cov_f2r1_multi, out_f2r1_multi, logger)
            cov_f1r2_total = {}
            cov_f2r1_total = {}
            for chrom in cov_f1r2_unique:
                cov_f1r2_total[chrom] = ChunkedCoverage(cov_f1r2_unique[chrom].total_length)
                cov_f2r1_total[chrom] = ChunkedCoverage(cov_f2r1_unique[chrom].total_length)
                num_chunks = min(len(cov_f1r2_total[chrom].chunks),
                                 len(cov_f1r2_unique[chrom].chunks),
                                 len(cov_f1r2_multi[chrom].chunks))
                for i in range(num_chunks):
                    target = cov_f1r2_total[chrom].chunks[i]
                    np.copyto(target, cov_f1r2_unique[chrom].chunks[i])
                    np.add(target, cov_f1r2_multi[chrom].chunks[i], out=target)
                    target = cov_f2r1_total[chrom].chunks[i]
                    np.copyto(target, cov_f2r1_unique[chrom].chunks[i])
                    np.add(target, cov_f2r1_multi[chrom].chunks[i], out=target)
            out_f1r2_total = os.path.join(output_dir, f"{prefix}_F1R2_Total.bedGraph")
            out_f2r1_total = os.path.join(output_dir, f"{prefix}_F2R1_Total.bedGraph")
            write_bedgraph(cov_f1r2_total, out_f1r2_total, logger)
            write_bedgraph(cov_f2r1_total, out_f2r1_total, logger)
            logger.info(f"Done (multi, stranded) => Stranded bedGraph files generated.")
        elapsed = time.time() - start_time
        logger.info(f"Total run time: {elapsed:.1f}s")
    except Exception as e:
        logger.exception(f"Error in main: {e}")
        sys.exit(1)
    finally:
        # Consolidate cleanup; free all coverage dictionaries if they exist.
        for cov_dict in [unique_coverage_dict, coverage_all, cov_f1r2_total, cov_f2r1_total]:
            if cov_dict:
                for cov in cov_dict.values():
                    cov.free_memory()
        for cov_dict in [cov_f1r2_unique, cov_f2r1_unique, cov_f1r2_multi, cov_f2r1_multi]:
            if cov_dict:
                for cov in cov_dict.values():
                    cov.free_memory()

if __name__ == "__main__":
    main()

