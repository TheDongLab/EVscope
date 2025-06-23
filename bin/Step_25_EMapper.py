#!/usr/bin/env python3
"""
EMapper.py: High-performance, multi-process EM-based read coverage estimation.

This industrial-grade version is memory-safe, portable, and robust. It uses a
sparse-data architecture, a highly optimized and lightweight EM algorithm,
and includes numerous enhancements for stability and error reporting.

Version: 1.0.0
"""

import os
import sys
import time
import logging
import subprocess
import numpy as np
import pysam
import argparse
import psutil
import gc
import pyBigWig
import shutil
from collections import defaultdict
from typing import List, Tuple, Dict, Any
from numba import njit
# For robust error tracking.
import traceback
# Use cross-platform contextlib.redirect_stderr for portability.
from contextlib import redirect_stderr
import multiprocessing as mp

# --- Script version ---
__version__ = "1.2.1"

# --- Global logger ---
global_logger: logging.Logger = None

###############################################################################
# Numba Accelerated Functions
###############################################################################
@njit(boundscheck=False, nogil=True)
def add_fractional_coverage_numba(array: np.ndarray, start: int, end: int, fraction_value: float) -> None:
    """Add fraction_value to a slice of a numpy array in-place."""
    start_idx = max(0, start)
    end_idx = min(len(array), end)
    if end_idx > start_idx:
        array[start_idx:end_idx] += fraction_value

@njit
def fast_em_update_numba(weight_array: np.ndarray) -> np.ndarray:
    """Compute normalized fractional assignments for EM update."""
    n = weight_array.shape[0]
    if n == 0:
        return np.empty(0, dtype=weight_array.dtype)
    total_weight = np.sum(weight_array) + 1e-20 # Epsilon for stability
    if total_weight > 0:
        return weight_array / total_weight
    return np.full(n, 1.0 / n, dtype=weight_array.dtype)

@njit
def generate_bigwig_intervals_numba(array_data: np.ndarray) -> np.ndarray:
    """Generate intervals for BigWig from a numpy array."""
    n = array_data.shape[0]
    if n == 0:
        return np.empty((0, 3), dtype=np.float64)
    diff = np.diff(array_data)
    idx = np.where(diff != 0)[0]
    starts = np.concatenate((np.array([0]), idx + 1))
    ends = np.concatenate((idx + 1, np.array([n])))
    values = array_data[starts]
    nonzero_mask = (values != 0) & (ends > starts)
    if not np.any(nonzero_mask):
        return np.empty((0, 3), dtype=np.float64)
    return np.vstack((starts[nonzero_mask], ends[nonzero_mask], values[nonzero_mask])).T

@njit
def get_sum_from_sparse_numba(starts, ends, values, query_start, query_end):
    """Calculates sum over a query interval from sorted, non-overlapping intervals."""
    total_sum = 0.0
    idx = np.searchsorted(ends, query_start, side='right')
    for i in range(idx, len(starts)):
        interval_start, interval_end, value = starts[i], ends[i], values[i]
        if interval_start >= query_end: break
        overlap_start = max(query_start, interval_start)
        overlap_end = min(query_end, interval_end)
        if overlap_end > overlap_start:
            total_sum += value * (overlap_end - overlap_start)
    return total_sum

@njit
def merge_intervals_numba(intervals: np.ndarray) -> np.ndarray:
    if intervals.shape[0] <= 1: return intervals
    intervals = intervals[intervals[:, 0].argsort()]
    merged_list = []
    if intervals.shape[0] == 0: return np.empty((0, 2), dtype=np.int64)
    current_start, current_end = intervals[0]
    for i in range(1, intervals.shape[0]):
        next_start, next_end = intervals[i]
        if next_start < current_end:
            current_end = max(current_end, next_end)
        else:
            merged_list.append((current_start, current_end))
            current_start, current_end = next_start, next_end
    merged_list.append((current_start, current_end))
    return np.array(merged_list, dtype=np.int64)

###############################################################################
# Helper Functions
###############################################################################
def merge_intervals(interval_list: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not interval_list: return []
    return [tuple(row) for row in merge_intervals_numba(np.array(interval_list, dtype=np.int64))]

def _get_read_blocks(read: pysam.AlignedSegment) -> List[Tuple[int, int]]:
    if read.is_unmapped or not read.cigar: return []
    blocks = []
    pos = read.reference_start
    for op, length in read.cigar:
        if op == 0: blocks.append((pos, pos + length))
        if op in {0, 2, 3}: pos += length
    return blocks

def _is_unique_mapping(read: pysam.AlignedSegment, args: argparse.Namespace, logger: logging.Logger) -> bool:
    """Determines if a read is uniquely mapped, prioritizing NH tag."""
    if read.is_secondary or read.is_supplementary: return False
    if read.has_tag("NH"): return read.get_tag("NH") == 1

    if not hasattr(logger, '_warned_about_nh_tag'):
        logger.warning(f"NH tag not found. Falling back to MAPQ >= {args.uniq_mapping_mapq} for uniqueness. This warning will not be shown again.")
        setattr(logger, '_warned_about_nh_tag', True)

    return read.mapping_quality >= args.uniq_mapping_mapq

def _get_alignment_score(read: pysam.AlignedSegment, logger: logging.Logger) -> float:
    """Gets a quality score for an alignment. Prefers AS tag, falls back to a default."""
    if read.has_tag("AS"):
        return float(read.get_tag("AS"))

    if not hasattr(logger, '_warned_about_as_tag'):
        logger.warning("AS tag not found. Using a default alignment score of 1.0 for such reads. This warning will not be shown again.")
        setattr(logger, '_warned_about_as_tag', True)

    return 1.0

def get_mapping_based_strand(read: pysam.AlignedSegment) -> str:
    """Determines strand (F1R2/F2R1) based purely on mapping orientation."""
    if read.is_read2:
        return "F2R1" if not read.is_reverse else "F1R2"
    else: # R1 or true Single-End
        return "F1R2" if not read.is_reverse else "F2R1"

###############################################################################
# Core Processing Logic
###############################################################################

def _create_fragments_truly_robust(read_group, args, logger):
    """Implements the robust, hierarchical pairing strategy."""
    fragments = []
    r1_list = [r for r in read_group if r.is_read1]
    r2_list = [r for r in read_group if r.is_read2]

    processed_r1_mask = [False] * len(r1_list)
    processed_r2_mask = [False] * len(r2_list)

    if r1_list and r2_list:
        for i, r1 in enumerate(r1_list):
            for j, r2 in enumerate(r2_list):
                if r1.reference_id == r2.reference_id and \
                   r1.next_reference_start == r2.reference_start and \
                   r2.next_reference_start == r1.reference_start:

                    template_len = abs(r1.template_length)
                    if template_len == 0 or template_len > args.max_insert_size:
                        continue

                    blocks1, blocks2 = _get_read_blocks(r1), _get_read_blocks(r2)
                    if not (blocks1 and blocks2): continue

                    merged_blocks = merge_intervals(blocks1 + blocks2)
                    if merged_blocks:
                        fragments.append({
                            'blocks': merged_blocks, 'chrom': r1.reference_name,
                            'score': min(_get_alignment_score(r1, logger), _get_alignment_score(r2, logger)),
                            'strand': get_mapping_based_strand(r1),
                            'is_unique': _is_unique_mapping(r1, args, logger) and _is_unique_mapping(r2, args, logger)
                        })
                        processed_r1_mask[i] = True
                        processed_r2_mask[j] = True
                        break

    for mask, read_list in [(processed_r1_mask, r1_list), (processed_r2_mask, r2_list)]:
        for i, r in enumerate(read_list):
            if not mask[i]:
                blocks = _get_read_blocks(r)
                if blocks:
                    fragments.append({
                        'blocks': blocks, 'chrom': r.reference_name, 'score': _get_alignment_score(r, logger),
                        'strand': get_mapping_based_strand(r), 'is_unique': _is_unique_mapping(r, args, logger)
                    })

    for r in (r for r in read_group if not r.is_paired):
        blocks = _get_read_blocks(r)
        if blocks:
            fragments.append({
                'blocks': blocks, 'chrom': r.reference_name, 'score': _get_alignment_score(r, logger),
                'strand': get_mapping_based_strand(r), 'is_unique': _is_unique_mapping(r, args, logger)
            })
    return fragments

def _process_read_group_worker(read_group, unique_f1r2_sparse, unique_f2r1_sparse, multi_mapping_dict, args, logger):
    if not read_group: return

    primary_reads = [r for r in read_group if not r.is_secondary and not r.is_supplementary]
    if not primary_reads: return

    fragments = _create_fragments_truly_robust(primary_reads, args, logger)
    if not fragments: return

    is_group_unique = len(fragments) == 1 and fragments[0]['is_unique']

    if is_group_unique:
        frag = fragments[0]
        target_dict = unique_f1r2_sparse if frag['strand'] == "F1R2" else unique_f2r1_sparse
        for start, end in frag['blocks']:
            target_dict[frag['chrom']].append((start, end, 1.0))
    else:
        all_mappable_reads = [r for r in read_group if not r.is_supplementary]
        em_fragments = _create_fragments_truly_robust(all_mappable_reads, args, logger)

        qname = read_group[0].query_name
        for frag in em_fragments:
            multi_mapping_dict[qname].append(
                (frag['chrom'], frag['blocks'], frag['score'], frag['strand'])
            )

def process_bam_chunk(args_tuple: Tuple) -> Dict:
    try:
        (bam_path, chunk_start, chunk_end, args, log_level) = args_tuple

        worker_logger = logging.getLogger(f'EMapper_worker_{os.getpid()}')
        if not worker_logger.handlers:
                worker_logger.setLevel(log_level)

        unique_f1r2_sparse, unique_f2r1_sparse, multi_mapping_dict = defaultdict(list), defaultdict(list), defaultdict(list)

        with redirect_stderr(open(os.devnull, 'w')):
            with pysam.AlignmentFile(bam_path, "rb", check_sq=False, require_index=False) as bam_file:
                if chunk_start > 0: bam_file.seek(chunk_start)

                try:
                    first_read = next(bam_file)
                    # After a seek, the first read group might be incomplete or belong to the previous chunk.
                    # We find its name and ensure we don't process this partial group if we are not chunk 0.
                    last_qname_in_prev_chunk = first_read.query_name
                    initial_group = [first_read] if chunk_start == 0 else []
                except StopIteration:
                    return {} # Chunk is empty

                current_query_name, read_group = last_qname_in_prev_chunk, initial_group
                for read in bam_file:
                    if chunk_end != -1 and read.offset >= chunk_end and read.query_name != current_query_name:
                        break

                    if read.is_unmapped or read.is_qcfail or read.is_duplicate: continue

                    query_name = read.query_name
                    if query_name != current_query_name:
                        if read_group:
                            # Skip the first group if it's the boundary-spanning one from the previous chunk
                            if not (chunk_start > 0 and current_query_name == last_qname_in_prev_chunk):
                                _process_read_group_worker(read_group, unique_f1r2_sparse, unique_f2r1_sparse, multi_mapping_dict, args, worker_logger)
                        read_group = []

                    read_group.append(read)
                    current_query_name = query_name

                if read_group and not (chunk_start > 0 and current_query_name == last_qname_in_prev_chunk):
                    _process_read_group_worker(read_group, unique_f1r2_sparse, unique_f2r1_sparse, multi_mapping_dict, args, worker_logger)
        return {
            "unique_f1r2_sparse": dict(unique_f1r2_sparse),
            "unique_f2r1_sparse": dict(unique_f2r1_sparse),
            "multi_map_dict": dict(multi_mapping_dict)
        }
    except Exception as e:
        return {"error": repr(e), "traceback": traceback.format_exc()}

def run_em_algorithm_lightweight(multi_mapping_dict, targeted_sparse_unique, max_iterations, tolerance, reference_lengths, logger):
    logger.info("Starting lightweight EM algorithm...")

    fraction_assignments = {key: np.full(len(aligns), 1.0/len(aligns), dtype=np.float64) for key, aligns in multi_mapping_dict.items()}

    for iteration in range(1, max_iterations + 1):
        start_time = time.time()
        logger.info(f"[EM iter {iteration}] Memory: {psutil.Process().memory_info().rss / 1024**2:.1f} MB")

        multi_coverage_dense = {chrom: np.zeros(length, dtype=np.float32) for chrom, length in reference_lengths.items()}
        for qname, fracs in fraction_assignments.items():
            for i, frac in enumerate(fracs):
                chrom, blocks, _, _ = multi_mapping_dict[qname][i]
                if chrom in multi_coverage_dense:
                    for start, end in blocks:
                        add_fractional_coverage_numba(multi_coverage_dense[chrom], start, end, frac)

        new_fractions = {}
        total_difference, total_base = 0.0, 0.0
        for qname, alignments in multi_mapping_dict.items():
            weights = np.zeros(len(alignments), dtype=np.float64)
            for i, (chrom, blocks, score, _) in enumerate(alignments):
                read_coverage = 0.0
                if chrom in targeted_sparse_unique:
                    s, e, v = targeted_sparse_unique[chrom]
                    for start, end in blocks:
                        read_coverage += get_sum_from_sparse_numba(s, e, v, start, end)

                if chrom in multi_coverage_dense:
                    for start, end in blocks:
                        read_coverage += np.sum(multi_coverage_dense[chrom][start:end])
                weights[i] = read_coverage * score

            new_fracs = fast_em_update_numba(weights)
            total_difference += np.sum(np.abs(new_fracs - fraction_assignments[qname]))
            total_base += np.sum(new_fracs + fraction_assignments[qname])
            new_fractions[qname] = new_fracs

        fraction_assignments = new_fractions
        fraction_change = total_difference / (total_base + 1e-9)
        logger.info(f"[EM iter {iteration}] Time: {time.time() - start_time:.1f}s, Frac change: {fraction_change:.4e}")
        if fraction_change < tolerance:
            logger.info(f"EM converged after {iteration} iterations.")
            break

    return fraction_assignments

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=f"EMapper v{__version__}", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input_bam", required=True, help="Input BAM file (must be name-sorted).")
    parser.add_argument("--sample_name", default=None, help="Sample name prefix for output files.")
    parser.add_argument("--output_dir", default=None, help="Output directory. Defaults to <sample_name>_out.")
    parser.add_argument("--num_threads", type=int, default=max(1, mp.cpu_count() // 2), help="Number of threads.")
    parser.add_argument("--split_by_strand", choices=["yes", "no"], default="yes", help="Generate separate BigWigs for F1R2 and F2R1 strands.")
    parser.add_argument("--mode", choices=["uniq", "multi"], default="multi", help="'uniq' only; 'multi' runs EM.")
    parser.add_argument("--max_iter", type=int, default=100, help="Max iterations for EM.")
    parser.add_argument("--tol", type=float, default=1e-3, help="Convergence tolerance for EM.")
    parser.add_argument("--max_insert_size", type=int, default=2000, help="Max insert size for a proper pair.")
    parser.add_argument("--uniq_mapping_mapq", type=int, default=255, help="MAPQ threshold to use for unique mapping if NH tag is absent.")
    parser.add_argument("--keep_temp_files", action="store_true", help="Keep the intermediate name-sorted BAM file.")
    return parser.parse_args()

def write_bigwig_from_sparse(sparse_data, multi_coverage_dense, output_path, reference_lengths, logger):
    logger.info(f"Writing BigWig to {os.path.basename(output_path)}")
    active_chroms = set(sparse_data.keys())
    if multi_coverage_dense: active_chroms.update(multi_coverage_dense.keys())

    with pyBigWig.open(output_path, "w") as bw:
        bw.addHeader(sorted(list(reference_lengths.items())))
        for chrom in sorted(list(active_chroms)):
            length = reference_lengths.get(chrom)
            if not length: continue

            coverage_array = np.zeros(length, dtype=np.float32)
            if chrom in sparse_data:
                for start, end, value in sparse_data[chrom]: add_fractional_coverage_numba(coverage_array, start, end, value)
            if multi_coverage_dense and chrom in multi_coverage_dense:
                coverage_array += multi_coverage_dense[chrom]

            intervals = generate_bigwig_intervals_numba(coverage_array)
            if intervals.shape[0] > 0:
                starts, ends, values = intervals[:, 0].astype(np.int32), intervals[:, 1].astype(np.int32), intervals[:, 2].astype(np.float32)
                bw.addEntries([chrom] * len(starts), starts.tolist(), ends=ends.tolist(), values=values.tolist())

def find_bam_chunks(bam_path: str, num_chunks: int, logger: logging.Logger) -> List[Tuple[int, int]]:
    """
    Handles `OSError: [Errno 0] Success` by wrapping the file operation. This catches
    the rare error during `close()` on a file handle that previously had a read error.
    """
    file_size = os.path.getsize(bam_path)
    SAFETY_MARGIN = 128 * 1024

    if num_chunks <= 1 or file_size < SAFETY_MARGIN * 2:
        logger.info("Processing file in a single chunk (file is small or num_chunks=1).")
        return [(0, -1)]

    chunk_size = file_size // num_chunks
    seek_points = [i * chunk_size for i in range(1, num_chunks)]

    chunks = []
    start_vo = 0

    try:
        with redirect_stderr(open(os.devnull, 'w')):
            with pysam.AlignmentFile(bam_path, "rb", check_sq=False, require_index=False) as bam_file:
                for point in seek_points:
                    if point >= file_size - SAFETY_MARGIN:
                        logger.info(f"Seek point {point} is too close to EOF ({file_size}). Halting chunk creation.")
                        break

                    bam_file.seek(point)
                    try:
                        read = next(bam_file)
                        last_qname = read.query_name
                        for read in bam_file:
                            if read.query_name != last_qname:
                                end_vo = read.offset
                                chunks.append((start_vo, end_vo))
                                start_vo = end_vo
                                break
                        else:
                            break
                    except (OSError, StopIteration) as e:
                        logger.warning(f"Could not find next read record after seeking to byte {point} (likely at EOF). Details: {e}")
                        break
    except OSError as e:
        if e.errno == 0:
            logger.warning(f"Ignoring non-fatal error during BAM file close: {e}. This is a known issue and does not affect results.")
        else:
            raise

    chunks.append((start_vo, -1))
    logger.info(f"Defined {len(chunks)} chunks for parallel processing using virtual offsets.")
    return chunks

def setup_logging(log_path: str) -> logging.Logger:
    logger = logging.getLogger('EMapper')
    logger.handlers.clear()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('[%(asctime)s] [%(levelname)s] [%(processName)s] %(message)s', "%Y-%m-%d %H:%M:%S")
    file_handler = logging.FileHandler(log_path, mode='w')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    return logger

def check_and_prepare_bam(input_bam, output_dir, prefix, num_threads, logger):
    logger.info("--- BAM Preparation Step ---")
    try:
        with redirect_stderr(open(os.devnull, 'w')):
            with pysam.AlignmentFile(input_bam, "rb", check_sq=False, require_index=False) as bam_file:
                if bam_file.header.get('HD', {}).get('SO') == 'queryname':
                    logger.info("BAM header indicates sorting by queryname. Using original file.")
                    return input_bam
    except (ValueError, IOError) as e:
        logger.warning(f"Could not read BAM header from {os.path.basename(input_bam)}: {e}. Assuming it needs sorting.")

    logger.warning("Input BAM is not confirmed to be sorted by queryname. A new name-sorted file will be created.")
    sorted_bam_path = os.path.join(output_dir, f"{prefix}_namesorted.bam")
    logger.info(f"Sorting BAM {os.path.basename(input_bam)} by queryname to {os.path.basename(sorted_bam_path)}")
    command = ["samtools", "sort", "-n", "-@", str(num_threads), "-o", sorted_bam_path, input_bam]
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
    except FileNotFoundError:
        logger.error("samtools not found. Please ensure samtools is installed and in your PATH.")
        raise
    except subprocess.CalledProcessError as e:
        logger.error(f"Samtools sort command failed with exit code {e.returncode}: {' '.join(command)}\nSTDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}")
        raise

    if not os.path.exists(sorted_bam_path) or os.path.getsize(sorted_bam_path) == 0:
        raise IOError(f"samtools sort failed to create a valid output file: {sorted_bam_path}")

    logger.info("Successfully created name-sorted BAM file.")
    return sorted_bam_path

def main():
    global global_logger
    args = parse_arguments()

    prefix = args.sample_name if args.sample_name else os.path.basename(args.input_bam).rsplit('.', 1)[0]
    output_dir = args.output_dir if args.output_dir else f"{prefix}_out"
    os.makedirs(output_dir, exist_ok=True)
    log_path = os.path.join(output_dir, f"{prefix}_EMapper.log")
    global_logger = setup_logging(log_path)
    logger = global_logger

    bam_to_process = None
    is_temp_bam = False

    try:
        start_time = time.time()
        logger.info(f"EMapper v{__version__}")
        logger.info(f"Parameters: {vars(args)}")

        logger.info("Warming up JIT-compiled functions...")
        add_fractional_coverage_numba(np.array([0.0], dtype=np.float32), 0, 1, 1.0)
        fast_em_update_numba(np.array([1.0], dtype=np.float64))
        generate_bigwig_intervals_numba(np.array([0.0, 1.0], dtype=np.float64))
        merge_intervals_numba(np.array([[0, 2]], dtype=np.int64))
        get_sum_from_sparse_numba(np.array([0], dtype=np.int32), np.array([10], dtype=np.int32), np.array([1.0], dtype=np.float32), 0, 1)
        logger.info("Warm-up complete.")

        if not os.path.isfile(args.input_bam):
            logger.error(f"FATAL: Input BAM file not found: {args.input_bam}")
            sys.exit(1)

        bam_to_process = check_and_prepare_bam(args.input_bam, output_dir, prefix, args.num_threads, logger)
        is_temp_bam = (bam_to_process != args.input_bam)

        with redirect_stderr(open(os.devnull, 'w')):
            with pysam.AlignmentFile(bam_to_process, "rb", require_index=False) as bam:
                reference_lengths = {ref: length for ref, length in zip(bam.references, bam.lengths)}

        chunks = find_bam_chunks(bam_to_process, args.num_threads, logger)
        pool_args = [(bam_to_process, start, end, args, logger.level) for start, end in chunks]

        final_unique_f1r2_sparse, final_unique_f2r1_sparse, final_multi_map_dict = defaultdict(list), defaultdict(list), defaultdict(list)

        num_processes = min(args.num_threads, len(chunks))
        logger.info(f"Starting parallel processing with {num_processes} processes.")
        with mp.Pool(processes=num_processes) as pool:
            try:
                results_iterator = pool.imap_unordered(process_bam_chunk, pool_args)
                for i, res in enumerate(results_iterator):
                    if not res: continue
                    if res.get("error"):
                        logger.error(f"FATAL error in worker process: {res['error']}\n{res['traceback']}")
                        raise RuntimeError(f"Worker process failed: {res['error']}")

                    for chrom, data in res["unique_f1r2_sparse"].items(): final_unique_f1r2_sparse[chrom].extend(data)
                    for chrom, data in res["unique_f2r1_sparse"].items(): final_unique_f2r1_sparse[chrom].extend(data)
                    for qname, data in res["multi_map_dict"].items(): final_multi_map_dict[qname].extend(data)

                    if (i + 1) % 10 == 0 or (i + 1) == len(chunks):
                        logger.info(f"Aggregated {i + 1}/{len(chunks)} chunks...")
            except Exception as e:
                logger.error(f"A critical error occurred during parallel processing: {e}")
                pool.terminate()
                pool.join()
                raise

        logger.info("Finished processing all BAM chunks.")
        gc.collect()

        write_bigwig_from_sparse(final_unique_f1r2_sparse, None, os.path.join(output_dir, f"{prefix}_unique_F1R2.bw"), reference_lengths, logger)
        write_bigwig_from_sparse(final_unique_f2r1_sparse, None, os.path.join(output_dir, f"{prefix}_unique_F2R1.bw"), reference_lengths, logger)

        if args.mode == 'multi' and final_multi_map_dict:
            logger.info(f"Found {len(final_multi_map_dict)} read names with multi-mappings. Preparing for lightweight EM.")

            # BUG FIX: Create a new dictionary for combined unique reads to avoid modifying the original
            # 'final_unique_f1r2_sparse' by creating a deep copy. This prevents data contamination.
            logger.info("Creating a clean, combined unique coverage map for EM...")
            unique_total_sparse = defaultdict(list)
            for chrom, data in final_unique_f1r2_sparse.items():
                unique_total_sparse[chrom].extend(data)
            for chrom, data in final_unique_f2r1_sparse.items():
                unique_total_sparse[chrom].extend(data)

            logger.info("Building targeted unique coverage index for EM...")
            targeted_sparse_unique = {}
            for chrom, length in reference_lengths.items():
                if chrom in unique_total_sparse:
                    temp_array = np.zeros(length, dtype=np.float32)
                    for start, end, val in unique_total_sparse[chrom]:
                        add_fractional_coverage_numba(temp_array, start, end, val)

                    intervals = generate_bigwig_intervals_numba(temp_array)
                    if intervals.shape[0] > 0:
                        targeted_sparse_unique[chrom] = (
                            intervals[:, 0].astype(np.int32), intervals[:, 1].astype(np.int32), intervals[:, 2].astype(np.float32)
                        )
                    del temp_array

            del unique_total_sparse
            gc.collect()

            final_fractions = run_em_algorithm_lightweight(final_multi_map_dict, targeted_sparse_unique, args.max_iter, args.tol, reference_lengths, logger)
            del targeted_sparse_unique
            gc.collect()

            multi_f1r2_dense, multi_f2r1_dense = {}, {}
            active_chroms = {chrom for alignments in final_multi_map_dict.values() for chrom, _, _, _ in alignments}
            for chrom in active_chroms:
                if chrom in reference_lengths:
                    multi_f1r2_dense[chrom] = np.zeros(reference_lengths[chrom], dtype=np.float32)
                    multi_f2r1_dense[chrom] = np.zeros(reference_lengths[chrom], dtype=np.float32)

            for qname, fracs in final_fractions.items():
                for i, frac in enumerate(fracs):
                    if frac > 0:
                        chrom, blocks, _, strand = final_multi_map_dict[qname][i]
                        if chrom in reference_lengths:
                            target_array = multi_f1r2_dense[chrom] if strand == "F1R2" else multi_f2r1_dense[chrom]
                            for start, end in blocks: add_fractional_coverage_numba(target_array, start, end, frac)

            del final_multi_map_dict, final_fractions
            gc.collect()

            write_bigwig_from_sparse({}, multi_f1r2_dense, os.path.join(output_dir, f"{prefix}_multi_F1R2.bw"), reference_lengths, logger)
            write_bigwig_from_sparse({}, multi_f2r1_dense, os.path.join(output_dir, f"{prefix}_multi_F2R1.bw"), reference_lengths, logger)

            if args.split_by_strand == 'yes':
                write_bigwig_from_sparse(final_unique_f1r2_sparse, multi_f1r2_dense, os.path.join(output_dir, f"{prefix}_final_F1R2.bw"), reference_lengths, logger)
                write_bigwig_from_sparse(final_unique_f2r1_sparse, multi_f2r1_dense, os.path.join(output_dir, f"{prefix}_final_F2R1.bw"), reference_lengths, logger)

            unstranded_unique = defaultdict(list)
            for chrom, data in final_unique_f1r2_sparse.items(): unstranded_unique[chrom].extend(data)
            for chrom, data in final_unique_f2r1_sparse.items(): unstranded_unique[chrom].extend(data)
            final_unstranded_multi = {c: multi_f1r2_dense.get(c, 0) + multi_f2r1_dense.get(c, 0) for c in reference_lengths if c in multi_f1r2_dense or c in multi_f2r1_dense}
            write_bigwig_from_sparse(unstranded_unique, final_unstranded_multi, os.path.join(output_dir, f"{prefix}_final_unstranded.bw"), reference_lengths, logger)
        else:
            logger.info("No multi-mapping reads found or mode is 'uniq'. Writing final files from unique coverage only.")
            if args.split_by_strand == 'yes' and os.path.exists(os.path.join(output_dir, f"{prefix}_unique_F1R2.bw")):
                 shutil.copyfile(os.path.join(output_dir, f"{prefix}_unique_F1R2.bw"), os.path.join(output_dir, f"{prefix}_final_F1R2.bw"))
                 shutil.copyfile(os.path.join(output_dir, f"{prefix}_unique_F2R1.bw"), os.path.join(output_dir, f"{prefix}_final_F2R1.bw"))

            unstranded_unique = defaultdict(list)
            for chrom, data in final_unique_f1r2_sparse.items(): unstranded_unique[chrom].extend(data)
            for chrom, data in final_unique_f2r1_sparse.items(): unstranded_unique[chrom].extend(data)
            write_bigwig_from_sparse(unstranded_unique, None, os.path.join(output_dir, f"{prefix}_final_unstranded.bw"), reference_lengths, logger)

        logger.info(f"Total execution time: {time.time() - start_time:.2f} seconds")
        logger.info("EMapper run completed successfully.")

    except Exception as e:
        if global_logger:
            global_logger.error(f"An unhandled error occurred in the main process: {e}", exc_info=True)
        else:
            print(f"An unhandled error occurred: {e}", file=sys.stderr)
        sys.exit(1)
    finally:
        if is_temp_bam and not args.keep_temp_files:
            if bam_to_process and os.path.exists(bam_to_process):
                try:
                    logger.info(f"Removing temporary name-sorted BAM: {os.path.basename(bam_to_process)}")
                    os.remove(bam_to_process)
                    bai_path = bam_to_process + ".bai"
                    if os.path.exists(bai_path):
                        os.remove(bai_path)
                except OSError as e:
                    logger.warning(f"Failed to remove temporary file {bam_to_process}: {e}")

if __name__ == "__main__":
    mp.freeze_support()
    main()

