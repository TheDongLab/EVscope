"""
EMapper.py: EMapper: High-performance EM-based RNA-seq quantification tool.
===============================================================================
This script provides robust, memory-safe quantification of RNA-seq data from
BAM files using a two-round Expectation-Maximization (EM) algorithm:
  Round 1 (Positional EM): Resolves positional ambiguity for multi-mapping reads
                           by weighting positions based on coverage evidence.
  Round 2 (Gene EM): Resolves gene-level ambiguity for reads overlapping multiple
                     genes using abundance-based probability weights.
Key Features:
  - Splice junction detection from CIGAR N operations and GTF annotations
  - Polyadenylation site (PAS) detection with strict strand orientation
RPM Normalization:
  - Denominator = total mapped templates (all read groups in BAM)
  - Provides consistent library-size normalization across samples
Count Classification:
  - Counts_Uniq: Reads with unique genomic position (NH=1 or single alignment)
  - Counts_Multi: Reads with multiple genomic positions (multi-mappers)
  - This classification is based on GENOMIC position, not gene assignment
Version: 1.0.0
===============================================================================
"""
import os
import sys
import time
import logging
import re
import subprocess
import pickle
import gc
import atexit
import numpy as np
import argparse
from collections import defaultdict
from typing import List, Tuple, Dict, Any, Generator, Set, Optional, Union
import multiprocessing as mp
import logging.handlers
from functools import partial
import itertools
__version__ = "1.0.0"
# =============================================================================
# Global module-level variables for optional dependencies
# =============================================================================
pysam = None
pyBigWig = None
njit = None
_dependencies_imported = False
def _import_dependencies():
    """Import heavy dependencies lazily. Safe to call multiple times."""
    global pysam, pyBigWig, njit, _dependencies_imported
    if _dependencies_imported:
        return True
    try:
        import pysam as _pysam
        pysam = _pysam
    except ImportError:
        sys.stderr.write("ERROR: pysam not found. Install with: pip install pysam\n")
        return False
    try:
        from numba import njit as _njit
        njit = _njit
    except ImportError:
        sys.stderr.write("ERROR: numba not found. Install with: pip install numba\n")
        return False
    try:
        import pyBigWig as _pyBigWig
        pyBigWig = _pyBigWig
    except ImportError:
        sys.stderr.write("ERROR: pyBigWig not found. Install with: pip install pyBigWig\n")
        return False
    _dependencies_imported = True
    return True
# =============================================================================
# Constants
# =============================================================================
EM_STABILITY_EPSILON = 1e-9
BAM_SCAN_LIMIT = 50_000_000
POLYADENOSINE_MIN_LENGTH = 8
INTERNAL_PRIMING_WINDOW = 20
INTERNAL_PRIMING_A_FRACTION = 0.8
MAX_PAIRS_PER_GROUP = 1000
DEFAULT_MIN_INTRON_LEN = 20
MAX_GENE_OVERLAPS_PER_FRAGMENT = 1000
GENE_EM_PSEUDOCOUNT = 0.1
FRAGMENT_LENGTH_SAMPLE_SIZE = 100000
# =============================================================================
# Global Worker Variables (Process-level, not thread-level)
# =============================================================================
_worker_gene_structs = None
_worker_gene_map_rev = None
_worker_gene_info = None
_worker_log_queue = None
_worker_splice_junctions = None
_worker_single_exon_genes = None
_worker_fasta_path = None
_worker_fasta_handle = None
_bin_array_lock = None
_fasta_handles_to_close = []
def _cleanup_fasta_handles():
    """Cleanup function registered with atexit to close all FASTA handles."""
    global _fasta_handles_to_close, _worker_fasta_handle
    for handle in _fasta_handles_to_close:
        try:
            if handle is not None:
                handle.close()
        except Exception:
            pass
    _fasta_handles_to_close.clear()
    if _worker_fasta_handle is not None:
        try:
            _worker_fasta_handle.close()
        except Exception:
            pass
        _worker_fasta_handle = None
atexit.register(_cleanup_fasta_handles)
def _get_worker_fasta_handle():
    """Get process-local FASTA handle, creating if needed."""
    global _worker_fasta_handle, _fasta_handles_to_close
    if _worker_fasta_path is None:
        return None
    if _worker_fasta_handle is None:
        try:
            _worker_fasta_handle = pysam.FastaFile(_worker_fasta_path)
            _fasta_handles_to_close.append(_worker_fasta_handle)
        except (IOError, ValueError, OSError) as e:
            if _worker_log_queue:
                try:
                    log_record = logging.LogRecord(
                        'EMapper', logging.WARNING, '', 0,
                        f"Failed to open FASTA: {e}", None, None)
                    _worker_log_queue.put_nowait(log_record)
                except Exception:
                    pass
            _worker_fasta_handle = None
    return _worker_fasta_handle
def _validate_fasta_file(fasta_path: str, logger: logging.Logger) -> bool:
    """Validate FASTA file is readable and properly indexed."""
    try:
        with pysam.FastaFile(fasta_path) as fh:
            refs = fh.references
            if not refs:
                logger.error(f"FASTA has no sequences: {fasta_path}")
                return False
            test_ref = refs[0]
            ref_len = fh.get_reference_length(test_ref)
            fh.fetch(test_ref, 0, min(10, ref_len))
        return True
    except (IOError, ValueError, OSError) as e:
        logger.error(f"FASTA validation failed: {e}")
        return False
def _check_fasta_bam_compatibility(fasta_path: str, bam_references: Set[str],
                                   logger: logging.Logger) -> Tuple[bool, Set[str]]:
    """Check chromosome name compatibility between FASTA and BAM.
    Returns:
        Tuple of (is_compatible, matching_chroms)
        is_compatible: True if at least some chromosomes match
        matching_chroms: Set of chromosome names present in both
    """
    try:
        with pysam.FastaFile(fasta_path) as fh:
            fasta_refs = set(fh.references)
    except Exception as e:
        logger.error(f"Failed to read FASTA references: {e}")
        return False, set()
    matching_chroms = bam_references.intersection(fasta_refs)
    if not matching_chroms:
        logger.error("FATAL: No chromosome names match between FASTA and BAM!")
        logger.error(f"  BAM chromosomes (first 10): {sorted(list(bam_references))[:10]}")
        logger.error(f"  FASTA chromosomes (first 10): {sorted(list(fasta_refs))[:10]}")
        logger.error("  Common naming mismatches: 'chr1' vs '1', 'chrM' vs 'MT'")
        return False, set()
    missing_in_fasta = bam_references - fasta_refs
    if missing_in_fasta:
        logger.warning(f"Some BAM chromosomes not in FASTA: {len(missing_in_fasta)} missing")
        if len(missing_in_fasta) <= 10:
            logger.warning(f"  Missing: {sorted(missing_in_fasta)}")
        else:
            logger.warning(f"  Missing (first 10): {sorted(list(missing_in_fasta))[:10]}")
    logger.info(f"FASTA/BAM compatibility: {len(matching_chroms)}/{len(bam_references)} "
                f"BAM chromosomes found in FASTA")
    return True, matching_chroms
# =============================================================================
# Fragment Length Estimation
# =============================================================================
def estimate_fragment_length_from_bam(bam_path: str, sample_size: int = FRAGMENT_LENGTH_SAMPLE_SIZE) -> Tuple[Optional[float], Optional[float], bool]:
    """Estimate mean fragment length from paired-end BAM file.
    Args:
        bam_path: Path to BAM file
        sample_size: Number of proper pairs to sample
    Returns:
        Tuple of (mean_length, std_length, is_paired_end)
        Returns (None, None, False) for single-end data
    """
    sizes = []
    try:
        with pysam.AlignmentFile(bam_path, "rb", check_sq=False, require_index=False) as bam:
            for read in bam:
                if read.is_proper_pair and read.template_length > 0 and not read.is_secondary and not read.is_supplementary:
                    sizes.append(abs(read.template_length))
                if len(sizes) >= sample_size:
                    break
    except Exception:
        return None, None, False
    if len(sizes) >= 100:
        sizes_arr = np.array(sizes, dtype=np.float64)
        q1, q3 = np.percentile(sizes_arr, [25, 75])
        iqr = q3 - q1
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr
        filtered_sizes = sizes_arr[(sizes_arr >= lower_bound) & (sizes_arr <= upper_bound)]
        if len(filtered_sizes) >= 50:
            return float(np.mean(filtered_sizes)), float(np.std(filtered_sizes)), True
        return float(np.mean(sizes_arr)), float(np.std(sizes_arr)), True
    return None, None, False
# =============================================================================
# Numba Function Definitions (Module-level for consistent compilation)
# =============================================================================
_numba_compiled = False
_compiled_fast_em_update = None
_compiled_get_sum_sparse = None
_compiled_sum_blocks_bins = None
_compiled_find_overlaps = None
def _compile_numba_functions():
    """Compile Numba functions at module level. Must be called after _import_dependencies()."""
    global _numba_compiled
    global _compiled_fast_em_update, _compiled_get_sum_sparse
    global _compiled_sum_blocks_bins, _compiled_find_overlaps
    if _numba_compiled:
        return
    if njit is None:
        raise RuntimeError("Numba not imported. Call _import_dependencies() first.")
    @njit(cache=True)
    def _fast_em_update_impl(weight_array, epsilon):
        """Normalize weight array to probability distribution."""
        n = weight_array.shape[0]
        if n == 0:
            return np.empty(0, dtype=weight_array.dtype)
        weights_positive = np.maximum(weight_array, 0.0)
        total_weight = np.sum(weights_positive)
        if total_weight <= epsilon:
            return np.full(n, 1.0 / n, dtype=weight_array.dtype)
        return weights_positive / total_weight
    @njit(cache=True)
    def _get_sum_sparse_impl(starts, ends, values, query_start, query_end):
        """Query sparse coverage data for total value over a range."""
        total_sum = 0.0
        idx = np.searchsorted(ends, query_start, side='right')
        for i in range(idx, len(starts)):
            interval_start = starts[i]
            interval_end = ends[i]
            value = values[i]
            if interval_start >= query_end:
                break
            overlap_start = max(query_start, interval_start)
            overlap_end = min(query_end, interval_end)
            if overlap_end > overlap_start:
                total_sum += max(0.0, value) * (overlap_end - overlap_start)
        return total_sum
    @njit(cache=True)
    def _sum_blocks_bins_impl(bin_values, bin_size, block_starts, block_ends):
        """Sum binned coverage over blocks."""
        total = 0.0
        for i in range(len(block_starts)):
            s = block_starts[i]
            e = block_ends[i]
            if e <= s:
                continue
            b0 = s // bin_size
            b1 = (e - 1) // bin_size
            if b1 >= len(bin_values):
                b1 = len(bin_values) - 1
            if b0 > b1 or b0 < 0:
                continue
            if b0 == b1:
                total += max(0.0, bin_values[b0]) * (e - s)
            else:
                left_edge = (b0 + 1) * bin_size
                total += max(0.0, bin_values[b0]) * (left_edge - s)
                for b in range(b0 + 1, b1):
                    total += max(0.0, bin_values[b]) * bin_size
                right_edge = b1 * bin_size
                total += max(0.0, bin_values[b1]) * (e - right_edge)
        return total
    @njit(nogil=True, cache=True)
    def _find_overlaps_impl(starts, ends, prefix_max_ends, gene_indices,
                            block_starts, block_ends, out_gene_idx, out_overlap_len):
        """Find gene overlaps for read blocks using augmented interval tree logic."""
        if len(starts) == 0:
            return 0
        n_found = 0
        max_out = len(out_gene_idx)
        capacity_exceeded = False
        for i in range(len(block_starts)):
            b_start = block_starts[i]
            b_end = block_ends[i]
            insertion_idx = np.searchsorted(starts, b_start, side='right')
            for j in range(insertion_idx - 1, -1, -1):
                if prefix_max_ends[j] <= b_start:
                    break
                g_start = starts[j]
                g_end = ends[j]
                if g_end > b_start:
                    overlap_start = max(b_start, g_start)
                    overlap_end = min(b_end, g_end)
                    if overlap_end > overlap_start:
                        gene_idx = gene_indices[j]
                        overlap_len = overlap_end - overlap_start
                        found_existing = False
                        for k in range(n_found):
                            if out_gene_idx[k] == gene_idx:
                                out_overlap_len[k] += overlap_len
                                found_existing = True
                                break
                        if not found_existing:
                            if n_found < max_out:
                                out_gene_idx[n_found] = gene_idx
                                out_overlap_len[n_found] = overlap_len
                                n_found += 1
                            else:
                                capacity_exceeded = True
            for j in range(insertion_idx, len(starts)):
                g_start = starts[j]
                g_end = ends[j]
                if g_start >= b_end:
                    break
                overlap_start = max(b_start, g_start)
                overlap_end = min(b_end, g_end)
                if overlap_end > overlap_start:
                    gene_idx = gene_indices[j]
                    overlap_len = overlap_end - overlap_start
                    found_existing = False
                    for k in range(n_found):
                        if out_gene_idx[k] == gene_idx:
                            out_overlap_len[k] += overlap_len
                            found_existing = True
                            break
                    if not found_existing:
                        if n_found < max_out:
                            out_gene_idx[n_found] = gene_idx
                            out_overlap_len[n_found] = overlap_len
                            n_found += 1
                        else:
                            capacity_exceeded = True
        if capacity_exceeded:
            return -n_found
        return n_found
    _compiled_fast_em_update = _fast_em_update_impl
    _compiled_get_sum_sparse = _get_sum_sparse_impl
    _compiled_sum_blocks_bins = _sum_blocks_bins_impl
    _compiled_find_overlaps = _find_overlaps_impl
    _numba_compiled = True
def _ensure_numba_compiled():
    """Ensure Numba functions are compiled. Safe to call multiple times."""
    if not _numba_compiled:
        _compile_numba_functions()
def fast_em_update_numba(weight_array, epsilon):
    """Normalize weights to probability distribution."""
    return _compiled_fast_em_update(weight_array, epsilon)
def get_sum_from_sparse_numba(starts, ends, values, query_start, query_end):
    """Query sparse coverage data."""
    return _compiled_get_sum_sparse(starts, ends, values, query_start, query_end)
def _sum_blocks_on_bins(bin_values, bin_size, block_starts, block_ends):
    """Sum binned coverage over blocks."""
    return _compiled_sum_blocks_bins(bin_values, bin_size, block_starts, block_ends)
def _find_gene_overlaps_numba(starts, ends, prefix_max_ends, gene_indices,
                               block_starts, block_ends, out_gene_idx, out_overlap_len):
    """Find gene overlaps for read blocks."""
    return _compiled_find_overlaps(starts, ends, prefix_max_ends, gene_indices,
                                    block_starts, block_ends, out_gene_idx, out_overlap_len)
def _warmup_numba_functions():
    """Pre-compile Numba functions with dummy data."""
    fast_em_update_numba(np.array([1.0, 2.0], dtype=np.float64), 1e-9)
    _find_gene_overlaps_numba(
        np.array([10], dtype=np.int32), np.array([20], dtype=np.int32),
        np.array([20], dtype=np.int32), np.array([0], dtype=np.int64),
        np.array([12], dtype=np.int32), np.array([18], dtype=np.int32),
        np.zeros(MAX_GENE_OVERLAPS_PER_FRAGMENT, dtype=np.int64),
        np.zeros(MAX_GENE_OVERLAPS_PER_FRAGMENT, dtype=np.int64))
    _sum_blocks_on_bins(np.zeros(10, dtype=np.float64), 100,
                        np.array([50], dtype=np.int32), np.array([150], dtype=np.int32))
    get_sum_from_sparse_numba(np.array([0], dtype=np.int32), np.array([10], dtype=np.int32),
                              np.array([1.0], dtype=np.float64), 0, 10)
# =============================================================================
# Process-local overlap arrays
# =============================================================================
_process_overlap_gene_idx = None
_process_overlap_lengths = None
def _get_overlap_arrays():
    """Get process-local pre-allocated arrays for overlap detection."""
    global _process_overlap_gene_idx, _process_overlap_lengths
    if _process_overlap_gene_idx is None:
        _process_overlap_gene_idx = np.zeros(MAX_GENE_OVERLAPS_PER_FRAGMENT, dtype=np.int64)
        _process_overlap_lengths = np.zeros(MAX_GENE_OVERLAPS_PER_FRAGMENT, dtype=np.int64)
    return _process_overlap_gene_idx, _process_overlap_lengths
# =============================================================================
# Interval Utilities
# =============================================================================
def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge overlapping intervals into non-overlapping sorted list."""
    if not intervals:
        return []
    intervals.sort(key=lambda x: x[0])
    merged = []
    current_start, current_end = intervals[0]
    for i in range(1, len(intervals)):
        next_start, next_end = intervals[i]
        if next_start <= current_end:
            current_end = max(current_end, next_end)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = next_start, next_end
    merged.append((current_start, current_end))
    return merged
# =============================================================================
# GTF Parsing
# =============================================================================
def _parse_gtf_attributes(attr_string: str) -> Dict[str, str]:
    """Parse GTF attribute string."""
    attributes = {}
    pattern = re.compile(r'(\w+)\s+"([^"]*)"')
    for match in pattern.finditer(attr_string):
        attributes[match.group(1)] = match.group(2)
    if not attributes:
        for field in attr_string.split(';'):
            field = field.strip()
            if not field:
                continue
            match = re.match(r'([a-zA-Z0-9_.-]+)\s+(?:"([^"]*)"|(\S+))', field)
            if match:
                key = match.group(1)
                value = match.group(2) if match.group(2) is not None else match.group(3)
                attributes[key] = value
    return attributes
def read_and_process_gtf(gtf_path: str, bam_references: Set[str], logger: logging.Logger,
                         output_dir: str, prefix: str) -> Tuple[Dict, Dict, Dict, np.ndarray, Set, Dict, Set]:
    """Read GTF and build NumPy-based gene structures for fast lookup."""
    logger.info(f"Reading GTF: {os.path.basename(gtf_path)}")
    logger.info("NOTE: Quantification based on 'exon' features.")
    logger.info("COORDINATE SYSTEM: Converting GTF 1-based to internal 0-based.")
    gene_exon_groups = defaultdict(list)
    gtf_chroms = set()
    splice_junctions_by_chrom = defaultdict(set)
    try:
        with open(gtf_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9 or fields[2] != 'exon':
                    continue
                try:
                    attrs = _parse_gtf_attributes(fields[8])
                    gene_id = attrs.get('gene_id')
                    if not gene_id:
                        continue
                    chrom = fields[0]
                    start = int(fields[3]) - 1
                    end = int(fields[4])
                    strand = fields[6]
                    if start < 0 or end <= start or strand not in ['+', '-']:
                        continue
                    gtf_chroms.add(chrom)
                    gene_exon_groups[(gene_id, chrom, strand)].append((start, end))
                except (ValueError, IndexError):
                    continue
    except Exception as e:
        logger.error(f"Failed to read GTF: {e}", exc_info=True)
        raise
    if not gene_exon_groups:
        logger.error(f"FATAL: No 'exon' features in GTF: {gtf_path}")
        sys.exit(1)
    valid_gtf_chroms = gtf_chroms.intersection(bam_references)
    if not valid_gtf_chroms:
        logger.error("FATAL: No chromosome overlap between GTF and BAM.")
        sys.exit(1)
    logger.info(f"Processing genes on {len(valid_gtf_chroms)} matching chromosomes.")
    gene_loci = defaultdict(set)
    for (gene_id, chrom, strand), _ in gene_exon_groups.items():
        if chrom in valid_gtf_chroms:
            gene_loci[gene_id].add((chrom, strand))
    conflicting_gene_ids = {gid for gid, loci in gene_loci.items() if len(loci) > 1}
    if conflicting_gene_ids:
        logger.warning(f"DISCARDING {len(conflicting_gene_ids)} genes on multiple loci.")
        conflict_file = os.path.join(output_dir, f"{prefix}_gtf_conflicts.txt")
        try:
            with open(conflict_file, 'w') as f:
                f.write("# Genes discarded due to multiple chromosomes/strands.\n")
                for gid in sorted(conflicting_gene_ids):
                    loci_str = "; ".join(f"{c}:{s}" for c, s in sorted(gene_loci[gid]))
                    f.write(f"{gid}\t{loci_str}\n")
        except IOError:
            pass
    gene_data_by_chrom = defaultdict(list)
    gene_lengths = defaultdict(int)
    gene_map = {}
    gene_idx_counter = 0
    gene_properties = {}
    gene_exon_counts = defaultdict(int)
    for (gene_id, chrom, strand), exons in sorted(gene_exon_groups.items()):
        if chrom not in valid_gtf_chroms or gene_id in conflicting_gene_ids:
            continue
        if gene_id not in gene_map:
            gene_map[gene_id] = gene_idx_counter
            gene_idx_counter += 1
            gene_properties[gene_id] = (chrom, strand)
        gene_idx = gene_map[gene_id]
        merged_exons = merge_intervals(exons)
        gene_exon_counts[gene_id] = len(merged_exons)
        for i in range(len(merged_exons) - 1):
            donor_pos = merged_exons[i][1]
            acceptor_pos = merged_exons[i + 1][0]
            if acceptor_pos > donor_pos:
                splice_junctions_by_chrom[chrom].add((donor_pos, acceptor_pos))
        for start, end in merged_exons:
            gene_data_by_chrom[chrom].append((start, end, gene_idx))
            gene_lengths[gene_id] += (end - start)
    single_exon_genes = {gid for gid, count in gene_exon_counts.items() if count == 1}
    logger.info(f"Identified {len(single_exon_genes)} single-exon genes.")
    lightweight_gene_structs = {}
    for chrom, data in gene_data_by_chrom.items():
        if not data:
            continue
        n = len(data)
        starts = np.array([d[0] for d in data], dtype=np.int32)
        ends = np.array([d[1] for d in data], dtype=np.int32)
        gene_indices = np.array([d[2] for d in data], dtype=np.int64)
        sort_idx = np.argsort(starts)
        starts = starts[sort_idx]
        ends = ends[sort_idx]
        gene_indices = gene_indices[sort_idx]
        prefix_max_ends = np.maximum.accumulate(ends)
        lightweight_gene_structs[chrom] = {
            'starts': starts, 'ends': ends,
            'prefix_max_ends': prefix_max_ends, 'gene_indices': gene_indices}
    inverse_gene_map = {v: k for k, v in gene_map.items()}
    gene_info_array = np.zeros(len(gene_map), dtype=np.int8)
    for gene_id, (chrom, strand) in gene_properties.items():
        if gene_id in gene_map:
            gene_info_array[gene_map[gene_id]] = 1 if strand == '+' else -1
    valid_gene_ids = {gid for gid, length in gene_lengths.items() if length > 0 and gid in gene_map}
    logger.info(f"Built structures for {len(valid_gene_ids)} valid genes.")
    logger.info(f"Extracted {sum(len(v) for v in splice_junctions_by_chrom.values())} splice junctions.")
    return (lightweight_gene_structs, gene_lengths, inverse_gene_map, gene_info_array,
            valid_gene_ids, splice_junctions_by_chrom, single_exon_genes)
# =============================================================================
# LightweightRead Class
# =============================================================================
class LightweightRead:
    """Pickle-safe lightweight container for AlignedSegment attributes."""
    __slots__ = (
        'query_name', 'is_read1', 'is_read2', 'is_paired', 'is_reverse',
        'is_secondary', 'is_supplementary', 'reference_id', 'reference_name',
        'reference_start', 'reference_end', 'mapping_quality', 'cigar',
        'mate_is_unmapped', 'template_length', 'query_sequence', 'query_qualities_tuple',
        'tag_nh', 'tag_as')
    def __init__(self, read=None, keep_sequence: bool = True):
        """Initialize from pysam AlignedSegment or with defaults."""
        if read is None:
            for name in self.__slots__:
                setattr(self, name, None)
            self.cigar = ()
            return
        self.query_name = read.query_name
        self.is_read1 = read.is_read1
        self.is_read2 = read.is_read2
        self.is_paired = read.is_paired
        self.is_reverse = read.is_reverse
        self.is_secondary = read.is_secondary
        self.is_supplementary = read.is_supplementary
        self.reference_id = read.reference_id
        self.reference_name = read.reference_name
        self.reference_start = read.reference_start
        self.reference_end = read.reference_end
        self.mapping_quality = read.mapping_quality
        cigar_tuples = read.cigartuples
        self.cigar = tuple(cigar_tuples) if cigar_tuples else ()
        self.template_length = read.template_length
        self.mate_is_unmapped = read.mate_is_unmapped
        self.query_sequence = read.query_sequence if keep_sequence else None
        quals = read.query_qualities
        self.query_qualities_tuple = tuple(quals) if (keep_sequence and quals is not None) else None
        try:
            self.tag_nh = read.get_tag("NH")
        except KeyError:
            self.tag_nh = None
        try:
            self.tag_as = read.get_tag("AS")
        except KeyError:
            self.tag_as = None
    def __getstate__(self):
        """Pickle serialization with None safety."""
        return tuple(getattr(self, name, None) for name in self.__slots__)
    def __setstate__(self, state):
        """Pickle deserialization with None safety."""
        for name, value in zip(self.__slots__, state):
            setattr(self, name, value)
        if getattr(self, 'cigar', None) is None:
            self.cigar = ()
    def has_tag(self, tag: str) -> bool:
        """Check if a tag is present (only NH and AS cached)."""
        if tag == "NH":
            return self.tag_nh is not None
        if tag == "AS":
            return self.tag_as is not None
        return False
    def get_tag(self, tag: str) -> Any:
        """Get tag value (only NH and AS cached)."""
        if tag == "NH":
            return self.tag_nh
        if tag == "AS":
            return self.tag_as
        return None
# =============================================================================
# Fragment and Read Utilities
# =============================================================================
def get_fragment_locus_signature(frag: Dict) -> Tuple:
    """Generate unique signature for fragment based on genomic location."""
    blocks_r1 = frag.get('blocks_r1', [])
    blocks_r2 = frag.get('blocks_r2', [])
    return (frag['chrom'], tuple(blocks_r1) if blocks_r1 else (),
            tuple(blocks_r2) if blocks_r2 else (),
            frag.get('is_reverse_r1'), frag.get('is_reverse_r2'))
def deduplicate_fragments_by_locus(fragments: List[Dict]) -> List[Dict]:
    """Deduplicate fragments, keeping highest score per locus."""
    if len(fragments) <= 1:
        return fragments
    seen_loci = {}
    for frag in fragments:
        signature = get_fragment_locus_signature(frag)
        current_score = frag.get('score', 0)
        if signature not in seen_loci:
            seen_loci[signature] = frag
        elif current_score > seen_loci[signature].get('score', 0):
            seen_loci[signature] = frag
    return list(seen_loci.values())
def get_read_alignment_key(r: LightweightRead) -> tuple:
    """Create hashable identifier for a specific alignment."""
    return (r.reference_name, r.reference_start, r.reference_end, r.cigar, r.is_read1)
def _get_read_blocks(read: LightweightRead) -> List[Tuple[int, int]]:
    """Extract aligned blocks from CIGAR."""
    if not read.cigar or read.reference_start is None:
        return []
    blocks = []
    pos = read.reference_start
    current_block_start = None
    for op, length in read.cigar:
        if op in (0, 7, 8):
            if current_block_start is None:
                current_block_start = pos
            pos += length
        elif op == 2:
            pos += length
        elif op == 3:
            if current_block_start is not None:
                blocks.append((current_block_start, pos))
                current_block_start = None
            pos += length
    if current_block_start is not None:
        blocks.append((current_block_start, pos))
    return blocks
def _get_alignment_span(read: LightweightRead) -> int:
    """Get alignment span as fallback when blocks not available."""
    if read.reference_start is not None and read.reference_end is not None:
        return max(1, read.reference_end - read.reference_start)
    return 1
def _is_unique_mapping(read: LightweightRead, effective_mapq: int) -> bool:
    """Check if read is uniquely mapped using NH tag or MAPQ fallback."""
    if read.tag_nh is not None:
        return read.tag_nh == 1
    return read.mapping_quality >= effective_mapq
def _get_alignment_score(read: LightweightRead) -> float:
    """Get alignment score (AS tag), defaulting to 1.0."""
    return float(read.tag_as) if read.tag_as is not None else 1.0
def get_bw_bucket(is_read2: bool, is_reverse: bool) -> str:
    """Determine BigWig bucket for strand-specific coverage assignment.
    Strand assignment logic for paired-end reads:
      - F1R2: Forward R1 + Reverse R2 orientation
      - F2R1: Reverse R1 + Forward R2 orientation
    XOR logic corrected for proper strand assignment.
    """
    return "F2R1" if is_read2 ^ is_reverse else "F1R2"
def _strand_matches_gene(read_is_reverse: Optional[bool], is_read1: bool,
                         gene_strand_val: int, strandness: str) -> bool:
    """Check if single read strand is compatible with gene strand.
    For reverse (dUTP) protocol:
      - R1 should align antisense to gene (R1 strand != gene strand)
      - R2 should align sense to gene (R2 strand == gene strand)
    For forward protocol:
      - R1 should align sense to gene (R1 strand == gene strand)
      - R2 should align antisense to gene (R2 strand != gene strand)
    """
    if strandness == 'unstrand' or read_is_reverse is None:
        return True
    gene_strand = '+' if gene_strand_val == 1 else '-'
    read_strand = '-' if read_is_reverse else '+'
    if strandness == 'reverse':
        return (read_strand != gene_strand) if is_read1 else (read_strand == gene_strand)
    elif strandness == 'forward':
        return (read_strand == gene_strand) if is_read1 else (read_strand != gene_strand)
    return False
def _check_fragment_strand_compatibility(frag: Dict, gene_strand_val: int, strandness: str) -> bool:
    """Check if fragment passes strand compatibility for gene counting.
    Strict paired-end logic:
      - If both R1 and R2 present: BOTH must have correct orientation
      - If only R1 present: R1 must have correct orientation
      - If only R2 present: R2 must have correct orientation
    This filters out:
      - Discordant pairs (one correct, one wrong)
      - Completely wrong pairs (both wrong)
    """
    if strandness == 'unstrand':
        return True
    has_r1 = frag.get('is_reverse_r1') is not None
    has_r2 = frag.get('is_reverse_r2') is not None
    if not has_r1 and not has_r2:
        return False
    match_r1 = _strand_matches_gene(frag.get('is_reverse_r1'), True, gene_strand_val, strandness) if has_r1 else None
    match_r2 = _strand_matches_gene(frag.get('is_reverse_r2'), False, gene_strand_val, strandness) if has_r2 else None
    if has_r1 and has_r2:
        return match_r1 and match_r2
    elif has_r1:
        return match_r1
    elif has_r2:
        return match_r2
    return False
def _read_is_strand_consistent(read: LightweightRead, is_read1: bool, strandness: str) -> bool:
    """Check if read is sense strand for splice/PAS detection."""
    if strandness == 'unstrand':
        return True
    if strandness == 'reverse':
        return not is_read1
    elif strandness == 'forward':
        return is_read1
    return True
def _is_proper_pair(r1: LightweightRead, r2: LightweightRead, max_insert: int, min_insert: int = 10) -> bool:
    """Check for proper paired-end alignment (FR orientation, valid insert size)."""
    if r1.reference_id != r2.reference_id or r1.reference_id == -1:
        return False
    if r1.is_reverse == r2.is_reverse or r1.mate_is_unmapped or r2.mate_is_unmapped:
        return False
    start1, end1 = r1.reference_start, r1.reference_end
    start2, end2 = r2.reference_start, r2.reference_end
    if None in (start1, end1, start2, end2):
        return False
    if r1.is_reverse:
        if start1 < start2:
            return False
    else:
        if start2 < start1:
            return False
    span = max(end1, end2) - min(start1, start2)
    return min_insert <= span <= max_insert
def _get_read_mapped_length(blocks: List[Tuple[int, int]]) -> int:
    """Calculate total mapped length from blocks."""
    return sum(end - start for start, end in blocks)
def _detect_splice_junction_from_cigar(cigar: Tuple, min_intron_len: int) -> bool:
    """Detect splice junction from CIGAR N operations."""
    if not cigar:
        return False
    for op, length in cigar:
        if op == 3 and length >= min_intron_len:
            return True
    return False
def _detect_splice_junction(blocks: List[Tuple[int, int]], known_junctions: Set[Tuple[int, int]]) -> bool:
    """Detect if read blocks span a known splice junction from GTF."""
    if len(blocks) < 2:
        return False
    for i in range(len(blocks) - 1):
        if (blocks[i][1], blocks[i + 1][0]) in known_junctions:
            return True
    return False
def get_effective_frag_len(frag: Dict) -> int:
    """Get effective fragment length for overlap fraction calculation."""
    total = _get_read_mapped_length(frag.get('blocks_r1', [])) + \
            _get_read_mapped_length(frag.get('blocks_r2', []))
    if total <= 0:
        total = max(1, frag.get('r1_span', 0) + frag.get('r2_span', 0))
    return total
# =============================================================================
# Polyadenylation Detection
# =============================================================================
def _get_transcript_3prime_position(read: LightweightRead) -> Optional[int]:
    """Get precise 3' end position of transcript from read alignment."""
    if not read.cigar or read.reference_start is None:
        return None
    ref_pos = read.reference_start
    for op, length in read.cigar:
        if op in (0, 2, 3, 7, 8):
            ref_pos += length
    return read.reference_start if read.is_reverse else ref_pos
def _check_internal_priming(chrom: str, genomic_pos: int, is_reverse: bool, fasta_handle,
                            window: int = INTERNAL_PRIMING_WINDOW,
                            a_fraction: float = INTERNAL_PRIMING_A_FRACTION) -> bool:
    """Check if genomic region is A/T-rich (internal priming artifact)."""
    if fasta_handle is None:
        return False
    try:
        if is_reverse:
            start = max(0, genomic_pos - window)
            end = genomic_pos
        else:
            start = genomic_pos
            end = genomic_pos + window
        seq = fasta_handle.fetch(chrom, start, end).upper()
        if len(seq) == 0:
            return False
        a_frac = seq.count('A') / len(seq)
        t_frac = seq.count('T') / len(seq)
        return max(a_frac, t_frac) >= a_fraction
    except (KeyError, ValueError) as e:
        return False
    except Exception:
        return False
def _find_longest_consecutive_run(seq: str, base: str) -> Tuple[int, int, int]:
    """Find longest consecutive run of base. Returns (length, start, end)."""
    max_run, max_start, max_end = 0, 0, 0
    current_run, current_start = 0, 0
    for i, c in enumerate(seq):
        if c == base:
            if current_run == 0:
                current_start = i
            current_run += 1
            if current_run > max_run:
                max_run = current_run
                max_start = current_start
                max_end = i + 1
        else:
            current_run = 0
    return max_run, max_start, max_end
def _detect_polyadenylation(read: LightweightRead, fasta_handle, strandness: str,
                            min_a_length: int = POLYADENOSINE_MIN_LENGTH) -> bool:
    """Detect poly(A) tail in soft-clipped region with quality filtering."""
    if not read.query_sequence or not read.cigar:
        return False
    is_read1 = read.is_read1 or not read.is_paired
    if not _read_is_strand_consistent(read, is_read1, strandness):
        return False
    seq = read.query_sequence
    quals = read.query_qualities_tuple
    cigar = read.cigar
    if read.is_reverse:
        if cigar and cigar[0][0] == 4:
            clip_len = cigar[0][1]
            clip_seq = seq[:clip_len]
            clip_quals = quals[:clip_len] if quals else None
            if _is_polya_stretch(clip_seq, clip_quals, 'T', min_a_length):
                genomic_pos = _get_transcript_3prime_position(read)
                if genomic_pos is not None and _check_internal_priming(
                        read.reference_name, genomic_pos, True, fasta_handle):
                    return False
                return True
    else:
        if cigar and cigar[-1][0] == 4:
            clip_len = cigar[-1][1]
            clip_seq = seq[-clip_len:]
            clip_quals = quals[-clip_len:] if quals else None
            if _is_polya_stretch(clip_seq, clip_quals, 'A', min_a_length):
                genomic_pos = _get_transcript_3prime_position(read)
                if genomic_pos is not None and _check_internal_priming(
                        read.reference_name, genomic_pos, False, fasta_handle):
                    return False
                return True
    return False
def _is_polya_stretch(seq: str, quals: Optional[Tuple[int, ...]], base: str,
                      min_len: int, min_qual: int = 20) -> bool:
    """Check for consecutive poly(A/T) with quality filtering on the run."""
    if len(seq) < min_len:
        return False
    longest_run, run_start, run_end = _find_longest_consecutive_run(seq, base)
    if longest_run < min_len:
        return False
    if quals is not None and len(quals) >= run_end:
        run_quals = quals[run_start:run_end]
        if sum(1 for q in run_quals if q >= min_qual) < min_len:
            return False
    return True
# =============================================================================
# Bin Array Utilities
# =============================================================================
def _ensure_bin_array(chrom_bins: Dict[str, np.ndarray], chrom: str, chrom_len: int, bin_size: int):
    """Bin array initialization (no lock needed in single-process context per worker)."""
    if chrom not in chrom_bins:
        n_bins = (chrom_len + bin_size - 1) // bin_size
        chrom_bins[chrom] = np.zeros(n_bins, dtype=np.float64)
def _accumulate_blocks_into_bins(chrom_bins: Dict[str, np.ndarray], chrom_len_map: Dict[str, int],
                                 chrom: str, blocks: List[Tuple[int, int]], value: float, bin_size: int):
    """Accumulate coverage into bins (stores density per bp)."""
    if not blocks or chrom not in chrom_len_map:
        return
    chrom_len = chrom_len_map[chrom]
    _ensure_bin_array(chrom_bins, chrom, chrom_len, bin_size)
    arr = chrom_bins[chrom]
    n_bins = len(arr)
    for s, e in blocks:
        if e <= s:
            continue
        s = max(0, s)
        e = min(chrom_len, e)
        if e <= s:
            continue
        b0 = s // bin_size
        b1 = (e - 1) // bin_size
        if b1 >= n_bins:
            b1 = n_bins - 1
        if b0 > b1 or b0 < 0:
            continue
        if b0 == b1:
            arr[b0] += value * (e - s) / bin_size
        else:
            arr[b0] += value * ((b0 + 1) * bin_size - s) / bin_size
            if b1 > b0 + 1:
                arr[b0 + 1:b1] += value
            arr[b1] += value * (e - b1 * bin_size) / bin_size
# =============================================================================
# Gene Overlap Detection
# =============================================================================
def _get_gene_overlaps(blocks: List[Tuple[int, int]], chrom_gene_struct: Dict) -> Tuple[Dict[int, int], bool]:
    """Find gene overlaps using pre-allocated arrays. Returns (overlaps, overflow)."""
    if not blocks:
        return {}, False
    block_starts = np.array([s for s, e in blocks], dtype=np.int32)
    block_ends = np.array([e for s, e in blocks], dtype=np.int32)
    out_gene_idx, out_overlap_len = _get_overlap_arrays()
    out_gene_idx.fill(0)
    out_overlap_len.fill(0)
    n_found = _find_gene_overlaps_numba(
        chrom_gene_struct['starts'], chrom_gene_struct['ends'],
        chrom_gene_struct['prefix_max_ends'], chrom_gene_struct['gene_indices'],
        block_starts, block_ends, out_gene_idx, out_overlap_len)
    capacity_exceeded = n_found < 0
    if capacity_exceeded:
        n_found = -n_found
    return {int(out_gene_idx[i]): int(out_overlap_len[i]) for i in range(n_found)}, capacity_exceeded
def format_count(value: float) -> str:
    """Format count value, removing unnecessary trailing zeros and decimal points."""
    if value == 0:
        return '0'
    if value == int(value):
        return str(int(value))
    formatted = f"{value:.2f}"
    if '.' in formatted:
        formatted = formatted.rstrip('0').rstrip('.')
    return formatted if formatted else '0'
# =============================================================================
# Coverage and Count Accumulation
# =============================================================================
def _accumulate_frag_for_coverage(frag: Dict, args, chrom_len_map: Dict[str, int],
                                  unique_f1r2: Dict, unique_f2r1: Dict,
                                  unique_bins_f1r2: Dict, unique_bins_f2r1: Dict,
                                  weight: float = 1.0):
    """Accumulate fragment coverage into BigWig structures."""
    per_end_weight = weight
    for blocks, is_r2, is_rev_key in [
        (frag.get('blocks_r1'), False, 'is_reverse_r1'),
        (frag.get('blocks_r2'), True, 'is_reverse_r2')
    ]:
        if not blocks:
            continue
        is_rev = frag.get(is_rev_key)
        if is_rev is None:
            continue
        bucket = get_bw_bucket(is_r2, is_rev)
        if args.bin_size == 1:
            target = unique_f1r2 if bucket == "F1R2" else unique_f2r1
            target[frag['chrom']].extend((s, e, per_end_weight) for s, e in blocks)
        else:
            target_bins = unique_bins_f1r2 if bucket == "F1R2" else unique_bins_f2r1
            _accumulate_blocks_into_bins(target_bins, chrom_len_map, frag['chrom'],
                                         blocks, per_end_weight, args.bin_size)
# =============================================================================
# Fragment Creation
# =============================================================================
def _create_fragments_from_read_group(read_group: List[LightweightRead], args) -> List[Dict]:
    """Create fragments from read group with deduplication."""
    fragments = []
    r1_list = [r for r in read_group if r.is_read1]
    r2_list = [r for r in read_group if r.is_read2]
    paired_read_keys = set()
    min_intron_len = getattr(args, 'min_intron_len', DEFAULT_MIN_INTRON_LEN)
    strandness = getattr(args, 'strandness', 'reverse')
    if r1_list and r2_list:
        potential_pairs = []
        r1s = sorted(r1_list, key=_get_alignment_score, reverse=True)[:args.max_alignments_per_mate]
        r2s = sorted(r2_list, key=_get_alignment_score, reverse=True)[:args.max_alignments_per_mate]
        pair_count = 0
        for r1 in r1s:
            for r2 in r2s:
                if pair_count >= MAX_PAIRS_PER_GROUP:
                    break
                if _is_proper_pair(r1, r2, args.max_insert_size, args.min_insert_size):
                    score = _get_alignment_score(r1) + _get_alignment_score(r2)
                    potential_pairs.append((score, r1, r2))
                    pair_count += 1
            if pair_count >= MAX_PAIRS_PER_GROUP:
                break
        if args.positional_assignment == 'EM':
            seen_pairs = set()
            for score, r1, r2 in potential_pairs:
                pair_key = tuple(sorted((get_read_alignment_key(r1), get_read_alignment_key(r2))))
                if pair_key in seen_pairs:
                    continue
                seen_pairs.add(pair_key)
                has_splice = False
                if _read_is_strand_consistent(r1, True, strandness):
                    has_splice |= _detect_splice_junction_from_cigar(r1.cigar, min_intron_len)
                if _read_is_strand_consistent(r2, False, strandness):
                    has_splice |= _detect_splice_junction_from_cigar(r2.cigar, min_intron_len)
                fragments.append({
                    'blocks_r1': _get_read_blocks(r1), 'blocks_r2': _get_read_blocks(r2),
                    'chrom': r1.reference_name, 'score': score, 'is_paired': True,
                    'is_reverse_r1': r1.is_reverse, 'is_reverse_r2': r2.is_reverse,
                    'is_unique': _is_unique_mapping(r1, args.unique_mapping_threshold) and
                                 _is_unique_mapping(r2, args.unique_mapping_threshold),
                    'read1_obj': r1, 'read2_obj': r2, 'has_splice_cigar': has_splice,
                    'r1_span': _get_alignment_span(r1), 'r2_span': _get_alignment_span(r2)})
                paired_read_keys.add(get_read_alignment_key(r1))
                paired_read_keys.add(get_read_alignment_key(r2))
        else:
            potential_pairs.sort(key=lambda x: x[0], reverse=True)
            if potential_pairs:
                score, r1, r2 = potential_pairs[0]
                has_splice = False
                if _read_is_strand_consistent(r1, True, strandness):
                    has_splice |= _detect_splice_junction_from_cigar(r1.cigar, min_intron_len)
                if _read_is_strand_consistent(r2, False, strandness):
                    has_splice |= _detect_splice_junction_from_cigar(r2.cigar, min_intron_len)
                fragments.append({
                    'blocks_r1': _get_read_blocks(r1), 'blocks_r2': _get_read_blocks(r2),
                    'chrom': r1.reference_name, 'score': score, 'is_paired': True,
                    'is_reverse_r1': r1.is_reverse, 'is_reverse_r2': r2.is_reverse,
                    'is_unique': _is_unique_mapping(r1, args.unique_mapping_threshold) and
                                 _is_unique_mapping(r2, args.unique_mapping_threshold),
                    'read1_obj': r1, 'read2_obj': r2, 'has_splice_cigar': has_splice,
                    'r1_span': _get_alignment_span(r1), 'r2_span': _get_alignment_span(r2)})
                paired_read_keys.add(get_read_alignment_key(r1))
                paired_read_keys.add(get_read_alignment_key(r2))
    unpaired_reads = [r for r in read_group if get_read_alignment_key(r) not in paired_read_keys]
    seen_single = set()
    for r in unpaired_reads:
        alignment_key = get_read_alignment_key(r)
        if alignment_key in seen_single:
            continue
        seen_single.add(alignment_key)
        blocks = _get_read_blocks(r)
        if not blocks:
            continue
        is_r1_type = r.is_read1 or not r.is_paired
        is_sense = _read_is_strand_consistent(r, is_r1_type, strandness)
        has_splice = _detect_splice_junction_from_cigar(r.cigar, min_intron_len) if is_sense else False
        fragments.append({
            'blocks_r1': blocks if is_r1_type else [],
            'blocks_r2': blocks if not is_r1_type else [],
            'chrom': r.reference_name, 'score': _get_alignment_score(r), 'is_paired': False,
            'is_reverse_r1': r.is_reverse if is_r1_type else None,
            'is_reverse_r2': r.is_reverse if not is_r1_type else None,
            'is_unique': _is_unique_mapping(r, args.unique_mapping_threshold),
            'read1_obj': r if is_r1_type else None, 'read2_obj': r if not is_r1_type else None,
            'has_splice_cigar': has_splice,
            'r1_span': _get_alignment_span(r) if is_r1_type else 0,
            'r2_span': _get_alignment_span(r) if not is_r1_type else 0})
    return deduplicate_fragments_by_locus(fragments)
# =============================================================================
# Worker Processing Function
# =============================================================================
def _process_read_group_worker(read_group: List[LightweightRead], args_dict: Dict) -> Dict:
    """Worker function to process a read group."""
    args = argparse.Namespace(**args_dict)
    qname = read_group[0].query_name if read_group else "N/A"
    try:
        if not read_group:
            return {"qname": qname, "has_valid_fragments": False}
        rpm_count_alignments = sum(1 for r in read_group if not r.is_secondary)
        if args.positional_assignment == 'unique':
            for r in read_group:
                if not _is_unique_mapping(r, args.unique_mapping_threshold):
                    return {"qname": qname, "rpm_count_alignments": rpm_count_alignments,
                            "skipped_multi": True, "has_valid_fragments": False}
        fragments = _create_fragments_from_read_group(read_group, args)
        if not fragments:
            return {"qname": qname, "has_valid_fragments": False}
        fasta_handle = _get_worker_fasta_handle() if args_dict.get('polyA_tail_detection', False) else None
        strandness = args_dict.get('strandness', 'unstrand')
        gene_overlap_overflow_count = 0
        for frag in fragments:
            chrom = frag['chrom']
            frag['has_splice'] = frag.get('has_splice_cigar', False)
            if _worker_splice_junctions and chrom in _worker_splice_junctions:
                known_junctions = _worker_splice_junctions[chrom]
                if frag.get('blocks_r1'):
                    frag['has_splice'] |= _detect_splice_junction(frag['blocks_r1'], known_junctions)
                if frag.get('blocks_r2'):
                    frag['has_splice'] |= _detect_splice_junction(frag['blocks_r2'], known_junctions)
            has_pas = False
            if args_dict.get('polyA_tail_detection', False):
                min_polya = args_dict.get('min_polya_length', POLYADENOSINE_MIN_LENGTH)
                if frag.get('read1_obj'):
                    has_pas |= _detect_polyadenylation(frag['read1_obj'], fasta_handle, strandness, min_polya)
                if frag.get('read2_obj'):
                    has_pas |= _detect_polyadenylation(frag['read2_obj'], fasta_handle, strandness, min_polya)
            frag['has_pas'] = has_pas
        if _worker_gene_structs:
            for frag in fragments:
                chrom = frag['chrom']
                if chrom in _worker_gene_structs:
                    chrom_struct = _worker_gene_structs[chrom]
                    overlaps_r1, overflow_r1 = _get_gene_overlaps(frag.get('blocks_r1', []), chrom_struct)
                    overlaps_r2, overflow_r2 = _get_gene_overlaps(frag.get('blocks_r2', []), chrom_struct)
                    if overflow_r1 or overflow_r2:
                        gene_overlap_overflow_count += 1
                    all_overlaps_idx = defaultdict(int)
                    for gene_idx, length in overlaps_r1.items():
                        all_overlaps_idx[gene_idx] += length
                    for gene_idx, length in overlaps_r2.items():
                        all_overlaps_idx[gene_idx] += length
                    filtered_overlaps = {}
                    if _worker_gene_info is not None and args.strandness != 'unstrand':
                        for gene_idx, overlap_len in all_overlaps_idx.items():
                            gene_strand_val = _worker_gene_info[gene_idx]
                            if _check_fragment_strand_compatibility(frag, gene_strand_val, args.strandness):
                                gene_id = _worker_gene_map_rev.get(gene_idx)
                                if gene_id:
                                    filtered_overlaps[gene_id] = overlap_len
                    else:
                        for gene_idx, overlap_len in all_overlaps_idx.items():
                            gene_id = _worker_gene_map_rev.get(gene_idx)
                            if gene_id:
                                filtered_overlaps[gene_id] = overlap_len
                    frag['gene_overlaps'] = filtered_overlaps
                else:
                    frag['gene_overlaps'] = {}
        is_pos_unique = len(fragments) == 1 and fragments[0].get('is_unique', False)
        result = {"fragments": fragments, "rpm_count_alignments": rpm_count_alignments,
                  "is_pos_unique": is_pos_unique, "qname": qname, "has_valid_fragments": True}
        if gene_overlap_overflow_count > 0:
            result["gene_overlap_overflow"] = gene_overlap_overflow_count
        return result
    except Exception as e:
        return {"error": str(e), "qname": qname, "has_valid_fragments": False}
# =============================================================================
# BAM Reading and EM Algorithms
# =============================================================================
def read_group_generator(bam_path: str, keep_duplicates: bool, logger: logging.Logger) -> Generator[List, None, None]:
    """Read name-sorted BAM and yield read groups."""
    logger.info(f"Reading BAM: {os.path.basename(bam_path)}")
    read_group = []
    current_qname = None
    reads_processed = 0
    try:
        with pysam.AlignmentFile(bam_path, "rb", check_sq=False, require_index=False) as bam:
            filter_flags = pysam.FUNMAP | pysam.FQCFAIL | pysam.FSECONDARY | pysam.FSUPPLEMENTARY
            if not keep_duplicates:
                filter_flags |= pysam.FDUP
            for read in (r for r in bam if not (r.flag & filter_flags)):
                reads_processed += 1
                if reads_processed % 10_000_000 == 0:
                    logger.info(f"Read ~{reads_processed // 1_000_000}M alignments...")
                if read.query_name != current_qname:
                    if read_group:
                        yield read_group
                    read_group = [read]
                    current_qname = read.query_name
                else:
                    read_group.append(read)
            if read_group:
                yield read_group
    except Exception as e:
        logger.error(f"FATAL error reading BAM: {e}", exc_info=True)
        raise
    finally:
        logger.info(f"Finished reading BAM. Total alignments: {reads_processed:,}")
def run_positional_em(multi_mapping_dict, targeted_sparse_unique, max_iterations, tolerance, logger,
                      chrom_len_map, bin_size=1, unique_bins=None):
    """Round 1 EM: Resolve multi-mapping reads using coverage integral."""
    logger.info("Starting Round 1: Positional EM algorithm...")
    fraction_assignments = {}
    for key, aligns in multi_mapping_dict.items():
        scores = np.array([a['score'] for a in aligns], dtype=np.float64)
        delta = np.clip(scores - np.max(scores), -700, 700)
        fraction_assignments[key] = fast_em_update_numba(np.exp(delta), EM_STABILITY_EPSILON)
    for iteration in range(1, max_iterations + 1):
        start_time = time.time()
        merged_multi_coverage = {}
        merged_multi_bins = {}
        if bin_size == 1:
            multi_coverage_sparse = defaultdict(list)
            for qname, fracs in fraction_assignments.items():
                for i, frac in enumerate(fracs):
                    if frac < EM_STABILITY_EPSILON:
                        continue
                    frag = multi_mapping_dict[qname][i]
                    all_blocks = frag.get('blocks_r1', []) + frag.get('blocks_r2', [])
                    for start, end in all_blocks:
                        if end > start:
                            multi_coverage_sparse[frag['chrom']].append((start, end, frac))
            for chrom, intervals in multi_coverage_sparse.items():
                if not intervals:
                    continue
                scanline = defaultdict(float)
                for start, end, value in intervals:
                    scanline[start] += value
                    scanline[end] -= value
                sorted_pos = sorted(scanline.keys())
                if not sorted_pos:
                    continue
                merged = []
                level = 0.0
                last_pos = sorted_pos[0]
                for pos in sorted_pos:
                    if pos > last_pos and level > EM_STABILITY_EPSILON:
                        merged.append((last_pos, pos, level))
                    level += scanline[pos]
                    last_pos = pos
                if merged:
                    starts = np.array([m[0] for m in merged], dtype=np.int32)
                    ends = np.array([m[1] for m in merged], dtype=np.int32)
                    values = np.array([max(0.0, m[2]) for m in merged], dtype=np.float64)
                    merged_multi_coverage[chrom] = (starts, ends, values)
        else:
            for qname, fracs in fraction_assignments.items():
                for i, frac in enumerate(fracs):
                    if frac < EM_STABILITY_EPSILON:
                        continue
                    frag = multi_mapping_dict[qname][i]
                    all_blocks = frag.get('blocks_r1', []) + frag.get('blocks_r2', [])
                    _accumulate_blocks_into_bins(merged_multi_bins, chrom_len_map, frag['chrom'],
                                                 all_blocks, frac, bin_size)
        new_fractions = {}
        max_alpha_diff = 0.0
        for qname, alignments in multi_mapping_dict.items():
            n_aligns = len(alignments)
            weights = np.zeros(n_aligns, dtype=np.float64)
            scores = np.array([align['score'] for align in alignments], dtype=np.float64)
            delta = np.clip(scores - np.max(scores), -700, 700)
            score_weights = np.exp(delta)
            for i, frag in enumerate(alignments):
                all_blocks = frag.get('blocks_r1', []) + frag.get('blocks_r2', [])
                chrom = frag['chrom']
                unique_cov = 0.0
                multi_cov = 0.0
                if bin_size == 1:
                    if targeted_sparse_unique and chrom in targeted_sparse_unique:
                        starts, ends, values = targeted_sparse_unique[chrom]
                        for st, en in all_blocks:
                            unique_cov += get_sum_from_sparse_numba(starts, ends, values, st, en)
                    if chrom in merged_multi_coverage:
                        starts, ends, values = merged_multi_coverage[chrom]
                        for st, en in all_blocks:
                            multi_cov += get_sum_from_sparse_numba(starts, ends, values, st, en)
                else:
                    block_starts = np.array([s for s, _ in all_blocks], dtype=np.int32)
                    block_ends = np.array([e for _, e in all_blocks], dtype=np.int32)
                    if unique_bins and chrom in unique_bins:
                        unique_cov = _sum_blocks_on_bins(unique_bins[chrom], bin_size, block_starts, block_ends)
                    if chrom in merged_multi_bins:
                        multi_cov = _sum_blocks_on_bins(merged_multi_bins[chrom], bin_size, block_starts, block_ends)
                total_cov = max(0.0, unique_cov + multi_cov)
                weights[i] = (total_cov + EM_STABILITY_EPSILON) * score_weights[i]
            new_fracs = fast_em_update_numba(weights, EM_STABILITY_EPSILON)
            current_max_diff = np.max(np.abs(new_fracs - fraction_assignments[qname]))
            if current_max_diff > max_alpha_diff:
                max_alpha_diff = current_max_diff
            new_fractions[qname] = new_fracs
        fraction_assignments = new_fractions
        logger.info(f"[Position EM iter {iteration}] Time: {time.time() - start_time:.1f}s, Max Alpha Diff: {max_alpha_diff:.4e}")
        if max_alpha_diff < tolerance:
            logger.info(f"Positional EM converged after {iteration} iterations.")
            break
    return fraction_assignments
def run_gene_em_algorithm_scientific(ambiguous_reads_list, initial_uniq_counts, initial_multi_counts,
                                     gene_lengths, valid_gene_ids, single_exon_genes, args, logger,
                                     mean_frag_len: float = 250.0):
    """Round 2 EM: Scientific gene-level EM using effective length and abundance weights.
    Returns:
      Tuple of 6 dictionaries + convergence flag:
        (em_uniq_counts, em_multi_counts, em_splice_uniq, em_splice_multi, em_pas_uniq, em_pas_multi, converged)
    """
    logger.info("Starting Round 2: Scientific Gene-level EM algorithm...")
    logger.info(f"Using mean fragment length: {mean_frag_len:.1f} bp")
    empty_result = (defaultdict(float), defaultdict(float), defaultdict(float),
                    defaultdict(float), defaultdict(float), defaultdict(float), True)
    if not ambiguous_reads_list:
        logger.warning("No gene-ambiguous reads for Gene EM.")
        return empty_result
    gene_eff_lengths = {}
    for gene_id in valid_gene_ids:
        gene_len = gene_lengths.get(gene_id, 0)
        if gene_len > 0:
            gene_eff_lengths[gene_id] = max(1.0, gene_len - mean_frag_len + 1.0)
        else:
            gene_eff_lengths[gene_id] = 1.0
    initial_counts = defaultdict(float)
    for gene_id, count in initial_uniq_counts.items():
        initial_counts[gene_id] += count
    for gene_id, count in initial_multi_counts.items():
        initial_counts[gene_id] += count
    def calculate_abundance(counts_dict):
        """Calculate normalized abundance alpha for each gene."""
        alpha = {}
        total_alpha = 0.0
        for gene_id in valid_gene_ids:
            count = counts_dict.get(gene_id, 0.0) + GENE_EM_PSEUDOCOUNT
            eff_len = gene_eff_lengths.get(gene_id, 1.0)
            a = count / eff_len
            alpha[gene_id] = a
            total_alpha += a
        if total_alpha > 0:
            for gene_id in alpha:
                alpha[gene_id] /= total_alpha
        return alpha
    current_alpha = calculate_abundance(initial_counts)
    all_alpha_zero = all(current_alpha.get(g, 0.0) <= EM_STABILITY_EPSILON for g in valid_gene_ids)
    if all_alpha_zero:
        logger.warning("All initial abundances are zero. Using uniform distribution.")
        n_genes = len(valid_gene_ids)
        for gene_id in valid_gene_ids:
            current_alpha[gene_id] = 1.0 / n_genes
    final_em_uniq = defaultdict(float)
    final_em_multi = defaultdict(float)
    final_splice_uniq = defaultdict(float)
    final_splice_multi = defaultdict(float)
    final_pas_uniq = defaultdict(float)
    final_pas_multi = defaultdict(float)
    em_converged = False
    for iteration in range(1, args.max_iter + 1):
        start_time = time.time()
        iter_em_uniq = defaultdict(float)
        iter_em_multi = defaultdict(float)
        iter_splice_uniq = defaultdict(float)
        iter_splice_multi = defaultdict(float)
        iter_pas_uniq = defaultdict(float)
        iter_pas_multi = defaultdict(float)
        for compatible_genes, pos_prob, is_genomic_unique, frag in ambiguous_reads_list:
            if pos_prob < EM_STABILITY_EPSILON:
                continue
            has_splice = frag.get('has_splice', False)
            has_pas = frag.get('has_pas', False)
            if has_splice:
                splice_compatible = {g for g in compatible_genes if g not in single_exon_genes}
                effective_set = splice_compatible if splice_compatible else compatible_genes
            else:
                effective_set = compatible_genes
            gene_list = list(effective_set)
            n_genes = len(gene_list)
            if n_genes == 0:
                continue
            alpha_values = np.array([current_alpha.get(g, 0.0) for g in gene_list], dtype=np.float64)
            alpha_sum = np.sum(alpha_values)
            if alpha_sum <= EM_STABILITY_EPSILON:
                probs = np.full(n_genes, 1.0 / n_genes, dtype=np.float64)
            else:
                probs = alpha_values / alpha_sum
            for i, gene_id in enumerate(gene_list):
                val = pos_prob * probs[i]
                if val < EM_STABILITY_EPSILON:
                    continue
                if is_genomic_unique:
                    iter_em_uniq[gene_id] += val
                    if has_splice and gene_id not in single_exon_genes:
                        iter_splice_uniq[gene_id] += val
                    if has_pas:
                        iter_pas_uniq[gene_id] += val
                else:
                    iter_em_multi[gene_id] += val
                    if has_splice and gene_id not in single_exon_genes:
                        iter_splice_multi[gene_id] += val
                    if has_pas:
                        iter_pas_multi[gene_id] += val
        new_total_counts = defaultdict(float)
        for gene_id in valid_gene_ids:
            new_total_counts[gene_id] = (initial_counts.get(gene_id, 0.0) +
                                         iter_em_uniq.get(gene_id, 0.0) +
                                         iter_em_multi.get(gene_id, 0.0))
        new_alpha = calculate_abundance(new_total_counts)
        max_alpha_diff = 0.0
        for gene_id in valid_gene_ids:
            diff = abs(new_alpha.get(gene_id, 0.0) - current_alpha.get(gene_id, 0.0))
            if diff > max_alpha_diff:
                max_alpha_diff = diff
        current_alpha = new_alpha
        final_em_uniq = iter_em_uniq
        final_em_multi = iter_em_multi
        final_splice_uniq = iter_splice_uniq
        final_splice_multi = iter_splice_multi
        final_pas_uniq = iter_pas_uniq
        final_pas_multi = iter_pas_multi
        logger.info(f"[Gene EM iter {iteration}] Time: {time.time() - start_time:.1f}s, "
                   f"Max Alpha Diff: {max_alpha_diff:.4e}")
        if max_alpha_diff < args.tol:
            logger.info(f"Gene-level EM converged after {iteration} iterations.")
            em_converged = True
            break
    return (final_em_uniq, final_em_multi, final_splice_uniq, final_splice_multi,
            final_pas_uniq, final_pas_multi, em_converged)
# =============================================================================
# Argument Parsing
# =============================================================================
def parse_arguments() -> argparse.Namespace:
    """Parse and validate command-line arguments."""
    default_threads = 10
    description_text = f"""
================================================================================
EMapper v{__version__}: High-Performance EM-based RNA-seq Quantification Tool
================================================================================
A robust, memory-safe quantification tool for RNA-seq data that uses a two-round
Expectation-Maximization (EM) algorithm to resolve multi-mapping reads and
gene-level ambiguity.

Key Features:
  - Two-round EM: Positional EM for multi-mappers, Gene EM for overlapping genes
  - Splice junction detection from CIGAR and GTF annotations
  - Poly(A) site detection with internal priming filtering
  - Strand-specific quantification support (forward/reverse/unstranded)
  - BigWig coverage track generation with RPM normalization

Output Files:
  - *_readcounts.txt: Gene-level count matrix with TPM values
  - *_unstranded.bw: Combined coverage BigWig track
  - *_F1R2.bw / *_F2R1.bw: Strand-specific BigWig tracks (if stranded)
  - *_EMapper.log: Detailed processing log

Example Usage:
  python EMapper.py --input_bam sample.bam --gtf genes.gtf --strandness reverse
  python EMapper.py --input_bam sample.bam --gtf genes.gtf --polyA_tail_detection --reference_fasta genome.fa
================================================================================
"""
    parser = argparse.ArgumentParser(
        description=description_text,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    req = parser.add_argument_group('REQUIRED ARGUMENTS')
    req.add_argument("--input_bam", required=True, metavar="FILE",
                     help="Input BAM file path (name-sorted recommended; will auto-sort if coordinate-sorted)")
    req.add_argument("--gtf", required=True, metavar="FILE",
                     help="GTF/GFF annotation file with exon features for gene quantification")
    core = parser.add_argument_group('CORE PROCESSING OPTIONS')
    core.add_argument("--strandness", choices=["forward", "reverse", "unstrand"], default="unstrand", metavar="MODE",
                      help="Library strand-specificity: 'forward' (sense R1), 'reverse' (dUTP, antisense R1), 'unstrand' (unstranded). [Default: unstrand]")
    core.add_argument("--positional_assignment", choices=["unique", "primary", "EM"], default="EM", metavar="MODE",
                      help="Multi-mapping read handling: 'unique' (discard multi-mappers), 'primary' (use primary alignment), 'EM' (probabilistic assignment). [Default: EM]")
    core.add_argument("--gene_disambiguation", choices=["unique", "EM"], default="EM", metavar="MODE",
                      help="Gene overlap resolution: 'unique' (discard ambiguous), 'EM' (probabilistic assignment based on abundance). [Default: EM]")
    core.add_argument("--keep_duplicates", action="store_true",
                      help="Include PCR/optical duplicates in quantification (reads with 0x400 flag). [Default: exclude duplicates]")
    out = parser.add_argument_group('OUTPUT OPTIONS')
    out.add_argument("--output_dir", metavar="DIR",
                     help="Output directory for all result files. [Default: <bam_basename>_EMapper_output]")
    out.add_argument("--prefix", metavar="STR",
                     help="Prefix for output file names. [Default: input BAM basename]")
    out.add_argument("--normalize_bigwig", choices=['RPM', 'none'], default='RPM', metavar="METHOD",
                     help="BigWig normalization method: 'RPM' (reads per million), 'none' (raw counts). [Default: RPM]")
    out.add_argument("--bin_size", type=int, default=1, metavar="INT",
                     help="BigWig resolution in base pairs; 1 = single-base resolution, higher values reduce file size. [Default: 1]")
    out.add_argument("--no-cleanup", action="store_true", dest="no_cleanup",
                     help="Retain intermediate files for debugging: name-sorted BAM (*_namesorted.bam) and "
                          "diagnostic BigWig files separating unique/multi coverage (*_uniq_only.bw, *_multi_only.bw). "
                          "By default, these files are not generated/retained.")
    polya = parser.add_argument_group('POLY(A) DETECTION OPTIONS')
    polya.add_argument("--polyA_tail_detection", action="store_true",
                       help="Enable poly(A) tail detection from soft-clipped read sequences; requires --reference_fasta")
    polya.add_argument("--reference_fasta", metavar="FILE",
                       help="Reference genome FASTA file (indexed with samtools faidx) for internal priming filter")
    polya.add_argument("--min_polya_length", type=int, default=POLYADENOSINE_MIN_LENGTH, metavar="INT",
                       help=f"Minimum consecutive A/T bases to call poly(A) tail. [Default: {POLYADENOSINE_MIN_LENGTH}]")
    perf = parser.add_argument_group('PERFORMANCE OPTIONS')
    perf.add_argument("--num_threads", type=int, default=default_threads, metavar="INT",
                      help=f"Number of parallel worker threads for read processing. [Default: {default_threads}]")
    perf.add_argument("--chunksize", type=int, default=2000, metavar="INT",
                      help="Number of read groups per worker batch; higher values improve throughput but increase memory. [Default: 2000]")
    perf.add_argument("--max_ambiguous_groups", type=int, default=5_000_000, metavar="INT",
                      help="Maximum number of ambiguous read groups to store in memory; excess groups are discarded. [Default: 5000000]")
    em = parser.add_argument_group('EM ALGORITHM TUNING')
    em.add_argument("--max_iter", type=int, default=100, metavar="INT",
                    help="Maximum EM iterations before forced convergence. [Default: 100]")
    em.add_argument("--tol", type=float, default=1e-3, metavar="FLOAT",
                    help="Convergence tolerance; EM stops when max probability change < tol. [Default: 0.001]")
    filt = parser.add_argument_group('READ FILTERING OPTIONS')
    filt.add_argument("--min_gene_overlap_bp", type=int, default=10, metavar="INT",
                      help="Minimum overlap in base pairs between read and gene exons for assignment. [Default: 10]")
    filt.add_argument("--min_gene_overlap_frac", type=float, default=0.01, metavar="FLOAT",
                      help="Minimum overlap as fraction of read length (0.0-1.0) for gene assignment. [Default: 0.01]")
    filt.add_argument("--min_insert_size", type=int, default=10, metavar="INT",
                      help="Minimum fragment insert size for proper pair detection. [Default: 10]")
    filt.add_argument("--max_insert_size", type=int, default=2000000, metavar="INT",
                      help="Maximum fragment insert size for proper pair detection. [Default: 2000000]")
    filt.add_argument("--min_intron_len", type=int, default=DEFAULT_MIN_INTRON_LEN, metavar="INT",
                      help=f"Minimum CIGAR N operation length to call splice junction. [Default: {DEFAULT_MIN_INTRON_LEN}]")
    filt.add_argument("--unique_mapping_threshold", type=int, default=None, metavar="INT",
                      help="MAPQ threshold for unique mapping when NH tag absent; auto-detected from BAM if not specified. [Default: auto]")
    adv = parser.add_argument_group('ADVANCED OPTIONS')
    adv.add_argument("--max_alignments_per_mate", type=int, default=1000, metavar="INT",
                     help="Maximum alignments to consider per read mate (R1/R2) when pairing; limits memory for highly multi-mapped reads. [Default: 1000]")
    adv.add_argument("--bin_value", choices=["mean", "sum"], default="sum", metavar="MODE",
                     help="Bin aggregation method for binned BigWig: 'sum' (total coverage), 'mean' (average depth). [Default: sum]")
    adv.add_argument("--mean_frag_len", type=float, default=None, metavar="FLOAT", help="Mean fragment length (bp) for TPM calculation. Auto-detected for paired-end; defaults to 250 for single-end.")
    args = parser.parse_args()
    if not os.path.exists(args.input_bam):
        parser.error(f"Input BAM file not found: {args.input_bam}")
    if not os.path.exists(args.gtf):
        parser.error(f"GTF file not found: {args.gtf}")
    if args.reference_fasta and not os.path.exists(args.reference_fasta):
        parser.error(f"Reference FASTA not found: {args.reference_fasta}")
    if args.bin_size <= 0:
        parser.error("--bin_size must be positive.")
    if not (0.0 <= args.min_gene_overlap_frac <= 1.0):
        parser.error("--min_gene_overlap_frac must be between 0.0 and 1.0.")
    if args.polyA_tail_detection and not args.reference_fasta:
        parser.error("--polyA_tail_detection requires --reference_fasta.")
    if args.num_threads <= 0:
        parser.error("--num_threads must be positive.")
    if args.max_iter <= 0:
        parser.error("--max_iter must be positive.")
    if args.tol <= 0:
        parser.error("--tol must be positive.")
    if args.mean_frag_len is not None and args.mean_frag_len <= 0:
        parser.error("--mean_frag_len must be positive.")
    return args
# =============================================================================
# BigWig Writing
# =============================================================================
def write_bigwig_memory_safe(coverage_data: Dict, output_path: str, ref_lengths: List[Tuple[str, int]],
                             logger: logging.Logger, scale_factor: float = 1.0):
    """Write coverage to BigWig using scanline algorithm."""
    logger.info(f"Writing BigWig: {os.path.basename(output_path)}")
    try:
        with pyBigWig.open(output_path, "w") as bw:
            bw.addHeader(ref_lengths)
            for chrom, length in ref_lengths:
                if chrom not in coverage_data:
                    continue
                scanline = defaultdict(float)
                for start, end, value in coverage_data[chrom]:
                    start = max(0, start)
                    end = min(length, end)
                    if start >= end:
                        continue
                    scanline[start] += value * scale_factor
                    scanline[end] -= value * scale_factor
                if not scanline:
                    continue
                sorted_pos = sorted(scanline.keys())
                starts, ends, values = [], [], []
                level = 0.0
                last_pos = sorted_pos[0]
                for pos in sorted_pos:
                    if pos > last_pos and level > EM_STABILITY_EPSILON:
                        starts.append(int(last_pos))
                        ends.append(int(pos))
                        values.append(max(0.0, level))
                    level += scanline[pos]
                    last_pos = pos
                if starts:
                    bw.addEntries([chrom] * len(starts), starts, ends=ends, values=values)
    except Exception as e:
        logger.error(f"Failed to write BigWig {output_path}: {e}", exc_info=True)
        raise
def _write_binned_bigwig(bin_dict: Dict[str, np.ndarray], out_path: str, ref_lengths: List[Tuple[str, int]],
                         logger, scale_factor: float, bin_size: int, bin_value_mode: str):
    """Write binned coverage to BigWig."""
    logger.info(f"Writing binned BigWig: {os.path.basename(out_path)}")
    try:
        with pyBigWig.open(out_path, "w") as bw:
            bw.addHeader(ref_lengths)
            for chrom, chrom_len in ref_lengths:
                arr = bin_dict.get(chrom)
                if arr is None or len(arr) == 0:
                    continue
                values = arr * bin_size if bin_value_mode == "sum" else arr.copy()
                if abs(scale_factor - 1.0) > 1e-12:
                    values = values * scale_factor
                n_bins = len(values)
                starts = np.arange(n_bins, dtype=np.int32) * bin_size
                ends = np.minimum(starts + bin_size, chrom_len)
                mask = (ends > starts) & (values > EM_STABILITY_EPSILON)
                if np.any(mask):
                    bw.addEntries([chrom] * int(np.sum(mask)), starts[mask].tolist(),
                                  ends=ends[mask].tolist(), values=values[mask].tolist())
    except Exception as e:
        logger.error(f"Failed to write binned BigWig: {e}", exc_info=True)
        raise
def setup_logging(log_path: str, log_queue: mp.Queue) -> Tuple[logging.Logger, logging.handlers.QueueListener]:
    """Set up logging with queue handler for multiprocessing."""
    logger = logging.getLogger('EMapper')
    logger.handlers.clear()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('[%(asctime)s] [%(levelname)s] [%(processName)s] %(message)s', "%Y-%m-%d %H:%M:%S")
    file_handler = logging.FileHandler(log_path, mode='w')
    file_handler.setFormatter(formatter)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    listener = logging.handlers.QueueListener(log_queue, file_handler, console_handler, respect_handler_level=True)
    queue_handler = logging.handlers.QueueHandler(log_queue)
    logger.addHandler(queue_handler)
    logger.propagate = False
    return logger, listener
def check_and_prepare_bam(input_bam: str, output_dir: str, prefix: str, num_threads: int, logger: logging.Logger) -> str:
    """Check if BAM is name-sorted, sort if needed."""
    logger.info("--- BAM Preparation Step ---")
    try:
        with pysam.AlignmentFile(input_bam, "rb", check_sq=False, require_index=False) as bam_file:
            if bam_file.header.get('HD', {}).get('SO') == 'queryname':
                logger.info("BAM is name-sorted.")
                return input_bam
    except Exception:
        pass
    sorted_bam_path = os.path.join(output_dir, f"{prefix}_namesorted.bam")
    logger.info(f"Sorting BAM by queryname to {os.path.basename(sorted_bam_path)}")
    try:
        subprocess.run(["samtools", "sort", "-n", "-@", str(num_threads), "-o", sorted_bam_path, input_bam],
                       check=True, text=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    except FileNotFoundError:
        logger.error("samtools not found in PATH.")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        logger.error(f"Samtools sort failed: {e.stderr}")
        sys.exit(1)
    logger.info("Successfully created name-sorted BAM.")
    return sorted_bam_path
def calculate_tpm_values(counts_dict: Dict[str, float], gene_lengths: Dict[str, int],
                         mean_frag_len: float = 250.0, total_rpk_override: float = None) -> Dict[str, float]:
    """Calculate TPM values using effective length."""
    rpk = {}
    total_rpk = 0.0
    for gene_id, count in counts_dict.items():
        gene_len = gene_lengths.get(gene_id, 0)
        if gene_len <= 0 or count <= 0:
            continue
        eff_len = max(1.0, gene_len - mean_frag_len + 1.0)
        current_rpk = count / (eff_len / 1000.0)
        rpk[gene_id] = current_rpk
        total_rpk += current_rpk
    tpm = {}
    normalization_base = total_rpk_override if total_rpk_override is not None else total_rpk
    if normalization_base > 0:
        scale = 1e6 / normalization_base
        for gene_id, current_rpk in rpk.items():
            tpm[gene_id] = current_rpk * scale
    return tpm
def write_count_files(output_dir: str, prefix: str, uniq_counts: Dict, multi_counts: Dict,
                      gene_lengths: Dict, valid_gene_ids: Set, args, logger: logging.Logger,
                      splice_counts: Dict, pas_counts: Dict, mean_frag_len: float):
    """Write gene count file."""
    logger.info("Writing gene count file...")
    output_path = os.path.join(output_dir, f"{prefix}_readcounts.txt")
    total_counts = {}
    for gene_id in valid_gene_ids:
        total_counts[gene_id] = uniq_counts.get(gene_id, 0.0) + multi_counts.get(gene_id, 0.0)
    total_rpk = 0.0
    for gene_id, count in total_counts.items():
        gene_len = gene_lengths.get(gene_id, 0)
        if gene_len <= 0 or count <= 0:
            continue
        eff_len = max(1.0, gene_len - mean_frag_len + 1.0)
        total_rpk += count / (eff_len / 1000.0)
    tpm_total = calculate_tpm_values(total_counts, gene_lengths, mean_frag_len, total_rpk_override=total_rpk)
    tpm_uniq = calculate_tpm_values(uniq_counts, gene_lengths, mean_frag_len, total_rpk_override=total_rpk)
    all_gene_ids = sorted((set(uniq_counts.keys()) | set(multi_counts.keys()) |
                           set(splice_counts.keys()) | set(pas_counts.keys())) & valid_gene_ids)
    with open(output_path, "w") as f:
        f.write("GeneID\tLength\tCounts_Uniq\tTPM_Uniq\tCounts_Multi\tCounts_Total\tTPM\t"
                "Counts_Splice\tCounts_PAS\n")
        for gene_id in all_gene_ids:
            uniq = uniq_counts.get(gene_id, 0.0)
            multi = multi_counts.get(gene_id, 0.0)
            total = uniq + multi
            gene_len = gene_lengths.get(gene_id, 0)
            tpm_u = tpm_uniq.get(gene_id, 0.0)
            tpm_t = tpm_total.get(gene_id, 0.0)
            splice = splice_counts.get(gene_id, 0.0)
            pas = pas_counts.get(gene_id, 0.0)
            tpm_u_str = format_count(tpm_u) if tpm_u == int(tpm_u) else f"{tpm_u:.2f}".rstrip('0').rstrip('.')
            tpm_t_str = format_count(tpm_t) if tpm_t == int(tpm_t) else f"{tpm_t:.2f}".rstrip('0').rstrip('.')
            f.write(f"{gene_id}\t{gene_len}\t{format_count(uniq)}\t{tpm_u_str}\t"
                    f"{format_count(multi)}\t{format_count(total)}\t{tpm_t_str}\t"
                    f"{format_count(splice)}\t{format_count(pas)}\n")
    logger.info(f"Gene counts written to {os.path.basename(output_path)}")
def init_worker(gene_structs_data, gene_map_rev_data, gene_info_data, log_queue: mp.Queue,
                splice_junctions_data, single_exon_genes_data, fasta_path):
    """Initialize worker process with all required data."""
    global _worker_gene_structs, _worker_gene_map_rev, _worker_gene_info, _worker_log_queue
    global _worker_splice_junctions, _worker_single_exon_genes, _worker_fasta_path
    global _worker_fasta_handle, _process_overlap_gene_idx, _process_overlap_lengths
    if not _import_dependencies():
        raise RuntimeError("Failed to import dependencies in worker process.")
    _ensure_numba_compiled()
    _worker_gene_structs = gene_structs_data
    _worker_gene_map_rev = gene_map_rev_data
    _worker_gene_info = gene_info_data
    _worker_splice_junctions = splice_junctions_data
    _worker_single_exon_genes = single_exon_genes_data
    _worker_log_queue = log_queue
    _worker_fasta_path = fasta_path
    _worker_fasta_handle = None
    _process_overlap_gene_idx = None
    _process_overlap_lengths = None
    worker_logger = logging.getLogger('EMapper')
    worker_logger.handlers.clear()
    worker_logger.addHandler(logging.handlers.QueueHandler(_worker_log_queue))
    worker_logger.setLevel(logging.WARNING)
    worker_logger.propagate = False
def _aggregate_worker_results(results_iter, args, valid_gene_ids, chrom_len_map, logger: logging.Logger,
                              single_exon_genes: Set[str]):
    """Aggregate worker results with genomic-position-based count classification."""
    final_unique_f1r2 = defaultdict(list)
    final_unique_f2r1 = defaultdict(list)
    final_unique_bins_f1r2 = {}
    final_unique_bins_f2r1 = {}
    final_pos_multi_map_dict = {}
    genomic_uniq_counts = defaultdict(float)
    splice_counts_uniq = defaultdict(float)
    pas_counts_uniq = defaultdict(float)
    pos_ambig_groups_count = 0
    total_read_groups = 0
    uniq_multigene_discard_count = 0
    skipped_multi_count = 0
    gene_overlap_overflow_total = 0
    num_processed_groups = 0
    error_count = 0
    final_gene_ambig_list = []
    gene_ambig_count = 0
    for num_processed_groups, result in enumerate(results_iter, 1):
        if not result:
            continue
        if "error" in result:
            error_count += 1
            if error_count <= 10:
                logger.warning(f"Worker error for {result.get('qname', 'unknown')}: {result['error']}")
            continue
        if num_processed_groups % 50000 == 0:
            logger.info(f"Aggregating ~{num_processed_groups:,} read groups...")
        total_read_groups += 1
        if result.get("skipped_multi", False):
            skipped_multi_count += 1
            continue
        if not result.get("has_valid_fragments", False):
            continue
        fragments = result.get("fragments")
        if not fragments:
            continue
        if result.get("gene_overlap_overflow", 0) > 0:
            gene_overlap_overflow_total += result["gene_overlap_overflow"]
        qname = result.get("qname", f"unnamed_{num_processed_groups}")
        is_pos_unique = result.get("is_pos_unique", False)
        if is_pos_unique:
            frag = fragments[0]
            _accumulate_frag_for_coverage(frag, args, chrom_len_map, final_unique_f1r2, final_unique_f2r1,
                                          final_unique_bins_f1r2, final_unique_bins_f2r1)
            if args.gtf:
                multi_genes = _process_frag_gene_counts_simple(frag, args, valid_gene_ids)
                if multi_genes is None:
                    pass
                elif len(multi_genes) == 1:
                    gene_id = next(iter(multi_genes))
                    genomic_uniq_counts[gene_id] += 1.0
                    if frag.get('has_splice', False) and gene_id not in single_exon_genes:
                        splice_counts_uniq[gene_id] += 1.0
                    if frag.get('has_pas', False):
                        pas_counts_uniq[gene_id] += 1.0
                else:
                    if args.gene_disambiguation == 'EM':
                        entry = (multi_genes, 1.0, True, frag)
                        if len(final_gene_ambig_list) >= args.max_ambiguous_groups:
                            logger.error(f"FATAL: Ambiguous groups exceed limit.")
                            sys.exit(1)
                        final_gene_ambig_list.append(entry)
                        gene_ambig_count += 1
                    else:
                        uniq_multigene_discard_count += 1
        elif args.positional_assignment == 'EM':
            pos_ambig_groups_count += 1
            if len(final_pos_multi_map_dict) >= args.max_ambiguous_groups:
                logger.error(f"FATAL: Positional ambiguous groups exceed limit.")
                sys.exit(1)
            final_pos_multi_map_dict[qname] = fragments
        else:
            for frag in fragments:
                _accumulate_frag_for_coverage(frag, args, chrom_len_map, final_unique_f1r2, final_unique_f2r1,
                                              final_unique_bins_f1r2, final_unique_bins_f2r1)
    if error_count > 0:
        logger.warning(f"Total worker errors: {error_count}")
    if skipped_multi_count > 0:
        logger.info(f"Skipped {skipped_multi_count:,} multi-mapping read groups (unique mode).")
    if gene_overlap_overflow_total > 0:
        logger.warning(f"Gene overlap capacity exceeded for {gene_overlap_overflow_total:,} fragments.")
    return (final_unique_f1r2, final_unique_f2r1, final_unique_bins_f1r2, final_unique_bins_f2r1,
            final_pos_multi_map_dict, final_gene_ambig_list, genomic_uniq_counts, pos_ambig_groups_count,
            total_read_groups, uniq_multigene_discard_count, num_processed_groups, total_read_groups,
            splice_counts_uniq, pas_counts_uniq, gene_ambig_count)
def _process_frag_gene_counts_simple(frag: Dict, args, valid_gene_ids: Set[str]) -> Optional[Set[str]]:
    """Get compatible genes for fragment."""
    total_len = get_effective_frag_len(frag)
    gene_overlaps = frag.get('gene_overlaps', {})
    compatible_genes = set()
    for gene_id, overlap_len in gene_overlaps.items():
        if (overlap_len >= args.min_gene_overlap_bp and
            (overlap_len / total_len) >= args.min_gene_overlap_frac and
            gene_id in valid_gene_ids):
            compatible_genes.add(gene_id)
    return compatible_genes if compatible_genes else None
def _resolve_positional_ambiguity(final_pos_multi_map_dict, final_unique_f1r2, final_unique_f2r1,
                                  final_unique_bins_f1r2, final_unique_bins_f2r1, args, chrom_len_map,
                                  valid_gene_ids, logger, single_exon_genes, max_ambiguous_groups):
    """Resolve positional ambiguity."""
    targeted_sparse_unique = None
    unique_bins = None
    if args.bin_size == 1:
        targeted_sparse_unique = {}
        all_chroms = set(final_unique_f1r2) | set(final_unique_f2r1)
        for chrom in all_chroms:
            scanline = defaultdict(float)
            for s, e, v in itertools.chain(final_unique_f1r2.get(chrom, []), final_unique_f2r1.get(chrom, [])):
                scanline[s] += v
                scanline[e] -= v
            if not scanline:
                continue
            sorted_pos = sorted(scanline.keys())
            intervals = []
            level = 0.0
            last_pos = sorted_pos[0]
            for pos in sorted_pos:
                if pos > last_pos and level > EM_STABILITY_EPSILON:
                    intervals.append((last_pos, pos, level))
                level += scanline[pos]
                last_pos = pos
            if intervals:
                starts = np.array([i[0] for i in intervals], dtype=np.int32)
                ends = np.array([i[1] for i in intervals], dtype=np.int32)
                values = np.array([max(0.0, i[2]) for i in intervals], dtype=np.float64)
                targeted_sparse_unique[chrom] = (starts, ends, values)
    else:
        unique_bins = {}
        for chrom in set(final_unique_bins_f1r2) | set(final_unique_bins_f2r1):
            arr1 = final_unique_bins_f1r2.get(chrom)
            arr2 = final_unique_bins_f2r1.get(chrom)
            if arr1 is None and arr2 is None:
                continue
            if arr1 is None:
                unique_bins[chrom] = np.maximum(0.0, arr2.copy())
            elif arr2 is None:
                unique_bins[chrom] = np.maximum(0.0, arr1.copy())
            else:
                n = max(len(arr1), len(arr2))
                a1 = np.pad(arr1, (0, n - len(arr1))) if len(arr1) < n else arr1
                a2 = np.pad(arr2, (0, n - len(arr2))) if len(arr2) < n else arr2
                unique_bins[chrom] = np.maximum(0.0, a1 + a2)
    final_fractions = run_positional_em(final_pos_multi_map_dict, targeted_sparse_unique, args.max_iter,
                                        args.tol, logger, chrom_len_map, bin_size=args.bin_size, unique_bins=unique_bins)
    multi_f1r2 = defaultdict(list)
    multi_f2r1 = defaultdict(list)
    multi_bins_f1r2 = {}
    multi_bins_f2r1 = {}
    genomic_multi_gene_single = defaultdict(float)
    splice_counts_multi = defaultdict(float)
    pas_counts_multi = defaultdict(float)
    all_ambiguous_for_gene_em = []
    discarded_count = 0
    for qname, fracs in final_fractions.items():
        frags = final_pos_multi_map_dict[qname]
        for i, frac in enumerate(fracs):
            if frac <= EM_STABILITY_EPSILON:
                continue
            frag = frags[i]
            _accumulate_frag_for_coverage(frag, args, chrom_len_map, multi_f1r2, multi_f2r1,
                                          multi_bins_f1r2, multi_bins_f2r1, weight=frac)
            if args.gtf:
                compatible = _process_frag_gene_counts_simple(frag, args, valid_gene_ids)
                if compatible is None:
                    continue
                elif len(compatible) == 1:
                    gene_id = next(iter(compatible))
                    genomic_multi_gene_single[gene_id] += frac
                    if frag.get('has_splice', False) and gene_id not in single_exon_genes:
                        splice_counts_multi[gene_id] += frac
                    if frag.get('has_pas', False):
                        pas_counts_multi[gene_id] += frac
                else:
                    if args.gene_disambiguation == 'EM':
                        if len(all_ambiguous_for_gene_em) >= max_ambiguous_groups:
                            discarded_count += 1
                            continue
                        all_ambiguous_for_gene_em.append((compatible, frac, False, frag))
    if discarded_count > 0:
        logger.warning(f"Post-EM gene-ambiguous groups exceeded limit. Discarded {discarded_count:,}.")
    return (multi_f1r2, multi_f2r1, multi_bins_f1r2, multi_bins_f2r1,
            genomic_multi_gene_single, splice_counts_multi, pas_counts_multi, all_ambiguous_for_gene_em)
def _merge_coverage_data(coverage_dicts: List[Dict], is_binned: bool) -> Dict:
    """Merge multiple coverage dictionaries."""
    if not any(coverage_dicts):
        return {}
    all_chroms = set.union(*(set(d.keys()) for d in coverage_dicts if d))
    if is_binned:
        merged = {}
        for chrom in all_chroms:
            arrays = [d.get(chrom) for d in coverage_dicts if d.get(chrom) is not None]
            if not arrays:
                continue
            max_len = max(len(arr) for arr in arrays)
            summed = np.zeros(max_len, dtype=np.float64)
            for arr in arrays:
                summed[:len(arr)] += arr
            merged[chrom] = summed
    else:
        merged = defaultdict(list)
        for d in coverage_dicts:
            if d:
                for chrom, data in d.items():
                    merged[chrom].extend(data)
    return merged
def _write_final_outputs(prefix, output_dir, ref_lengths_list, args, scale_factor, logger,
                         final_unique_f1r2, final_unique_f2r1, multi_f1r2, multi_f2r1,
                         final_unique_bins_f1r2, final_unique_bins_f2r1, multi_bins_f1r2, multi_bins_f2r1):
    """Write final BigWig output files."""
    bw_prefix = os.path.join(output_dir, prefix)
    if args.bin_size == 1:
        data = _merge_coverage_data([final_unique_f1r2, final_unique_f2r1, multi_f1r2, multi_f2r1], False)
        write_bigwig_memory_safe(data, f"{bw_prefix}_unstranded.bw", ref_lengths_list, logger, scale_factor)
    else:
        bins = _merge_coverage_data([final_unique_bins_f1r2, final_unique_bins_f2r1,
                                     multi_bins_f1r2, multi_bins_f2r1], True)
        _write_binned_bigwig(bins, f"{bw_prefix}_unstranded.bw", ref_lengths_list, logger,
                            scale_factor, args.bin_size, args.bin_value)
    if args.strandness != 'unstrand':
        if args.bin_size == 1:
            f1r2_data = _merge_coverage_data([final_unique_f1r2, multi_f1r2], False)
            write_bigwig_memory_safe(f1r2_data, f"{bw_prefix}_F1R2.bw", ref_lengths_list, logger, scale_factor)
            f2r1_data = _merge_coverage_data([final_unique_f2r1, multi_f2r1], False)
            write_bigwig_memory_safe(f2r1_data, f"{bw_prefix}_F2R1.bw", ref_lengths_list, logger, scale_factor)
        else:
            f1r2_bins = _merge_coverage_data([final_unique_bins_f1r2, multi_bins_f1r2], True)
            _write_binned_bigwig(f1r2_bins, f"{bw_prefix}_F1R2.bw", ref_lengths_list, logger,
                                scale_factor, args.bin_size, args.bin_value)
            f2r1_bins = _merge_coverage_data([final_unique_bins_f2r1, multi_bins_f2r1], True)
            _write_binned_bigwig(f2r1_bins, f"{bw_prefix}_F2R1.bw", ref_lengths_list, logger,
                                scale_factor, args.bin_size, args.bin_value)
    if args.no_cleanup:
        logger.info("Writing diagnostic BigWig files (unique/multi separated)...")
        if args.bin_size == 1:
            uniq_unstranded = _merge_coverage_data([final_unique_f1r2, final_unique_f2r1], False)
            write_bigwig_memory_safe(uniq_unstranded, f"{bw_prefix}_unstranded_uniq_only.bw", ref_lengths_list, logger, scale_factor)
            multi_unstranded = _merge_coverage_data([multi_f1r2, multi_f2r1], False)
            write_bigwig_memory_safe(multi_unstranded, f"{bw_prefix}_unstranded_multi_only.bw", ref_lengths_list, logger, scale_factor)
            if args.strandness != 'unstrand':
                write_bigwig_memory_safe(final_unique_f1r2, f"{bw_prefix}_F1R2_uniq_only.bw", ref_lengths_list, logger, scale_factor)
                write_bigwig_memory_safe(final_unique_f2r1, f"{bw_prefix}_F2R1_uniq_only.bw", ref_lengths_list, logger, scale_factor)
                write_bigwig_memory_safe(multi_f1r2, f"{bw_prefix}_F1R2_multi_only.bw", ref_lengths_list, logger, scale_factor)
                write_bigwig_memory_safe(multi_f2r1, f"{bw_prefix}_F2R1_multi_only.bw", ref_lengths_list, logger, scale_factor)
        else:
            uniq_bins_unstranded = _merge_coverage_data([final_unique_bins_f1r2, final_unique_bins_f2r1], True)
            _write_binned_bigwig(uniq_bins_unstranded, f"{bw_prefix}_unstranded_uniq_only.bw", ref_lengths_list, logger, scale_factor, args.bin_size, args.bin_value)
            multi_bins_unstranded = _merge_coverage_data([multi_bins_f1r2, multi_bins_f2r1], True)
            _write_binned_bigwig(multi_bins_unstranded, f"{bw_prefix}_unstranded_multi_only.bw", ref_lengths_list, logger, scale_factor, args.bin_size, args.bin_value)
            if args.strandness != 'unstrand':
                _write_binned_bigwig(final_unique_bins_f1r2, f"{bw_prefix}_F1R2_uniq_only.bw", ref_lengths_list, logger, scale_factor, args.bin_size, args.bin_value)
                _write_binned_bigwig(final_unique_bins_f2r1, f"{bw_prefix}_F2R1_uniq_only.bw", ref_lengths_list, logger, scale_factor, args.bin_size, args.bin_value)
                _write_binned_bigwig(multi_bins_f1r2, f"{bw_prefix}_F1R2_multi_only.bw", ref_lengths_list, logger, scale_factor, args.bin_size, args.bin_value)
                _write_binned_bigwig(multi_bins_f2r1, f"{bw_prefix}_F2R1_multi_only.bw", ref_lengths_list, logger, scale_factor, args.bin_size, args.bin_value)
# =============================================================================
# Main Function
# =============================================================================
def main():
    """Main execution function."""
    args = None
    logger = None
    log_listener = None
    pool_context = None
    is_temp_bam = False
    bam_to_process = None
    run_completed = False
    start_time = time.time()
    try:
        args = parse_arguments()
        if not _import_dependencies():
            sys.exit(1)
        args_dict = vars(args)
        bam_basename = os.path.basename(args.input_bam).rsplit('.', 1)[0]
        prefix = args.prefix if args.prefix else bam_basename
        output_dir = args.output_dir if args.output_dir else f"{prefix}_EMapper_output"
        os.makedirs(output_dir, exist_ok=True)
        log_queue = mp.Queue(-1)
        log_path = os.path.join(output_dir, f"{prefix}_EMapper.log")
        logger, log_listener = setup_logging(log_path, log_queue)
        log_listener.start()
        logger.info(f"--- EMapper v{__version__} ---")
        logger.info("Compiling Numba functions...")
        _ensure_numba_compiled()
        _warmup_numba_functions()
        if not os.path.isfile(args.input_bam):
            logger.error(f"FATAL: BAM not found: {args.input_bam}")
            sys.exit(1)
        if args.strandness == 'unstrand':
            logger.warning("NOTE: Running in UNSTRANDED mode.")
        bam_to_process = check_and_prepare_bam(args.input_bam, output_dir, prefix, args.num_threads, logger)
        is_temp_bam = (bam_to_process != args.input_bam)
        logger.info("--- Fragment Length Detection ---")
        if args.mean_frag_len is not None:
            logger.info(f"Using user-specified fragment length: {args.mean_frag_len} bp")
        else:
            mean_len, std_len, is_paired_end = estimate_fragment_length_from_bam(bam_to_process)
            if mean_len is not None and is_paired_end:
                args.mean_frag_len = mean_len
                logger.info(f"Auto-detected fragment length from paired-end data: {mean_len:.1f} ± {std_len:.1f} bp")
            else:
                args.mean_frag_len = 250.0
                logger.info(f"Single-end data or insufficient proper pairs; using default fragment length: 250.0 bp")
        args_dict['mean_frag_len'] = args.mean_frag_len
        with pysam.AlignmentFile(bam_to_process, "rb", check_sq=False, require_index=False) as bam:
            ref_lengths_list = [(r, l) for r, l in zip(bam.references, bam.lengths) if l > 0]
            bam_references = {r for r, l in ref_lengths_list}
        if args.reference_fasta:
            logger.info("Validating FASTA...")
            if not _validate_fasta_file(args.reference_fasta, logger):
                logger.error("FATAL: FASTA validation failed.")
                sys.exit(1)
            is_compatible, matching_chroms = _check_fasta_bam_compatibility(
                args.reference_fasta, bam_references, logger)
            if not is_compatible:
                logger.error("FATAL: FASTA and BAM chromosome names are incompatible.")
                logger.error("PolyA detection will fail. Please ensure chromosome names match.")
                sys.exit(1)
        with pysam.AlignmentFile(bam_to_process, "rb", check_sq=False, require_index=False) as bam:
            logger.info("Scanning BAM for aligner properties...")
            has_nh, has_as, max_mapq = False, False, 0
            scan_flags = pysam.FUNMAP | pysam.FQCFAIL | pysam.FDUP | pysam.FSECONDARY | pysam.FSUPPLEMENTARY
            count = 0
            for count, r in enumerate(itertools.islice(
                    (r for r in bam if not (r.flag & scan_flags)), BAM_SCAN_LIMIT), 1):
                if not has_nh and r.has_tag("NH"):
                    has_nh = True
                if not has_as and r.has_tag("AS"):
                    has_as = True
                if r.mapping_quality > max_mapq:
                    max_mapq = r.mapping_quality
                if has_nh and has_as and max_mapq >= 255:
                    break
            logger.info(f"Scanned {count:,} alignments. max MAPQ={max_mapq}, NH={has_nh}, AS={has_as}")
        auto_thr = None
        if has_nh:
            logger.info("NH tag detected. Using NH==1 for unique mapping.")
        else:
            logger.warning("NH tag not found. Using MAPQ threshold.")
            auto_thr = 255 if max_mapq >= 255 else (60 if max_mapq >= 60 else (40 if max_mapq >= 40 else 30))
            logger.info(f"Auto-detected unique_mapping_threshold = {auto_thr}")
        if args.unique_mapping_threshold is None:
            args.unique_mapping_threshold = auto_thr if auto_thr else 255
        args_dict['unique_mapping_threshold'] = args.unique_mapping_threshold
        args_dict['min_polya_length'] = args.min_polya_length
        args_dict['max_alignments_per_mate'] = args.max_alignments_per_mate
        logger.info("--- Runtime Parameters ---")
        for key, value in sorted(args_dict.items()):
            logger.info(f"  {key:<30}: {value}")
        gene_structs = None
        gene_lengths = None
        gene_map_rev = None
        gene_info_array = None
        valid_gene_ids = set()
        splice_junctions = {}
        single_exon_genes = set()
        if args.gtf:
            (gene_structs, gene_lengths, gene_map_rev, gene_info_array,
             valid_gene_ids, splice_junctions, single_exon_genes) = read_and_process_gtf(
                 args.gtf, bam_references, logger, output_dir, prefix)
        fasta_path = args.reference_fasta if args.polyA_tail_detection else None
        chrom_len_map = {r: l for r, l in ref_lengths_list}
        pool_init_args = (gene_structs, gene_map_rev, gene_info_array, log_queue,
                          splice_junctions, single_exon_genes, fasta_path)
        pool_context = mp.Pool(processes=args.num_threads, initializer=init_worker, initargs=pool_init_args)
        logger.info(f"Starting parallel processing with {args.num_threads} threads.")
        raw_read_group_iter = read_group_generator(bam_to_process, args.keep_duplicates, logger)
        keep_sequence = args.polyA_tail_detection
        lightweight_iter = ([LightweightRead(r, keep_sequence) for r in group] for group in raw_read_group_iter)
        results_iter = pool_context.imap_unordered(partial(_process_read_group_worker, args_dict=args_dict),
                                                   lightweight_iter, chunksize=args.chunksize)
        (final_unique_f1r2, final_unique_f2r1, final_unique_bins_f1r2, final_unique_bins_f2r1,
         final_pos_multi_map_dict, final_gene_ambig_list, genomic_uniq_counts, pos_ambig_groups_count,
         total_rpm_denominator, uniq_multigene_discard_count, num_processed_groups, total_read_groups,
         splice_counts_uniq, pas_counts_uniq, gene_ambig_count) = _aggregate_worker_results(
             results_iter, args, valid_gene_ids, chrom_len_map, logger, single_exon_genes)
        pool_context.close()
        pool_context.join()
        pool_context = None
        logger.info("Finished processing all BAM data.")
        if num_processed_groups == 0:
            logger.error("No read groups processed.")
            sys.exit(1)
        gc.collect()
        all_ambiguous_for_gene_em = final_gene_ambig_list if final_gene_ambig_list else []
        multi_f1r2 = defaultdict(list)
        multi_f2r1 = defaultdict(list)
        multi_bins_f1r2 = {}
        multi_bins_f2r1 = {}
        genomic_multi_counts = defaultdict(float)
        splice_counts_multi = defaultdict(float)
        pas_counts_multi = defaultdict(float)
        if args.positional_assignment == 'EM' and final_pos_multi_map_dict:
            (multi_f1r2, multi_f2r1, multi_bins_f1r2, multi_bins_f2r1,
             genomic_multi_gene_single, splice_multi, pas_multi, ambig_from_pos_em) = _resolve_positional_ambiguity(
                final_pos_multi_map_dict, final_unique_f1r2, final_unique_f2r1,
                final_unique_bins_f1r2, final_unique_bins_f2r1, args, chrom_len_map,
                valid_gene_ids, logger, single_exon_genes, args.max_ambiguous_groups)
            for gene_id, count in genomic_multi_gene_single.items():
                genomic_multi_counts[gene_id] += count
            for gene_id, count in splice_multi.items():
                splice_counts_multi[gene_id] += count
            for gene_id, count in pas_multi.items():
                pas_counts_multi[gene_id] += count
            all_ambiguous_for_gene_em.extend(ambig_from_pos_em)
        em_uniq_counts = defaultdict(float)
        em_multi_counts = defaultdict(float)
        em_splice_uniq = defaultdict(float)
        em_splice_multi = defaultdict(float)
        em_pas_uniq = defaultdict(float)
        em_pas_multi = defaultdict(float)
        gene_em_converged = True
        if args.gtf and all_ambiguous_for_gene_em and args.gene_disambiguation == 'EM':
            (em_uniq_counts, em_multi_counts, em_splice_uniq, em_splice_multi,
             em_pas_uniq, em_pas_multi, gene_em_converged) = run_gene_em_algorithm_scientific(
                all_ambiguous_for_gene_em, genomic_uniq_counts, genomic_multi_counts,
                gene_lengths, valid_gene_ids, single_exon_genes, args, logger,
                mean_frag_len=args.mean_frag_len)
        for gene_id, count in em_uniq_counts.items():
            genomic_uniq_counts[gene_id] += count
        for gene_id, count in em_multi_counts.items():
            genomic_multi_counts[gene_id] += count
        for gene_id, count in em_splice_uniq.items():
            splice_counts_uniq[gene_id] += count
        for gene_id, count in em_splice_multi.items():
            splice_counts_multi[gene_id] += count
        for gene_id, count in em_pas_uniq.items():
            pas_counts_uniq[gene_id] += count
        for gene_id, count in em_pas_multi.items():
            pas_counts_multi[gene_id] += count
        total_splice = defaultdict(float)
        total_pas = defaultdict(float)
        for gene_id in valid_gene_ids:
            total_splice[gene_id] = splice_counts_uniq.get(gene_id, 0.0) + splice_counts_multi.get(gene_id, 0.0)
            total_pas[gene_id] = pas_counts_uniq.get(gene_id, 0.0) + pas_counts_multi.get(gene_id, 0.0)
        if args.gtf:
            write_count_files(output_dir, prefix, genomic_uniq_counts, genomic_multi_counts,
                              gene_lengths, valid_gene_ids, args, logger,
                              total_splice, total_pas, args.mean_frag_len)
        scale_factor = 1.0
        if args.normalize_bigwig == 'RPM':
            if total_rpm_denominator > 0:
                scale_factor = 1_000_000.0 / total_rpm_denominator
                logger.info(f"RPM normalization: scale={scale_factor:.6f}, total_mapped_read_groups={total_rpm_denominator:,.0f}")
            else:
                logger.warning("RPM denominator is zero. BigWigs not scaled.")
        _write_final_outputs(prefix, output_dir, ref_lengths_list, args, scale_factor, logger,
                             final_unique_f1r2, final_unique_f2r1, multi_f1r2, multi_f2r1,
                             final_unique_bins_f1r2, final_unique_bins_f2r1, multi_bins_f1r2, multi_bins_f2r1)
        logger.info("--- Run Summary ---")
        logger.info(f"Total read groups (library size for RPM): {total_read_groups:,}")
        logger.info(f"RPM normalization denominator: {total_rpm_denominator:,.0f}")
        if args.positional_assignment == 'EM':
            logger.info(f"Positional ambiguous groups: {pos_ambig_groups_count:,}")
        if args.gtf:
            uniq_genes = len({g for g in genomic_uniq_counts if genomic_uniq_counts[g] > 0})
            multi_genes = len({g for g in genomic_multi_counts if genomic_multi_counts[g] > 0})
            logger.info(f"Genes with genomic-unique counts: {uniq_genes:,}")
            logger.info(f"Genes with genomic-multi counts: {multi_genes:,}")
        run_completed = True
    except KeyboardInterrupt:
        if logger:
            logger.error("Interrupted by user.")
    except SystemExit as e:
        run_completed = (e.code == 0)
    except Exception:
        if logger:
            logger.error("Unhandled error:", exc_info=True)
        else:
            import traceback
            traceback.print_exc()
    finally:
        if pool_context:
            pool_context.terminate()
            pool_context.join()
        if log_listener:
            log_listener.stop()
        if args and is_temp_bam and bam_to_process and os.path.exists(bam_to_process):
            if run_completed and not args.no_cleanup:
                try:
                    os.remove(bam_to_process)
                    if logger:
                        logger.info(f"Cleaned up temporary BAM: {os.path.basename(bam_to_process)}")
                except OSError:
                    pass
        if logger:
            logger.info(f"Total time: {time.time() - start_time:.2f}s")
            logger.info("EMapper completed." if run_completed else "EMapper failed.")
        if not run_completed:
            sys.exit(1)
# =============================================================================
# Entry Point
# =============================================================================
if __name__ == "__main__":
    try:
        mp.set_start_method("spawn")
    except RuntimeError:
        pass
    mp.freeze_support()
    main()
