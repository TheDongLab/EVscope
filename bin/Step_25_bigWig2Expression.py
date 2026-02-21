#!/usr/bin/env python3
"""
bigWig2Expression.py: Calculate CPM and TPM values for genes from BigWig files.

This script computes both CPM (Counts Per Million) and TPM (Transcripts Per Million)
for genes using strand-specific or unstranded BigWig coverage files. It uses an
exon-based model that merges overlapping exons to create non-redundant gene
representations, similar to featureCounts.

IMPORTANT NOTES ON NORMALIZATION:
    - CPM: Signal Per Million (SPM) - normalized by total signal, NOT fragment count
      This is NOT traditional CPM which uses fragment counts.
    - TPM: All TPM values (unstranded, forward, reverse) use the SAME denominator
      (total_spk_unstranded) to ensure comparability.

Output Format (7 columns, similar to STAR ReadsPerGene.out.tab):
    1. gene_id: Gene identifier
    2. unstranded_CPM: CPM from both strands combined
    3. forward_CPM: CPM from forward strand (F1R2 for + genes, F2R1 for - genes)
    4. reverse_CPM: CPM from reverse strand (F2R1 for + genes, F1R2 for - genes)
    5. unstranded_TPM: TPM from both strands combined
    6. forward_TPM: TPM from forward strand
    7. reverse_TPM: TPM from reverse strand

Coordinate Systems:
    - GTF: 1-based closed intervals [start, end] → converted to 0-based internally
    - BED: 0-based half-open intervals [start, end) → used as-is
    - BigWig: 0-based half-open intervals [start, end)

Strand Specificity:
    Library type interpretation (matches RNA direction):
    - forward: F1R2 (Read2) matches RNA strand direction
    - reverse: F2R1 (Read1) matches RNA strand direction (Illumina dUTP protocol)

Version: 1.0.0
"""

import os
import time
import sys
import logging
import argparse
import pyBigWig
import re
import tempfile
import gzip
import math  # Added for nan checking
from typing import Dict, List, Tuple, Optional, TextIO
from collections import defaultdict

__version__ = "1.0.0"
global_logger: Optional[logging.Logger] = None

###############################################################################
# Utility Functions
###############################################################################

def setup_logging(log_path: str) -> logging.Logger:
    """
    Configure logging to both file and console.

    Args:
        log_path: Path to log file

    Returns:
        Configured logger instance
    """
    logger = logging.getLogger('bigWig2Expression')
    logger.handlers.clear()
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter(
        '[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # File handler
    file_handler = logging.FileHandler(log_path, mode='w')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    return logger

def parse_arguments() -> argparse.Namespace:
    """
    Parse and validate command-line arguments.

    Returns:
        Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description=f'bigWig2Expression v{__version__}: Calculate CPM and TPM from BigWig files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
DESCRIPTION:
  This tool calculates both CPM and TPM values from BigWig coverage files,
  similar to STAR's ReadsPerGene.out.tab output format. It uses merged exons
  to avoid double-counting overlapping regions.

INPUT FILE FORMATS:

  BigWig Files (strand-specific or unstranded):
    - F1R2: Read2/Forward strand coverage
    - F2R1: Read1/Reverse strand coverage
    - Combined: Unstranded coverage (both strands merged)

  Annotation Files:

    GTF format:
      - Must contain 'exon' features with gene_id attribute
      - Coordinates are 1-based (will be converted to 0-based internally)
      - Example: chr1  HAVANA  exon  100  200  .  +  .  gene_id "ENSG00000001"

    BED format (BED3 or BED6 only):
      - BED3 (3 columns): chr, start, end
        * Treated as UNSTRANDED
        * gene_id format: chr_start:end (e.g., chr1_100000:200000)

      - BED6 (6 columns): chr, start, end, name, score, strand
        * Strand column (field 6) values:
          + or - : Treated as STRANDED
          . or other: Treated as UNSTRANDED
        * gene_id format:
          - If name is valid (not "."): chr_start:end:name
          - If name is ".": chr_start:end
        * Coordinates are 0-based

      - Example BED3:
        chr1    100000  200000

      - Example BED6:
        chr1    100000  200000  GENE1   0   +
        chr2    300000  400000  .       0   .

OUTPUT FORMAT (7 columns):
  gene_id  unstranded_CPM  forward_CPM  reverse_CPM  unstranded_TPM  forward_TPM  reverse_TPM

STRAND CONVENTIONS:
  - forward library: F1R2 (Read2) matches RNA strand
  - reverse library: F2R1 (Read1) matches RNA strand (Illumina dUTP)

  For a gene on + strand:
    - forward_CPM/TPM: from F1R2
    - reverse_CPM/TPM: from F2R1

  For a gene on - strand:
    - forward_CPM/TPM: from F2R1
    - reverse_CPM/TPM: from F1R2

NORMALIZATION NOTES:
  - CPM is actually "Signal Per Million" (total signal, not fragment count)
  - All TPM values use the same denominator for comparability
  - For unstranded data, forward/reverse CPM/TPM are set to 0

EXAMPLES:

  1. Stranded RNA-seq data:
      %(prog)s \\
        --input_F1R2_bw sample.F1R2.bw \\
        --input_F2R1_bw sample.F2R1.bw \\
        --gtf genes.gtf \\
        --output sample_expression.tsv

  2. Unstranded RNA-seq data:
      %(prog)s \\
        --input_combined_bw sample.bw \\
        --gtf genes.gtf \\
        --output sample_expression.tsv

  3. Using BED6 annotation:
      %(prog)s \\
        --input_F1R2_bw sample.F1R2.bw \\
        --input_F2R1_bw sample.F2R1.bw \\
        --bed genes.bed \\
        --output sample_expression.tsv

NOTES:
  - Overlapping exons are automatically merged to prevent double-counting
  - Genes spanning multiple chromosomes are processed separately per chromosome
  - Uses bw.values() for accurate signal extraction (not stats with binning)
  - All coordinates are internally converted to 0-based half-open intervals
        """
    )

    # Input files
    input_group = parser.add_argument_group('Input BigWig Files')
    input_group.add_argument(
        '--input_F1R2_bw',
        metavar='FILE',
        help='BigWig file for F1R2 strand (Read2/forward). Required for stranded data.'
    )
    input_group.add_argument(
        '--input_F2R1_bw',
        metavar='FILE',
        help='BigWig file for F2R1 strand (Read1/reverse). Required for stranded data.'
    )
    input_group.add_argument(
        '--input_combined_bw',
        metavar='FILE',
        help='BigWig file for combined unstranded coverage. Use for unstranded libraries.'
    )

    # Annotation files
    annot_group = parser.add_argument_group('Annotation Files')
    annot_group.add_argument(
        '--gtf',
        metavar='FILE',
        help='GTF annotation file (.gtf or .gtf.gz). Uses exon features only.'
    )
    annot_group.add_argument(
        '--bed',
        metavar='FILE',
        help='BED annotation file (.bed or .bed.gz). Supports BED3 or BED6 format only. '
              'BED3 is treated as unstranded. BED6 with strand +/- is stranded, '
              'with strand . or other is unstranded.'
    )

    # Output
    output_group = parser.add_argument_group('Output Options')
    output_group.add_argument(
        '--output',
        required=True,
        metavar='FILE',
        help='Output TSV file with 7 columns (gene_id, CPMs, TPMs)'
    )

    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')

    args = parser.parse_args()

    # Validate BigWig file combinations
    if args.input_combined_bw:
        if args.input_F1R2_bw or args.input_F2R1_bw:
            parser.error(
                'Cannot use --input_combined_bw with --input_F1R2_bw/--input_F2R1_bw. '
                'Choose one mode: combined (unstranded) OR F1R2+F2R1 (stranded).'
            )
    elif not (args.input_F1R2_bw and args.input_F2R1_bw):
        parser.error(
            'Must provide either --input_combined_bw OR both --input_F1R2_bw and --input_F2R1_bw'
        )

    # Validate annotation file
    if not args.gtf and not args.bed:
        parser.error('Must provide either --gtf or --bed annotation file')
    if args.gtf and args.bed:
        parser.error('Cannot use both --gtf and --bed simultaneously')

    return args

def validate_file(file_path: str, logger: logging.Logger, file_type: str) -> None:
    """
    Validate that file exists and is not empty.

    Args:
        file_path: Path to file
        logger: Logger instance
        file_type: File type description for error messages

    Raises:
        FileNotFoundError: If file does not exist
        ValueError: If file is empty
    """
    if not os.path.exists(file_path):
        logger.error(f'{file_type} file not found: {file_path}')
        raise FileNotFoundError(f'{file_type} file not found: {file_path}')

    if os.path.getsize(file_path) == 0:
        logger.error(f'{file_type} file is empty: {file_path}')
        raise ValueError(f'{file_type} file is empty: {file_path}')

def open_file(file_path: str) -> TextIO:
    """
    Open file with automatic gzip decompression if needed.

    Args:
        file_path: Path to file (can be .gz)

    Returns:
        File handle
    """
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    return open(file_path, 'r')

def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Merge overlapping intervals to create non-redundant regions.

    This prevents double-counting of overlapping exons for the same gene.

    Args:
        intervals: List of (start, end) tuples in 0-based coordinates

    Returns:
        List of merged non-overlapping intervals, sorted by start position

    Example:
        >>> merge_intervals([(100, 200), (150, 250), (300, 400)])
        [(100, 250), (300, 400)]
    """
    if not intervals:
        return []

    # Sort by start position
    intervals.sort(key=lambda x: x[0])

    merged = [intervals[0]]
    for current_start, current_end in intervals[1:]:
        last_start, last_end = merged[-1]

        if current_start <= last_end:
            # Overlapping or adjacent intervals, merge them
            merged[-1] = (last_start, max(last_end, current_end))
        else:
            # No overlap, add as new interval
            merged.append((current_start, current_end))

    return merged

def get_signal_sum(bw: pyBigWig.pyBigWig, chrom: str, start: int, end: int) -> float:
    """
    Extract signal sum from BigWig using values() method for accuracy.

    This method is more accurate than stats(nBins=1) which uses binning
    and may drop partial values, especially for large regions.

    FIX: Handles nan values returned by pyBigWig for regions with no coverage.

    Args:
        bw: BigWig file handle
        chrom: Chromosome name
        start: Start position (0-based)
        end: End position (0-based, exclusive)

    Returns:
        Sum of signal values in the region
    """
    try:
        values = bw.values(chrom, start, end)
        # Sum all non-None and non-nan values
        # pyBigWig returns nan for uncovered regions in some versions/files
        return sum(v for v in values if v is not None and not math.isnan(v))
    except RuntimeError:
        # Chromosome not in BigWig or other error
        return 0.0

###############################################################################
# Annotation Parsing Functions
###############################################################################

def read_and_merge_exons_from_gtf(
    gtf_path: str,
    logger: logging.Logger
) -> Dict[str, List[Tuple[str, str, List[Tuple[int, int]]]]]:
    """
    Read exon features from GTF file and merge overlapping exons per gene.

    GTF coordinates are 1-based closed intervals [start, end].
    They are converted to 0-based half-open intervals [start, end) internally.

    Args:
        gtf_path: Path to GTF file
        logger: Logger instance

    Returns:
        Dictionary: {gene_id: [(chrom, strand, merged_exons), ...]}
        where merged_exons is [(start, end), ...] in 0-based coordinates

    Note:
        Only processes lines with feature type 'exon'.
        Requires gene_id attribute in column 9.
        Genes spanning multiple chromosomes are kept separate per chromosome.
    """
    logger.info(f'Reading exons from GTF: {os.path.basename(gtf_path)}')

    gene_exon_groups = defaultdict(list)
    line_count = 0
    exon_count = 0

    with open_file(gtf_path) as f:
        for line in f:
            line_count += 1

            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')

            if len(fields) < 9:
                continue

            # Only process exon features
            if fields[2] != 'exon':
                continue

            try:
                chrom = fields[0]
                start = int(fields[3]) - 1  # Convert GTF 1-based to 0-based
                end = int(fields[4])         # GTF end is inclusive, now half-open [start, end)
                strand = fields[6]
                attributes = fields[8]

                # Validate coordinates
                if start < 0 or end <= start:
                    logger.warning(f'Line {line_count}: Invalid coordinates {chrom}:{start}-{end}')
                    continue

                # Validate strand
                if strand not in ['+', '-']:
                    logger.warning(f'Line {line_count}: Invalid strand "{strand}"')
                    continue

                # Extract gene_id from attributes
                gene_id_match = re.search(r'gene_id\s+"([^"]+)"', attributes)
                if not gene_id_match:
                    gene_id_match = re.search(r'gene_id\s+(\S+)', attributes)

                if not gene_id_match:
                    logger.warning(f'Line {line_count}: No gene_id found')
                    continue

                gene_id = gene_id_match.group(1)

                # Store exon grouped by (gene_id, chrom, strand)
                # This ensures genes on different chromosomes are kept separate
                gene_exon_groups[(gene_id, chrom, strand)].append((start, end))
                exon_count += 1

            except (ValueError, IndexError) as e:
                logger.warning(f'Line {line_count}: Parsing error - {e}')
                continue

    logger.info(f'Parsed {exon_count} exons from {line_count} GTF lines')

    # Merge overlapping exons for each gene
    processed_genes = defaultdict(list)
    for (gene_id, chrom, strand), exons in gene_exon_groups.items():
        merged = merge_intervals(exons)
        processed_genes[gene_id].append((chrom, strand, merged))

    if not processed_genes:
        logger.error('No valid exon entries found in GTF')
        raise ValueError('No valid exon entries in GTF file')

    # Check for genes spanning multiple chromosomes
    multi_chrom_genes = sum(1 for parts in processed_genes.values() if len(parts) > 1)
    if multi_chrom_genes > 0:
        logger.warning(
            f'{multi_chrom_genes} genes span multiple chromosomes/strands. '
            'Each chromosome will be processed separately.'
        )

    logger.info(f'Processed {len(processed_genes)} genes with merged exons')

    return processed_genes

def read_and_merge_exons_from_bed(
    bed_path: str,
    logger: logging.Logger
) -> Dict[str, List[Tuple[str, str, List[Tuple[int, int]]]]]:
    """
    Read exon features from BED file and merge overlapping exons.

    Supports BED3 and BED6 formats:
    - BED3 (chr, start, end): Treated as UNSTRANDED, strand='.'
      gene_id format: chr_start:end

    - BED6 (chr, start, end, name, score, strand):
      * If strand is + or -: STRANDED
      * If strand is . or other: UNSTRANDED (strand='.')
      gene_id format: chr_start:end:name (if name != '.'), otherwise chr_start:end

    BED coordinates are 0-based half-open [start, end), used as-is.

    Args:
        bed_path: Path to BED file
        logger: Logger instance

    Returns:
        Dictionary: {gene_id: [(chrom, strand, merged_exons), ...]}

    Note:
        Each BED line is treated as an exon. Lines with the same gene_id
        (and chrom, strand) will be merged.
        Genes spanning multiple chromosomes are kept separate per chromosome.
    """
    logger.info(f'Reading exons from BED: {os.path.basename(bed_path)}')

    gene_exon_groups = defaultdict(list)
    line_count = 0
    exon_count = 0
    format_stats = {'BED3': 0, 'BED6_stranded': 0, 'BED6_unstranded': 0}

    with open_file(bed_path) as f:
        for line in f:
            line_count += 1

            # Skip comments and track lines
            if line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                continue

            if not line.strip():
                continue

            fields = line.strip().split('\t')
            num_fields = len(fields)

            # Only support BED3 or BED6
            if num_fields not in [3, 6]:
                logger.warning(
                    f'Line {line_count}: BED file must be BED3 (3 columns) or BED6 (6 columns), '
                    f'found {num_fields} columns. Skipping.'
                )
                continue

            try:
                chrom = fields[0]
                start = int(fields[1])  # BED is 0-based
                end = int(fields[2])

                # Validate coordinates
                if start < 0 or end <= start:
                    logger.warning(f'Line {line_count}: Invalid coordinates {chrom}:{start}-{end}')
                    continue

                if num_fields == 3:
                    # BED3 format: chr_start:end
                    gene_id = f'{chrom}_{start}:{end}'
                    strand = '.'  # Unstranded
                    format_stats['BED3'] += 1

                else:  # num_fields == 6
                    # BED6 format
                    name = fields[3].strip()
                    strand_field = fields[5].strip()

                    # Determine strand
                    if strand_field in ['+', '-']:
                        strand = strand_field
                        format_stats['BED6_stranded'] += 1
                    else:
                        # Strand is '.', 'N/A', or other non-standard value
                        strand = '.'
                        format_stats['BED6_unstranded'] += 1

                    # Construct gene_id
                    if name and name != '.':
                        gene_id = f'{chrom}_{start}:{end}:{name}'
                    else:
                        gene_id = f'{chrom}_{start}:{end}'

                # Store exon grouped by (gene_id, chrom, strand)
                gene_exon_groups[(gene_id, chrom, strand)].append((start, end))
                exon_count += 1

            except (ValueError, IndexError) as e:
                logger.warning(f'Line {line_count}: Parsing error - {e}')
                continue

    # Log format statistics
    logger.info(f'Parsed {exon_count} exons from {line_count} BED lines')
    logger.info(f'Format breakdown: {format_stats}')

    # Merge overlapping exons for each gene
    processed_genes = defaultdict(list)
    for (gene_id, chrom, strand), exons in gene_exon_groups.items():
        merged = merge_intervals(exons)
        processed_genes[gene_id].append((chrom, strand, merged))

    if not processed_genes:
        logger.error('No valid entries found in BED file')
        raise ValueError('No valid entries in BED file')

    # Log strand statistics
    stranded_genes = sum(
        1 for parts in processed_genes.values()
        for _, strand, _ in parts if strand in ['+', '-']
    )
    unstranded_genes = len(processed_genes) - stranded_genes

    logger.info(
        f'Processed {len(processed_genes)} genes with merged exons '
        f'(stranded: {stranded_genes}, unstranded: {unstranded_genes})'
    )

    return processed_genes

###############################################################################
# Expression Calculation Functions
###############################################################################

def compute_total_signal(
    bw_f1r2: Optional[pyBigWig.pyBigWig],
    bw_f2r1: Optional[pyBigWig.pyBigWig],
    bw_combined: Optional[pyBigWig.pyBigWig],
    logger: logging.Logger
) -> Tuple[float, float, float]:
    """
    Compute total signal from BigWig files for normalization.

    Args:
        bw_f1r2: F1R2 BigWig handle (or None)
        bw_f2r1: F2R1 BigWig handle (or None)
        bw_combined: Combined BigWig handle (or None)
        logger: Logger instance

    Returns:
        Tuple of (total_unstranded, total_f1r2, total_f2r1)

    Note:
        - For unstranded mode:
          total_unstranded = signal from combined file
          total_f1r2 = total_unstranded (for CPM calculation consistency)
          total_f2r1 = total_unstranded (for CPM calculation consistency)
        - For stranded mode:
          total_unstranded = total_f1r2 + total_f2r1
    """
    logger.info('Computing total signal from BigWig files')

    if bw_combined:
        # Unstranded mode: single BigWig file
        total_signal = 0.0
        for chrom, length in bw_combined.chroms().items():
            if length > 0:
                signal = get_signal_sum(bw_combined, chrom, 0, length)
                total_signal += signal

        # In unstranded mode, set all three totals to the same value
        # This ensures correct CPM calculation for forward/reverse
        total_unstranded = total_signal
        total_f1r2 = total_signal
        total_f2r1 = total_signal

    else:
        # Stranded mode: two BigWig files
        total_f1r2 = 0.0
        total_f2r1 = 0.0

        chroms = set(bw_f1r2.chroms().keys()) | set(bw_f2r1.chroms().keys())

        for chrom in chroms:
            length = max(
                bw_f1r2.chroms().get(chrom, 0),
                bw_f2r1.chroms().get(chrom, 0)
            )

            if length == 0:
                continue

            f1r2_sig = get_signal_sum(bw_f1r2, chrom, 0, length)
            f2r1_sig = get_signal_sum(bw_f2r1, chrom, 0, length)

            total_f1r2 += f1r2_sig
            total_f2r1 += f2r1_sig

        total_unstranded = total_f1r2 + total_f2r1

    logger.info(
        f'Total signal - Unstranded: {total_unstranded:.2f}, '
        f'F1R2: {total_f1r2:.2f}, F2R1: {total_f2r1:.2f}'
    )

    if total_unstranded <= 0:
        logger.error('Total signal is zero, negative, or nan')
        raise ValueError('Total signal must be positive')

    return total_unstranded, total_f1r2, total_f2r1

def calculate_gene_expression(
    processed_genes: Dict[str, List[Tuple[str, str, List[Tuple[int, int]]]]],
    bw_f1r2: Optional[pyBigWig.pyBigWig],
    bw_f2r1: Optional[pyBigWig.pyBigWig],
    bw_combined: Optional[pyBigWig.pyBigWig],
    total_unstranded: float,
    total_f1r2: float,
    total_f2r1: float,
    logger: logging.Logger
) -> Dict[str, Tuple[float, float, float, float, float, float]]:
    """
    Calculate CPM and TPM values for all genes using a two-pass algorithm.

    Pass 1: Calculate Signal Per Kilobase (SPK) for each gene
    Pass 2: Normalize SPK to TPM using UNIFIED denominator (total_spk_unstranded)

    CRITICAL CORRECTIONS:
    1. CPM uses strand-specific totals (total_f1r2/total_f2r1) as denominators
    2. All TPM values use the SAME denominator (total_spk_unstranded) for comparability
    3. Unstranded genes: forward=f1r2, reverse=f2r1 (NO double counting)
    4. Unstranded mode (bw_combined): forward/reverse signals set to 0
    5. Uses get_signal_sum() for accurate signal extraction

    Args:
        processed_genes: Gene structure from read_and_merge_exons_*
        bw_f1r2: F1R2 BigWig handle
        bw_f2r1: F2R1 BigWig handle
        bw_combined: Combined BigWig handle
        total_unstranded: Total signal from both strands
        total_f1r2: Total F1R2 signal
        total_f2r1: Total F2R1 signal
        logger: Logger instance

    Returns:
        Dictionary: {gene_id: (unstranded_CPM, forward_CPM, reverse_CPM,
                               unstranded_TPM, forward_TPM, reverse_TPM)}

    Strand Interpretation for STRANDED data:
        For genes on + strand:
            - forward = F1R2 signal
            - reverse = F2R1 signal
        For genes on - strand:
            - forward = F2R1 signal (reverse library: F2R1 matches RNA)
            - reverse = F1R2 signal
        For unstranded genes (strand='.'):
            - unstranded = F1R2 + F2R1
            - forward = F1R2 (NOT doubled)
            - reverse = F2R1 (NOT doubled)

    For UNSTRANDED data (bw_combined):
        - unstranded = signal from combined file
        - forward = 0 (cannot determine strand)
        - reverse = 0 (cannot determine strand)
    """
    logger.info('Calculating expression values (CPM and TPM)')
    logger.info('Pass 1/2: Calculating Signal Per Kilobase (SPK)')

    gene_data = {}

    total_spk_unstranded = 0.0

    genes_processed = 0
    genes_with_signal = 0
    genes_skipped = 0

    for gene_id, parts in processed_genes.items():
        genes_processed += 1

        # Accumulate signal and length across all exons for this gene
        signal_unstranded = 0.0
        signal_f1r2 = 0.0  # Raw F1R2 signal
        signal_f2r1 = 0.0  # Raw F2R1 signal
        total_length = 0

        # Track strand consistency for this gene
        gene_strands = set()

        for chrom, strand, merged_exons in parts:
            gene_strands.add(strand)

            for start, end in merged_exons:
                exon_length = end - start
                total_length += exon_length

                if bw_combined:
                    # Unstranded mode: single BigWig
                    sig = get_signal_sum(bw_combined, chrom, start, end)
                    signal_unstranded += sig
                    # In unstranded mode, we cannot determine strand-specific signal
                    # So F1R2 and F2R1 remain 0

                else:
                    # Stranded mode: two BigWig files
                    f1r2_sig = get_signal_sum(bw_f1r2, chrom, start, end)
                    f2r1_sig = get_signal_sum(bw_f2r1, chrom, start, end)

                    signal_unstranded += f1r2_sig + f2r1_sig
                    signal_f1r2 += f1r2_sig
                    signal_f2r1 += f2r1_sig

        if total_length == 0:
            logger.warning(f'Gene {gene_id} has zero length, skipping')
            genes_skipped += 1
            continue

        if signal_unstranded > 0:
            genes_with_signal += 1

        # Calculate SPK (Signal Per Kilobase)
        length_kb = total_length / 1000.0
        spk_unstranded = signal_unstranded / length_kb

        # Now determine forward/reverse signals based on strand
        # CRITICAL: This must not double-count for unstranded genes

        if bw_combined:
            # Unstranded data mode: cannot determine strand-specific signal
            signal_forward = 0.0
            signal_reverse = 0.0
            spk_forward = 0.0
            spk_reverse = 0.0

        else:
            # Stranded data mode
            # Check if gene has consistent strand annotation
            if len(gene_strands) > 1:
                logger.warning(
                    f'Gene {gene_id} has inconsistent strand annotations: {gene_strands}. '
                    'Treating as unstranded.'
                )
                # Treat as unstranded
                signal_forward = signal_f1r2
                signal_reverse = signal_f2r1
            elif '.' in gene_strands or len(gene_strands & {'+', '-'}) == 0:
                # Unstranded gene annotation
                # CRITICAL FIX: Do NOT double count
                # forward gets F1R2, reverse gets F2R1
                signal_forward = signal_f1r2
                signal_reverse = signal_f2r1
            elif '+' in gene_strands:
                # Positive strand gene
                # forward = F1R2 (Read2 matches RNA)
                # reverse = F2R1
                signal_forward = signal_f1r2
                signal_reverse = signal_f2r1
            else:  # '-' in gene_strands
                # Negative strand gene
                # forward = F2R1 (for reverse library, Read1 matches RNA)
                # reverse = F1R2
                signal_forward = signal_f2r1
                signal_reverse = signal_f1r2

            # Calculate strand-specific SPK
            spk_forward = signal_forward / length_kb
            spk_reverse = signal_reverse / length_kb

        # Calculate CPM (Counts Per Million)
        # CRITICAL: Use strand-specific totals as denominators
        # Check against >0 and ensure not nan
        cpm_unstranded = (signal_unstranded / total_unstranded) * 1e6 if total_unstranded > 0 else 0.0

        if bw_combined:
            # Unstranded mode: cannot calculate strand-specific CPM
            cpm_forward = 0.0
            cpm_reverse = 0.0
        else:
            # Stranded mode: use appropriate denominators
            cpm_forward = (signal_forward / total_f1r2) * 1e6 if total_f1r2 > 0 else 0.0
            cpm_reverse = (signal_reverse / total_f2r1) * 1e6 if total_f2r1 > 0 else 0.0

        # Store intermediate data
        gene_data[gene_id] = {
            'spk_unstranded': spk_unstranded,
            'spk_forward': spk_forward,
            'spk_reverse': spk_reverse,
            'cpm_unstranded': cpm_unstranded,
            'cpm_forward': cpm_forward,
            'cpm_reverse': cpm_reverse
        }

        # Accumulate SPK for TPM calculation
        # Only accumulate from genes with valid length
        total_spk_unstranded += spk_unstranded

    logger.info(
        f'Processed {genes_processed} genes: '
        f'{genes_with_signal} with signal > 0, '
        f'{genes_skipped} skipped (zero length)'
    )
    logger.info(f'Total SPK (unstranded) for TPM normalization: {total_spk_unstranded:.2f}')

    # Pass 2: Calculate TPM from SPK
    # CRITICAL: All TPM values use the SAME denominator (total_spk_unstranded)
    logger.info('Pass 2/2: Calculating TPM values using unified denominator')

    gene_expression = {}

    for gene_id, data in gene_data.items():
        # TPM = (SPK / total_spk_unstranded) × 1,000,000
        # CRITICAL: All three TPM values use the SAME denominator
        if total_spk_unstranded > 0:
            tpm_unstranded = (data['spk_unstranded'] / total_spk_unstranded) * 1e6
            tpm_forward = (data['spk_forward'] / total_spk_unstranded) * 1e6
            tpm_reverse = (data['spk_reverse'] / total_spk_unstranded) * 1e6
        else:
            tpm_unstranded = 0.0
            tpm_forward = 0.0
            tpm_reverse = 0.0

        gene_expression[gene_id] = (
            data['cpm_unstranded'],
            data['cpm_forward'],
            data['cpm_reverse'],
            tpm_unstranded,
            tpm_forward,
            tpm_reverse
        )

    # Verify TPM sums (for debugging)
    total_tpm_check = sum(vals[3] for vals in gene_expression.values())
    logger.info(f'TPM sum verification (should be ~1e6): {total_tpm_check:.2f}')

    return gene_expression

###############################################################################
# Output Functions
###############################################################################

def write_output(
    gene_expression: Dict[str, Tuple[float, float, float, float, float, float]],
    output_path: str,
    logger: logging.Logger
) -> None:
    """
    Write gene expression values to TSV file (7 columns).

    Output format (similar to STAR ReadsPerGene.out.tab):
        gene_id  unstranded_CPM  forward_CPM  reverse_CPM
                 unstranded_TPM  forward_TPM  reverse_TPM

    Args:
        gene_expression: Dictionary mapping gene_id to 6-tuple of values
        output_path: Path to output TSV file
        logger: Logger instance
    """
    logger.info(f'Writing expression values to {os.path.basename(output_path)}')

    try:
        output_dir = os.path.dirname(output_path) or '.'

        with tempfile.NamedTemporaryFile(
            mode='w',
            delete=False,
            dir=output_dir,
            suffix='.tmp'
        ) as temp_file:

            # Write header
            temp_file.write(
                'gene_id\t'
                'unstranded_CPM\tforward_CPM\treverse_CPM\t'
                'unstranded_TPM\tforward_TPM\treverse_TPM\n'
            )

            # Write data (sorted by gene_id)
            for gene_id in sorted(gene_expression.keys()):
                values = gene_expression[gene_id]
                temp_file.write(
                    f'{gene_id}\t'
                    f'{values[0]:.6f}\t{values[1]:.6f}\t{values[2]:.6f}\t'
                    f'{values[3]:.6f}\t{values[4]:.6f}\t{values[5]:.6f}\n'
                )

        # Atomic file replacement
        os.replace(temp_file.name, output_path)

        logger.info(f'Successfully wrote {len(gene_expression)} genes')

    except IOError as e:
        logger.error(f'Failed to write output: {e}')
        raise

###############################################################################
# Main Function
###############################################################################

def main():
    """Main execution function."""
    global global_logger

    try:
        # Parse arguments
        args = parse_arguments()

        # Setup logging
        output_directory = os.path.dirname(os.path.abspath(args.output))
        os.makedirs(output_directory, exist_ok=True)

        log_filename = os.path.splitext(os.path.basename(args.output))[0] + '.log'
        log_path = os.path.join(output_directory, log_filename)

        global_logger = setup_logging(log_path)
        logger = global_logger

        # Log run information
        logger.info(f'bigWig2Expression v{__version__}')
        logger.info(f'Command: {" ".join(sys.argv)}')
        logger.info(f'Output file: {args.output}')
        logger.info(f'Log file: {log_path}')

        start_time = time.time()

        # Validate input files
        logger.info('Validating input files...')

        if args.input_combined_bw:
            validate_file(args.input_combined_bw, logger, 'Combined BigWig')
            logger.info('Mode: Unstranded (single BigWig file)')
        else:
            validate_file(args.input_F1R2_bw, logger, 'F1R2 BigWig')
            validate_file(args.input_F2R1_bw, logger, 'F2R1 BigWig')
            logger.info('Mode: Stranded (F1R2 + F2R1 BigWig files)')

        # Read and process annotation
        if args.gtf:
            validate_file(args.gtf, logger, 'GTF')
            processed_genes = read_and_merge_exons_from_gtf(args.gtf, logger)
        else:
            validate_file(args.bed, logger, 'BED')
            processed_genes = read_and_merge_exons_from_bed(args.bed, logger)

        # Open BigWig files
        logger.info('Opening BigWig files...')
        bw_f1r2, bw_f2r1, bw_combined = None, None, None

        try:
            if args.input_combined_bw:
                bw_combined = pyBigWig.open(args.input_combined_bw)
            else:
                bw_f1r2 = pyBigWig.open(args.input_F1R2_bw)
                bw_f2r1 = pyBigWig.open(args.input_F2R1_bw)

            # Compute total signal for normalization
            total_unstranded, total_f1r2, total_f2r1 = compute_total_signal(
                bw_f1r2, bw_f2r1, bw_combined, logger
            )

            # Calculate expression values
            gene_expression = calculate_gene_expression(
                processed_genes,
                bw_f1r2,
                bw_f2r1,
                bw_combined,
                total_unstranded,
                total_f1r2,
                total_f2r1,
                logger
            )

            # Write output
            write_output(gene_expression, args.output, logger)

        finally:
            # Close BigWig files
            if bw_f1r2:
                bw_f1r2.close()
            if bw_f2r1:
                bw_f2r1.close()
            if bw_combined:
                bw_combined.close()

        # Log completion
        elapsed_time = time.time() - start_time
        logger.info(f'Total runtime: {elapsed_time:.2f} seconds')
        logger.info('Analysis completed successfully')

    except Exception as e:
        if global_logger:
            global_logger.exception(f'Fatal error: {e}')
        else:
            print(f'Fatal error: {e}', file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
