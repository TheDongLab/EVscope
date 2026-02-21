#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# Usage example:
#   python Step_22_run_RNA_deconvolution_ARIC.py \
#       --input_expr_file sample_expr.tsv.TPM \
#       --ref_expr_file references/deconvolution_HG38/LM22.csv \
#       --output_dir ./results \
#       --male_specific_tissue "Prostate" "Testis" \
#       --female_specific_tissue "Breast" "Ovary" "Uterus" \
#       --sex Male

Generalized RNA deconvolution via ARIC.

Key features:
- Input file (--input_expr_file) for a single sample (bulk, pseudobulk, etc.). Must be 2 columns: [GeneID, SampleData].
- **Duplicate GeneIDs in the input file will cause an error.**
- Reference file (--ref_expr_file) must be multi-column with a header.
- Auto-detects delimiter (tab/comma) for both files.
- Auto-detects header for --input_expr_file.
- Isoform versions (e.g., .1, .2) are stripped *before* ID type detection.
- Auto-detects gene ID type (Ensembl vs Gene Symbol) based on >10000 count threshold.
- Robust sex-inference finds Symbols (XIST) or Ensembl IDs (ENSG00000229807).
- Auto-skips sex filtering for brain-related references.
- User-definable lists for sex-specific tissues via CLI (case-insensitive matching).
- 'Unknown' or 'Mixed' sex will filter *all* user-defined sex-specific tissues.
- Input matrix and Reference matrix are strictly aligned:
  - The Reference gene list is used as the master list.
  - The Input matrix is re-indexed to match the Reference gene list; missing genes are filled with 0.
  - Both matrices are sorted to have the identical gene order.
- Outputs:
    * Deconvolution result TSV (tab-separated), rounded to 4 decimals, 0s removed.
    * Sex annotation TSV (if --sex Auto is used).
    * Publication-ready barplot (PNG/PDF) with vector fonts, dynamic height, and sorted high-to-low.
    * A clean, user-facing run log written to the output directory.
- Intermediate files are auto-deleted.
"""

import argparse
import csv
import io
import logging
import math
import os
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Set, Union

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # Use Agg backend for headless environments
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from matplotlib import font_manager
mpl.rcParams['pdf.fonttype'] = 42 # Embed fonts in PDF for vector output
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

try:
    from ARIC import ARIC  # ARIC deconvolution function
except ImportError:
    sys.stderr.write("Error: ARIC library not found. Please ensure 'ARIC.py' is in your Python path or the same directory.\n")
    sys.exit(1)


# --------------------------- Logging ---------------------------------

def init_logger(outdir: Union[str, Path], input_expr_path: str) -> logging.Logger:
    """Initializes the logger to output to both console and file."""
    out_path = Path(outdir)
    out_path.mkdir(parents=True, exist_ok=True)

    # Get base name for log file
    base_expr_name = Path(input_expr_path).name.split('.')[0]
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_path = out_path / f"Deconvolution_{timestamp}_{base_expr_name}.log"

    logger = logging.getLogger("ARICDeconv")
    logger.setLevel(logging.INFO)

    if logger.hasHandlers():
        logger.handlers.clear()

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    fh = logging.FileHandler(log_path, encoding="utf-8")
    fh.setLevel(logging.INFO)

    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    ch.setFormatter(fmt)
    fh.setFormatter(fmt)

    logger.addHandler(ch)
    logger.addHandler(fh)
    logger.info(f"Log file: {log_path.name}")
    return logger


# --------------------------- IO helpers ---------------------------------

_DELIM_GUESSERS = [",", "\t", ";", "|"]

def _sniff_delimiter(sample: str) -> str:
    """Tries to auto-detect delimiter via csv.Sniffer, falling back to a count-based method."""
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t;|")
        return dialect.delimiter
    except Exception:
        line1 = sample.splitlines()[0]
        best = max(_DELIM_GUESSERS, key=lambda d: line1.count(d))
        return best

def read_single_sample_expr_file(path: Union[str, Path], logger: logging.Logger) -> pd.DataFrame:
    """
    Reads a single-sample expression matrix (2 columns: GeneID, SampleData).
    Auto-detects delimiter and header.
    If no header, sample name is inferred from filename.
    """
    f_path = Path(path)
    if not f_path.exists():
        logger.error(f"Input file not found: {path}")
        raise FileNotFoundError(f"Input file not found: {path}")

    with open(f_path, "r", encoding="utf-8", errors="ignore") as f:
        buf = f.read(4096)
        if not buf:
            raise ValueError(f"Empty input file: {path}")
        delim = _sniff_delimiter(buf)

    logger.info(f"Detected delimiter for {f_path.name}: {delim!r}")

    # Read CSV, warning on (and skipping) bad lines
    df_try = pd.read_csv(f_path, sep=delim, engine="python", dtype=str, on_bad_lines='warn')

    if df_try.empty:
        raise ValueError(f"Input file is empty or could not be parsed: {path}")

    # Enforce 2-column rule for single sample
    if df_try.shape[1] > 2:
        logger.error(f"Input file {f_path.name} has {df_try.shape[1]} columns.")
        raise ValueError("Input expression file (--input_expr_file) must contain exactly 2 columns: [GeneID, SampleData].")
    elif df_try.shape[1] < 2:
        raise ValueError("Input file must contain at least 2 columns.")

    # FORCE: treat input as NO HEADER
    logger.info("Forced no-header mode: using filename as sample ID.")

    df = pd.read_csv(f_path, sep=delim, engine="python", header=None, dtype=str, on_bad_lines='warn')

    if df.shape[1] != 2:
        raise ValueError(f"Input file has {df.shape[1]} columns. Must be exactly 2: [GeneID, Expression].")

    filename_stem = f_path.stem
    df.columns = ["GeneID", filename_stem]
    logger.info(f"Assigned sample name from filename: '{filename_stem}'")

    df = df.dropna(how="all")
    if df.shape[1] != 2:
        raise ValueError("After cleaning, matrix does not have 2 columns.")
    return df

def read_reference_matrix(path: Union[str, Path], logger: logging.Logger) -> pd.DataFrame:
    """
    Reads a multi-column reference matrix.
    ASSUMES a header is always present (header=0).
    Does NOT rename columns.
    """
    f_path = Path(path)
    if not f_path.exists():
        logger.error(f"Reference file not found: {path}")
        raise FileNotFoundError(f"Reference file not found: {path}")

    with open(f_path, "r", encoding="utf-8", errors="ignore") as f:
        buf = f.read(4096)
        if not buf:
            raise ValueError(f"Empty reference file: {path}")
        delim = _sniff_delimiter(buf)

    logger.info(f"Detected delimiter for {f_path.name}: {delim!r}")

    df = pd.read_csv(f_path, sep=delim, engine="python", dtype=str, header=0, on_bad_lines='warn')

    if df.empty:
        raise ValueError(f"Reference file is empty or could not be parsed: {path}")
    if df.shape[1] < 2:
        raise ValueError("Reference file must contain at least 2 columns (gene ID + samples).")

    df = df.dropna(how="all")
    if df.shape[1] < 2:
        raise ValueError("After cleaning, reference matrix has <2 columns.")
    return df


# --------------------------- Gene ID utilities ---------------------------------

ENSEMBL_RE_FLEX = re.compile(r"ENSG\d{9,}", re.IGNORECASE)

def strip_ensembl_version(gene_id: str) -> str:
    """Strips the version suffix (e.g., '.12') from an Ensembl ID."""
    return re.sub(r"\.\d+$", "", gene_id)

def contains_ensembl(gene_id: str) -> bool:
    """Checks if an Ensembl ID is anywhere in the string."""
    return bool(ENSEMBL_RE_FLEX.search(gene_id.strip()))

def extract_symbol_from_mixed(id_str: str) -> Optional[str]:
    """
    Tries to extract a SYMBOL from mixed ID strings like:
    'ENSG...|XIST', 'ENSG...;XIST', 'XIST|ENSG...', 'ENSG... XIST'
    """
    tokens = re.split(r"[|;,\s/]+", id_str.strip())
    for t in tokens:
        if t and not t.upper().startswith("ENSG") and re.match(r"^[A-Za-z][A-Za-z0-9._-]*$", t):
            return t.upper()
    return None

def determine_id_type(stripped_gene_ids: pd.Series, logger: logging.Logger) -> str:
    """
    Decides whether IDs are predominantly Ensembl or Symbols based on a
    hard count threshold (> 10000) of ENSG matches from *pre-stripped* IDs.
    """
    s = stripped_gene_ids.astype(str)
    ensg_count = s.apply(contains_ensembl).sum()

    if ensg_count > 10000:
        id_type = "Ensembl_ID"
    else:
        id_type = "Gene_symbol"

    logger.info(f"Auto-detected input gene ID type (count > 10000 rule): {id_type} (ENSG count={ensg_count})")
    return id_type


# --------------------------- Sex inference ---------------------------------

SEX_SYMBOLS_F: List[str] = ["XIST"]
SEX_SYMBOLS_M: List[str] = ["RPS4Y1", "RPS4Y2", "KDM5D", "DDX3Y", "USP9Y"]

# Map for both Ensembl (no version) and Symbol
SEX_GENE_MAP = {
    "ENSG00000229807": "XIST",
    "ENSG00000129824": "RPS4Y1",
    "ENSG00000280969": "RPS4Y2",
    "ENSG00000012817": "KDM5D",
    "ENSG00000067048": "DDX3Y",
    "ENSG00000114374": "USP9Y",
    "XIST": "XIST",
    "RPS4Y1": "RPS4Y1",
    "RPS4Y2": "RPS4Y2",
    "KDM5D": "KDM5D",
    "DDX3Y": "DDX3Y",
    "USP9Y": "USP9Y",
}

def _collect_sex_gene_rows(df: pd.DataFrame, id_type: str) -> Dict[str, List[int]]:
    """
    Return row indices for each target sex gene symbol.
    This function is robust to Symbols, Ensembl IDs, and Mixed IDs.
    """
    idx_map: Dict[str, List[int]] = {g: [] for g in SEX_SYMBOLS_F + SEX_SYMBOLS_M}

    for i, raw in enumerate(df["GeneID"].astype(str).values):

        gid_stripped = strip_ensembl_version(raw.strip())
        gid_upper = gid_stripped.upper()

        found_symbol = None

        if gid_upper in SEX_GENE_MAP:
            found_symbol = SEX_GENE_MAP[gid_upper]

        if not found_symbol:
            sym_from_mixed = extract_symbol_from_mixed(gid_stripped)
            if sym_from_mixed and sym_from_mixed in idx_map:
                found_symbol = sym_from_mixed

        if found_symbol:
            idx_map[found_symbol].append(i)

    return idx_map

def infer_sample_sex(df_expr: pd.DataFrame, logger: logging.Logger) -> Tuple[str, Dict[str, Dict[str, Any]]]:
    """
    Infer sex for the single sample.
    Returns:
        overall_label: "Male" | "Female" | "Unknown" | "Mixed" (from logic)
        details: {sample_name: {"X_median_expr": x, "Y_median_expr": y, "X_score": x_s, "Y_score": y_s, "decision": label}}
    """
    temp_stripped_ids = df_expr["GeneID"].astype(str).apply(strip_ensembl_version)
    id_type = determine_id_type(temp_stripped_ids, logger)

    idx_map = _collect_sex_gene_rows(df_expr, id_type)
    have_x = len(idx_map["XIST"]) > 0
    have_y = sum(len(idx_map[g]) for g in SEX_SYMBOLS_M) >= 1 # Use 1 for robustness

    sample_name = df_expr.columns[1]
    details: Dict[str, Dict[str, Any]] = {}

    if not (have_x or have_y):
        logger.warning("Sex inference: none of target genes found; returning 'Unknown'.")
        details[sample_name] = {
            "X_median_expr": 0.0, "Y_median_expr": 0.0,
            "X_score": 0.0, "Y_score": 0.0, "decision": "Unknown"
        }
        return "Unknown", details

    x_vals = []
    y_vals = []

    # Values are already numeric floats from main()
    if have_x:
        for r in idx_map["XIST"]:
            val = df_expr.iloc[r][sample_name]
            x_vals.append(val if not pd.isna(val) else 0.0)
    for yg in SEX_SYMBOLS_M:
        for r in idx_map[yg]:
            val = df_expr.iloc[r][sample_name]
            y_vals.append(val if not pd.isna(val) else 0.0)

    # Store median expression *before* log transformation
    x_median_expr = float(np.median(x_vals) if x_vals else 0.0)
    y_median_expr = float(np.median(y_vals) if y_vals else 0.0)

    X_score = float(np.log2(1.0 + x_median_expr))
    Y_score = float(np.log2(1.0 + y_median_expr))

    decision = "Unknown"
    if have_x and have_y:
        if (Y_score - X_score) >= 1.0:
            decision = "Male"
        elif (X_score - Y_score) >= 1.0:
            decision = "Female"
    elif have_y and not have_x:
        decision = "Male" if Y_score >= 1.0 else "Unknown"
    elif have_x and not have_y:
        decision = "Female" if X_score >= 1.0 else "Unknown"

    details[sample_name] = {
        "X_median_expr": x_median_expr, "Y_median_expr": y_median_expr,
        "X_score": X_score, "Y_score": Y_score, "decision": decision
    }
    overall = decision

    logger.info(f"Sex inference details: {details}")
    logger.info(f"Collapsed sex label: {overall}")
    return overall, details


# --------------------------- Reference processing ---------------------------------

def process_reference(ref_path: Union[str, Path],
                      id_type_needed: str,
                      outdir: Union[str, Path],
                      pid: int,
                      logger: logging.Logger) -> Tuple[Path, str]:
    """
    Standardizes the reference matrix:
    1. Intelligently selects the correct gene column (e.g., 'Ensembl_ID' or 'Gene_symbol').
    2. *** Automatically discards the *other* gene ID column ***
    3. Strips isoform versions from the selected gene column.
    4. Coerces expression columns to numeric.
    5. Keeps all reference genes (filtering of NaN/0 rows is disabled).
    6. Saves as CSV (with GeneID as index) and returns the file path AND the final ID type used.
    """
    ref_df = read_reference_matrix(ref_path, logger)

    gene_col_to_use = None
    original_id_type = id_type_needed
    other_type = "Gene_symbol" if id_type_needed == "Ensembl_ID" else "Ensembl_ID"

    # Standardize column names for searching
    col_map_lower = {col.lower().replace("_", ""): col for col in ref_df.columns}
    
    # Define standardized search names
    target_names_std = {id_type_needed.lower().replace("_", "")}
    other_names_std = {other_type.lower().replace("_", "")}
    
    other_gene_col_name = None # This will store the *original* column name to be dropped

    for std_name, original_col_name in col_map_lower.items():
        if std_name in target_names_std:
            gene_col_to_use = original_col_name
        elif std_name in other_names_std:
            other_gene_col_name = original_col_name # Found the *other* ID col

    if not gene_col_to_use:
        # Fallback: Try to find the *other* type if the target wasn't found
        for std_name, original_col_name in col_map_lower.items():
             if std_name in other_names_std:
                gene_col_to_use = original_col_name
                logger.warning(f"Input type is {original_id_type}, but reference only seems to have {gene_col_to_use}. Using it.")
                id_type_needed = other_type
                # We are using the 'other' type, so the 'target' type (if present) becomes the one to drop
                for std_name_T, original_col_name_T in col_map_lower.items():
                    if std_name_T in target_names_std:
                        other_gene_col_name = original_col_name_T
                        break
                break

    if not gene_col_to_use:
        # Critical Fallback: Could not find any standard ID column, use first column
        gene_col_to_use = ref_df.columns[0]
        logger.warning(f"Could not find 'Ensembl_ID' or 'Gene_symbol' in reference. Using first column '{gene_col_to_use}' as GeneID.")
        # Try to guess the *other* column to drop
        if len(ref_df.columns) > 1 and ("symbol" in ref_df.columns[1].lower() or "ensembl" in ref_df.columns[1].lower()):
            other_gene_col_name = ref_df.columns[1]
            logger.warning(f"Guessed that '{other_gene_col_name}' is the other ID column and will be dropped.")
        
        # Update id_type_needed based on this guess
        if "ensembl" in gene_col_to_use.lower():
            id_type_needed = "Ensembl_ID"
        elif "symbol" in gene_col_to_use.lower() or "gene" in gene_col_to_use.lower():
             id_type_needed = "Gene_symbol"

    logger.info(f"Using reference column '{gene_col_to_use}' as the Gene ID (processing as {id_type_needed}).")

    # Get *all* other columns as potential expression columns
    expr_cols = [c for c in ref_df.columns if c != gene_col_to_use]

    # ** CRITICAL STEP **
    # Remove the *other* gene ID column, if it was found
    if other_gene_col_name and other_gene_col_name in expr_cols:
        expr_cols.remove(other_gene_col_name)
        logger.info(f"Dropping non-primary gene ID column: '{other_gene_col_name}'")
    
    if not expr_cols:
        raise ValueError(f"Reference matrix file {ref_path} has no expression columns after processing.")

    # Select only the chosen GeneID column and the determined expression columns
    ref_df_processed = ref_df[[gene_col_to_use] + expr_cols].copy()
    ref_df_processed = ref_df_processed.rename(columns={gene_col_to_use: "GeneID"})

    ref_df_processed["GeneID"] = ref_df_processed["GeneID"].astype(str).apply(strip_ensembl_version)

    if id_type_needed == "Gene_symbol":
         ref_df_processed["GeneID"] = ref_df_processed["GeneID"].astype(str).str.upper()
    
    # Handle duplicate GeneIDs in the REFERENCE file
    if ref_df_processed['GeneID'].duplicated().any():
        logger.warning("Duplicate GeneIDs found in reference file. Aggregating by taking the mean.")
        # Need to separate numeric and non-numeric... oh wait, all expr_cols are numeric now
        # We can just groupby and mean.
        numeric_cols = expr_cols
        ref_df_processed[numeric_cols] = ref_df_processed[numeric_cols].apply(pd.to_numeric, errors='coerce')
        ref_df_processed = ref_df_processed.groupby('GeneID').mean()
        # After groupby, 'GeneID' becomes the index. We need to reset it to be a column
        ref_df_processed = ref_df_processed.reset_index()


    out_path = Path(outdir)
    out_path.mkdir(parents=True, exist_ok=True)

    # Append PID to create a unique intermediate filename (prevents race conditions)
    out_file = out_path / f"processed_{Path(ref_path).stem}_{pid}.csv"

    # Save as CSV with GeneID as index (for ARIC compatibility)
    ref_df_processed.set_index("GeneID").to_csv(out_file, sep=",")

    return out_file, id_type_needed

def maybe_filter_reference_by_sex(ref_csv_path: Path,
                                  sex_mode: str,
                                  inferred_sex: Optional[str],
                                  male_tissues: List[str],
                                  female_tissues: List[str],
                                  pid: int,
                                  logger: logging.Logger) -> Path:
    """
    Filters sex-specific tissues from the reference based on sex_mode and user-provided lists.
    'Unknown' or 'Mixed' will filter *all* sex-specific tissues.
    Matching is case-insensitive.
    """
    if sex_mode is None or sex_mode == "None":
        logger.info("Sex filtering disabled (sex=None).")
        return ref_csv_path

    # Auto-skip filtering for brain references
    if re.search(r"brain|cortex|neuron|cerebellum", ref_csv_path.name.lower()):
        logger.info("Detected brain-related reference. Skipping sex-based tissue filtering.")
        return ref_csv_path

    # Read CSV, assuming GeneID is the index
    df = pd.read_csv(ref_csv_path, sep=",", index_col=0)
    cols = list(df.columns)


    # Use case-insensitive matching
    SEX_SPECIFIC_FEMALE = set(t.lower() for t in female_tissues)
    SEX_SPECIFIC_MALE = set(t.lower() for t in male_tissues)


    present_f = [c for c in cols if any(key in c.lower() for key in SEX_SPECIFIC_FEMALE)]
    present_m = [c for c in cols if any(key in c.lower() for key in SEX_SPECIFIC_MALE)]

    if not present_f and not present_m:
        logger.info("Reference has no known sex-specific tissue columns; skipping sex filtering.")
        return ref_csv_path

    target_sex = None
    if sex_mode in ("Male", "Female"):
        target_sex = sex_mode
    elif sex_mode == "Auto":
        target_sex = inferred_sex

    keep_cols = [] # GeneID (index) is kept automatically
    dropped_count = 0

    if target_sex == "Male":
        keep_cols += [c for c in cols if c not in present_f]
        dropped_count = len(present_f)
        logger.info(f"Target sex is {target_sex}; dropping {dropped_count} female-specific tissues.")
    elif target_sex == "Female":
        keep_cols += [c for c in cols if c not in present_m]
        dropped_count = len(present_m)
        logger.info(f"Target sex is {target_sex}; dropping {dropped_count} male-specific tissues.")
    elif target_sex in ("Unknown", "Mixed"):
        # Filter both male and female tissues
        keep_cols += [c for c in cols if c not in present_f and c not in present_m]
        dropped_count = len(present_f) + len(present_m)
        logger.info(f"Inferred sex is {target_sex}; dropping ALL {dropped_count} sex-specific tissues.")
    else:
        logger.info("No sex filtering applied.")
        return ref_csv_path

    filtered_df = df[keep_cols].copy()

    # Append PID to create a unique intermediate filename (prevents race conditions)
    out_file_path = ref_csv_path.with_suffix(f".sexFiltered_{target_sex}_{pid}.csv")

    # Save filtered CSV, keeping GeneID as index
    filtered_df.to_csv(out_file_path, sep=",")

    return out_file_path


# --------------------------- Overlap & plotting ---------------------------------

def compute_overlap_ratio(expr_df: pd.DataFrame, ref_df_indexed: pd.DataFrame) -> float:
    """Computes the gene ID overlap ratio (Input & Ref) / (Input Total)."""
    s_expr = set(expr_df["GeneID"].astype(str).unique())
    s_ref = set(ref_df_indexed.index.astype(str).unique())
    if not s_expr:
        return 0.0
    return len(s_expr & s_ref) / float(len(s_expr))

def plot_fractions(tsv_path: Union[str, Path], out_prefix: str, logger: logging.Logger) -> None:
    """Generates a publication-quality horizontal bar plot using Seaborn."""

    # Safe font selection with fallback
    available_fonts = [f.name.lower() for f in font_manager.fontManager.ttflist]
    if "arial" in available_fonts:
        font_family = 'Arial'
    else:
        font_family = 'DejaVu Sans'
        if 'dejavu sans' not in available_fonts:
            font_family = 'sans-serif' # Absolute fallback
        logger.warning(f"Arial font not found. Falling back to '{font_family}'.")
    plt.rcParams['font.family'] = font_family

    # Set Seaborn theme
    sns.set_theme(style="whitegrid", font=font_family, rc={
        "axes.linewidth": 0.6,
        "axes.edgecolor": "0.4",
        "grid.linestyle": "--",
        "grid.linewidth": 0.3,
        "xtick.bottom": True,
        "ytick.left": True,
    })

    try:
        # Read the FINAL TSV file
        df = pd.read_csv(tsv_path, sep="\t")

        if df.empty or df.shape[1] < 2:
            logger.warning("Deconvolution result is empty (all fractions were 0); skip plotting.")
            return

        ct_col, frac_col = df.columns[0], df.columns[1]

        # Convert to percentage if it's in 0-1 range
        if (df[frac_col].max() <= 1.0 + 1e-6) and (df[frac_col].min() >= 0.0 - 1e-6):
            df[frac_col] = df[frac_col] * 100.0

        # Sort high-to-low for plotting
        df = df.sort_values(frac_col, ascending=False)

        # Dynamic figure height
        n_bars = df.shape[0]
        bar_height = 0.2
        spacing_factor = 0.2
        fig_height = max(2.0, n_bars * spacing_factor)

        plt.figure(figsize=(4.8, fig_height))

        # Added hue=ct_col and legend=False to fix seaborn FutureWarning
        ax = sns.barplot(
            x=frac_col,
            y=ct_col,
            data=df,
            palette="crest",
            edgecolor="none",
            hue=ct_col,
            legend=False
        )

        plt.xlabel("Percentage contribution (%)", fontsize=9)
        plt.ylabel("")
        plt.title("RNA source composition", fontsize=10, fontweight="bold", pad=8)
        plt.tick_params(labelsize=8)

        # Set x-axis ticks
        max_x = max(100, df[frac_col].max() * 1.05)
        plt.xlim(0, max_x)
        ax.xaxis.set_major_locator(MultipleLocator(10))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))

        # Style borders
        sns.despine(left=True, bottom=False)

        # Add text labels to bars
        for i, v in enumerate(df[frac_col].values):
            if v > 1:
                plt.text(v + 0.5, i, f"{v:.1f}", va="center", fontsize=7, color="dimgray")

        plt.tight_layout()

        # Save figures
        png_path = Path(f"{out_prefix}.png")
        pdf_path = Path(f"{out_prefix}.pdf")

        plt.savefig(png_path, dpi=300, bbox_inches="tight", transparent=True)
        plt.savefig(pdf_path, dpi=300, bbox_inches="tight", transparent=True)
        logger.info(f"Saved plots: {png_path.name} , {pdf_path.name}")

    except Exception as e:
        logger.error(f"Failed to generate plots: {e}")
    finally:
        plt.close()


# --------------------------- Main ---------------------------------

def main() -> None:

    epilog_text = (
        "Example:\n"
        "  python3 %(prog)s --input_expr_file my_sample.tsv --ref_expr_file gtex_ref.tsv --sex Auto --output_dir ./results \\\n"
        "    --male_specific_tissue \"Prostate\" \"Testis\" \\\n"
        "    --female_specific_tissue \"Ovary\" \"Uterus\" \"Vagina\" \"Cervix uteri\" \"Fallopian tube\"\n"
    )

    parser = argparse.ArgumentParser(
        description="Generalized RNA deconvolution via ARIC. This script processes a single-sample expression file against a reference matrix, handles gene ID and sex-based filtering, and outputs deconvolution fractions and a plot.",
        epilog=epilog_text,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("--input_expr_file", required=True,
                        help="Path to the input expression matrix for a single sample (e.g., bulk, pseudobulk). "
                             "Must contain exactly two columns: [GeneID, SampleData]. "
                             "Auto-detects header and delimiter.")
    parser.add_argument("--ref_expr_file", required=True,
                        help="Path to the reference expression matrix (e.g., TSV/CSV). "
                             "MUST have a header (e.g., 'Ensembl_ID', 'Gene_symbol', ...). "
                             "The script will auto-select the correct gene column.")
    parser.add_argument("--gene_id_type", choices=["Gene_symbol", "Ensembl_ID"], default=None,
                        help="Optionally, force the gene ID type ('Gene_symbol' or 'Ensembl_ID'). "
                             "If omitted, the script will auto-detect based on the input gene IDs.")
    parser.add_argument("--output_dir", required=True,
                        help="Path to the output directory. All results (logs, plots, TSV) "
                             "will be saved here. Will be created if it does not exist.")
    parser.add_argument("--sex", choices=["None", "Male", "Female", "Auto"], default="None",
                        help="Sex-specific tissue filtering for GTEx-like references. "
                             "'None': No filtering (use for brain references). "
                             "'Auto': Infer sex from expression and filter. 'Unknown'/'Mixed' sex results will filter *all* sex-specific tissues.")

    parser.add_argument("--male_specific_tissue", nargs='+',
                        default=["Prostate", "Testis"],
                        help="List of male-specific tissues to filter from reference. (Default: \"Prostate\" \"Testis\")")
    parser.add_argument("--female_specific_tissue", nargs='+',
                        default=["Ovary", "Uterus", "Vagina", "Cervix uteri", "Fallopian tube"],
                        help="List of female-specific tissues to filter from reference. (Default: \"Ovary\" \"Uterus\" ...)")

    parser.add_argument("--iter_num", type=int, default=10,
                        help="Number of iterations for ARIC. (Default: 10)")
    parser.add_argument("--confidence", type=float, default=0.75,
                        help="Confidence threshold for ARIC. (Default: 0.75)")
    parser.add_argument("--scale", type=float, default=0.1,
                        help="Scale parameter for ARIC. (Default: 0.1)")

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger = init_logger(output_dir, args.input_expr_file)
    logger.info(f"Args: {vars(args)}")

    # Initialize intermediate file paths for cleanup
    ref_processed_path = None
    ref_for_run_path = None
    input_csv_path = None
    ref_csv_path = None

    try:
        # Append PID to create a unique intermediate filename (prevents race conditions)
        pid = os.getpid()

        # 1. Read and process input file
        df_expr_raw = read_single_sample_expr_file(args.input_expr_file, logger)

        # Strip versions *before* detection
        stripped_ids_for_detection = df_expr_raw["GeneID"].astype(str).apply(strip_ensembl_version)

        in_id_type = args.gene_id_type or determine_id_type(stripped_ids_for_detection, logger)
        logger.info(f"Using input gene ID type: {in_id_type}")

        if in_id_type == "Ensembl_ID":
            df_expr_raw["GeneID"] = stripped_ids_for_detection
        else:
            df_expr_raw["GeneID"] = stripped_ids_for_detection.apply(
                lambda x: (extract_symbol_from_mixed(x) or x).upper()
            )

        sample_col_name = df_expr_raw.columns[1]
        df_expr_raw[sample_col_name] = pd.to_numeric(df_expr_raw[sample_col_name], errors="coerce")
        df_expr_raw = df_expr_raw.dropna(subset=[sample_col_name], how="all")

        # Handle duplicate GeneIDs before sex inference or indexing
        if df_expr_raw['GeneID'].duplicated().any():
            duplicates = df_expr_raw[df_expr_raw['GeneID'].duplicated()]['GeneID'].unique()
            dup_str = ", ".join(map(str, duplicates))
            logger.error(f"Duplicate GeneIDs found in input file: {dup_str}")
            logger.error("Please clean your input file to ensure all GeneIDs are unique before proceeding.")
            raise ValueError("Duplicate GeneIDs detected in input file. Aborting.")

        # 2. Infer Sex (if Auto)
        # Sex inference is run on the full (and now unique) expression data
        inferred_sex_label = None
        sex_detail = {}
        if args.sex == "Auto":
            inferred_sex_label, sex_detail = infer_sample_sex(df_expr_raw, logger)
            logger.info(f"Sex inference details: {sex_detail}")

        # 3. Process Reference
        # This function now correctly handles dropping the *other* gene ID column
        ref_processed_path, actual_ref_id_type = process_reference(
            args.ref_expr_file,
            id_type_needed=in_id_type,
            outdir=output_dir,
            pid=pid, # Pass PID for unique filename
            logger=logger
        )
        if in_id_type != actual_ref_id_type:
            logger.warning(f"Input ID type was {in_id_type}, but reference processing adjusted to {actual_ref_id_type}.")
            in_id_type = actual_ref_id_type # Update in_id_type to match what was *actually* used

        # 4. Filter Reference by Sex
        ref_for_run_path = maybe_filter_reference_by_sex(
            ref_processed_path,
            sex_mode=args.sex,
            inferred_sex=inferred_sex_label,
            male_tissues=args.male_specific_tissue,
            female_tissues=args.female_specific_tissue,
            pid=pid, # Pass PID for unique filename
            logger=logger
        )

        # 5. Load reference and create indexed input
        ref_df_for_overlap = pd.read_csv(ref_for_run_path, sep=",", index_col=0)

        if ref_df_for_overlap.index.empty:
            logger.error("Reference gene list is empty after processing. Aborting.")
            sys.exit(1)

        # This input_df_indexed now contains ALL genes from the input file (including 0s)
        # Duplicates were handled earlier, so set_index is safe.
        input_df_indexed = df_expr_raw.set_index("GeneID")

        # 6. Define gene list based on the REFERENCE matrix
        # The reference matrix defines the complete set of genes required for deconvolution.
        ref_genes_set = set(ref_df_for_overlap.index.astype(str))
        
        # Use a sorted list for consistent row order in both matrices
        # This is the MASTER list of genes.
        common_genes_list = sorted(list(ref_genes_set))
        common_genes_count = len(common_genes_list) # This is the count for filenames

        if not common_genes_list:
            logger.error("No genes found in the reference matrix. Aborting.")
            sys.exit(1)

        # Check how many of these reference genes are actually in the input file
        input_genes_set = set(input_df_indexed.index.astype(str))
        found_in_input_set = input_genes_set.intersection(ref_genes_set)
        found_in_input_count = len(found_in_input_set)

        ov = found_in_input_count / float(common_genes_count) if common_genes_count else 0.0
        
        logger.info(f"Aligning to the {common_genes_count} genes defined in the reference matrix.")
        logger.info(f"{found_in_input_count} of these genes were found in the input file.")
        logger.info(f"Gene overlap ratio (Found in Input / Total in Ref): {ov*100:.2f}%")
        logger.info("Genes present in Reference but missing from Input will be assigned 0 expression.")

        if ov < 0.5:
            logger.warning(f"Low gene overlap ratio ({ov*100:.2f}%). Results may be unreliable.")

        # 7. Align both input and reference to the common, sorted gene list
        
        logger.info(f"Aligning both input and reference matrices to {len(common_genes_list)} common genes.")

        # Align input matrix to the MASTER list.
        # .reindex() selects genes from common_genes_list,
        # and fills any missing genes (that are in ref but not input) with NaN.
        # .fillna(0.0) converts these NaNs to 0, as requested.
        df_expr_final = input_df_indexed.reindex(common_genes_list).fillna(0.0)
        # This ensures df_expr_final has exactly the genes in common_genes_list, in that order.

        # Align reference matrix to the *same* common list and order
        # .loc is used here because we know all genes in common_genes_list
        # are in ref_df_for_overlap (by definition).
        ref_df_final = ref_df_for_overlap.loc[common_genes_list].copy()

        # Ensure no NaNs in the final input (e.g., from weird dtype issues)
        df_expr_final = df_expr_final.fillna(0.0)

        # 8. Define output paths
        base_expr_name = Path(args.input_expr_file).name.split('.')[0]
        base_ref_name = Path(args.ref_expr_file).stem
        sex_tag = (inferred_sex_label if args.sex == "Auto" else (args.sex if args.sex != "None" else "NoSex"))

        # 9. Save the filtered matrices as intermediate CSVs (with index)
        # Corrected 'filterd' to 'filtered'
        input_csv_path = output_dir / f"Deconvolution_filtered_input_expr_file_{common_genes_count}_genes_{base_expr_name}.csv"
        df_expr_final.to_csv(input_csv_path, sep=",")
        logger.info(f"Processed and filtered input saved: {input_csv_path.name} (Rows: {len(df_expr_final)})")

        ref_csv_path = output_dir / f"Deconvolution_filtered_ref_expr_file_{common_genes_count}_genes_{base_ref_name}_{base_expr_name}.csv"
        
        ref_df_final.to_csv(ref_csv_path, sep=",")
        logger.info(f"Processed and filtered reference saved: {ref_csv_path.name} (Rows: {len(ref_df_final)})")


        # 10. Define final output paths
        # Using a clearer name for the output prefix
        prefix_str = str(output_dir / f"Deconvolution_fraction_{common_genes_count}_genes_from_ref_{base_ref_name}_sex-{sex_tag}_{base_expr_name}")


        out_tsv_path = Path(f"{prefix_str}.tsv") # Final output is TSV
        tmp_csv_path = out_tsv_path.with_suffix(".csv.tmp") # ARIC's output is CSV

        # 11. Save Sex Annotation TSV (if applicable)
        if args.sex == "Auto" and sex_detail:
            sex_tsv_path = output_dir / f"Deconvolution_sex_annotation_{base_expr_name}.tsv"
            df_sex = pd.DataFrame.from_dict(sex_detail, orient='index')
            df_sex = df_sex.reset_index().rename(columns={'index': 'Sample'})
            cols = ['Sample', 'X_median_expr', 'Y_median_expr', 'X_score', 'Y_score', 'decision']
            df_sex = df_sex[cols]
            df_sex.to_csv(sex_tsv_path, sep="\t", index=False, float_format="%.4f")
            logger.info(f"Sex annotation table saved (TSV): {sex_tsv_path.name}")

        # 12. Run ARIC
        logger.info("Running ARIC...")
        ARIC(
            str(input_csv_path),
            str(ref_csv_path), # Use the new, unique, common-gene-filtered ref file
            save_path=str(tmp_csv_path),
            marker_path=None,
            selected_marker=False,
            scale=args.scale,
            iter_num=args.iter_num,
            confidence=args.confidence,
            delcol_factor=10,
            w_thresh=10,
            unknown=False,
            is_methylation=False
        )

        if not tmp_csv_path.exists() or tmp_csv_path.stat().st_size == 0:
            logger.error(f"ARIC did not produce expected file or file is empty: {tmp_csv_path}")
            sys.exit(1)

        # 13. Format and save final TSV
        df_out = pd.read_csv(tmp_csv_path)
        if df_out.shape[1] >= 2:
            df_out.columns = ["Label", "Fraction"] + list(df_out.columns[2:])

        df_out["Fraction"] = pd.to_numeric(df_out["Fraction"], errors='coerce').fillna(0.0)
        df_out = df_out[df_out["Fraction"] > 0]
        df_out["Fraction"] = df_out["Fraction"].round(4)



        df_out.to_csv(out_tsv_path, sep="\t", index=False)
        logger.info(f"Deconvolution result saved (TSV): {out_tsv_path.name}")
        tmp_csv_path.unlink() # Delete temp file

        # 14. Plot
        plot_fractions(out_tsv_path, prefix_str, logger)

        logger.info("All done.")

    except Exception as e:
        logger.error(f"An error occurred: {e}", exc_info=True)
        sys.exit(1)

    finally:
        # 15. Cleanup intermediate files (silently)
        if input_csv_path and input_csv_path.exists():
            input_csv_path.unlink(missing_ok=True)

        if ref_csv_path and ref_csv_path.exists():
            ref_csv_path.unlink(missing_ok=True)

        if ref_for_run_path and ref_for_run_path.exists():
            ref_for_run_path.unlink(missing_ok=True)

        if (ref_processed_path and
            ref_processed_path.exists() and
            ref_processed_path != ref_for_run_path): # Don't delete twice
            ref_processed_path.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
