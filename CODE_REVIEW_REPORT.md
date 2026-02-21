# EVscope Code Review Report

**Date:** 2025
**Version:** EVscope v1.0.0
**Scope:** All Python scripts in `bin/`, shell scripts, configuration files, and documentation

---

## Executive Summary

EVscope is a comprehensive 27-step modular RNA-seq pipeline for extracellular vesicle analysis, covering circRNA detection, tissue deconvolution, and extensive quality control. The codebase demonstrates high quality with robust error handling, modular design, and thorough documentation. This review identified **5 bugs**, **1 security issue**, and **several code quality improvements**, all of which have been resolved in v1.0.0.

---

## Issues Identified and Resolved

### Bugs

| ID | File | Description | Severity | Status |
|----|------|-------------|----------|--------|
| BUG-1 | `Step_03_UMIAdapterTrimR1.py` | Typo "readthrouth" in variable names and TSV headers | Medium | Fixed |
| BUG-2 | `Step_24_generate_QC_matrix.py` | Typo "Fagments)" in QC metric label | Medium | Fixed |
| BUG-3 | `Step_10_circRNA_merge.py` | Venn diagram only saved as PNG; CLI described PDF output | Low | Fixed (saves both PNG and PDF) |
| BUG-4 | `Step_03_UMIAdapterTrimR1.py` | `logging.info()` called before logger configured | Low | Fixed |
| BUG-5 | `Step_02_plot_fastq2UMI_motif.py` | Temp file not cleaned up on error | Low | Fixed (try/finally) |

### Security

| ID | File | Description | Severity | Status |
|----|------|-------------|----------|--------|
| SEC-1 | `Step_07_bam2strand.py` | Shell injection via `subprocess.run(shell=True)` with f-string | High | Fixed (list args, `shell=False`) |

### Code Quality

| ID | File | Description | Status |
|----|------|-------------|--------|
| CQ-1 | Multiple files | Inconsistent shebang (`python` vs `python3`) | Standardized to `python3` |
| CQ-2 | `Step_09_convert_CIRI2CPM.py` | Unused `sys` import | Removed |
| CQ-3 | `Step_02_plot_fastq2UMI_motif.py` | Hardcoded log file path in working directory | Redirected to output directory |
| CQ-4 | `Step_25_bigWig2Expression.py` | Type annotation mismatch (`Logger` vs `None`) | Fixed to `Optional[logging.Logger]` |
| CQ-5 | `EVscope.conf` | Hardcoded absolute paths for Kraken2 | Changed to `${EVscope_PATH}` variables |
| CQ-6 | `Step_22_run_RNA_deconvolution_ARIC.py` | Hardcoded absolute path in usage example | Replaced with relative path |
| CQ-7 | `Step_15_plot_top_expressed_genes.py` | Missing shebang and module docstring | Added |

---

## Codebase Strengths

1. **Robust error handling** — strict mode (`set -eEuo pipefail`), trap handlers, and comprehensive input validation in `EVscope.sh`
2. **Performance optimization** — Numba JIT compilation in EMapper and UMI trimming for computationally intensive operations
3. **Modular architecture** — each of the 27 steps is self-contained with clear CLI interfaces via `argparse`
4. **Publication-quality visualization** — consistent use of Arial font, PDF Type 42 embedding, and 300 DPI across all plotting scripts
5. **Robust deconvolution** — Step_22 (ARIC) features auto-detection of delimiters, gene ID types, and sex-specific tissue filtering
6. **Lazy dependency loading** — EMapper prevents import failures for optional dependencies (e.g., `pyBigWig`)
7. **Race condition prevention** — PID-based temporary filenames in parallel-safe scripts
