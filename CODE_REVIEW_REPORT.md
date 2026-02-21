# EVscope Code Review Report

**Date:** 2026-02-21
**Reviewer:** Claude (Automated Code Review)
**Version:** EVscope v3.0.0
**Scope:** All Python scripts in `bin/`, configuration file `EVscope.conf`

---

## Executive Summary

EVscope is a well-engineered, production-grade RNA-seq pipeline with 27 steps covering circRNA detection, tissue deconvolution, and comprehensive QC. The codebase is generally high quality with good error handling, modular design, and thorough documentation. This review identified **5 bugs**, **2 security issues**, and **several minor code quality improvements**.

---

## Critical Issues (Bugs)

### BUG-1: Typo "readthrouth" in Step_03_UMIAdapterTrimR1.py (Lines 303, 307, 326-327)
- **Severity:** Medium (affects output file headers and variable names)
- **Description:** The word "readthrough" is consistently misspelled as "readthrouth" in variable names and output TSV headers. Downstream scripts or users parsing these headers will fail if they expect the correct spelling.
- **Fix:** Rename all instances of `readthrouth` → `readthrough`.

### BUG-2: Typo in QC metric name in Step_24_generate_QC_matrix.py (Line 265)
- **Severity:** Medium (corrupted QC metric label)
- **Description:** `"Percentage of UMI-dedup Fagments)"` has two issues:
  1. Missing opening parenthesis — should be `"(Fragments)"`
  2. "Fagments" should be "Fragments"
- **Fix:** Change to `"Percentage of UMI-dedup Fragments"`.

### BUG-3: Venn diagram saves as PNG but description says PDF in Step_10_circRNA_merge.py (Line 138)
- **Severity:** Low (output format mismatch with CLI description)
- **Description:** The `--out_venn` argument help says "Output PDF file for Venn diagram" but `plt.savefig()` at line 138 uses `format='png'`. The file extension may also mismatch.
- **Fix:** Save both PNG and PDF formats.

### BUG-4: `logging.info()` called before logger is configured in Step_03 (Line 258)
- **Severity:** Low (message silently dropped)
- **Description:** `logging.info("Starting the paired-end UMI trimming process...")` is called at line 258 but no logging handler is configured. The message is never seen.
- **Fix:** Configure basic logging or use `print()` instead.

### BUG-5: Temp file not cleaned up on error in Step_02_plot_fastq2UMI_motif.py (Lines 104-121)
- **Severity:** Low (temp files accumulate on failures)
- **Description:** If `parse_fastq_and_count()` or plotting fails, the temporary FASTQ file created by `seqtk` at line 104 is not cleaned up.
- **Fix:** Use `try/finally` to ensure cleanup.

---

## Security Issues

### SEC-1: Shell injection vulnerability in Step_07_bam2strand.py (Line 42-43)
- **Severity:** High
- **Description:** The command is constructed using an f-string and executed with `shell=True`:
  ```python
  cmd = f"infer_experiment.py -i {bam_file} -r {refgene_bed} -s {test_read_num}"
  subprocess.run(cmd, shell=True, ...)
  ```
  If `bam_file` or `refgene_bed` contain special shell characters (spaces, semicolons, etc.), this could cause command injection or unexpected behavior.
- **Fix:** Use a list of arguments with `shell=False`.

### SEC-2: File handle leak in Step_02_calculate_ACC_motif_fraction.py (Lines 53-55)
- **Severity:** Low
- **Description:** The file is opened outside the `with` block. If the `with` statement itself fails, the handle may leak. However, the current code structure makes this unlikely in practice.
- **Fix:** Already handled by the `with f:` pattern — no change needed.

---

## Code Quality Issues

### CQ-1: Inconsistent shebang lines
- **Files:** `Step_03_UMIAdapterTrimR1.py`, `Step_07_bam2strand.py`, `Step_18_plot_reads_mapping_stats.py`
- **Issue:** Use `#!/usr/bin/env python` instead of `#!/usr/bin/env python3`. This may invoke Python 2 on some systems.
- **Fix:** Standardize to `#!/usr/bin/env python3`.

### CQ-2: Unused imports
- **Step_09_convert_CIRI2CPM.py:** `sys` imported but never used.
- **Step_15_plot_top_expressed_genes.py:** `num_types` variable assigned but never used (line 77).

### CQ-3: Hardcoded log file path in Step_02_plot_fastq2UMI_motif.py (Line 35)
- **Description:** `logging.FileHandler('processing.log')` writes to the current working directory, not the output directory. This could fail or pollute unexpected locations.
- **Fix:** Write log to the output directory.

### CQ-4: `global_logger` type annotation in Step_25_bigWig2Expression.py (Line 52)
- **Description:** `global_logger: logging.Logger = None` — type annotation says `Logger` but initial value is `None`.
- **Fix:** Use `Optional[logging.Logger]`.

### CQ-5: Empty README.md
- **Description:** The project README is empty, making it hard for new users to understand the pipeline.

---

## Summary of Changes Made

| ID | File | Change | Severity |
|----|------|--------|----------|
| BUG-1 | Step_03_UMIAdapterTrimR1.py | Fixed "readthrouth" → "readthrough" | Medium |
| BUG-2 | Step_24_generate_QC_matrix.py | Fixed "Fagments)" → "Fragments" | Medium |
| BUG-3 | Step_10_circRNA_merge.py | Added PDF output alongside PNG | Low |
| BUG-4 | Step_03_UMIAdapterTrimR1.py | Added logging.basicConfig() | Low |
| BUG-5 | Step_02_plot_fastq2UMI_motif.py | Added try/finally for temp file cleanup | Low |
| SEC-1 | Step_07_bam2strand.py | Changed shell=True to shell=False with list args | High |
| CQ-1 | Multiple files | Standardized shebang to python3 | Low |
| CQ-2 | Step_09_convert_CIRI2CPM.py | Removed unused `sys` import | Low |
| CQ-3 | Step_02_plot_fastq2UMI_motif.py | Log file written to output directory | Low |
| CQ-4 | Step_25_bigWig2Expression.py | Fixed type annotation to Optional | Low |

---

## Positive Observations

1. **Excellent error handling** throughout, especially in `EVscope.sh` (strict mode, trap handlers)
2. **Performance optimization** with Numba JIT in EMapper and UMI trimming
3. **Robust design** in Step_22 (ARIC deconvolution) with auto-detection of delimiters, ID types, and sex inference
4. **Good modularity** — each step is self-contained with clear CLI interfaces
5. **Publication-quality visualization** settings across all plotting scripts (Arial font, PDF type 42)
6. **Lazy dependency loading** in EMapper prevents import failures for optional tools
7. **Race condition prevention** using PID-based temp filenames in Step_22
