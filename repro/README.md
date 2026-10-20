# EVscope reproducibility index

This directory records release-level reproducibility metadata for EVscope.

Large reference files, SRA sequencing data, controlled-access AMP-PD data, and some full benchmark matrices are not stored directly in Git. Use the DOI/accessions and notes in `manifest.tsv`, and synchronize final public-release assets with Zenodo/GitHub before journal upload when redistribution is permitted.

Minimum local release-candidate checks:

```bash
bash tests/smoke/run_smoke.sh
sha256sum -c repro/checksums.sha256
```

The smoke test covers syntax, kit-aware CLI validation, SMARTer smRNA structure trimming/no-UMI alignment with synthetic reads, STAR argument injection from config, cross-kit resume protection, standalone concise-summary report logic, read-length plotting when matplotlib is available, and a tiny synthetic EMapper/BigWig/bigWig2Expression run when `pysam`, `pyBigWig`, `numba`, and `numpy` are available.

`repro/checksums.sha256` covers tracked release-critical files except itself. Text files use LF line endings through `.gitattributes` so hashes match across platforms. Git commit/tag identity remains the authoritative source-tree integrity record.
