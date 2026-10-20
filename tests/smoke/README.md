# EVscope smoke tests

These tests validate repository installability, kit-aware command-line behavior, and selected script-level transformations without requiring the full human reference bundle. They are intended for CI and release checks.

The v1.1 smoke coverage includes:

- Bash and Python syntax plus CLI help/version checks.
- The canonical `--kit` values `smarter_total_rna` and `smarter_smrna`, including rejection of unsupported values.
- CLI help for the backward-compatible `smarter_total_rna` default.
- SMARTer smRNA paired-end validation, with R1 treated as sense and an implicit `forward` strand when `--strand` is not supplied.
- Kit-aware Steps 2, 3, and 6 behavior: the Total RNA UMI/ACC workflow remains available, while the smRNA profile applies the R1 `-u 3`/poly(A) trim, the R2 14-nt poly(T)/3' adapter trim, Q20 filtering of both mates, and paired no-UMI/no-deduplication alignment.
- STAR common/profile parameters and global/kit-specific extra parameters are injected from the configuration arrays without unsafe word splitting.
- SMARTer smRNA STAR profiles retain chimeric-junction output so circRNA Steps 8-10 can run without UMI deduplication.
- Step metadata records the normalized kit and effective strand, and cross-kit `--resume` is rejected.
- Selected toy FASTQ, read-length plotting, and EMapper/BigWig transformations when their optional dependencies are installed.
- Standalone Step 27 concise-summary generation, including human/contaminant mapping, genomic-region, strand-orientation, and top-gene reporting logic.

The SMARTer smRNA pipeline checks use deterministic test doubles for cutadapt, STAR, and samtools, plus placeholder plots inside that mocked branch. They verify the synthetic trimming fixtures, command arguments, output routing, and metadata; they do not validate production alignment or circRNA accuracy. The separate optional plotting and EMapper checks run the real Python modules.

The smoke suite does not run the complete 27-step workflow or render the comprehensive R Markdown report. Full end-to-end tests require the Zenodo reference bundle and R environment documented in the main README. All committed sequencing fixtures are synthetic; generated logs, reports, metadata, and temporary test environments stay under the ignored `work/` directory.

Run:

```bash
bash tests/smoke/run_smoke.sh
```

Temporary outputs are written under the ignored `tests/smoke/work/` directory. Optional checks report that they were skipped when their lightweight dependencies are unavailable.

For a release candidate, also validate the tracked release files from the repository root:

```bash
sha256sum -c repro/checksums.sha256
```
