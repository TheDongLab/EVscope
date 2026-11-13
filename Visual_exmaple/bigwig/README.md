# Chromosome 1 strand tracks

These BigWig files are chromosome-1 subsets of the same example analysis shown in the HTML report. They contain the combined unique-read and EM-weighted multi-mapping coverage for each original EMapper orientation group.

| Download | EMapper alignment-orientation group | Size |
| --- | --- | --- |
| [Example_001_chr1_F1R2.bw](Example_001_chr1_F1R2.bw) | Forward R1 + reverse R2 | 1.43 MB |
| [Example_001_chr1_F2R1.bw](Example_001_chr1_F2R1.bw) | Reverse R1 + forward R2 | 1.49 MB |

Each file contains only `chr1`, with chromosome length 248,956,422 bp, matching the source GRCh38/hg38 reference. Load both files in a genome browser such as IGV with the matching reference, then navigate to chromosome 1.

The F1R2/F2R1 labels describe paired-read alignment orientation; RNA sense/antisense interpretation depends on the library protocol. Original coordinate positions, signal values, sign, and scale are preserved. No smoothing, binning, renormalization, or sign inversion was applied.

Every retained interval and its float32 value was compared with the original full-genome track: 151,217 intervals in F1R2 and 160,893 intervals in F2R1. [validation.json](validation.json) records these checks and the output SHA-256 values. The source files were not modified.

Sample identifiers have been replaced by `Example_001`. The files retain real coverage data. Full-genome tracks and the separate unique-only/multi-only diagnostic tracks are excluded from this small example.

Use GitHub's **Download raw file** button for each track, or obtain the files by cloning the repository. File integrity is also covered by [the repository checksum list](../../repro/checksums.sha256).
