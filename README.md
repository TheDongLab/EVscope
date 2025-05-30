# EVscope Pipeline

**EVscope** is a comprehensive bioinformatics pipeline for processing extracellular vesicle (EV) RNA sequencing data from raw paired-end FASTQ files. It performs quality control, UMI extraction, adapter trimming, alignment, deduplication, circular RNA detection, expression quantification, taxonomic classification, rRNA detection, and deconvolution. Optimized for data generated using the SMARTer-seq protocol with 14-bp UMIs in Read2, EVscope integrates tools like FastQC, STAR, CIRCexplorer2, CIRI2, RSEM, and Kraken2 to deliver robust analyses.

![EVscope Pipeline Overview](figures/EVscope_pipeline.png)

## Table of Contents
- [Features](#features)
- [Directory Structure](#directory-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Input Data Format](#input-data-format)
- [Pipeline Steps](#pipeline-steps)
- [Output Structure](#output-structure)
- [Troubleshooting](#troubleshooting)
- [FAQ](#faq)
- [Contributing](#contributing)
- [Feedback](#feedback)
- [License](#license)
- [Contact](#contact)

## Features
- **Quality Control**: Evaluates raw and trimmed FASTQ files using FastQC.
- **UMI Processing**: Extracts and deduplicates 14-bp UMIs from Read2.
- **Alignment**: Maps reads to the human genome (hg38) with STAR and BWA.
- **Circular RNA Detection**: Identifies circRNAs using CIRCexplorer2 and CIRI2, with merged results.
- **Expression Quantification**: Generates TPM/CPM matrices with featureCounts and RSEM.
- **Contamination Detection**: Detects bacterial (BBSplit) and microbial (Kraken2) reads.
- **rRNA Detection**: Filters ribosomal RNA with RiboDetector.
- **Deconvolution**: Estimates cell type proportions using GTEx and single-cell references.
- **Visualization**: Produces bigWig tracks, density plots, and an HTML report.
- **Flexible Workflow**: Supports running specific steps or the full pipeline.

## Directory Structure
The repository is organized as follows:
```
EVscope_pipeline/
├── EVscope.def           # Pipeline definition file
├── EVscope.sh            # Main pipeline script
├── EVscope_v1.sh         # Alternative pipeline script
├── README.md             # This documentation
├── bin/                  # Custom scripts
│   ├── Step_02_calculate_ACC_motif_ratio.py
│   ├── Step_02_plot_fastq2UMI_motif.py
│   ├── Step_03_UMIAdapterTrimR1.py
│   └── ... (other Python and shell scripts)
├── config/               # Configuration files (currently empty)
├── references/           # Reference files
│   ├── annotations_HG38/ # HG38 annotations (GTF, SAF, BED)
│   ├── deconvolution_HG38/ # GTEx and single-cell references
│   ├── genome/           # hg38, Mycoplasma, E. coli FASTA
│   └── index/            # STAR, BWA, RSEM indexes
├── soft/                 # External tools
│   ├── bbmap/            # BBMap for BBSplit
│   ├── CIRI_v2.0.6/      # CIRI2 for circRNA detection
│   ├── kraken2/          # Kraken2 for taxonomic classification
│   └── RSEM_v1.3.3/      # RSEM for expression quantification
├── TEMP/                 # Temporary files
└── test_data/            # Test datasets and logs
    ├── DNAnexus_CTRL_chr21_fastq/ # Sample FASTQ files
    ├── small_data/                # Small test data
    ├── test/                      # Test output
    ├── processing.log             # Processing logs
    └── run.sh                     # Test run script
```

## Requirements

### Software
- **Operating System**: Linux (e.g., Ubuntu) or macOS with a terminal.
- **Bash**: Version 4.0 or higher.
- **Conda**: Miniconda or Anaconda for environment management.
- **Tools** (install via Conda where possible):
  - FastQC
  - MultiQC
  - umi_tools
  - cutadapt
  - trim_galore
  - ribodetector_cpu
  - seqtk
  - kraken2
  - KronaTools (includes `kreport2krona.py`)
  - STAR (v2.7.10b)
  - samtools
  - deepTools (for bamCoverage, computeMatrix, plotProfile)
  - featureCounts (subread package)
  - CIRCexplorer2
  - CIRI2 (v2.0.6, included in `soft/`)
  - BBMap (included in `soft/`)
  - Picard
  - RSEM (v1.3.3, included in `soft/`)
  - R with `rmarkdown` package
  - Python 3.8+ with libraries: `pandas`, `numpy`, `matplotlib`, `seaborn`

### Reference Files
Ensure the following are in the `references/` directory:
- **Genomes**:
  - `genome/hg38/hg38.p14.whole.genome.fa`
  - `genome/mycoplasma/mycoplasma_ge.fa`
  - `genome/ecoli/E.coli_ge.fa`
- **Indexes**:
  - `index/RSEM_bowtie2_index/RSEM_REF_HG38_3659642_RNNAs/`
  - `index/Index_STAR-2.7.10b_sjdbOverhang99_GeneCode45/`
  - `index/Index_BWA_V0.7.18/hg38.p14.whole.genome.fa`
- **Annotations**:
  - `annotations_HG38/HG38_3659642_combined_RNAs_with_gold_standard_piRNAs.gtf`
  - `annotations_HG38/gencode.v45.chr_patch_hapl_scaff.annotation_UCSC.gtf`
  - SAF/BED files: `HG38_3UTR_noOverlap.saf`, `HG38_5UTR_noOverlap.saf`, etc.
  - `annotations_HG38/GeneCodeV45_combined_Nucl_Mt_rRNA_589.interval_list`
- **Deconvolution**:
  - `deconvolution_HG38/GTEx_v10_SMTS_AVG_TPM.csv`
  - `deconvolution_HG38/GTEx_v10_SMTSD_AVG_TPM.csv`
  - `deconvolution_HG38/human_brain_single_cell_atlasV1_31superclusters_AVP_CPM.csv`
- **Kraken2 Database**:
  - `soft/kraken2/standard_krakendb/`

### Hardware
- **CPU**: Multi-core processor (recommended: 20+ threads).
- **RAM**: Minimum 64 GB (Picard requires up to 250 GB).
- **Storage**: 500 GB+ for input, reference, and output files.

## Installation

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/yourusername/EVscope.git
   cd EVscope
   ```

2. **Install Conda** (if not installed):
   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   source ~/.bashrc
   ```

3. **Create Conda Environments**:
   Create environments for the pipeline and dependencies:
   ```bash
   conda env create -f environments/evscope_env.yml
   conda env create -f environments/picard_env.yml
   conda env create -f environments/kraken2_env.yml
   ```
   Example `environments/evscope_env.yml`:
   ```yaml
   name: evscope_env
   channels:
     - bioconda
     - conda-forge
   dependencies:
     - fastqc
     - multiqc
     - umi_tools
     - cutadapt
     - trim_galore
     - ribodetector
     - seqtk
     - star=2.7.10b
     - samtools
     - deeptools
     - subread
     - circos
     - python=3.8
     - pandas
     - numpy
     - matplotlib
     - seaborn
     - r-base
     - r-rmarkdown
   ```
   Example `environments/picard_env.yml`:
   ```yaml
   name: picard_env
   channels:
     - bioconda
     - conda-forge
   dependencies:
     - picard
     - imagemagick
   ```
   Example `environments/kraken2_env.yml`:
   ```yaml
   name: kraken2_env
   channels:
     - bioconda
     - conda-forge
   dependencies:
     - kraken2
     - kronatools
   ```
   > **Note**: Create these YAML files in an `environments/` directory or install dependencies manually.

4. **Install CIRCexplorer2**:
   ```bash
   conda activate evscope_env
   pip install CIRCexplorer2
   ```

5. **Verify Included Tools**:
   The following tools are bundled in `soft/`:
   - BBMap: `soft/bbmap/`
   - CIRI2: `soft/CIRI_v2.0.6/`
   - RSEM: `soft/RSEM_v1.3.3/`
   Ensure they are accessible:
   ```bash
   ls soft/bbmap/bbsplit.sh
   ls soft/CIRI_v2.0.6/CIRI2.pl
   ls soft/RSEM_v1.3.3/rsem-calculate-expression
   ```

6. **Prepare Reference Files**:
   Download and place reference files in `references/` as shown in [Requirements](#requirements). Example:
   ```bash
   mkdir -p references/genome/hg38
   wget -O references/genome/hg38/hg38.p14.whole.genome.fa <hg38_download_link>
   ```
   > **Note**: Obtain references from trusted sources (e.g., UCSC, Gencode, GTEx).

7. **Set Up Kraken2 Database**:
   ```bash
   mkdir -p soft/kraken2/standard_krakendb
   conda activate kraken2_env
   kraken2-build --standard --db soft/kraken2/standard_krakendb --threads 20
   ```

8. **Test Installation**:
   ```bash
   conda activate evscope_env
   fastqc --version
   STAR --version
   samtools --version
   python --version
   CIRCexplorer2 --version
   ```

## Usage

### Command Syntax
```bash
bash EVscope.sh <EVscope_path> -t <threads> -o <output_dir> --step <step_list> --input_fastq <R1_fastq.gz> <R2_fastq.gz>
```

- `<EVscope_path>`: Path to the EVscope directory (e.g., `/home/yz2474/yiyong_2023/EVscope_pipeline`).
- `-t <threads>`: Number of CPU threads (default: 1, recommended: 10-20).
- `-o <output_dir>`: Output directory.
- `--step <step_list>`: Steps to run, e.g., `[1,2,3]` (default: `[1-18]`).
- `--input_fastq <R1_fastq.gz> <R2_fastq.gz>`: Paired-end FASTQ files.

### Example
Run steps 1-3 with 10 threads using test data:
```bash
bash EVscope.sh /home/yz2474/yiyong_2023/EVscope_pipeline -t 10 -o test_output --step [1,2,3] --input_fastq test_data/DNAnexus_CTRL_chr21_fastq/chr21_R1_001.fastq.gz test_data/DNAnexus_CTRL_chr21_fastq/chr21_R2_001.fastq.gz
```

### Test Run
Verify setup with the included test data:
```bash
cd test_data
bash run.sh
```
Or manually:
```bash
bash ../EVscope.sh $(pwd)/.. -t 4 -o test --step [1,2] --input_fastq small_data/test_R1.fastq.gz small_data/test_R2.fastq.gz
```

### Notes
- **Sample Name**: Extracted from R1 FASTQ filename (e.g., `chr21` from `chr21_R1_001.fastq.gz`).
- **Paths with Spaces**: Quote paths with spaces, e.g., `"/path with spaces/fastq.gz"`.
- **Step Selection**: Use `--step` to run specific steps. Steps with `step.done` are skipped.
- **Logs**: Check step-specific logs (e.g., `Step_03_cutadapt/UMI_extract.log`) for debugging.

## Input Data Format
- **FASTQ Files**: Paired-end, gzipped FASTQ files (`R1.fastq.gz`, `R2.fastq.gz`).
- **Protocol**: Designed for EV RNA-seq data from the SMARTer-seq protocol with 14-bp UMIs in Read2.
- **Naming**: R1 and R2 should share a prefix (e.g., `Sample_007_R1.fastq.gz`, `Sample_007_R2.fastq.gz`).
- **Quality**: High-quality reads with sufficient coverage.

> **Warning**: Non-SMARTer-seq data may require modifying UMI extraction in `bin/Step_02_*.py`.

## Pipeline Steps
EVscope includes 18 steps, each generating specific outputs:
1. **Raw FASTQ QC**: FastQC and 14-bp UMI motif analysis (`Step_02_plot_fastq2UMI_motif.py`).
2. **UMI Extraction & Trimming**: UMI extraction (`umi_tools`) and adapter trimming (`cutadapt`).
3. **Post-Trimming QC**: FastQC on trimmed FASTQ.
4. **Bacterial Detection**: Identifies E. coli/Mycoplasma reads (BBSplit).
5. **STAR Alignment & Deduplication**: STAR alignment, UMI deduplication, and re-alignment.
6. **Strand Detection**: RNA strandness analysis (`Step_07_bam2strand.py`).
7. **circRNA Detection (CIRCexplorer2)**: Detects circRNAs from STAR chimeric junctions.
8. **circRNA Detection (BWA & CIRI2)**: Detects circRNAs with BWA and CIRI2, merges results.
9. **Picard Metrics**: RNA-seq and insert size metrics.
10. **Expression Quantification**: Read counting with featureCounts and RSEM.
11. **Expression Matrix & Plots**: TPM/CPM matrices and RNA distribution plots.
12. **Meta-Feature Mapping**: Quantifies reads in genomic regions (e.g., 3'UTR, exons).
13. **Taxonomic Classification**: Downsamples and classifies reads (Kraken2).
14. **Deconvolution**: Cell type proportion estimation (GTEx, single-cell).
15. **rRNA Detection**: Identifies rRNA (RiboDetector).
16. **QC Matrix**: Compiles quality metrics (`Step_16_generate_QC_matrix.py`).
17. **bigWig Generation**: Normalized coverage tracks (bamCoverage).
18. **HTML Report**: Summary report (`Step_18_html_report.Rmd`).

> **Note**: Steps 13 and 14 are optional and fail gracefully.

## Output Structure
Results are stored in:
```
<output_dir>/<sample_name>_EVscope_output/
├── Step_01_raw_fastqc/
│   └── <sample_name>_R1_fastqc.html
├── Step_02_UMI_motif/
│   └── <sample_name>_UMI_motif_plot.png
├── Step_03_cutadapt/
│   └── <sample_name>_R1_clean.fq.gz
├── Step_04_cutadapt_fastqc/
├── Step_05_BBSplit/
├── Step_06_STAR/
│   ├── 1st_STAR_UMI_tools/
│   └── 2nd_STAR/
├── Step_07_strand_detection/
├── Step_08_circRNA_detection/
│   ├── CIRCexplorer2/
│   ├── BWA_CIRI2/
│   └── merge_CIRI2_CIRCexplorer2/
├── Step_09_picard_metrics/
├── Step_10_RNA_expr_quantification/
│   ├── featureCounts/
│   └── RSEM/
├── Step_11_RNA_distribution/
├── Step_12_reads_mapping_stats/
├── Step_13_Kraken/
├── Step_14_bulk_RNA_deconvolution/
├── Step_15_RiboDetector/
├── Step_16_QC_matrix_generation/
│   └── <sample_name>_QC_matrix.tsv
├── Step_17_bigWig/
│   └── bamCoverage/
├── Step_18_HTML_report/
│   └── <sample_name>_final_report.html
└── TEMP/
```

Key outputs:
- **QC**: FastQC reports, UMI motif plots.
- **FASTQ**: Trimmed reads.
- **BAM**: Aligned, deduplicated BAM files.
- **circRNA**: TSV files, Venn diagrams.
- **Expression**: TPM/CPM matrices, RNA plots.
- **Deconvolution**: Cell type estimates.
- **QC Matrix**: Metrics in `Step_16_QC_matrix_generation/`.
- **Coverage**: bigWig files, density plots.
- **Report**: HTML summary in `Step_18_HTML_report/`.

## Troubleshooting
- **"Unknown option" Error**:
  - Check syntax: `<EVscope_path>` is first, followed by options.
  - Example: `bash EVscope.sh /path/to/EVscope -t 10 -o output --step [1] --input_fastq R1.fastq.gz R2.fastq.gz`
- **Dependency Issues**:
  - Verify tools: `conda list -n evscope_env`.
  - Activate environment: `conda activate evscope_env`.
- **Reference File Errors**:
  - Ensure files are in `references/` with correct paths.
  - Check permissions: `chmod -R u+rw references/`.
- **Memory Errors**:
  - Picard needs 250 GB RAM. Reduce threads or use a high-memory server.
- **Kraken2 Failure**:
  - Confirm `soft/kraken2/standard_krakendb/` exists.
- **Step Skipped**:
  - Delete `step.done`: `rm <output_dir>/<sample_name>_EVscope_output/Step_XX/step.done`.

## FAQ
**Q: Can I use non-SMARTer-seq data?**  
A: Yes, but adjust UMI extraction in `bin/Step_02_*.py`.

**Q: How do I run specific steps?**  
A: Use `--step`, e.g., `--step [1,2,3]`.

**Q: Why did a step fail?**  
A: Check logs in the step directory (e.g., `Step_03_cutadapt/`). Steps 13/14 are optional.

**Q: How do I view the HTML report?**  
A: Open `Step_18_HTML_report/<sample_name>_final_report.html` in a browser.

**Q: Can I use a different genome?**  
A: Requires updating `references/` with new genome files and indexes.

## Contributing
Contributions are welcome! To contribute:
1. Fork the repository.
2. Create a branch: `git checkout -b feature/YourFeature`.
3. Commit changes: `git commit -m 'Add YourFeature'`.
4. Push: `git push origin feature/YourFeature`.
5. Open a pull request.

Report issues or suggest features via [Issues](https://github.com/yourusername/EVscope/issues).

## Feedback
We value your input to improve EVscope! Share feedback, report bugs, or ask questions:
- **GitHub Issues**: [https://github.com/yourusername/EVscope/issues](https://github.com/yourusername/EVscope/issues)
- **Email**: your.email@example.com

## License
EVscope is licensed under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/). You are free to share and adapt the material, provided you give appropriate credit, provide a link to the license, and indicate if changes were made. See the [LICENSE](LICENSE) file for details.

## Contact
For support or inquiries:
- **Email**: your.email@example.com
- **GitHub**: [https://github.com/yourusername/EVscope](https://github.com/yourusername/EVscope)
