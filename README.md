---

# EVscope: A Comprehensive RNA-Seq Pipeline for Extracellular Vesicle Analysis

**EVscope** is an open-source, modular bioinformatics pipeline designed for the analysis of extracellular vesicle (EV) RNA sequencing data. Tailored to handle the unique challenges of EV RNA-seq—such as low RNA yield, fragmentation, diverse RNA biotypes, and contamination—EVscope processes paired-end or single-end FASTQ files through a robust workflow. It includes quality control, UMI-based deduplication, two-pass STAR alignment, circular RNA detection, expression quantification, contamination screening, tissue deconvolution, and comprehensive reporting. Optimized for the SMARTer Stranded Total RNA-Seq Kit v3 (Pico Input), EVscope introduces a novel expectation-maximization (EM) algorithm for multi-mapping read assignment and a unique read-through detection method for Read1 trimming.

![EVscope Pipeline Overview](figures/EVscope_pipeline.png)

---

## Table of Contents

- [Key Features](#key-features)
- [Motivation](#motivation)
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
- [Credits](#credits)
- [License](#license)
- [Contact](#contact)

---

## Key Features

- **Novel Read-Through Detection**: Trims UMI-derived adapter sequences from Read1 using reverse-complemented Read2 UMIs (`bin/Step_03_UMIAdapterTrimR1.py`).
- **EM-Based Multi-Mapping Resolution**: Assigns multi-mapped reads at single-base resolution using a genome-wide expectation-maximization algorithm (`bin/Step_25_EMapper.py`).
- **Comprehensive RNA Annotation**: Supports 3,659,642 RNAs across 20 biotypes (e.g., protein-coding, lncRNAs, miRNAs, piRNAs, retrotransposons) from GENCODE v45, piRBase v3.0, and RepeatMasker.
- **Dual circRNA Detection**: Integrates CIRCexplorer2 and CIRI2 for robust circular RNA identification, with merged results for enhanced sensitivity (`bin/Step_10_circRNA_merge.py`).
- **Tissue Deconvolution**: Infers EV RNA cellular origins using GTEx v10 and Human Brain Cell Atlas v1.0 references (`bin/Step_22_run_RNA_deconvolution_ARIC.py`).
- **Contamination Screening**: Filters bacterial (BBSplit) and microbial (Kraken2) contamination, with optional genomic DNA correction via strand-specific subtraction.
- **Extensive Quality Control**: Validates raw and trimmed FASTQs (FastQC), UMI motifs, and alignment metrics (`bin/Step_24_generate_QC_matrix.py`).
- **Expression Quantification**: Produces TPM/CPM matrices using featureCounts and RSEM, with RNA distribution visualizations (`bin/Step_15_plot_RNA_distribution_*.py`).
- **Interactive Reporting**: Generates bigWig tracks, density plots, and a comprehensive HTML report via R Markdown (`bin/Step_27_html_report.Rmd`).
- **Reproducibility**: Single-command Bash script with Conda environments, containerization support, and detailed logging for robust execution.

---

## Motivation

Extracellular vesicles (EVs) are critical mediators of intercellular communication, carrying diverse RNAs that serve as potential biomarkers for diseases like cancer and neurodegeneration. However, EV RNA sequencing is challenging due to low RNA abundance, fragmented transcripts, contamination (genomic DNA, bacterial RNA), and the presence of non-polyadenylated RNAs (e.g., miRNAs, lncRNAs). Standard RNA-seq pipelines are ill-suited for EV data, often failing to handle multi-mapping reads, diverse RNA biotypes, or contamination effectively.

EVscope addresses these challenges with a specialized, end-to-end pipeline optimized for EV total RNA-seq. It combines innovative algorithms, comprehensive RNA annotations, and robust quality control to deliver accurate and reproducible results, enabling biomarker discovery and advancing EV research.

---

## Directory Structure

```
EVscope/
├── EVscope.conf          # Configuration file for paths and tools
├── EVscope.sh            # Main pipeline script
├── README.md             # This documentation
├── bin/                  # Custom Python/R/shell scripts
│   ├── Step_02_calculate_ACC_motif_fraction.py
│   ├── Step_03_UMIAdapterTrimR1.py
│   ├── Step_07_bam2strand.py
│   ├── Step_10_circRNA_merge.py
│   ├── Step_15_plot_RNA_distribution_*.py
│   ├── Step_25_EMapper.py
│   ├── Step_27_html_report.Rmd
│   └── ... (other pipeline scripts)
├── example_data/         # Test datasets and run script
│   ├── Example_Data/
│   ├── Example_Data_EVscope_output/
│   └── run_EVscope.sh
├── figures/              # Pipeline overview image
│   └── EVscope_pipeline.png
├── references/           # Reference genomes, annotations, and indices
│   ├── annotations_HG38/
│   ├── deconvolution_HG38/
│   ├── genome/
│   └── index/
├── Resource/             # Additional resources and older pipeline versions
├── soft/                 # Bundled external tools
│   ├── bbmap/
│   ├── CIRI_v2.0.6/
│   ├── kraken2/
│   ├── KrakenTools/
│   └── RSEM_v1.3.3/
```

---

## Requirements

### Software
- **Operating System**: Linux (e.g., Ubuntu 20.04+) or macOS.
- **Bash**: Version 4.0 or higher.
- **Conda**: Miniconda or Anaconda for environment management.
- **Core Tools**:
  - FastQC (v0.12.1)
  - umi_tools (v1.1.5)
  - cutadapt (v4.9)
  - STAR (v2.7.11b)
  - samtools (v1.21)
  - featureCounts (v2.0.6)
  - CIRCexplorer2 (v2.3.8)
  - CIRI2 (v2.0.6, in `soft/`)
  - RSEM (v1.3.3, in `soft/`)
  - BBMap (v39.15, in `soft/`)
  - Kraken2 (in `soft/`)
  - KronaTools
  - ribodetector (v0.3.1)
  - seqtk (v1.4)
  - BWA (v0.7.18)
  - Picard (v3.3.0)
  - deepTools (v3.5.5)
  - R (v4.3.1) with packages: `rmarkdown`, `DT`, `kableExtra`, `bookdown`, `ggplot2`, `dplyr`
  - Python (v3.10.0) with packages: `pandas`, `numpy`, `matplotlib`, `biopython`, `numba`, `pyBigWig`, `pysam`

### Reference Files
- **Genomes**:
  - Human: `references/genome/hg38/hg38.p14.whole.genome.fa`
  - Mycoplasma: `references/genome/mycoplasma/mycoplasma_ge.fa`
  - E. coli: `references/genome/ecoli/E.coli_ge.fa`
- **Indices**:
  - STAR: `references/index/Index_STAR-2.7.10b_sjdbOverhang99_GeneCode45/`
  - BWA: `references/index/Index_BWA_V0.7.18/hg38.p14.whole.genome.fa`
  - RSEM: `references/index/RSEM_bowtie2_index/RSEM_REF_HG38_3659642_RNNAs/`
- **Annotations**:
  - GENCODE v45: `references/annotations_HG38/gencode.v45.chr_patch_hapl_scaff.annotation_UCSC.gtf`
  - Comprehensive RNA: `references/annotations_HG38/HG38_3659642_combined_RNAs_with_gold_standard_piRNAs.gtf`
  - SAF/BED files for genomic regions (e.g., `HG38_3UTR_noOverlap.saf`)
  - rRNA intervals: `references/annotations_HG38/GeneCodeV45_combined_Nucl_Mt_rRNA_589.interval_list`
- **Deconvolution References**:
  - GTEx v10: `references/deconvolution_HG38/GTEx_v10_SMTS_AVG_TPM.csv`
  - Brain Atlas: `references/deconvolution_HG38/human_brain_single_cell_atlasV1_31superclusters_AVP_CPM.csv`
- **Kraken2 Database**: `soft/kraken2/standard_krakendb/`

### Hardware
- **CPU**: 20+ threads recommended for optimal performance.
- **RAM**: Minimum 64 GB; 250 GB recommended for Picard tools.
- **Storage**: 500 GB+ for input data, references, and outputs.

---

## Installation

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/TheDongLab/EVscope.git
   cd EVscope
   ```

2. **Install Conda** (if not already installed):
   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   source ~/.bashrc
   ```

3. **Create Conda Environments**:
   Create three environments for EVscope, Picard, and Kraken2:
   ```bash
   conda env create -f environments/evscope_env.yml
   conda env create -f environments/picard_env.yml
   conda env create -f environments/kraken2_env.yml
   ```
   Example `evscope_env.yml`:
   ```yaml
   name: evscope_env
   channels:
     - bioconda
     - conda-forge
   dependencies:
     - fastqc=0.12.1
     - umi_tools=1.1.5
     - cutadapt=4.9
     - star=2.7.11b
     - samtools=1.21
     - subread=2.0.6
     - circos
     - ribodetector=0.3.1
     - seqtk=1.4
     - bwa=0.7.18
     - deeptools=3.5.5
     - python=3.10.0
     - pandas=2.2.3
     - numpy
     - matplotlib=3.9.1
     - biopython=1.78
     - numba=0.60.0
     - pybigwig=0.3.22
     - pysam
     - r-base=4.3.1
     - r-rmarkdown=2.29
     - r-dt=0.33
     - r-kableextra=1.4.0
     - r-bookdown=0.42
     - r-ggplot2=3.5.1
     - r-dplyr=1.1.4
   ```
   Create similar YAML files for `picard_env` (Picard, ImageMagick) and `kraken2_env` (Kraken2, KronaTools).

4. **Install CIRCexplorer2**:
   ```bash
   conda activate evscope_env
   pip install CIRCexplorer2==2.3.8
   ```

5. **Verify Bundled Tools**:
   Ensure bundled tools are present in `soft/`:
   ```bash
   ls soft/bbmap/bbsplit.sh
   ls soft/CIRI_v2.0.6/CIRI2.pl
   ls soft/RSEM_v1.3.3/rsem-calculate-expression
   ls soft/kraken2/krakendb
   ```

6. **Download Reference Files**:
   Organize genomes, indices, and annotations in `references/`:
   ```bash
   mkdir -p references/{genome,annotations_HG38,deconvolution_HG38,index}
   # Example: Download hg38 genome
   # wget -O references/genome/hg38/hg38.p14.whole.genome.fa <URL>
   ```
   Reference files are available at [EVscope GitHub Releases](https://github.com/TheDongLab/EVscope/releases).

7. **Set Up Kraken2 Database**:
   ```bash
   mkdir -p soft/kraken2/standard_krakendb
   conda activate kraken2_env
   kraken2-build --standard --db soft/kraken2/standard_krakendb --threads 20
   ```

8. **Test Installation**:
   Verify key tools:
   ```bash
   conda activate evscope_env
   fastqc --version
   STAR --version
   samtools --version
   python --version
   CIRCexplorer2 --version
   ```

---

## Usage

### Command Syntax
```bash
bash EVscope.sh --sample_name <name> --input_fastqs <files> [options]
```

**Required Arguments**:
- `--sample_name <name>`: Unique sample identifier (used for output files).
- `--input_fastqs <files>`: Comma-separated FASTQ file paths (e.g., `R1.fq.gz,R2.fq.gz` for paired-end).

**Optional Arguments**:
- `--threads <int>`: Number of CPU threads (default: 1, recommended: 20).
- `--run_steps <list>`: Steps to run (e.g., `1,3,5-8`, `all`; default: `all`).
- `--skip_steps <list>`: Steps to skip (e.g., `2,4`).
- `--circ_tool <tool>`: circRNA detection tool (`CIRCexplorer2`, `CIRI2`, `both`; default: `both`).
- `--read_count_mode <mode>`: Read counting strategy (`uniq` for featureCounts, `multi` for RSEM; default: `uniq`).
- `--gDNA_correction <yes|no>`: Apply genomic DNA correction (default: `no`).
- `--strandedness <strand>`: Library strandedness (`forward`, `reverse`, `unstrand`; default: `reverse`).
- `--config <path>`: Custom configuration file (default: `EVscope.conf`).
- `-V, --verbosity <level>`: Logging level (1=DEBUG, 2=INFO, 3=WARN, 4=ERROR; default: 2).
- `-h, --help`: Display help message.
- `-v, --version`: Show pipeline version.

### Example: Full Pipeline
Run all steps for a paired-end sample:
```bash
bash EVscope.sh --sample_name Example_Data \
    --input_fastqs example_data/Example_Data/chr21_2000_reads_R1_001.fastq.gz,example_data/Example_Data/chr21_2000_reads_R2_001.fastq.tgz \
    --threads 20 \
    --run_steps all \
    --gDNA_correction yes \
    --strandedness reverse
```

### Example: Specific Steps
Run steps 1-3 for quality control:
```bash
bash EVscope.sh --sample_name Test_Sample \
    --input_fastqs test_R1.fq.gz,test_R2.fq.gz \
    --threads 10 \
    --run_steps 1-3
```

### Test Run with Example Data
Use the provided test dataset:
```bash
cd example_data
bash run_EVscope.sh
```

### Configuration File
Edit `EVscope.conf` to specify paths to tools, genomes, indices, and annotations. Example:
```bash
EVscope_PATH="/path/to/EVscope"
HUMAN_GENOME_FASTA="${EVscope_PATH}/references/genome/hg38/hg38.p14.whole.genome.fa"
STAR_INDEX="${EVscope_PATH}/references/index/Index_STAR-2.7.10b_sjdbOverhang99_GeneCode45"
```

> **Note**: Ensure all paths in `EVscope.conf` are absolute and files exist. The pipeline validates essential references before execution.

---

## Input Data Format

- **FASTQ Files**: Gzipped, paired-end (`R1.fastq.gz`, `R2.fastq.gz`) or single-end (`R1.fastq.gz`).
- **Sequencing Protocol**: Optimized for SMARTer Stranded Total RNA-Seq Kit v3 (Pico Input) with 14-bp UMIs in Read2.
- **Naming Convention**: R1 and R2 files should share a common prefix (e.g., `Sample_001_R1.fastq.gz`, `Sample_001_R2.fastq.gz`).
- **Quality**: High-quality reads suitable for EV RNA-seq.

> **Warning**: For non-SMARTer-seq data, adjust UMI parameters in `bin/Step_02_*.py` and `bin/Step_03_UMIAdapterTrimR1.py`.

---

## Pipeline Steps

EVscope comprises 27 modular steps for comprehensive EV RNA-seq analysis:

1. **Raw FASTQ QC**: Generates FastQC reports for input FASTQs (`Step_01_Raw_QC`).
2. **UMI Motif Analysis**: Visualizes UMI motifs and calculates ACC motif fractions (`Step_02_UMI_Analysis`).
3. **UMI Extraction & Trimming**: Extracts UMIs, trims adapters, and removes read-through UMIs (`Step_03_UMI_Adaptor_Trim`).
4. **Trimmed FASTQ QC**: FastQC on trimmed FASTQs (`Step_04_Trimmed_QC`).
5. **Bacterial Contamination Screening**: Filters E. coli and Mycoplasma reads using BBSplit (`Step_05_Bacterial_Filter`).
6. **Two-Pass STAR Alignment**: Performs initial and refined alignments with UMI deduplication (`Step_06_Alignment_Initial`, `Step_06_Alignment_Refined`).
7. **Strandedness Detection**: Determines library strandedness and detects gDNA (`Step_07_Strand_Detection`).
8. **CIRCexplorer2 circRNA Detection**: Identifies circular RNAs from STAR junctions (`Step_08_CIRCexplorer2_circRNA`).
9. **CIRI2 circRNA Detection**: Detects circular RNAs using BWA alignments (`Step_09_CIRI2_circRNA`).
10. **circRNA Merging**: Combines CIRCexplorer2 and CIRI2 results (`Step_10_circRNA_Merge`).
11. **RNA-Seq Metrics**: Collects alignment and insert size metrics using Picard (`Step_11_RNA_Metrics`).
12. **featureCounts Quantification**: Quantifies unique-mapping reads (`Step_12_featureCounts_Quant`).
13. **gDNA-Corrected Quantification**: Applies genomic DNA correction to featureCounts (`Step_13_gDNA_Corrected_Quant`).
14. **RSEM Quantification**: Quantifies multi-mapping reads (`Step_14_RSEM_Quant`).
15. **featureCounts Expression Analysis**: Generates TPM matrices and RNA distribution plots (`Step_15_featureCounts_Expression`).
16. **gDNA-Corrected Expression Analysis**: Expression analysis with gDNA correction (`Step_16_gDNA_Corrected_Expression`).
17. **RSEM Expression Analysis**: Expression analysis using RSEM results (`Step_17_RSEM_Expression`).
18. **Genomic Region Mapping**: Quantifies reads in 3'UTR, 5'UTR, introns, etc. (`Step_18_Genomic_Regions`).
19. **Taxonomic Classification**: Classifies reads using Kraken2 (`Step_19_Taxonomy`).
20. **featureCounts Deconvolution**: Infers tissue origins for featureCounts results (`Step_20_featureCounts_Deconvolution`).
21. **gDNA-Corrected Deconvolution**: Deconvolution with gDNA correction (`Step_21_gDNA_Corrected_Deconvolution`).
22. **RSEM Deconvolution**: Deconvolution for RSEM results (`Step_22_RSEM_Deconvolution`).
23. **rRNA Detection**: Identifies rRNA reads using ribodetector (`Step_23_rRNA_Detection`).
24. **QC Summary**: Compiles comprehensive quality control metrics (`Step_24_QC_Summary`).
25. **Coverage Analysis**: Generates bigWig tracks using EMapper (`Step_25_EMapper_BigWig_Quantification`).
26. **Density Plots**: Visualizes read coverage across RNA types and meta-gene regions (`Step_26_BigWig_Density_Plot`).
27. **HTML Report**: Produces an interactive final report (`Step_27_HTML_Report`).

> **Note**: Steps 13, 16, 19, 20, 21, 22, and 23 are optional and fail gracefully if inputs are unavailable or conditions are not met (e.g., `--gDNA_correction=no` skips steps 13, 16, 21).

---

## Output Structure

The pipeline generates a structured output directory:

```
<sample_name>_EVscope_output/
├── EVscope_pipeline.log                # Main pipeline log
├── EVscope_pipeline.png                # Pipeline overview
├── Step_01_Raw_QC/                     # FastQC reports for raw FASTQs
├── Step_02_UMI_Analysis/               # UMI motif plots and ACC fraction
├── Step_03_UMI_Adaptor_Trim/           # Trimmed FASTQs and logs
├── Step_04_Trimmed_QC/                 # FastQC reports for trimmed FASTQs
├── Step_05_Bacterial_Filter/           # Contamination-filtered FASTQs
├── Step_06_Alignment_Initial/          # Initial STAR alignment
├── Step_06_Alignment_Refined/          # Final BAM files
├── Step_07_Strand_Detection/           # Strandedness results
├── Step_08_CIRCexplorer2_circRNA/      # CIRCexplorer2 circRNA results
├── Step_09_CIRI2_circRNA/              # CIRI2 circRNA results
├── Step_10_circRNA_Merge/              # Merged circRNA results
├── Step_11_RNA_Metrics/                # Picard metrics
├── Step_12_featureCounts_Quant/        # featureCounts results
├── Step_13_gDNA_Corrected_Quant/       # gDNA-corrected counts
├── Step_14_RSEM_Quant/                 # RSEM quantification
├── Step_15_featureCounts_Expression/   # TPM matrices and plots
├── Step_16_gDNA_Corrected_Expression/  # gDNA-corrected expression
├── Step_17_RSEM_Expression/            # RSEM expression matrices
├── Step_18_Genomic_Regions/            # Genomic region counts
├── Step_19_Taxonomy/                   # Kraken2 classification
├── Step_20_featureCounts_Deconvolution/ # Tissue deconvolution
├── Step_21_gDNA_Corrected_Deconvolution/ # gDNA-corrected deconvolution
├── Step_22_RSEM_Deconvolution/         # RSEM deconvolution
├── Step_23_rRNA_Detection/             # rRNA detection results
├── Step_24_QC_Summary/                 # QC metrics matrix
├── Step_25_EMapper_BigWig_Quantification/ # bigWig tracks
├── Step_26_BigWig_Density_Plot/        # Density plots
└── Step_27_HTML_Report/                # Final HTML report
```

**Key Outputs**:
- **Quality Control**: FastQC reports, UMI motif plots, QC matrix (`Step_24_QC_Summary/<sample_name>_QC_matrix.tsv`).
- **FASTQ**: Trimmed reads (`Step_03_UMI_Adaptor_Trim/`).
- **Alignments**: BAM files (`Step_06_Alignment_Refined/`).
- **circRNAs**: Merged CIRCexplorer2/CIRI2 results and Venn diagrams (`Step_10_circRNA_Merge/`).
- **Expression**: TPM/CPM matrices and RNA distribution plots (`Step_15_*/`, `Step_16_*/`, `Step_17_*/`).
- **Deconvolution**: Tissue/cell type estimates (`Step_20_*/`, `Step_21_*/`, `Step_22_*/`).
- **Coverage**: bigWig tracks and density plots (`Step_25_*/`, `Step_26_*/`).
- **Report**: Interactive HTML report (`Step_27_HTML_Report/<sample_name>_final_report.html`).

---

## Troubleshooting

- **Command Syntax Error**:
  - Ensure correct argument order: `--sample_name`, `--input_fastqs`, etc.
  - Quote paths with spaces: `"path with spaces/fastq.gz"`.
- **Dependency Not Found**:
  - Verify Conda environments: `conda list -n evscope_env`.
  - Activate environment: `conda activate evscope_env`.
- **Reference File Missing**:
  - Check `EVscope.conf` paths and file existence.
  - Set permissions: `chmod -R u+rw references/`.
- **Memory Issues**:
  - Picard requires up to 250 GB RAM. Reduce `--threads` or use a high-memory server.
- **Kraken2 Failure**:
  - Ensure `soft/kraken2/standard_krakendb/` exists and is built.
- **Step Skipped Unexpectedly**:
  - Check for `step.done` files: `rm Step_XX/step.done` to rerun.
- **Step Failure**:
  - Review logs: `<output_dir>/Step_XX/step.stderr.log` or `EVscope_pipeline.log`.
  - Steps 19, 20, 21, 22, and 23 are optional and may fail gracefully.

---

## FAQ

**Q: Can EVscope process non-SMARTer-seq data?**  
A: Yes, but you must modify UMI parameters in `bin/Step_02_*.py` and `bin/Step_03_UMIAdapterTrimR1.py`.

**Q: How do I run specific pipeline steps?**  
A: Use `--run_steps`, e.g., `--run_steps 1,3,5-8`.

**Q: Why did an optional step fail?**  
A: Steps like Kraken2 or deconvolution may fail if no relevant reads are detected or references are missing. Check logs for details.

**Q: How do I view the final report?**  
A: Open `Step_27_HTML_Report/<sample_name>_final_report.html` in a web browser.

**Q: Can I use a different genome?**  
A: Yes, update `references/` with new genomes, indices, and annotations, and modify `EVscope.conf`.

**Q: How do I increase pipeline speed?**  
A: Increase `--threads` (e.g., 20) and ensure sufficient RAM (250 GB for Picard).

---

## Contributing

We welcome contributions to improve EVscope! To contribute:
1. Fork the repository: [https://github.com/TheDongLab/EVscope](https://github.com/TheDongLab/EVscope).
2. Create a feature branch: `git checkout -b feature/YourFeature`.
3. Commit changes: `git commit -m 'Add YourFeature'`.
4. Push to your fork: `git push origin feature/YourFeature`.
5. Submit a pull request.

Please report bugs or suggest features via [GitHub Issues](https://github.com/TheDongLab/EVscope/issues).

---

## Feedback

Your feedback helps us improve EVscope! Share suggestions, issues, or questions:
- **GitHub Issues**: [https://github.com/TheDongLab/EVscope/issues](https://github.com/TheDongLab/EVscope/issues)
- **Email**: [xianjun.dong@yale.edu](mailto:xianjun.dong@yale.edu)

---

## Credits

**Authors**:
- **Yiyong Zhao**: Data curation, Formal analysis, Software, Visualization
- **Himanshu Chintalapudi**: Visualization
- **Ziqian Xu**: Resources
- **Weiqiang Liu**: Data curation
- **Yuxuan Hu**: Validation
- **Ewa Grassin, Minsun Song, SoonGweon Hong, Luke P. Lee**: Resources
- **Xianjun Dong**: Conceptualization, Methodology, Funding, Supervision

**Affiliations**:
1. Department of Neurology and Biomedical Informatics and Data Science, Yale School of Medicine, Yale University, New Haven, CT, USA
2. Aligning Science Across Parkinson’s (ASAP) Collaborative Research Network, Chevy Chase, MD, USA
3. Department of Medicine, Brigham and Women’s Hospital, Harvard Medical School, Boston, MA, USA
4. Department of Bioengineering, University of California at Berkeley, Berkeley, CA, USA

**Funding**:
- NIH Grants: 1R01NS124916, 1R24NS132738
- Aligning Science Across Parkinson’s: ASAP-000301, ASAP-000529 (via Michael J. Fox Foundation)
- Dong Lab computational resources

**Corresponding Author**: Xianjun Dong ([xianjun.dong@yale.edu](mailto:xianjun.dong@yale.edu))

**Citation**:
> Zhao, Y., et al. (2025). EVscope: A Comprehensive Pipeline for Extracellular Vesicle RNA-Seq Analysis. *In preparation*. DOI: [10.5281/zenodo.15577789](https://doi.org/10.5281/zenodo.15577789)

---

## License

EVscope is licensed under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/). You are free to share and adapt the material, provided you give appropriate credit, include a link to the license, and indicate any changes made. See [LICENSE](LICENSE) for details.

---

## Contact

For inquiries, contact:
- **Xianjun Dong**: [xianjun.dong@yale.edu](mailto:xianjun.dong@yale.edu)
- **GitHub**: [https://github.com/TheDongLab/EVscope](https://github.com/TheDongLab/EVscope)

---

### Notes on Updates
- **Pipeline Steps**: Updated to reflect the 27 steps in the new `EVscope.sh`, with detailed descriptions.
- **Output Structure**: Aligned with the new output directories and key files.
- **Usage**: Rewritten to match the new command-line interface, emphasizing `--sample_name` and `--input_fastqs`.
- **Configuration**: Highlighted the importance of `EVscope.conf` and provided an example.
- **Installation**: Streamlined instructions with specific tool versions and Conda environment setup.
- **Visuals**: Maintained the pipeline overview image and added clarity to section headings.
- **Professional Tone**: Enhanced readability with concise, clear language and consistent formatting.
- **Reproducibility**: Emphasized Conda environments, reference file setup, and test data usage.

This README is ready for GitHub, providing a professional and user-friendly guide to EVscope. Let me know if you need further refinements or additional sections!
