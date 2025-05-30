EVscope Pipeline
EVscope is a state-of-the-art, open-source bioinformatics pipeline for comprehensive analysis of extracellular vesicle (EV) RNA sequencing data. Designed to address the unique challenges of EV RNA-seq—low RNA abundance, fragmentation, and contamination—it processes paired-end FASTQ files through quality control, UMI-based deduplication, alignment, circular RNA detection, expression quantification, taxonomic classification, rRNA detection, and cellular origin inference. Optimized for the SMARTer Stranded Total RNA-Seq Kit v3 (Pico Input), EVscope introduces a novel genome-wide expectation-maximization (EM) algorithm for multi-mapping read assignment and a unique read-through detection method for Read1 trimming. It supports 20 RNA biotype annotations and integrates tools like FastQC, STAR, CIRCexplorer2, CIRI2, RSEM, and Kraken2, delivering robust results via an HTML report.

Table of Contents

Motivation
Features
Directory Structure
Requirements
Installation
Usage
Input Data Format
Pipeline Steps
Output Structure
Troubleshooting
FAQ
Contributing
Feedback
Credits
License
Contact

Motivation
Extracellular vesicles (EVs) are nanosized, membrane-bound structures that mediate intercellular communication through RNA transfer, making them promising biomarkers for diseases like cancer and neurodegeneration. However, EV RNA sequencing faces significant challenges: low RNA yield, fragmented transcripts, diverse tissue origins, and contamination from genomic DNA, bacterial RNA, and co-isolated proteins. Standard RNA-seq pipelines, designed for cellular RNA, fail to address these issues, often missing non-polyadenylated RNAs (e.g., miRNAs, lncRNAs) and producing unreliable results due to multi-mapping reads and contamination.
EVscope addresses these challenges with a specialized, end-to-end pipeline tailored for EV total RNA-seq. It introduces innovative methods like a genome-wide EM algorithm for precise multi-mapping read quantification, a novel read-through detection for Read1 UMI trimming, and comprehensive annotation of 20 RNA biotypes. By integrating robust quality control, contamination filtering, dual circRNA detection, and cellular origin inference, EVscope standardizes EV RNA-seq analysis, enhancing reproducibility and enabling biomarker discovery.
Features

Novel Read-Through Detection: Unique method trims UMI-derived adapter sequences from Read1 caused by short RNA inserts, using reverse-complemented Read2 UMIs (bin/Step_03_UMIAdapterTrimR1.py).
EM Algorithm for Multi-Mapping Reads: Assigns multi-mapped reads at single-base resolution with an expectation-maximization approach, improving expression quantification accuracy.
Comprehensive RNA Annotation: Supports 3,659,642 RNAs across 20 categories (e.g., protein-coding, lncRNAs, miRNAs, piRNAs, retrotransposons) from GENCODE v45, piRBase, and RepeatMasker.
Dual circRNA Detection: Combines CIRCexplorer2 and CIRI2 for sensitive and accurate circular RNA identification, with merged results for robustness.
Cellular Origin Inference: Deconvolves EV RNA sources using GTEx v10 and Human Brain Cell Atlas v1.0 references, estimating tissue/cell type contributions.
Contamination Filtering: Detects bacterial (BBSplit) and microbial (Kraken2) contamination, with gDNA correction via strand-specific subtraction.
Quality Control: Validates raw and trimmed FASTQ (FastQC) and UMI motifs (bin/Step_02_plot_fastq2UMI_motif.py).
Expression Quantification: Generates TPM/CPM matrices with featureCounts and RSEM, including RNA distribution plots.
Visualization: Produces bigWig tracks, density plots, and an interactive HTML report with R Markdown.
Reproducibility: Single-command Bash script, containerized with Singularity, and tested on diverse EV RNA-seq datasets.

Directory Structure
EVscope_pipeline/
├── EVscope.def           # Pipeline definition file
├── EVscope.sh            # Main pipeline script
├── EVscope_v1.sh         # Alternative pipeline script
├── README.md             # Documentation
├── bin/                  # Custom scripts
│   ├── Step_02_calculate_ACC_motif_ratio.py
│   ├── Step_02_plot_fastq2UMI_motif.py
│   ├── Step_03_UMIAdapterTrimR1.py
│   ├── Step_07_bam2strand.py
│   ├── Step_08_combined_CIRCexplorer2_CIRI2.py
│   ├── Step_08_convert_CIRCexplorer2CPM.py
│   ├── Step_08_convert_CIRI2CPM.py
│   ├── Step_11_combine_total_RNA_expr_matrix.py
│   ├── Step_11_featureCounts2TPM.py
│   ├── Step_11_plot_RNA_distribution_1subplot.py
│   ├── Step_11_plot_RNA_distribution_20subplots.py
│   ├── Step_11_plot_RNA_distribution_2subplots.py
│   ├── Step_11_RSEM2expr_matrix.py
│   ├── Step_12_plot_reads_mapping_stats.py
│   ├── Step_14_run_RNA_deconvolution_ARIC.py
│   ├── Step_16_generate_QC_matrix.py
│   ├── Step_17_density_plot_over_meta_gene.sh
│   ├── Step_17_density_plot_over_RNA_types.sh
│   ├── Step_17_EMapper.py
│   └── Step_18_html_report.Rmd
├── config/               # Configuration files (empty)
├── references/           # Reference files
│   ├── annotations_HG38/
│   ├── deconvolution_HG38/
│   ├── genome/
│   └── index/
├── soft/                 # External tools
│   ├── bbmap/
│   ├── CIRI_v2.0.6/
│   ├── kraken2/
│   └── RSEM_v1.3.3/
├── TEMP/                 # Temporary files (empty)
└── test_data/            # Test datasets
    ├── DNAnexus_CTRL_chr21_fastq/
    ├── small_data/
    ├── test/
    ├── processing.log
    └── run.sh

Requirements
Software

OS: Linux (e.g., Ubuntu) or macOS terminal.
Bash: v4.0+.
Conda: Miniconda/Anaconda.
Tools:
FastQC
MultiQC
umi_tools
cutadapt
trim_galore
ribodetector_cpu
seqtk
kraken2
KronaTools (kreport2krona.py)
STAR (v2.7.10b)
samtools
deepTools (bamCoverage, computeMatrix, plotProfile)
featureCounts (subread)
CIRCexplorer2
CIRI2 (v2.0.6, in soft/)
BBMap (in soft/)
Picard
RSEM (v1.3.3, in soft/)
R (with rmarkdown, DT, kableExtra, bookdown, here)
Python 3.8+ (with pandas, numpy, matplotlib, seaborn, biopython, numba, concurrent.futures)



Reference Files

Genomes:
genome/hg38/hg38.p14.whole.genome.fa
genome/mycoplasma/mycoplasma_ge.fa
genome/ecoli/E.coli_ge.fa


Indexes:
index/RSEM_bowtie2_index/RSEM_REF_HG38_3659642_RNNAs/
index/Index_STAR-2.7.10b_sjdbOverhang99_GeneCode45/
index/Index_BWA_V0.7.18/hg38.p14.whole.genome.fa


Annotations:
annotations_HG38/HG38_3659642_combined_RNAs_with_gold_standard_piRNAs.gtf
annotations_HG38/gencode.v45.chr_patch_hapl_scaff.annotation_UCSC.gtf
SAF/BED: HG38_3UTR_noOverlap.saf, HG38_5UTR_noOverlap.saf, etc.
annotations_HG38/GeneCodeV45_combined_Nucl_Mt_rRNA_589.interval_list


Deconvolution:
deconvolution_HG38/GTEx_v10_SMTS_AVG_TPM.csv
deconvolution_HG38/GTEx_v10_SMTSD_AVG_TPM.csv
deconvolution_HG38/human_brain_single_cell_atlasV1_31superclusters_AVP_CPM.csv


Kraken2 Database: soft/kraken2/standard_krakendb/
Bacterial References: 240 Mycoplasma strains, E. coli (GCF_000005845.2_ASM584) in references/.

Hardware

CPU: 20+ threads recommended.
RAM: 64 GB minimum (Picard requires up to 250 GB).
Storage: 500 GB+ for inputs, references, and outputs.

Installation

Clone Repository:
git clone https://github.com/TheDongLab/EVscope.git
cd EVscope


Install Conda:
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc


Create Conda Environments:
conda env create -f environments/evscope_env.yml
conda env create -f environments/picard_env.yml
conda env create -f environments/kraken2_env.yml

Example environments/evscope_env.yml:
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
  - biopython
  - numba
  - r-base
  - r-rmarkdown
  - r-dt
  - r-kableextra
  - r-bookdown
  - r-here

Create similar YAML files for picard_env (Picard, ImageMagick) and kraken2_env (Kraken2, KronaTools).

Install CIRCexplorer2:
conda activate evscope_env
pip install CIRCexplorer2


Verify Bundled Tools:Check soft/:
ls soft/bbmap/bbsplit.sh
ls soft/CIRI_v2.0.6/CIRI2.pl
ls soft/RSEM_v1.3.3/rsem-calculate-expression


Set Up References:Download and organize in references/:
mkdir -p references/genome/hg38
# Example: wget -O references/genome/hg38/hg38.p14.whole.genome.fa <url>

Download annotations and deconvolution files from GitHub.

Kraken2 Database:
mkdir -p soft/kraken2/standard_krakendb
conda activate kraken2_env
kraken2-build --standard --db soft/kraken2/standard_krakendb --threads 20


Test Installation:
conda activate evscope_env
fastqc --version
STAR --version
samtools --version
python --version
CIRCexplorer2 --version



Usage
Command Syntax
bash EVscope.sh <EVscope_path> -t <threads> -o <output_dir> --step <step_list> --input_fastq <R1_fastq.gz> <R2_fastq.gz>


<EVscope_path>: Path to EVscope directory (e.g., /home/yz2474/yiyong_2023/EVscope_pipeline).
-t <threads>: CPU threads (default: 1, recommended: 10-20).
-o <output_dir>: Output directory.
--step <step_list>: Steps to run, e.g., [1,2,3] (default: [1-18]).
--input_fastq <R1_fastq.gz> <R2_fastq.gz>: Paired-end FASTQ files.

Example
Run steps 1-3:
bash EVscope.sh /home/yz2474/yiyong_2023/EVscope_pipeline -t 10 -o test_output --step [1,2,3] --input_fastq test_data/DNAnexus_CTRL_chr21_fastq/chr21_R1_001.fastq.gz test_data/DNAnexus_CTRL_chr21_fastq/chr21_R2_001.fastq.gz

Test Run
Use test_data/:
cd test_data
bash run.sh

Or:
bash ../EVscope.sh $(pwd)/.. -t 4 -o test --step [1,2] --input_fastq small_data/test_R1.fastq.gz small_data/test_R2.fastq.gz

Notes

Sample Name: Derived from R1 FASTQ (e.g., chr21 from chr21_R1_001.fastq.gz).
Paths with Spaces: Quote, e.g., "/path with spaces/fastq.gz".
Steps: Skipped if step.done exists in step directory.
Tutorials: Available at GitHub.

Input Data Format

FASTQ: Paired-end, gzipped (R1.fastq.gz, R2.fastq.gz).
Protocol: SMARTer Stranded Total RNA-Seq Kit v3 (Pico Input) with 14-bp UMIs in Read2.
Naming: Shared prefix for R1/R2 (e.g., Sample_007_R1.fastq.gz).
Quality: High-quality reads for EV RNA-seq.


Warning: Non-SMARTer-seq data requires adjusting UMI parameters in bin/Step_02_*.py and bin/Step_03_UMIAdapterTrimR1.py.

Pipeline Steps
EVscope’s 18 steps process EV RNA-seq data comprehensively (see Figure 1a in manuscript):

Raw FASTQ QC: FastQC and UMI motif visualization (Step_02_plot_fastq2UMI_motif.py, Figure S3).
UMI Extraction & Trimming: Extracts 14-bp UMIs (umi_tools), trims adapters (cutadapt), and removes Read1 read-through UMIs (Step_03_UMIAdapterTrimR1.py).
Post-Trimming QC: FastQC on trimmed FASTQ.
Bacterial Detection: Filters E. coli/Mycoplasma reads (BBSplit, Figure 1b).
STAR Alignment & Deduplication: Two-pass STAR alignment, UMI deduplication, and re-alignment.
Strand Detection: Assesses RNA strandness, detects gDNA (Step_07_bam2strand.py, Figure 1c).
circRNA Detection (CIRCexplorer2): Identifies circRNAs from STAR junctions.
circRNA Detection (BWA & CIRI2): Detects circRNAs, merges with CIRCexplorer2 (Step_08_combined_CIRCexplorer2_CIRI2.py).
Picard Metrics: RNA-seq and insert size metrics.
Expression Quantification: Counts reads (featureCounts, RSEM) with gDNA correction.
Expression Matrix & Plots: TPM/CPM matrices, RNA distribution plots (Step_11_plot_RNA_distribution_*.py, Figure 1d).
Meta-Feature Mapping: Quantifies reads in genomic regions (e.g., 3'UTR).
Taxonomic Classification: Downsamples, classifies reads (Kraken2, Figure 1b).
Deconvolution: Infers cell/tissue origins (ARIC, GTEx, brain atlas, Figure 1h).
rRNA Detection: Identifies rRNA (RiboDetector).
QC Matrix: Compiles metrics (Step_16_generate_QC_matrix.py).
bigWig Generation: Coverage tracks with EM algorithm (Step_17_density_*.sh, Figure 1e).
HTML Report: Interactive summary (Step_18_html_report.Rmd).


Note: Steps 13/14 are optional and fail gracefully.

Output Structure
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

Key Outputs:

QC: FastQC reports, UMI motif plots.
FASTQ: Trimmed reads.
BAM: Aligned, deduplicated files.
circRNA: TSVs, Venn diagrams.
Expression: TPM/CPM matrices, RNA plots.
Deconvolution: Cell type estimates.
QC Matrix: Metrics in Step_16_QC_matrix_generation/.
Coverage: bigWig tracks, density plots.
Report: HTML in Step_18_HTML_report/.

Troubleshooting

Syntax Error:
Ensure correct order: <EVscope_path> -t <threads> -o <output> --step <list> --input_fastq <R1> <R2>.


Dependency Issues:
Check: conda list -n evscope_env.
Activate: conda activate evscope_env.


Reference Errors:
Verify references/ paths.
Set permissions: chmod -R u+rw references/.


Memory Errors:
Picard needs 250 GB. Reduce -t or use high-memory server.


Kraken2 Failure:
Ensure soft/kraken2/standard_krakendb/ exists.


Step Skipped:
Remove step.done: rm <output_dir>/<sample_name>_EVscope_output/Step_XX/step.done.



FAQ
Q: Can I use non-SMARTer-seq data?A: Yes, modify UMI parameters in bin/Step_02_*.py and bin/Step_03_UMIAdapterTrimR1.py.
Q: How do I run specific steps?A: Use --step, e.g., --step [1,2,3].
Q: Why did a step fail?A: Check logs (e.g., Step_03_cutadapt/). Steps 13/14 are optional.
Q: How to view the HTML report?A: Open Step_18_HTML_report/<sample_name>_final_report.html in a browser.
Q: Can I use another genome?A: Update references/ with new genome/indexes.
Contributing
Contributions welcome! To contribute:

Fork: https://github.com/TheDongLab/EVscope.
Branch: git checkout -b feature/YourFeature.
Commit: git commit -m 'Add YourFeature'.
Push: git push origin feature/YourFeature.
Pull request.

Report issues at Issues.
Feedback
We value feedback to enhance EVscope! Share suggestions or issues:

GitHub Issues: https://github.com/TheDongLab/EVscope/issues
Email: xianjun.dong@yale.edu

Credits
Authors:

Yiyong Zhao (Data curation, Formal analysis, Software, Visualization)
Himanshu Chintalapudi (Visualization)
Ziqian Xu (Resources)
Weiqiang Liu (Data curation)
Yuxuan Hu (Validation)
Ewa Grassin, Minsun Song, SoonGweon Hong, Luke P. Lee (Resources)
Xianjun Dong (Conceptualization, Methodology, Funding, Supervision)

Affiliations:

Department of Neurology and Biomedical Informatics and Data Science, Yale School of Medicine, Yale University, New Haven, CT, USA
Aligning Science Across Parkinson’s (ASAP) Collaborative Research Network, Chevy Chase, MD, USA
Department of Medicine, Brigham and Women’s Hospital, Harvard Medical School, Boston, MA, USA
Department of Bioengineering, University of California at Berkeley, Berkeley, CA, USA

Funding:

NIH grants 1R01NS124916, 1R24NS132738
Aligning Science Across Parkinson’s [ASAP-000301, ASAP-000529] via Michael J. Fox Foundation
Dong Lab computational resources

Corresponding Author: Xianjun Dong (xianjun.dong@yale.edu)
License
EVscope is licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0). You may share and adapt the material with appropriate credit, a license link, and indication of changes. See LICENSE.
