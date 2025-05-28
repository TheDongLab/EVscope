#!/usr/bin/env bash
# Usage: bash run_SMARTer_pipe_v3_ZYY.sh <SampleName> <R1.fq.gz> <R2.fq.gz>
# Example: bash run_SMARTer_pipe_v3_ZYY.sh B4_TDP43_KD_caRNA_rep3 /path/to/R1.fq.gz /path/to/R2.fq.gz

set -eo pipefail
export TMPDIR=/home/yz2474/tmp
# Input parameters
SampleName="$1"
R1_path="$2"
R2_path="$3"

# Define paths to reference files and indexes
genomeDir="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Index/Index_STAR-2.7.10b_sjdbOverhang100_Gene"
combined_gene_gtf="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Annotation/refs_EV-RNA-Profiler/UCSC_HG38/combine/11711764_HG38_gene.gtf"
combined_gene_refflat="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Annotation/refs_EV-RNA-Profiler/UCSC_HG38/combine/11711764_HG38_gene.refflat"
genecodeV45_gene_gtf="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Annotation/refs_EV-RNA-Profiler/UCSC_HG38/combine/gencode.v45.chr_patch_hapl_scaff.annotation_UCSC.gtf"
genome_fa="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Sequence/hg38p14/hg38.p14.whole.genome.fa"
bwa_index_ref="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Index/Index_BWA_V0.7.18/hg38.p14.whole.genome.fa"
rRNA_interval="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Annotation/refs_EV-RNA-Profiler/UCSC_HG38/combine/rRNA/GeneCodeV45_combined_Nucl_Mt_rRNA_589.interval_list"
ref_myco="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Annotation/refs_EV-RNA-Profiler/Bacterial/240_Mycoplasma_ge/Mycoplasma_240_ge.fa"
ref_ecoli="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Annotation/refs_EV-RNA-Profiler/Bacterial/E.coli/E.coli_GCF_000005845.2_ASM584v2_ge.fa"
# Define the refFlat file (adjust the path as needed)

# Function to run each step; if the step has already been completed (step.done exists), it skips it.
run_step() {
    local step_dir="$1"
    local step_cmds="$2"
    mkdir -p "$step_dir"
    if [ ! -f "$step_dir/step.done" ]; then
        start_time=$(date +%s)
        eval "$step_cmds"
        echo "======> sart run: $step_cmds"
        ret=$?
        end_time=$(date +%s)
        elapsed=$((end_time - start_time))
        if [ $ret -eq 0 ]; then
            touch "$step_dir/step.done"
            echo -e "======> Finished this step: $step_dir in $elapsed seconds"
        else
            echo "Error: Step failed in $step_dir after $elapsed seconds. Check logs for details."
            exit 1
        fi
    else
        echo "======> already Finished this step: $step_dir"
    fi
}

# Create the output directory (named after the sample)
outputDir="${SampleName}"
mkdir -p "${outputDir}"

########################################
# Step 1: Raw FastQC QC
########################################
run_step "${outputDir}/01_raw_fastqc" "
fastqc -o ${outputDir}/01_raw_fastqc -t 20 $R1_path $R2_path
"

########################################
# Step 2: Extract 14-bp UMI motif from Read2
########################################
run_step "${outputDir}/02_umi_motif" "
python /home/yz2474/scripts/donglab/get_SMARTer_Read2_14BP_motif/get_umi_seq_motif_for_individual_fq.py \
  -head ${SampleName} \
  -fq $R2_path \
  -n 14 -r 100000 \
  -o ${outputDir}/02_umi_motif \
  -p 20
"

########################################
# Step 3: UMI labeling for Read1/2 and adapter trimming
########################################
run_step "${outputDir}/03_trim_galore" "
umi_tools extract --bc-pattern='NNNNNNNN' \
  --stdin=$R2_path \
  --stdout=${outputDir}/03_trim_galore/${SampleName}_R2_umi.fq.gz \
  --read2-in=$R1_path \
  --read2-out=${outputDir}/03_trim_galore/${SampleName}_R1_umi.fq.gz \
  --log=${outputDir}/03_trim_galore/umi_extract.log \
  --umi-separator='_'
cutadapt -u 14 -o ${outputDir}/03_trim_galore/${SampleName}_R2_umi_trim.fq.gz \
  ${outputDir}/03_trim_galore/${SampleName}_R2_umi.fq.gz
trim_galore --phred33 --illumina --gzip --stringency 3 --cores 10 --paired --length 10 \
  ${outputDir}/03_trim_galore/${SampleName}_R1_umi.fq.gz \
  ${outputDir}/03_trim_galore/${SampleName}_R2_umi_trim.fq.gz \
  -o ${outputDir}/03_trim_galore
mv ${outputDir}/03_trim_galore/${SampleName}_R1_umi_val_1.fq.gz ${outputDir}/03_trim_galore/${SampleName}_R1_trimed.fq.gz
mv ${outputDir}/03_trim_galorem/${SampleName}_R2_umi_trim_val_2.fq.gz ${outputDir}/03_trim_galore/${SampleName}_R2_trimed.fq.gz
"

########################################
# Step 4: QC after trim_galore
########################################
run_step "${outputDir}/04_trim_galore_fastqc" "
fastqc -o ${outputDir}/04_trim_galore_fastqc -t 20 \
  ${outputDir}/04_trim_galore_fastqc/${SampleName}_R1_trimed.fq.gz \
  ${outputDir}/04_trim_galore_fastqc/${SampleName}_R2_trimed.fq.gz
"

########################################
# Step 5: BBSplit (Bacterial read detection)
########################################
run_step "${outputDir}/05_BBSplit" "
bash /home/yz2474/yiyong_2023/soft/bbmap/bbsplit.sh \
  build=1 \
  in1=${outputDir}/03_trim_galore/${SampleName}_R1_trimed.fq.gz \
  in2=${outputDir}/03_trim_galore/${SampleName}_R2_trimed.fq.gz \
  ref=$ref_ecoli,$ref_myco \
  basename=${outputDir}/05_BBSplit/${SampleName}_%_R#.fq.gz \
  ambiguous=best
"

########################################
# Step 6: Alignment Using BWA
########################################
run_step "${outputDir}/06_BWA" "
bwa mem -t 20 -T 19 $bwa_index_ref ${outputDir}/03_trim_galore/${SampleName}_R1_trimed.fq.gz ${outputDir}/03_trim_galore/${SampleName}_R2_trimed.fq.gz > ${outputDir}/06_BWA/${SampleName}.sam
samtools sort -@ 10 -o ${outputDir}/06_BWA/${SampleName}.sorted.bam ${outputDir}/07_BWA/${SampleName}.sam
samtools index -@ 10 ${outputDir}/06_BWA/${SampleName}.sorted.bam
#rm ${outputDir}/06_BWA/${SampleName}.sam
"

########################################
# Step 7: UMI Deduplication and Sorting
########################################
run_step "${outputDir}/07_umi_dedup_bam" "
umi_tools dedup \
  -I ${outputDir}/06_BWA/${SampleName}.sorted.bam \
  -S ${outputDir}/07_umi_dedup_bam/${SampleName}_umi_dedup.sorted.bam \
  --log=${outputDir}/07_umi_dedup_bam/${SampleName}_umi_dedup.log \
  --extract-umi-method=read_id --paired
samtools index -@ 10 ${outputDir}/07_umi_dedup_bam/${SampleName}_umi.dedup.sorted.bam
samtools flagstat -@ 10 ${outputDir}/07_umi_dedup_bam/${SampleName}_umi.dedup.sorted.bam > ${outputDir}/07_umi_dedup_bam/${SampleName}_umi_dedup.flagstat
python /home/yz2474/scripts/donglab/run_infer_experiment_ZYY.py ${outputDir}/07_umi_dedup_bam/${SampleName}_umi.dedup.sorted.bam
"

########################################
# Step 8: Circular RNA Detection with CIRI2
########################################
run_step "${outputDir}/08_CIRI2" "
samtools view -@ 10 -h -o ${outputDir}/08_CIRI2/${SampleName}_umi.dedup.sorted.sam ${outputDir}/07_umi_dedup_bam/${SampleName}_umi_dedup.sorted.bam

perl /home/yz2474/yiyong_2023/soft/CIRI_v2.0.6/CIRI2.pl -T 19 -I ${outputDir}/08_CIRI2/${SampleName}_umi.dedup.sorted.sam  -O ${outputDir}/08_CIRI2/${SampleName}_CIRI2_out.txt -F $genome_fa -A $combined_gene_gtf 
rm ${outputDir}/08_CIRI2/${SampleName}_umi.dedup.sorted.sam 
"

########################################
# Step 9: Picard RNA Metrics
########################################
run_step "${outputDir}/9_picard_metrics" "
source /mnt/data/projects/donglab/yiyong_2023/soft/Miniforge_3/etc/profile.d/conda.sh
conda activate picard_env
picard -Xmx250g CollectRnaSeqMetrics \
  I=${outputDir}/07_umi_dedup_bam/${SampleName}_umi.dedup.sorted.bam \
  O=${outputDir}/09_picard_metrics/${SampleName}_umi_dedup.picard_metrics.txt \
  REF_FLAT=$combined_gene_refflat \
  STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
  RIBOSOMAL_INTERVALS=$rRNA_interval
picard CollectInsertSizeMetrics \
  I=${outputDir}/07_umi_dedup_bam/${SampleName}_umi.dedup.sorted.bam \
  O=${outputDir}/09_picard_metrics/${SampleName}_umi_dedup.insert_size_metrics.txt \
  H=${outputDir}/09_picard_metrics/${SampleName}_umi_dedup.insert_size_histogram.pdf
conda deactivate
"

########################################
# Step 10: FeatureCounts for Read counting
########################################
run_step "${outputDir}/12_featurecounts" "
featureCounts \
  -a $combined_gene_gtf \
  -o ${outputDir}/12_featurecounts/${SampleName}_featureCounts.txt \
  -p -T 20 -s 2 -g gene_id -t exon -B -C -M -O --fraction --countReadPairs \
  ${outputDir}/08_umi_dedup_bam/${SampleName}_umi.dedup.sorted.bam
"

########################################
# Step 11: gene expr TPM calculation
########################################
run_step "${outputDir}/11_gene_TPM" "
python /home/yz2474/scripts/donglab/run_gene_TPM_calculate_ZYY.py \
  -i ${outputDir}/12_featurecounts/${SampleName}_featureCounts.txt \
  -o ${outputDir}/11_gene_TPM/${SampleName}_gene_TPM.txt
"


########################################
# Step 12: Kraken for Taxonomic Classification (Downsampling)
########################################
run_step "${outputDir}/13_kraken" "
seqtk sample -s100 ${outputDir}/04_UMI_trim/${SampleName}_R1_trimed.fq.gz 100000 | gzip > ${outputDir}/13_kraken/${SampleName}_R1_downsampled_100K.fq.gz
seqtk sample -s100 ${outputDir}/04_UMI_trim/${SampleName}_R2_trimed.fq.gz 100000 | gzip > ${outputDir}/13_kraken/${SampleName}_R2_downsampled_100K.fq.gz

source /mnt/data/projects/donglab/yiyong_2023/soft/Miniforge_3/etc/profile.d/conda.sh
conda activate kraken2_env

kraken2 --db /home/yz2474/yiyong_2023/soft/kraken2/standard_krakendb --threads 20 \
  --report ${outputDir}/13_kraken/${SampleName}_report.txt --report-minimizer-data \
  --paired --gzip-compressed \
  ${outputDir}/13_kraken/${SampleName}_R1_downsampled_100K.fq.gz \
  ${outputDir}/13_kraken/${SampleName}_R2_downsampled_100K.fq.gz

python /home/yz2474/yiyong_2023/soft/KrakenTools/kreport2krona.py \
  -r ${outputDir}/13_kraken/${SampleName}_report.txt \
  -o ${outputDir}/13_kraken/${SampleName}_krona_input.txt

ktImportText ${outputDir}/13_kraken/${SampleName}_krona_input.txt \
  -o ${outputDir}/13_kraken/${SampleName}_krona.html

conda deactivate
"

########################################
# Step 13: rRNA Detection
########################################
run_step "${outputDir}/14_ribodetector" "
ribodetector_cpu -t 20 \
  -l 100 \
  -i ${outputDir}/04_UMI_trim/${SampleName}_R1_trimed.fq.gz \
     ${outputDir}/04_UMI_trim/${SampleName}_R2_trimed.fq.gz \
  -e rrna \
  --chunk_size 256 \
  -r ${outputDir}/14_ribodetector/${SampleName}_rRNA_R1_trimed.fq.gz \
     ${outputDir}/14_ribodetector/${SampleName}_rRNA_R2_trimed.fq.gz \
  -o /dev/null /dev/null
"

########################################
# Step 14: BigWig Density Generation
########################################
run_step "${outputDir}/15_bw_density" "
bamCoverage -b ${outputDir}/08_umi_dedup_bam/${SampleName}_umi.dedup.sorted.bam \
  -o ${outputDir}/15_bw_density/${SampleName}_umi_dedup_unstranded.bw \
  --normalizeUsing BPM --numberOfProcessors 20
bash ~/scripts/donglab/run_deeptools_plotprfile_bed_stackd_202412_ZYY.sh \
  ${outputDir}/15_bw_density/${SampleName}_umi_dedup_unstranded.bw
"

echo "Pipeline completed successfully for sample: ${SampleName}"
