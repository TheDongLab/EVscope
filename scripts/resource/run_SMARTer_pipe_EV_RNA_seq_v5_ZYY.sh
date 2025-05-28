#!/usr/bin/env bash
# Usage: bash run_SMARTer_pipe_v3_ZYY.sh <SampleName> <R1.fq.gz> <R2.fq.gz>
# Example: bash run_SMARTer_pipe_v3_ZYY.sh B4_TDP43_KD_caRNA_rep3 /path/to/R1.fq.gz /path/to/R2.fq.gz
set -eo pipefail
Header=$1
R1_path=$2
R2_path=$3
genomeDir="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Index/Index_STAR-2.7.10b_sjdbOverhang100_Gene"
Gene_GTF="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Annotation/refs_EV-RNA-Profiler/UCSC_HG38/combine/15907154_HG38_gene.gtf"
genome_fa="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Sequence/hg38p14/hg38.p14.whole.genome.fa"
bwa_index_ref="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Index/Index_BWA_V0.7.18/hg38.p14.whole.genome.fa"
rRNA_bed="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Annotation/refs_EV-RNA-Profiler/UCSC_HG38/combine/GeneCodeV45_combined_Nucl_Mt_rRNA_589.bed"
ref_myco="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Annotation/refs_EV-RNA-Profiler/Bacterial/240_Mycoplasma_ge/Mycoplasma_240_ge.fa"
ref_ecoli="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Annotation/refs_EV-RNA-Profiler/Bacterial/E.coli/E.coli_GCF_000005845.2_ASM584v2_ge.fa"

run_step() {
    local step_dir="$1"
    local step_cmds="$2"
    echo "Running in: $step_dir"
    mkdir -p "$step_dir"
    if [ ! -f "$step_dir/step.done" ]; then
        eval "$step_cmds"
        # Only create step.done if all commands run successfully
        if [ $? -eq 0 ]; then
            touch "$step_dir/step.done"
        else
            echo "Error: Step failed in $step_dir. Check logs for details."
            exit 1
        fi
    fi
}

outputDir="${Header}"
mkdir -p "${outputDir}"
########################################
# Step 1: Raw QC
########################################
run_step "${outputDir}/01_raw_fastqc" "
fastqc -o ${outputDir}/01_raw_fastqc -t 20 $R1_path $R2_path
"

########################################
# Step 2: 14-bp motif for R2
########################################
run_step "${outputDir}/02_umi_motif" "
python /home/yz2474/scripts/donglab/get_SMARTer_Read2_14BP_motif/get_umi_seq_motif_for_individual_fq.py \
  -head ${Header} \
  -fq $R2_path \
  -n 14 -r 100000 \
  -o ${outputDir}/02_umi_motif \
  -p 20
"

########################################
# Step 3: UMI extraction
########################################
run_step "${outputDir}/03_umi_marker_for_fq" "
umi_tools extract --bc-pattern='NNNNNNNN' \
  --stdin=$R2_path \
  --stdout=${outputDir}/03_umi_marker_for_fq/${Header}_R2_umi.fq.gz \
  --read2-in=$R1_path \
  --read2-out=${outputDir}/03_umi_marker_for_fq/${Header}_R1_umi.fq.gz \
  --log=${outputDir}/03_umi_marker_for_fq/umi_extract.log \
  --umi-separator='_'
"

########################################
# Step 4: Cut UMI + adapter trimming
########################################
run_step "${outputDir}/04_UMI_trim" "
cutadapt -10 -u 14 -m 10 -o ${outputDir}/04_UMI_trim/${Header}_R2_umi_trim.fq.gz \
  ${outputDir}/03_umi_marker_for_fq/${Header}_R2_umi.fq.gz

trim_galore --phred33 --illumina --gzip --stringency 3 --cores 10 --paired --length 10 \
  ${outputDir}/03_umi_marker_for_fq/${Header}_R1_umi.fq.gz \
  ${outputDir}/04_UMI_trim/${Header}_R2_umi_trim.fq.gz \
  -o ${outputDir}/04_UMI_trim
"

########################################
# Step 5: QC after trimming
########################################
run_step "${outputDir}/05_UMI_trim_fastqc" "
fastqc -o ${outputDir}/05_UMI_trim_fastqc -t 20 \
  ${outputDir}/04_UMI_trim/${Header}_R1_umi_val_1.fq.gz \
  ${outputDir}/04_UMI_trim/${Header}_R2_umi_trim_val_2.fq.gz
"

########################################
# Step 6: BBSplit
########################################
run_step "${outputDir}/06_BBSplit" "
bash /home/yz2474/yiyong_2023/soft/bbmap/bbsplit.sh \
  build=1 \
  in1=$R1_path \
  in2=$R2_path \
  ref_ecoli= $ref_ecoli \
  ref_myco=$ref_ecoli \
  outx1=${outputDir}/06_BBSplit/${Header}_ecoli_R1.fq.gz \
  outx2=${outputDir}/06_BBSplit/${Header}_ecoli_R2.fq.gz \
  outy1=${outputDir}/06_BBSplit/${Header}_myco_R1.fq.gz \
  outy2=${outputDir}/06_BBSplit/${Header}_myco_R2.fq.gz \
  ambiguous=best

########################################
# Step 7: BWA
########################################
run_step "${outputDir}/07_BWA" "
bwa mem -t 20 -T 19 $bwa_index_ref ${outputDir}/04_UMI_trim/${Header}_R1_umi_val_1.fq.gz ${outputDir}/04_UMI_trim/${Header}_R2_umi_trim_val_2.fq.gz > ${outputDir}/07_BWA/${Header}.sam
samtools view -@ 10 -bS ${outputDir}/07_BWA/${Header}.sam > ${outputDir}/07_BWA/${Header}.bam
rm ${outputDir}/07_BWA/${Header}.sam
"

########################################
# Step 8: UMI deduplication
########################################
run_step "${outputDir}/08_umi_dedup_bam" "
umi_tools dedup \
  -I ${outputDir}/07_BWA/${Header}.bam \
  -S ${outputDir}/08_umi_dedup_bam/${Header}_umi_dedup.bam \
  --log=${outputDir}/08_umi_dedup_bam/${Header}_umi_dedup.log \
  --extract-umi-method=read_id --paired
samtools index -@ 10 ${outputDir}/08_umi_dedup_bam/${Header}_umi_dedup.bam
samtools sort -@ 10 -n ${outputDir}/08_umi_dedup_bam/${Header}_umi_dedup.bam > ${outputDir}/08_umi_dedup_bam/${Header}_Aligned.sortedByCoord.umi.dedup.bam
samtools flagstat -@ 10 ${outputDir}/08_umi_dedup_bam/${Header}_Aligned.sortedByCoord.umi.dedup.bam > ${outputDir}/08_umi_dedup_bam/${Header}_umi_dedup.flagstat
python /home/yz2474/scripts/donglab/run_infer_experiment_ZYY.py \
  ${outputDir}/08_umi_dedup_bam/${Header}_Aligned.sortedByCoord.umi.dedup.bam
"

########################################
# Step 9: CIRI2
########################################
run_step "${outputDir}/09_CIRI2" "
samtools view -@ 10 -h -o ${outputDir}/9_CIRI2/${Header}_umi_dedup.sam ${outputDir}/08_CIRI2/${Header}_umi_dedup.bam
perl /home/yz2474/yiyong_2023/soft/CIRI_v2.0.6 -T 19 -I ${outputDir}/09_CIRI2/${Header}_umi_dedup.sam  -O ${outputDir}/09_CIRI2/${Header}_CIRI2_out.txt -F $genome_fa -A $Gene_GTF 
rm ${outputDir}/08_umi_dedup_bam/${Header}_umi_dedup.bam ${outputDir}/09_CIRI2/${Header}_umi_dedup.sam
"

########################################
# Step 10: Picard_metrics
########################################
run_step "${outputDir}/9_picard_metrics" "
source /mnt/data/projects/donglab/yiyong_2023/soft/Miniforge_3/etc/profile.d/conda.sh
conda activate picard_env
picard CollectRnaSeqMetrics \
  I=${outputDir}/8_umi_dedup_bam/${Header}_Aligned.sortedByCoord.umi.dedup.bam \
  O=${outputDir}/9_picard_metrics/${Header}_umi_dedup.picard_metrics.txt \
  REF_FLAT=$refflat \
  STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
  RIBOSOMAL_INTERVALS=$rRNA_bed
picard CollectInsertSizeMetrix \
  I=${outputDir}/8_umi_dedup_bam/${Header}_Aligned.sortedByCoord.umi.dedup.bam \
  O=${outputDir}/9_picard_metrics/${Header}_umi_dedup.insert_size_metrics.txt \
  H=${outputDir}/9_picard_metrics/${Header}_umi_dedup.insert_size_histogram.pdf
conda deactivate
"


########################################
# Step 11: salmon for TPM
########################################
run_step "${outputDir}/11_salmon" "
conda activate salmon_env
salmon quant --no-version-check --index salmonIndex \
--libType A --output sample_output --threads 20 \
--numBootstraps 100 --validateMappings --seqBias \
--gcBias --dumpEq \
--alignments ${outputDir}/8_umi_dedup_bam/${Header}_Aligned.sortedByCoord.umi.dedup.bam  \
--geneMap $Gene_GTF
conda deactivate
"

########################################
# Step 12: featureCounts
########################################
run_step "${outputDir}/12_featurecounts" "
featureCounts \
  -a $Gene_GTF \
  -o ${outputDir}/12_featurecounts/${Header}_featureCounts.txt \
  -p -T 20 -s 2 -g gene_id -t exon -B -C -M --fraction --countReadPairs \
  ${outputDir}/08_umi_dedup_bam/${Header}_Aligned.sortedByCoord.umi.dedup.bam
"

########################################
# Step 13: Kraken (downsampling)
########################################
run_step "${outputDir}/13_kraken" "
seqtk sample -s100 ${outputDir}/4_UMI_trim/${Header}_R1_umi_val_1.fq.gz 100000 | gzip > ${outputDir}/13_kraken/${Header}_R1_downsampled.fq.gz
seqtk sample -s100 ${outputDir}/4_UMI_trim/${Header}_R2_umi_trim_val_2.fq.gz 100000 | gzip > ${outputDir}/13_kraken/${Header}_R2_downsampled.fq.gz

source /mnt/data/projects/donglab/yiyong_2023/soft/Miniforge_3/etc/profile.d/conda.sh
conda activate kraken2_env

kraken2 --db /home/yz2474/yiyong_2023/soft/kraken2/standard_krakendb --threads 20 \
  --report ${outputDir}/13_kraken/${Header}_report.txt --report-minimizer-data \
  --paired --gzip-compressed \
  ${outputDir}/13_kraken/${Header}_R1_downsampled.fq.gz \
  ${outputDir}/13_kraken/${Header}_R2_downsampled.fq.gz

python /home/yz2474/yiyong_2023/soft/KrakenTools/kreport2krona.py \
  -r ${outputDir}/13_kraken/${Header}_report.txt \
  -o ${outputDir}/13_kraken/${Header}_krona_input.txt

ktImportText ${outputDir}/13_kraken/${Header}_krona_input.txt \
  -o ${outputDir}/13_kraken/${Header}_krona.html

conda deactivate
"

########################################
# Step 14: rRNA detection
########################################
run_step "${outputDir}/14_ribodetector" "
ribodetector_cpu -t 20 \
  -l 100 \
  -i ${outputDir}/8_bam2fastq/${Header}_umi_deduped_R1.fq.gz \
     ${outputDir}/8_bam2fastq/${Header}_umi_deduped_R2.fq.gz \
  -e rrna \
  --chunk_size 256 \
  -r ${outputDir}/14_ribodetector/${Header}_umi_deduped_rRNA_R1.fq.gz \
     ${outputDir}/14_ribodetector/${Header}_umi_deduped_rRNA_R2.fq.gz \
  -o /dev/null /dev/null
"

########################################
# Step 15: bigWig density
########################################
run_step "${outputDir}/15_bw_density" "
bamCoverage -b ${outputDir}/6_STAR/${Header}_Aligned.sortedByCoord.out.bam \
  -o ${outputDir}/15_bw_density/${Header}_umi_dedup_unstrand.bw \
  --normalizeUsing BPM --numberOfProcessors 20

bash ~/scripts/donglab/run_deeptools_plotprfile_bed_stackd_202412_ZYY.sh \
  ${outputDir}/15_bw_density/${Header}_umi_dedup_unstrand.bw
"

echo "Pipeline completed successfully for sample: $Header"