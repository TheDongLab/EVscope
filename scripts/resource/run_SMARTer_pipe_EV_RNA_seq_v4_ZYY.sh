#!/usr/bin/env bash
# Usage: bash run_SMARTer_pipe_v3_ZYY.sh <SampleName> <R1.fq.gz> <R2.fq.gz>
# Example: bash run_SMARTer_pipe_v3_ZYY.sh B4_TDP43_KD_caRNA_rep3 /path/to/R1.fq.gz /path/to/R2.fq.gz

set -eo pipefail

Header=$1
R1_path=$2
R2_path=$3
outputDir="${Header}"

genomeDir="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Index/Index_STAR-2.7.10b_sjdbOverhang100_Gene"
Gene_GTF="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Annotation/1_Build_custom_TE_GTF/3_gencode.v45.chr_patch_hapl_scaff.annotation_UCSC.gtf"
refflat="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Annotation/1_Build_custom_TE_GTF/4_gencode.v45.chr_patch_hapl_scaff.annotation_UCSC_RepeatMasker_v4.0.7_Dfam2.0_rmskJoinedCurrent.bed12_converted.LINE.refflat"
genome_fa="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Sequence/hg38p14/hg38.p14.whole.genome.fa"

run_step() {
    local step_dir="$1"
    local step_cmds="$2"
    echo "Running in: $step_dir"
    mkdir -p "$step_dir"
    if [ ! -f "$step_dir/step.done" ]; then
        eval "$step_cmds"
        touch "$step_dir/step.done"
    fi
}

mkdir -p "${outputDir}"

########################################
# Step 1: Raw QC
########################################
run_step "${outputDir}/1_raw_fastqc" "
fastqc -o ${outputDir}/1_raw_fastqc -t 20 $R1_path $R2_path
"

########################################
# Step 2: 14-bp motif for R2
########################################
run_step "${outputDir}/2_umi_motif" "
python /home/yz2474/scripts/donglab/get_SMARTer_Read2_14BP_motif/get_umi_seq_motif_for_individual_fq.py \
  -head ${Header} \
  -fq $R2_path \
  -n 14 -r 100000 \
  -o ${outputDir}/2_umi_motif \
  -p 20
"

########################################
# Step 3: UMI extraction
########################################
run_step "${outputDir}/3_umi_marker_for_fq" "
umi_tools extract --bc-pattern='NNNNNNNN' \
  --stdin=$R2_path \
  --stdout=${outputDir}/3_umi_marker_for_fq/${Header}_R2_umi.fq.gz \
  --read2-in=$R1_path \
  --read2-out=${outputDir}/3_umi_marker_for_fq/${Header}_R1_umi.fq.gz \
  --log=${outputDir}/3_umi_marker_for_fq/umi_extract.log \
  --umi-separator='_'
"

########################################
# Step 4: Cut UMI + adapter trimming
########################################
run_step "${outputDir}/4_UMI_trim" "
cutadapt -u 14 -o ${outputDir}/4_UMI_trim/${Header}_R2_umi_trim.fq.gz \
  ${outputDir}/3_umi_marker_for_fq/${Header}_R2_umi.fq.gz

trim_galore --phred33 --illumina --gzip --stringency 3 --cores 20 --paired --length 10 \
  ${outputDir}/3_umi_marker_for_fq/${Header}_R1_umi.fq.gz \
  ${outputDir}/4_UMI_trim/${Header}_R2_umi_trim.fq.gz \
  -o ${outputDir}/4_UMI_trim
"

########################################
# Step 5: QC after trimming
########################################
run_step "${outputDir}/5_UMI_trim_fastqc" "
fastqc -o ${outputDir}/5_UMI_trim_fastqc -t 20 \
  ${outputDir}/4_UMI_trim/${Header}_R1_umi_val_1.fq.gz \
  ${outputDir}/4_UMI_trim/${Header}_R2_umi_trim_val_2.fq.gz
"

########################################
# Step 6: STAR alignment
########################################
run_step "${outputDir}/6_STAR" "
ulimit -n 65535
STAR --genomeDir $genomeDir \
  --readFilesIn ${outputDir}/4_UMI_trim/${Header}_R1_umi_val_1.fq.gz \
                ${outputDir}/4_UMI_trim/${Header}_R2_umi_trim_val_2.fq.gz \
  --outFileNamePrefix ${outputDir}/6_STAR/${Header}_ \
  --runThreadN 20 --twopassMode Basic --runMode alignReads --quantMode GeneCounts \
  --readFilesCommand zcat --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 \
  --outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimJunctionOverhangMin 10 \
  --chimScoreMin 1 --chimOutType Junctions WithinBAM --outBAMsortingThreadN 20

samtools index -@ 20 ${outputDir}/6_STAR/${Header}_Aligned.sortedByCoord.out.bam
"

########################################
# Step 7: UMI deduplication
########################################
run_step "${outputDir}/7_umi_dedup_bam" "
umi_tools dedup \
  -I ${outputDir}/6_STAR/${Header}_Aligned.sortedByCoord.out.bam \
  -S ${outputDir}/7_umi_dedup_bam/${Header}_Aligned.sortedByCoord.umi.dedup.bam \
  --log=${outputDir}/7_umi_dedup_bam/${Header}_umi_dedup.log \
  --extract-umi-method=read_id --paired

samtools index -@ 20 ${outputDir}/7_umi_dedup_bam/${Header}_Aligned.sortedByCoord.umi.dedup.bam
python /home/yz2474/scripts/donglab/run_infer_experiment_ZYY.py \
  ${outputDir}/7_umi_dedup_bam/${Header}_Aligned.sortedByCoord.umi.dedup.bam
"

########################################
# Step 8: Convert deduplicated BAM to FASTQ
########################################
run_step "${outputDir}/8_bam2fastq" "
samtools fastq -@ 20 \
  -1 ${outputDir}/8_bam2fastq/${Header}_R1.fq.gz \
  -2 ${outputDir}/8_bam2fastq/${Header}_R2.fq.gz \
  -0 ${outputDir}/8_bam2fastq/${Header}_unpaired_R1.fq.gz \
  -s ${outputDir}/8_bam2fastq/${Header}_unpaired_R2.fq.gz \
  -n \
  ${outputDir}/7_umi_dedup_bam/${Header}_Aligned.sortedByCoord.umi.dedup.bam

cat ${outputDir}/8_bam2fastq/${Header}_R1.fq.gz ${outputDir}/8_bam2fastq/${Header}_unpaired_R1.fq.gz > ${outputDir}/8_bam2fastq/${Header}_umi_deduped_R1.fq.gz
cat ${outputDir}/8_bam2fastq/${Header}_R2.fq.gz ${outputDir}/8_bam2fastq/${Header}_unpaired_R2.fq.gz > ${outputDir}/8_bam2fastq/${Header}_umi_deduped_R2.fq.gz
rm ${outputDir}/8_bam2fastq/${Header}_R1.fq.gz \
   ${outputDir}/8_bam2fastq/${Header}_R2.fq.gz \
   ${outputDir}/8_bam2fastq/${Header}_unpaired_R1.fq.gz \
   ${outputDir}/8_bam2fastq/${Header}_unpaired_R2.fq.gz
"

########################################
# Step 9: Re-run STAR with deduped FASTQ
########################################
run_step "${outputDir}/9_rerun_STAR" "
ulimit -n 65535
STAR --genomeDir $genomeDir \
  --readFilesIn ${outputDir}/8_bam2fastq/${Header}_umi_deduped_R1.fq.gz \
                ${outputDir}/8_bam2fastq/${Header}_umi_deduped_R2.fq.gz \
  --outFileNamePrefix ${outputDir}/9_rerun_STAR/${Header}_ \
  --runThreadN 20 --twopassMode Basic --runMode alignReads --quantMode GeneCounts \
  --readFilesCommand zcat --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 \
  --outSAMtype BAM SortedByCoordinate --chimSegmentMin 10 --chimJunctionOverhangMin 10 \
  --chimScoreMin 1 --chimOutType Junctions WithinBAM --outBAMsortingThreadN 20

samtools index -@ 20 ${outputDir}/9_rerun_STAR/${Header}_Aligned.sortedByCoord.out.bam
"

########################################
# Step 10: featureCounts
########################################
run_step "${outputDir}/10_featurecounts" "
featureCounts \
  -a $Gene_GTF \
  -o ${outputDir}/10_featurecounts/${Header}_featureCounts.txt \
  -p -T 20 -s 2 -g gene_id -t exon -B -C -M --fraction --countReadPairs \
  ${outputDir}/9_rerun_STAR/${Header}_Aligned.sortedByCoord.out.bam
"

########################################
# Step 11: CIRCexplorer2
########################################
run_step "${outputDir}/11_CIRCexplorer2" "
CIRCexplorer2 parse -t STAR \
  -b ${outputDir}/11_CIRCexplorer2/${Header}_back_spliced_junction.bed \
  ${outputDir}/9_rerun_STAR/${Header}_Chimeric.out.junction

CIRCexplorer2 annotate \
  -r $refflat \
  -g $genome_fa \
  -b ${outputDir}/11_CIRCexplorer2/${Header}_back_spliced_junction.bed \
  -o ${outputDir}/11_CIRCexplorer2/${Header}_circularRNA_known.txt
"

########################################
# Step 12: Kraken (downsampling)
########################################
run_step "${outputDir}/12_kraken" "
seqtk sample -s100 ${outputDir}/4_UMI_trim/${Header}_R1_umi_val_1.fq.gz 100000 | gzip > ${outputDir}/12_kraken/${Header}_R1_downsampled.fq.gz
seqtk sample -s100 ${outputDir}/4_UMI_trim/${Header}_R2_umi_trim_val_2.fq.gz 100000 | gzip > ${outputDir}/12_kraken/${Header}_R2_downsampled.fq.gz

source /mnt/data/projects/donglab/yiyong_2023/soft/Miniforge_3/etc/profile.d/conda.sh
conda activate kraken2_env

kraken2 --db /home/yz2474/yiyong_2023/soft/kraken2/standard_krakendb --threads 20 \
  --report ${outputDir}/12_kraken/${Header}_report.txt --report-minimizer-data \
  --paired --gzip-compressed \
  ${outputDir}/12_kraken/${Header}_R1_downsampled.fq.gz \
  ${outputDir}/12_kraken/${Header}_R2_downsampled.fq.gz

python /home/yz2474/yiyong_2023/soft/KrakenTools/kreport2krona.py \
  -r ${outputDir}/12_kraken/${Header}_report.txt \
  -o ${outputDir}/12_kraken/${Header}_krona_input.txt

ktImportText ${outputDir}/12_kraken/${Header}_krona_input.txt \
  -o ${outputDir}/12_kraken/${Header}_krona.html

conda deactivate
"

########################################
# Step 13: Count reads in bed regions
########################################
run_step "${outputDir}/13_readcounts_for_bed" "
BED_FILES=(
  \"/home/yz2474/yiyong_2023/DB/UCSC_hg38/Annotation/1_Build_custom_TE_GTF/11_get_meta_bed/inner_dup_merge_inter_non_overlap/HG38_exon_noOverlap.bed\"
  \"/home/yz2474/yiyong_2023/DB/UCSC_hg38/Annotation/1_Build_custom_TE_GTF/11_get_meta_bed/inner_dup_merge_inter_non_overlap/HG38_intron_noOverlap.bed\"
  \"/home/yz2474/yiyong_2023/DB/UCSC_hg38/Annotation/1_Build_custom_TE_GTF/11_get_meta_bed/inner_dup_merge_inter_non_overlap/HG38_intergenic_noOverlap.bed\"
  \"/home/yz2474/yiyong_2023/DB/UCSC_hg38/Annotation/1_Build_custom_TE_GTF/11_get_meta_bed/inner_dup_merge_inter_non_overlap/hg38-blacklist.v2.bed\"
)
counts=()
bamFile=\"${outputDir}/9_rerun_STAR/${Header}_Aligned.sortedByCoord.out.bam\"

for bed_file in \"\${BED_FILES[@]}\"; do
  counts+=(\$(bedtools coverage -a \"\$bed_file\" -b \"\$bamFile\" -counts | awk '{sum+=\$NF} END{print sum}'))
done

outSummary=\"${outputDir}/13_readcounts_for_bed/${Header}_readcounts_summary_for_bed.txt\"
if [ ! -f \"\$outSummary\" ]; then
  echo -e \"Sample\\tExon\\tIntron\\tIntergenic\\tBlacklist\" > \"\$outSummary\"
fi
echo -e \"${Header}\\t\${counts[0]}\\t\${counts[1]}\\t\${counts[2]}\\t\${counts[3]}\" >> \"\$outSummary\"
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
