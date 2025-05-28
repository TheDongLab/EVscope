#!/usr/bin/env bash
# Usage: bash run_SMARTer_pipe_v3_ZYY.sh <SampleName> <R1.fq.gz> <R2.fq.gz>
set -eo pipefail
export TMPDIR=/home/yz2474/tmp

# Input parameters
SampleName="$1"
R1_path="$2"
R2_path="$3"

# Define paths to reference files and indexes
genomeDir="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Index/Index_STAR-2.7.10b_sjdbOverhang100_Gene"
combined_gene_gtf="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/combine/11711764_HG38_gene.gtf"
combined_gene_refflat="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/combine/11711764_HG38_gene.refflat"
genecodeV45_gene_gtf="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/combine/gencode.v45.chr_patch_hapl_scaff.annotation_UCSC.gtf"
genecodeV45_gene_refflat="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/combine/gencode.v45.chr_patch_hapl_scaff.annotation_UCSC.refflat"
genome_fa="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Sequence/hg38p14/hg38.p14.whole.genome.fa"
bwa_index_ref="/home/yz2474/yiyong_2023/DB/UCSC_hg38/Index/Index_BWA_V0.7.18/hg38.p14.whole.genome.fa"
rRNA_interval="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/combine/rRNA/GeneCodeV45_combined_Nucl_Mt_rRNA_589.interval_list"
ref_myco="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/Bacterial/240_Mycoplasma_ge/Mycoplasma_240_ge.fa"
ref_ecoli="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/Bacterial/E.coli/E.coli_GCF_000005845.2_ASM584v2_ge.fa"
HG38_3UTR_bed="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_3UTR_noOverlap.bed"
HG38_5UTR_bed="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_5UTR_noOverlap.bed"
HG38_downstream_2kb_bed="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_downstream_2kb_noOverlap.bed"
HG38_exon_bed="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_exon_noOverlap.bed"
HG38_intergenic_bed="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_intergenic_noOverlap.bed"
HG38_intron_bed="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_intron_noOverlap.bed"
HG38_promoter_1500_500bp_bed="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_promoter_1500_500bp_noOverlap.bed"
HG38_blacklist="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/HG38_GeneCodeV45_meta_bed/HG38_ENCODE_blacklist_V2.bed"
GeneID_meta_table="/home/yz2474/scripts/donglab/EV-RNA-Profiler/genome_anno/UCSC_HG38/combine/11711764_HG38_geneID_Symbol_RNAtype.tsv"

# Function to run each step; only when all commands run successfully will the file step.done be created.
# If any command in the step fails, the script will exit immediately.
run_step() {
    local step_dir="$1"
    local step_cmds="$2"
    mkdir -p "$step_dir"
    if [ ! -f "$step_dir/step.done" ]; then
        echo "======> Starting step in $step_dir"
        start_time=$(date +%s)
        if eval "$step_cmds"; then
            end_time=$(date +%s)
            elapsed=$((end_time - start_time))
            touch "$step_dir/step.done"
            echo "======> Finished step: $step_dir in $elapsed seconds"
        else
            end_time=$(date +%s)
            elapsed=$((end_time - start_time))
            echo "Error: Step failed in $step_dir after $elapsed seconds. Check logs for details."
            exit 1
        fi
    else
        echo "======> Already finished step: $step_dir"
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
run_step "${outputDir}/02_UMI_motif" "
python /home/yz2474/scripts/donglab/get_SMARTer_Read2_14BP_motif/get_UMI_seq_motif_for_individual_fq.py \
  -head ${SampleName} \
  -fq $R2_path \
  -n 14 -r 100000 \
  -o ${outputDir}/02_UMI_motif \
  -p 20
"

########################################
# Step 3: UMI labeling for Read1/2 and adapter trimming
########################################
run_step "${outputDir}/03_trim_galore" "
umi_tools extract --bc-pattern='NNNNNNNN' \
  --stdin=$R2_path \
  --stdout=${outputDir}/03_trim_galore/${SampleName}_R2_UMI.fq.gz \
  --read2-in=$R1_path \
  --read2-out=${outputDir}/03_trim_galore/${SampleName}_R1_UMI.fq.gz \
  --log=${outputDir}/03_trim_galore/UMI_extract.log \
  --umi-separator='_'
cutadapt -u 14 -o ${outputDir}/03_trim_galore/${SampleName}_R2_UMI_trim.fq.gz \
  ${outputDir}/03_trim_galore/${SampleName}_R2_UMI.fq.gz
trim_galore --phred33 --illumina --gzip --stringency 3 --cores 10 --paired --length 10 \
  ${outputDir}/03_trim_galore/${SampleName}_R1_UMI.fq.gz \
  ${outputDir}/03_trim_galore/${SampleName}_R2_UMI_trim.fq.gz \
  -o ${outputDir}/03_trim_galore
mv ${outputDir}/03_trim_galore/${SampleName}_R1_UMI_val_1.fq.gz ${outputDir}/03_trim_galore/${SampleName}_R1_trimed.fq.gz
mv ${outputDir}/03_trim_galore/${SampleName}_R2_UMI_trim_val_2.fq.gz ${outputDir}/03_trim_galore/${SampleName}_R2_trimed.fq.gz
"

########################################
# Step 4: QC after trim_galore
########################################
run_step "${outputDir}/04_trim_galore_fastqc" "
fastqc -o ${outputDir}/04_trim_galore_fastqc -t 20 \
  ${outputDir}/03_trim_galore/${SampleName}_R1_trimed.fq.gz \
  ${outputDir}/03_trim_galore/${SampleName}_R2_trimed.fq.gz
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
  ambiguous=best \
  path=${outputDir}/05_BBSplit
"

########################################
# Step 6: Alignment Using BWA
########################################
run_step "${outputDir}/06_BWA" "
bwa mem -t 20 -T 19 $bwa_index_ref ${outputDir}/03_trim_galore/${SampleName}_R1_trimed.fq.gz ${outputDir}/03_trim_galore/${SampleName}_R2_trimed.fq.gz > ${outputDir}/06_BWA/${SampleName}.sam
samtools sort -@ 10 -o ${outputDir}/06_BWA/${SampleName}.sorted.bam ${outputDir}/06_BWA/${SampleName}.sam
samtools index -@ 10 ${outputDir}/06_BWA/${SampleName}.sorted.bam
"

########################################
# Step 7: UMI Deduplication and Sorting
########################################
run_step "${outputDir}/07_UMI_dedup_bam" "
umi_tools dedup \
  -I ${outputDir}/06_BWA/${SampleName}.sorted.bam \
  -S ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam \
  --log=${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.log \
  --extract-umi-method=read_id --paired
samtools index -@ 10 ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam
samtools flagstat -@ 10 ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam > ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.flagstat
python /home/yz2474/scripts/donglab/run_infer_experiment_ZYY.py ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam
"

########################################
# Step 8: Picard RNA Metrics
########################################
run_step "${outputDir}/08_picard_metrics" "
source /mnt/data/projects/donglab/yiyong_2023/soft/Miniforge_3/etc/profile.d/conda.sh
conda activate picard_env
picard -Xmx250g CollectRnaSeqMetrics \
  I=${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam \
  O=${outputDir}/08_picard_metrics/${SampleName}_UMI_dedup.picard_metrics.tsv \
  REF_FLAT=$genecodeV45_gene_refflat \
  STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
  RIBOSOMAL_INTERVALS=$rRNA_interval
picard CollectInsertSizeMetrics \
  I=${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam \
  O=${outputDir}/08_picard_metrics/${SampleName}_UMI_dedup.insert_size_metrics.tsv \
  H=${outputDir}/08_picard_metrics/${SampleName}_UMI_dedup.insert_size_histogram.pdf
conda deactivate
"

########################################
# Step 9: featureCounts for read counting
########################################
run_step "${outputDir}/09_featureCounts_TPM" "
#Fragment counts
featureCounts \
  -a $combined_gene_gtf \
  -o ${outputDir}/09_featureCounts_TPM/${SampleName}_fragment_level_featureCounts.tsv \
  -p -T 20 -s 2 -g gene_id -t exon -B -C -M -O --fraction --countReadPairs \
  ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam
#Read counts (each read counted separately)
featureCounts \
  -a $combined_gene_gtf \
  -o ${outputDir}/09_featureCounts_TPM/${SampleName}_read_level_featureCounts.tsv \
  -p -T 20 -s 2 -g gene_id -t exon -B -C -M -O --fraction \
  ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam
"

########################################
# Step 10: Circular RNA Detection with CIRI2
########################################
run_step "${outputDir}/10_CIRI2" "
curr_dir=\$(pwd)
cd ${outputDir}/10_CIRI2
perl /home/yz2474/yiyong_2023/soft/CIRI_v2.0.6/CIRI2.pl -T 19 -I ../06_BWA/${SampleName}.sam  -O ${SampleName}_CIRI2_out.tsv -F $genome_fa -A $genecodeV45_gene_gtf -G ${SampleName}_CIRI2.log
cd \$curr_dir
"

########################################
# Step 11: Expression matrix generation with read counts and normalized TPM/CPM
########################################
# Fragment-level TPM
run_step "${outputDir}/11_expression_matrix" "
# Generate linear expression file (calculate TPM from featureCounts)
awk 'BEGIN { FS="\t"; OFS="\t" }
     # Process annotation file (first file)
     FNR==NR {
         if (NR > 1) { 
             ann[$1] = $2 "\t" $3 
         }
         next
     }
     # Skip header in the second file
     FNR==1 { next }
     {
         gene = $1; 
         len = $6 + 0; 
         count = $7 + 0;
         if (len > 0) { 
             rpk = count / (len/1000); 
         } else { 
             rpk = 0 
         }
         total += rpk;
         raw[gene] = count;
         rpk_arr[gene] = rpk;
     }
     END {
         print "GeneID","GeneSymbol","GeneType","ReadCounts:gene-readcounts:circRNA:","Normalized_Expression-gene:TPM-circRNA:CPM";
         for (g in rpk_arr) {
             if (raw[g] == 0) next;  # Skip genes with zero raw count
             tpm = rpk_arr[g] / total * 1000000;
             if (g in ann)
                 print g, ann[g], raw[g], tpm;
             else
                 print g, "NA", "NA", raw[g], tpm;
         }
     }' "$GeneID_meta_table" "${outputDir}/11_expression_matrix/${SampleName}_fragment_level_featureCounts.tsv" > "${outputDir}/11_expression_matrix/${SampleName}_fragment_level_TPM_final.tsv"

# Step 3: Process circRNA file.
python /home/yz2474/scripts/donglab/EV-RNA-Profiler/filter_CIRI2_by_UMI_dedup_BAM.py  --CIRI2_result ${outputDir}/10_CIRI2/${SampleName}_CIRI2_out.tsv --input_bam ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam --output ${outputDir}/11_expression_matrix/${SampleName}_CIRI2_dedup_out.tsv
awk 'BEGIN { FS="\t"; OFS="\t" }
     NR==1 { print "GeneID","GeneSymbol","GeneType","ReadCounts","Normalized_Expression-gene:TPM-circRNA:CPM"; next }
     { print $1, $1, "circRNAs", $2, $3 }' ${outputDir}/11_expression_matrix/${SampleName}_CIRI2_dedup_out.tsv > ${outputDir}/11_expression_matrix/${SampleName}_CIRI2_final_dedup_out.tsv

# Step 4: Combine linear and circRNA files.
cat "${outputDir}/11_expression_matrix/${SampleName}_fragment_level_TPM_final.tsv" \
    <(tail -n +2 "${outputDir}/11_expression_matrix/${SampleName}_CIRI2_final_dedup_out.tsv") \
    > "${outputDir}/11_expression_matrix/${SampleName}_Gene_circRNA_readcounts_matrix.tsv"

#draw stacked bar plot for expressed RNA type composition
#python /home/yz2474/scripts/donglab/EV-RNA-Profiler/StackedRNATypeCompositionPlot.py $combined_gene_gtf ${outputDir}/11_expression_matrix/${SampleName}_fragment_level_featureCounts.tsv ${outputDir}/11_expression_matrix/${SampleName}_fragment_level_TPM.tsv ${SampleName}
"

########################################
# Step 12: bedtools for read counting in GeneCodeV45 meta-features
########################################
run_step "${outputDir}/12_bedtools_meta-features_readcounts" "
bedtools intersect -a $HG38_3UTR_bed -b ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam -c > ${outputDir}/12_bedtools_meta-features_readcounts/${SampleName}_HG38_3UTR_bedtools_readcounts.tsv
bedtools intersect -a $HG38_5UTR_bed -b ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam -c > ${outputDir}/12_bedtools_meta-features_readcounts/${SampleName}_HG38_5UTR_bedtools_readcounts.tsv
bedtools intersect -a $HG38_downstream_2kb_bed -b ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam -c > ${outputDir}/12_bedtools_meta-features_readcounts/${SampleName}_HG38_downstream_2kb_bedtools_readcounts.tsv
bedtools intersect -a $HG38_exon_bed -b ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam -c > ${outputDir}/12_bedtools_meta-features_readcounts/${SampleName}_HG38_exon_bedtools_readcounts.tsv
bedtools intersect -a $HG38_intergenic_bed -b ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam -c > ${outputDir}/12_bedtools_meta-features_readcounts/${SampleName}_HG38_intergenic_bedtools_readcounts.tsv
bedtools intersect -a $HG38_intron_bed -b ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam -c > ${outputDir}/12_bedtools_meta-features_readcounts/${SampleName}_HG38_intron_bedtools_readcounts.tsv
bedtools intersect -a $HG38_promoter_1500_500bp_bed -b ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam -c > ${outputDir}/12_bedtools_meta-features_readcounts/${SampleName}_HG38_promoter_1500_500bp_bedtools_readcounts.tsv
bedtools intersect -a $HG38_blacklist -b ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam -c > ${outputDir}/12_bedtools_meta-features_readcounts/${SampleName}_HG38_GENECODE_blacklist_bedtools_readcounts.tsv
bedtools coverage -a $HG38_3UTR_bed -b ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam -c > ${outputDir}/12_bedtools_meta-features_readcounts/${SampleName}_HG38_3UTR_bedtools_cov.tsv
bedtools coverage -a $HG38_5UTR_bed -b ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam -c > ${outputDir}/12_bedtools_meta-features_readcounts/${SampleName}_HG38_5UTR_bedtools_cov.tsv
bedtools coverage -a $HG38_downstream_2kb_bed -b ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam -c > ${outputDir}/12_bedtools_meta-features_readcounts/${SampleName}_HG38_downstream_2kb_bedtools_cov.tsv
bedtools coverage -a $HG38_exon_bed -b ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam -c > ${outputDir}/12_bedtools_meta-features_readcounts/${SampleName}_HG38_exon_bedtools_cov.tsv
bedtools coverage -a $HG38_intergenic_bed -b ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam -c > ${outputDir}/12_bedtools_meta-features_readcounts/${SampleName}_HG38_intergenic_bedtools_cov.tsv
bedtools coverage -a $HG38_intron_bed -b ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam -c > ${outputDir}/12_bedtools_meta-features_readcounts/${SampleName}_HG38_intron_bedtools_cov.tsv
bedtools coverage -a $HG38_promoter_1500_500bp_bed -b ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam -c > ${outputDir}/12_bedtools_meta-features_readcounts/${SampleName}_HG38_promoter_1500_500bp_bedtools_cov.tsv
bedtools coverage -a $HG38_blacklist -b ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam -c > ${outputDir}/12_bedtools_meta-features_readcounts/${SampleName}_HG38_GENECODE_blacklist_bedtools_cov.tsv
"

########################################
# Step 13: Kraken for Taxonomic Classification (Downsampling)
########################################
run_step "${outputDir}/12_Kraken" "
seqtk sample -s100 ${outputDir}/03_trim_galore/${SampleName}_R1_trimed.fq.gz 100000 | gzip > ${outputDir}/12_Kraken/${SampleName}_R1_downsampled_100K.fq.gz
seqtk sample -s100 ${outputDir}/03_trim_galore/${SampleName}_R2_trimed.fq.gz 100000 | gzip > ${outputDir}/12_Kraken/${SampleName}_R2_downsampled_100K.fq.gz

source /mnt/data/projects/donglab/yiyong_2023/soft/Miniforge_3/etc/profile.d/conda.sh
conda activate kraken2_env

kraken2 --db /home/yz2474/yiyong_2023/soft/kraken2/standard_krakendb --threads 20 \
  --report ${outputDir}/12_Kraken/${SampleName}_report.tsv --report-minimizer-data \
  --paired --gzip-compressed \
  ${outputDir}/12_Kraken/${SampleName}_R1_downsampled_100K.fq.gz \
  ${outputDir}/12_Kraken/${SampleName}_R2_downsampled_100K.fq.gz

python /home/yz2474/yiyong_2023/soft/KrakenTools/kreport2krona.py \
  -r ${outputDir}/12_Kraken/${SampleName}_report.tsv \
  -o ${outputDir}/12_Kraken/${SampleName}_krona_input.tsv

ktImportText ${outputDir}/12_Kraken/${SampleName}_krona_input.tsv \
  -o ${outputDir}/12_Kraken/${SampleName}_krona.html
conda deactivate
"

########################################
# Step 14: rRNA Detection
########################################
run_step "${outputDir}/13_RiboDetector" "
ribodetector_cpu -t 20 \
  -l 100 \
  -i ${outputDir}/03_trim_galore/${SampleName}_R1_trimed.fq.gz \
     ${outputDir}/03_trim_galore/${SampleName}_R2_trimed.fq.gz \
  -e rrna \
  --chunk_size 256 \
  -r ${outputDir}/13_RiboDetector/${SampleName}_rRNA_R1_trimed.fq.gz \
     ${outputDir}/13_RiboDetector/${SampleName}_rRNA_R2_trimed.fq.gz \
  -o /dev/null /dev/null
"

########################################
# Step 15: bigWig Density Generation
########################################
run_step "${outputDir}/14_bigWig" "
bamCoverage -b ${outputDir}/07_UMI_dedup_bam/${SampleName}_UMI_dedup.sorted.bam \
  -o ${outputDir}/14_bigWig/${SampleName}_UMI_dedup_unstranded.bw \
  --normalizeUsing BPM --numberOfProcessors 20
bash ~/scripts/donglab/run_deeptools_plotprfile_bed_stackd_202412_ZYY.sh \
  ${outputDir}/14_bigWig/${SampleName}_UMI_dedup_unstranded.bw
"


########################################
# Step 16: QC_matrix table generation
########################################

echo "Pipeline completed successfully for sample: ${SampleName}"

