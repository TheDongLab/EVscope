#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$repo_dir"

mkdir -p tests/smoke/work
rm -rf tests/smoke/work/emapper tests/smoke/work/smrna
chmod +x tests/smoke/mock_bin/*

tr -d '\r' < EVscope.sh | bash -n
for sh in bin/*.sh; do tr -d '\r' < "$sh" | bash -n; done
python3 -m py_compile bin/*.py
python3 tests/test_summary_report.py
help_output=$(bash EVscope.sh --help)
grep -q -- '--kit <kit>' <<< "$help_output"
grep -q 'smarter_total_rna' <<< "$help_output"
grep -q 'smarter_smrna' <<< "$help_output"
bash EVscope.sh --version >/dev/null

if bash EVscope.sh --kit unsupported_kit >/dev/null 2>tests/smoke/work/invalid_kit.stderr; then
  echo "unsupported --kit value was unexpectedly accepted" >&2
  exit 1
fi
grep -q -- "--kit must be 'smarter_total_rna' or 'smarter_smrna'" tests/smoke/work/invalid_kit.stderr

python3 bin/Step_02_calculate_ACC_motif_fraction.py \
  --input_fastq tests/smoke/data/toy_R2.fastq.gz \
  --positions 1-3 \
  --output_tsv tests/smoke/work/acc.tsv

test -s tests/smoke/work/acc.tsv

smrna_work="$repo_dir/tests/smoke/work/smrna"
mkdir -p "$smrna_work"
(
  cd "$smrna_work"
  PATH="$repo_dir/tests/smoke/mock_bin:$PATH" \
  MOCK_CUTADAPT_LOG="$smrna_work/cutadapt.calls" \
  MOCK_STAR_LOG="$smrna_work/star.calls" \
  MOCK_SAMTOOLS_LOG="$smrna_work/samtools.calls" \
    bash "$repo_dir/EVscope.sh" \
      --sample_name smarter_smoke \
      --kit smarter_smrna \
      --input_fastqs "$repo_dir/tests/smoke/data/smarter_smrna_R1.fastq" "$repo_dir/tests/smoke/data/smarter_smrna_R2.fastq" \
      --run_steps 2-3,6 \
      --config "$repo_dir/tests/smoke/minimal.conf"
)

smrna_output="$smrna_work/smarter_smoke_EVscope_output"
clean_fastq="$smrna_output/Step_03_UMI_Adaptor_Trim/smarter_smoke_R1_clean.fq.gz"
clean_fastq2="$smrna_output/Step_03_UMI_Adaptor_Trim/smarter_smoke_R2_clean.fq.gz"
profile_tsv="$smrna_output/Step_02_UMI_Analysis/smarter_smoke_library_profile.tsv"
step_meta="$smrna_output/Step_03_UMI_Adaptor_Trim/step.meta.json"
alignment_meta="$smrna_output/Step_06_Alignment_Refined/step.meta.json"
analysis_bam="$smrna_output/Step_06_Alignment_Refined/smarter_smoke_STAR_no_dedup_Aligned.sortedByCoord.out.bam"
gzip -dc "$clean_fastq" | diff -u tests/smoke/expected/smarter_smrna_clean.fastq -
gzip -dc "$clean_fastq2" | diff -u tests/smoke/expected/smarter_smrna_R2_clean.fastq -
awk -F '\t' '$1 == "UMI" && $2 == "No" {found=1} END {exit !found}' "$profile_tsv"
awk -F '\t' '$1 == "Input layout" && $2 == "Paired-end Read 1 / Read 2" {found=1} END {exit !found}' "$profile_tsv"
grep -q '"kit": "smarter_smrna"' "$step_meta"
grep -q '"strand": "forward"' "$step_meta"
grep -Fq -- '-u 3 -U 14 -a AAAAAAAAAA;min_overlap=3' "$smrna_work/cutadapt.calls"
grep -Fq -- '-A NNNAGATCGGAAGAGC;min_overlap=10' "$smrna_work/cutadapt.calls"
grep -Fq -- '-q 20 -Q 20 --minimum-length 15' "$smrna_work/cutadapt.calls"
test ! -e "$smrna_output/Step_03_UMI_Adaptor_Trim/smarter_smoke_R1_umi_tools.fq.gz"
test -s "$analysis_bam"
grep -q '"kit": "smarter_smrna"' "$alignment_meta"
grep -q '"strand": "forward"' "$alignment_meta"
test -s "$smrna_output/Step_06_Alignment_Initial/smarter_smoke_R1_no_dedup.clean.fq.gz"
test -s "$smrna_output/Step_06_Alignment_Initial/smarter_smoke_R2_no_dedup.clean.fq.gz"
test ! -e "$smrna_output/Step_06_Alignment_Initial/smarter_smoke_R1_umi_dedup.clean.fq.gz"
test "$(wc -l < "$smrna_work/star.calls")" -eq 2
test "$(grep -Fc -- '--twopassMode None' "$smrna_work/star.calls")" -eq 2
test "$(grep -Fc -- '--alignEndsType EndToEnd' "$smrna_work/star.calls")" -eq 2
test "$(grep -Fc -- '--alignIntronMax 1' "$smrna_work/star.calls")" -eq 2
test "$(grep -Fc -- '--outFilterMismatchNmax 1' "$smrna_work/star.calls")" -eq 2
test "$(grep -Fc -- '--outFilterScoreMinOverLread 0' "$smrna_work/star.calls")" -eq 2
test "$(grep -Fc -- '--outFilterMatchNminOverLread 0' "$smrna_work/star.calls")" -eq 2
test "$(grep -Fc -- '--outFilterScoreMin 24' "$smrna_work/star.calls")" -eq 2
test "$(grep -Fc -- '--outFilterMatchNmin 28' "$smrna_work/star.calls")" -eq 2
test "$(grep -Fc -- '--mockGlobalInjected yes' "$smrna_work/star.calls")" -eq 2
test "$(grep -Fc -- '--mockSmrnaInjected yes' "$smrna_work/star.calls")" -eq 2
test "$(grep -Fc -- '--chimSegmentMin 10' "$smrna_work/star.calls")" -eq 2
test "$(grep -Fc -- '--chimJunctionOverhangMin 10' "$smrna_work/star.calls")" -eq 2
test "$(grep -Fc -- '--chimOutType Junctions WithinBAM' "$smrna_work/star.calls")" -eq 2
test -s "$smrna_output/Step_06_Alignment_Refined/smarter_smoke_STAR_no_dedup_Chimeric.out.junction"

if (
  cd "$smrna_work"
  PATH="$repo_dir/tests/smoke/mock_bin:$PATH" \
    bash "$repo_dir/EVscope.sh" \
      --sample_name smarter_smoke \
      --kit smarter_total_rna \
      --input_fastqs "$repo_dir/tests/smoke/data/smarter_smrna_R1.fastq" "$repo_dir/tests/smoke/data/smarter_smrna_R2.fastq" \
      --run_steps 2 \
      --resume \
      --config "$repo_dir/tests/smoke/minimal.conf"
) >/dev/null 2>"$smrna_work/cross_kit_resume.stderr"
then
  echo "cross-kit --resume unexpectedly reused SMARTer smRNA results" >&2
  exit 1
fi
grep -q 'metadata does not match current inputs/config/parameters' "$smrna_work/cross_kit_resume.stderr"

if python3 - <<'PY'
import importlib.util, sys
sys.exit(0 if importlib.util.find_spec('matplotlib') else 1)
PY
then
  python3 bin/Step_03_plot_fastq_read_length_dist.py \
    --input_fastqs tests/smoke/data/toy_R2.fastq.gz tests/smoke/data/toy_R2.fastq.gz tests/smoke/data/toy_R2.fastq.gz \
    --titles toy_R2 toy_R2 toy_R2 \
    --output_pdf tests/smoke/work/read_length.pdf \
    --output_png tests/smoke/work/read_length.png
  test -s tests/smoke/work/read_length.pdf
  test -s tests/smoke/work/read_length.png
else
  echo "matplotlib not installed; skipping read-length plot smoke test"
fi

if python3 - <<'PY'
import importlib.util, sys
missing = [m for m in ('pysam', 'pyBigWig', 'numba', 'numpy') if importlib.util.find_spec(m) is None]
if missing:
    print('missing dependencies:', ','.join(missing))
    sys.exit(1)
PY
then
  mkdir -p tests/smoke/work/emapper
  python3 - <<'PY'
from pathlib import Path
import pysam
work = Path('tests/smoke/work/emapper')
bam_path = work / 'tiny.bam'
header = {'HD': {'VN': '1.6', 'SO': 'queryname'}, 'SQ': [{'SN': 'chrToy', 'LN': 1000}]}
with pysam.AlignmentFile(str(bam_path), 'wb', header=header) as out:
    for i, start in enumerate([100, 120, 200], 1):
        read = pysam.AlignedSegment()
        read.query_name = f'read{i}'
        read.query_sequence = 'A' * 50
        read.flag = 0
        read.reference_id = 0
        read.reference_start = start
        read.mapping_quality = 255
        read.cigar = [(0, 50)]
        read.query_qualities = pysam.qualitystring_to_array('I' * 50)
        read.set_tag('NH', 1)
        out.write(read)
with open(work / 'tiny.gtf', 'w') as fh:
    fh.write('chrToy\tsmoke\texon\t90\t260\t.\t+\t.\tgene_id "GENE1"; gene_name "GENE1"; transcript_id "TX1"; gene_type "protein_coding";\n')
PY
  python3 bin/Step_25_EMapper.py \
    --input_bam tests/smoke/work/emapper/tiny.bam \
    --gtf tests/smoke/work/emapper/tiny.gtf \
    --strandness unstrand \
    --output_dir tests/smoke/work/emapper/out \
    --prefix tiny \
    --num_threads 1 \
    --no-cleanup \
    --max_iter 5 \
    --tol 0.01
  test -s tests/smoke/work/emapper/out/tiny_readcounts.txt
  test -s tests/smoke/work/emapper/out/tiny_unstranded.bw
  test -s tests/smoke/work/emapper/out/tiny_unstranded_uniq_only.bw
  python3 bin/Step_25_bigWig2Expression.py \
    --input_combined_bw tests/smoke/work/emapper/out/tiny_unstranded.bw \
    --gtf tests/smoke/work/emapper/tiny.gtf \
    --output tests/smoke/work/emapper/tiny_bigwig_expression.tsv
  test -s tests/smoke/work/emapper/tiny_bigwig_expression.tsv
  grep -q '^GENE1' tests/smoke/work/emapper/out/tiny_readcounts.txt
  grep -q '^GENE1' tests/smoke/work/emapper/tiny_bigwig_expression.tsv
else
  echo "EMapper dependencies not installed; skipping EMapper smoke test"
fi

sha256sum tests/smoke/work/acc.tsv > tests/smoke/work/sha256sums.txt

echo "EVscope smoke tests passed"
