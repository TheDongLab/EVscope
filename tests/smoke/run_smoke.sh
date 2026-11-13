#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$repo_dir"

mkdir -p tests/smoke/work
rm -rf tests/smoke/work/emapper

bash -n EVscope.sh
for sh in bin/*.sh; do bash -n "$sh"; done
python3 -m py_compile bin/*.py
bash EVscope.sh --help >/dev/null
bash EVscope.sh --version >/dev/null

if python3 - <<'PY'
import importlib.util, sys
sys.exit(0 if all(importlib.util.find_spec(name) for name in ('seaborn', 'ARIC')) else 1)
PY
then
  python3 bin/Step_22_run_RNA_deconvolution_ARIC.py --help >/dev/null
else
  echo "seaborn or ARIC not installed; skipping deconvolution import smoke test"
fi

python3 bin/Step_02_calculate_ACC_motif_fraction.py \
  --input_fastq tests/smoke/data/toy_R2.fastq.gz \
  --positions 1-3 \
  --output_tsv tests/smoke/work/acc.tsv

test -s tests/smoke/work/acc.tsv

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
