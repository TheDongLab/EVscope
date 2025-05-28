#!/usr/bin/env python3
#usage: python calculate_ACC_motif_ratio.py A1_Yale_P2_A1_S1_L001_R2_001.fastq.gz

import sys
import gzip



def calculate_motif_ratio(fastq_path, motif="ACC", start=8, end=11):
    total_reads = 0
    motif_reads = 0

    # 判断是否为 gzipped 文件，使用不同的打开方式
    if fastq_path.endswith(".gz"):
        f = gzip.open(fastq_path, "rt")  # 文本模式读取gzip文件
    else:
        f = open(fastq_path, "r")

    with f:
        while True:
            header = f.readline()
            if not header:
                break  # EOF
            seq = f.readline().strip()
            plus = f.readline()
            qual = f.readline()

            total_reads += 1
            # 检查序列长度足够，并判断第9-11个碱基是否为 motif
            if len(seq) >= end:
                if seq[start:end] == motif:
                    motif_reads += 1

    ratio = motif_reads / total_reads if total_reads > 0 else 0
    return total_reads, motif_reads, ratio

def main():
    if len(sys.argv) != 2:
        print("Usage: python calculate_ACC_motif_ratio.py <fastq_file>")
        sys.exit(1)

    fastq_file = sys.argv[1]
    total, count, ratio = calculate_motif_ratio(fastq_file)
    print(f"Total reads: {total}")
    print(f"Reads with bases 9-11 == 'ACC': {count}")
    print(f"Ratio: {ratio:.4f} ({ratio*100:.2f}%)")

if __name__ == "__main__":
    main()

