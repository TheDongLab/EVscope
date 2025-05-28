#USAGE: python run_infer_experiment_ZYY.py test.bam
import subprocess
import sys
import os

def run_infer_experiment(bam_file):
    # Check if the input BAM file exists
    if not os.path.isfile(bam_file):
        print(f"Error: The file {bam_file} does not exist.")
        return

    # Specify the path to the reference gene model BED file
    refgene_bed = "/home/yz2474/scripts/donglab/EVscope/genome_anno/UCSC_HG38/3_gencode.v45.chr_patch_hapl_scaff.annotation_UCSC.geneID.bed6"

    # Extract the base name (file name without path) from the input BAM file
    base_name = os.path.basename(bam_file).replace('.bam', '')
    
    bam_dir = os.path.dirname(bam_file)

    # Derive the output file name with the same directory as the BAM file
    output_file = os.path.join(bam_dir, f"{base_name}_bam2strandness.tsv")


    # Open the output file in write mode
    with open(output_file, "w") as file:
        # Write the header
        file.write("BAM\tdata type\t'1++,1--,2+-,2-+ (forward)'\t'1+-,1-+,2++,2-- (reverse)'\t%reads failed to determine:\n")

        try:
            # Build the command string
            command = f"infer_experiment.py -i {bam_file} -r {refgene_bed} -s 100000000"

            # Run the subprocess and capture the output
            process = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            # Parse the library type and results
            result_lines = process.stdout.strip().split('\n')
            library_type = result_lines[0].strip() if result_lines else "Library type not found"
            fraction_failed = result_lines[1].split(": ")[1] if len(result_lines) > 1 else ""
            fraction_1 = result_lines[2].split(": ")[1] if len(result_lines) > 2 else ""
            fraction_2 = result_lines[3].split(": ")[1] if len(result_lines) > 3 else ""

            # Write each file's information to the output file
            file.write(f"{bam_file}\t{library_type}\t{fraction_1}\t{fraction_2}\t{fraction_failed}\n")

        except subprocess.CalledProcessError as e:
            # Handle subprocess error
            file.write(f"{bam_file}\tError\tExit code: {e.returncode}, Error message: {e.stderr}\n")
            print(f"Error processing {bam_file}: {e.stderr}")

    print(f"Results have been written to {output_file}")

if __name__ == "__main__":
    # Check if BAM file is provided as a command-line argument
    if len(sys.argv) != 2:
        print("Usage: python 02_bam2strandness_infer_experiment.py <bam_file>")
        sys.exit(1)

    bam_file = sys.argv[1]
    run_infer_experiment(bam_file)

