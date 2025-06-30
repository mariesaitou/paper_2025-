#!/bin/bash
#SBATCH --job-name=metaT_processing
#SBATCH --mem=60G
#SBATCH --cpus-per-task=6
#SBATCH --partition=hugemem-avx2
#SBATCH --output=metaT_processing.%A_%a.txt
#SBATCH --mail-type=END

# Activate the appropriate conda environment
eval "$(conda shell.bash hook)"
conda activate /mnt/users/mariesai/.conda/envs/kraken2bracken

set -e  # Stop script on any error

# Define input/output directories and index
MERGED_DIR="/net/fs-2/scale/OrionStore/Projects/MSLab/Akira/2025/salmon_down/merged_fastq"
DOWNSAMPLED_DIR="downsampled"
SALMON_INDEX="/net/fs-2/scale/OrionStore/Projects/MSLab/Akira/2025/unmapped/salmon_index"
SALMON_OUT="salmon_output"

# Create output directories if not already present
mkdir -p "$DOWNSAMPLED_DIR" "$SALMON_OUT"

# Define fixed number of reads for downsampling (based on the smallest sample)
min_reads=368940
head_lines=$(( min_reads * 4 ))  # FASTQ format: 4 lines per read
echo "Minimum read count per sample: $min_reads"
echo "Extracting first $head_lines lines from each file..."

# Process each interleaved FASTQ file
for merged_file in ${MERGED_DIR}/*_merged.fq.gz; do
    output_prefix="${DOWNSAMPLED_DIR}/$(basename "$merged_file" _merged.fq.gz)"
    output_r1="${output_prefix}_downsampled_R1.fq.gz"
    output_r2="${output_prefix}_downsampled_R2.fq.gz"

    echo "Splitting $merged_file into R1 and R2..."

    # Step 1: Split interleaved FASTQ into R1 and R2 using modulo logic
    zcat "$merged_file" | awk 'NR%8<4' | pigz -p 8 > "${output_prefix}_R1.fq.gz"
    zcat "$merged_file" | awk 'NR%8>=4' | pigz -p 8 > "${output_prefix}_R2.fq.gz"

    echo "Subsampling first $head_lines lines for R1 and R2..."

    # Step 2: Extract the first N reads (deterministic downsampling)
    zcat "${output_prefix}_R1.fq.gz" | head -n "$head_lines" | pigz -p 8 > "$output_r1"
    zcat "${output_prefix}_R2.fq.gz" | head -n "$head_lines" | pigz -p 8 > "$output_r2"

    # Step 3: Remove intermediate files to save space
    rm -f "${output_prefix}_R1.fq.gz" "${output_prefix}_R2.fq.gz"

    # Step 4: Run Salmon quantification on the downsampled data
    echo "Running Salmon quantification on $output_r1 and $output_r2..."
    salmon quant -i "$SALMON_INDEX" \
        --libType A \
        --validateMappings \
        -1 "$output_r1" \
        -2 "$output_r2" \
        -o "${SALMON_OUT}/$(basename "$output_prefix")"
done

echo "All samples processed successfully."
