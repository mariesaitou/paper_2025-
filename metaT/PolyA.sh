#!/bin/bash
#SBATCH --job-name=metaT_processing
#SBATCH --mem=20G
#SBATCH --partition=hugemem-avx2  
#SBATCH --output=metaT_processing.%A_%a.txt
#SBATCH --mail-type=END

# Directory where the FASTQ files are located
FASTQ_DIR="/net/fs-2/scale/OrionStore/Projects/MSLab/Akira/2025/mt_unmapped"

# Output CSV file
OUTPUT_CSV="polyA_tail_24.csv"

# Threshold for Poly(A) tail and Poly(T) tail (number of consecutive A or T)
THRESHOLD=24

# Write header to the CSV
echo "File,Total_Reads,PolyA_Tail_Reads,PolyT_Tail_Reads,Percentage_PolyA,Percentage_PolyT" > "$OUTPUT_CSV"

# Process all FASTQ files in the directory
for fq in "$FASTQ_DIR"/*_nonmt.fq.gz; do
    # FASTQ file name
    fq_name=$(basename "$fq")

    # Total number of reads (since FASTQ has 4 lines per read, divide the line count by 4)
    total_reads=$(zcat "$fq" | wc -l)
    total_reads=$((total_reads / 4))

    # If it's an R1 file, detect Poly(A) tail; if it's R2, detect Poly(T) tail
    if [[ "$fq_name" == *"_R1_"* ]]; then
        # Count reads with Poly(A) tail
        polyA_reads=$(zcat "$fq" | grep -E "A{$THRESHOLD,}" | wc -l)
        polyT_reads=0
    elif [[ "$fq_name" == *"_R2_"* ]]; then
        # Count reads with Poly(T) tail
        polyT_reads=$(zcat "$fq" | grep -E "T{$THRESHOLD,}" | wc -l)
        polyA_reads=0
    else
        # Ignore files that are neither R1 nor R2 (or handle as necessary)
        continue
    fi

    # Calculate the percentage of Poly(A) and Poly(T) tails
    if [[ $total_reads -gt 0 ]]; then
        percentage_polyA=$(awk -v a="$polyA_reads" -v b="$total_reads" 'BEGIN { printf "%.2f", (a/b)*100 }')
        percentage_polyT=$(awk -v a="$polyT_reads" -v b="$total_reads" 'BEGIN { printf "%.2f", (a/b)*100 }')
    else
        percentage_polyA=0
        percentage_polyT=0
    fi

    # Append results to the CSV
    echo "$fq_name,$total_reads,$polyA_reads,$polyT_reads,$percentage_polyA,$percentage_polyT" >> "$OUTPUT_CSV"
done
