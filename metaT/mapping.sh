#!/bin/bash
#SBATCH --job-name=unmap
#SBATCH --output=unmap.%A_%a.txt
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --partition=smallmem  # ← Specify the appropriate available partition

# Directory paths
ALL_READS_DIR="/net/fs-2/scale/OrionStore/Projects/MSLab/Akira/samples_june05"
UNMAPPED_READS_DIR="/net/fs-2/scale/OrionStore/Projects/MSLab/Akira/2025/unmapped"

# Output CSV file
OUTPUT_CSV="read_counts.csv"

# Write the header to the CSV
echo "Sample,All_Reads,Unmapped_Reads" > $OUTPUT_CSV

# Count reads for each sample
for file in ${ALL_READS_DIR}/*_1.fq.gz; do
    # Get the sample name
    basename=$(basename $file)
    sample_name=$(echo $basename | sed 's/_1.fq.gz//')

    # Get the read count for All Reads (R1 + R2)
    all_r1="${ALL_READS_DIR}/${sample_name}_1.fq.gz"
    all_r2="${ALL_READS_DIR}/${sample_name}_2.fq.gz"
    all_count=0
    if [[ -f $all_r1 ]]; then
        all_count=$((all_count + $(zgrep -c "^@" $all_r1)))
    fi
    if [[ -f $all_r2 ]]; then
        all_count=$((all_count + $(zgrep -c "^@" $all_r2)))
    fi

    # Get the read count for Unmapped Reads (R1 + R2)
    unmapped_r1="${UNMAPPED_READS_DIR}/${sample_name}_R1_unmapped.fq.gz"
    unmapped_r2="${UNMAPPED_READS_DIR}/${sample_name}_R2_unmapped.fq.gz"
    unmapped_count=0
    if [[ -f $unmapped_r1 ]]; then
        unmapped_count=$((unmapped_count + $(zgrep -c "^@" $unmapped_r1)))
    fi
    if [[ -f $unmapped_r2 ]]; then
        unmapped_count=$((unmapped_count + $(zgrep -c "^@" $unmapped_r2)))
    fi

    # Write the results to the CSV
    echo "$sample_name,$all_count,$unmapped_count" >> $OUTPUT_CSV
done

echo "Read counting completed. Results saved in $OUTPUT_CSV."
