#!/bin/bash
#SBATCH --job-name=kraken_bracken
#SBATCH --mem=50G  # Increase memory
#SBATCH --cpus-per-task=4    
#SBATCH --mail-user=marie.saitou@nmbu.no
#SBATCH --output=kraken_bracken.%A_%a.txt
#SBATCH --mail-type=END

# Activate Conda environment
eval "$(conda shell.bash hook)"
conda activate ~/conda_envs/metaGTools_custom

# Define directories
DATA_DIR="/net/fs-2/scale/OrionStore/Projects/MSLab/Akira"
INPUT_DIR="${DATA_DIR}/2025/unmapped/fastq_pair_output"
KRAKEN_DB="/net/fs-2/scale/OrionStore/Databases/kraken2/Standard"
KRAKEN_RESULTS="${DATA_DIR}/2025/kraken2"
BRACKEN_RESULTS="${DATA_DIR}/2025/kraken2"

# Load Kraken2 database once
KRAKEN2_COMMAND="kraken2 --db \"$KRAKEN_DB\" --threads 4 --gzip-compressed"

# Run Kraken2 and then Bracken for each pair of reads
for R1 in ${INPUT_DIR}/*_R1_unmapped.fq.paired.fq.gz; do
    R2="${R1/_R1_/_R2_}"
    
    # Get the base name
    base_name=$(basename "$R1" | sed 's/_R1_unmapped.fq.paired.fq.gz//')
    
    if [[ ! -f "$R2" ]]; then
        echo "Error: Paired file for $R1 not found. Skipping..."
        continue
    fi

    # Run Kraken2 to generate the kraken_report
    echo "Running Kraken2 on $base_name..."
    $KRAKEN2_COMMAND --fastq-input --paired "$R1" "$R2" --output "${KRAKEN_RESULTS}/kraken_report_${base_name}.txt"
    
    # Run Bracken
    echo "Running Bracken on $base_name..."
    bracken -d "$KRAKEN_DB" \
            -i "${KRAKEN_RESULTS}/kraken_report_${base_name}.txt" \
            -o "${BRACKEN_RESULTS}/bracken_${base_name}.txt"

done

echo "Processing completed."
