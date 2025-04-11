#!/bin/bash
#SBATCH --job-name=matam
#SBATCH --mem=40G  # Increase memory
#SBATCH --cpus-per-task=16  # Number of CPU cores per job
#SBATCH --partition=hugemem-avx2  
#SBATCH --output=matam_output.%A_%a.txt
#SBATCH --mail-type=END
#SBATCH --array=1-56%10  # Change based on the number of samples (10 samples in this case)

# Set up Conda environment
eval "$(conda shell.bash hook)"
conda activate /mnt/users/mariesai/.conda/envs/matam_env3

# Input data directory and output directory
READS_DIR="/net/fs-2/scale/OrionStore/Projects/MSLab/Akira/2025/mt_unmapped"
OUTPUT_DIR="/net/fs-2/scale/OrionStore/Projects/MSLab/Akira/2025/unmapped/matam"
DB_DIR="/net/fs-2/scale/OrionStore/Projects/MSLab/Akira/2025/unmapped/matam_db"

# Default reference database
MATAM_DB_PREFIX="SILVA_138.1_SSURef_NR99_tax_silva_NR95"

# Create output and database directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$DB_DIR"

# Get all FASTQ files and store in an array
FASTQ_FILES=("$READS_DIR"/*.fq.gz)
NUM_FILES=${#FASTQ_FILES[@]}

# Check if the array job index is within the range of available files
if [[ $SLURM_ARRAY_TASK_ID -ge $NUM_FILES ]]; then
    echo "Job ID $SLURM_ARRAY_TASK_ID is greater than the available number of files ($NUM_FILES), skipping..."
    exit 1
fi

# Get the corresponding FASTQ file
READS_GZ="${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}"
SAMPLE_NAME=$(basename "$READS_GZ" .fq.gz)
SAMPLE_OUT_DIR="$OUTPUT_DIR/$SAMPLE_NAME"
TEMP_FASTQ="$SAMPLE_OUT_DIR/$SAMPLE_NAME.fastq"

mkdir -p "$SAMPLE_OUT_DIR"

echo "Starting processing of sample $SAMPLE_NAME (Job ID: $SLURM_ARRAY_TASK_ID)..."

# Decompress the `.gz` file and convert to `.fastq`
gunzip -c "$READS_GZ" > "$TEMP_FASTQ"

# Use `trap` to delete the temporary file when the script exits
trap "rm -f $TEMP_FASTQ" EXIT

# Perform 16S/23S rRNA assembly
matam_assembly.py -d "$DB_DIR/$MATAM_DB_PREFIX" -i "$TEMP_FASTQ" --cpu 16 --max_memory 40000 -v --perform_taxonomic_assignment -o "$SAMPLE_OUT_DIR"

echo "MATAM analysis for sample $SAMPLE_NAME is complete. Output is saved in $SAMPLE_OUT_DIR."

# Completion message
echo "Processing for Job ID $SLURM_ARRAY_TASK_ID is complete."
