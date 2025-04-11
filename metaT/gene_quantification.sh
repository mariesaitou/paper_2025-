#!/bin/bash
#SBATCH --job-name=salmon
#SBATCH --mem=40G  # Allocate memory
#SBATCH --cpus-per-task=6
#SBATCH --partition=hugemem-avx2  
#SBATCH --output=salmon.%A_%a.txt
#SBATCH --error=salmon.%A_%a.err
#SBATCH --mail-type=END

# Uncomment if using Conda
#eval "$(conda shell.bash hook)"
#conda activate  /mnt/users/mariesai/.conda/envs/assembly_env

# Run salmon index command (not used in this script)
# salmon index -t metaT_nucleotides.fna -i metaT_salmon_gene_index

# Load the Salmon module
module load Salmon/1.1.0-gompi-2019b

# Working directories
INDIR="/net/fs-2/scale/OrionStore/Projects/MSLab/Akira/2025/unmapped/"
OUTDIR="salmon_gene_output"
INDEX_DIR="/net/fs-2/scale/OrionStore/Projects/MSLab/Akira/2025/unmapped/megahit_output_20250224_230453/metaT_salmon_gene_index"

# Automatically retrieve a list of FASTQ files (looking for *_unmapped.fq.gz)
FASTQ_FILES=($(find "$INDIR" -name "*_unmapped.fq.gz" | LC_ALL=C sort))

# Run Salmon for each pair of FASTQ files
for ((i=0; i<${#FASTQ_FILES[@]}; i+=2))
do
  SAMPLE_R1=${FASTQ_FILES[$i]}
  SAMPLE_R2=${FASTQ_FILES[$i+1]}
  
  # Get the sample name (extract plate name only)
  SAMPLE_NAME=$(basename "$SAMPLE_R1" | sed -E 's/_R[12]_unmapped.fq.gz//')

  SAMPLE_OUTDIR="$OUTDIR/$SAMPLE_NAME"

  mkdir -p "$SAMPLE_OUTDIR"
  
  echo "Processing $SAMPLE_R1 and $SAMPLE_R2 for $SAMPLE_NAME"
  
  # Run Salmon quantification
  salmon quant -i "$INDEX_DIR" \
    -l A \
    -1 "$SAMPLE_R1" \
    -2 "$SAMPLE_R2" \
    -p 8 \
    --validateMappings \
    -o "$SAMPLE_OUTDIR" || echo "Error processing $SAMPLE_NAME"
done

# Completion message
echo "Salmon quantification completed."
