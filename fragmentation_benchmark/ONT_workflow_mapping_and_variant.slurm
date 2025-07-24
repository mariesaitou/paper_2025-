#!/bin/bash
#SBATCH --job-name=ONT_mapping_variant
#SBATCH --array=1-10
#SBATCH --partition=hugemem-avx2
#SBATCH --mem=50G
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mail-type=END

### Environment
module load Miniconda3 && eval "$(conda shell.bash hook)"
conda activate /mnt/courses/BIO326/EUK/condaenvironments/condaEUK

### Path setup
BASE_DIR="/net/fs-2/scale/OrionStore/Projects/MSLab/Marie/Bio326"
FASTQ_DIR="${BASE_DIR}/Analysis/FilteredFastq"
OUT_DIR="${BASE_DIR}/Analysis/Mapping/ARS-UCD2.0"
REF="${BASE_DIR}/GCF_002263795.3_ARS-UCD2.0_genomic.mmi"  # minimap2 pre-indexed reference

FILE_LIST="${BASE_DIR}/Analysis/nanoplot_all_files.txt"
# /net/fs-2/scale/OrionStore/Projects/MSLab/Marie/Bio326/2024
# /net/fs-2/scale/OrionStore/Projects/MSLab/Marie/Bio326/2025

### Retrieve file name
FASTQ_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILE_LIST")
BASENAME=$(basename "$FASTQ_FILE" .fastq.gz)

IN_FASTQ="${FASTQ_DIR}/${FASTQ_FILE}"
OUT_SAM="${OUT_DIR}/${BASENAME}.sam"
OUT_BAM="${OUT_DIR}/${BASENAME}.bam"
OUT_VCF="${OUT_DIR}/${BASENAME}.vcf"

mkdir -p "$OUT_DIR"

### Logging
echo "===== START: ${BASENAME} ====="
echo "Input FASTQ: $IN_FASTQ"
echo "Reference: $REF"
echo "Start time: $(date)"

### Mapping
zcat "$IN_FASTQ" | minimap2 -a -t 12 -K 50m "$REF" - > "$OUT_SAM"

### BAM conversion
samtools view -@ 4 -b "$OUT_SAM" | samtools sort -@ 4 -o "$OUT_BAM"
samtools index "$OUT_BAM"

### Structural variant calling
sniffles --input "$OUT_BAM" --vcf "$OUT_VCF" --threads 12

### Cleanup
rm "$OUT_SAM"

echo "===== DONE: ${BASENAME} ====="
echo "End time: $(date)"
