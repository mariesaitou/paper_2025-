#!/bin/bash
#SBATCH --job-name=salmon
#SBATCH --mem=40G 
#SBATCH --cpus-per-task=6
#SBATCH --partition=hugemem-avx2  
#SBATCH --output=salmon.%A_%a.txt
#SBATCH --error=salmon.%A_%a.err
#SBATCH --mail-user=marie.saitou@nmbu.no
#SBATCH --mail-type=END

eval "$(conda shell.bash hook)"
conda activate  /mnt/users/mariesai/.conda/envs/assembly_env



minimap2 -d mito.mmi salmon_mt.fasta

minimap2 -ax map-ont mito.mmi final.contigs.fa | samtools view -bS - | samtools sort -o mito_mapped.bam
samtools index mito_mapped.bam


samtools view mito_mapped.bam | awk '{print $1}' | sort | uniq > mito_matches.txt

seqkit grep -v -f mito_matches.txt final.contigs.fa > filtered_contigs.fa


