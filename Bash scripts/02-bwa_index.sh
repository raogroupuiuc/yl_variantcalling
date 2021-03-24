#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --mem=10g
#SBATCH -p normal
#SBATCH --mail-user=deewan2@illinois.edu
#SBATCH --mail-type=ALL
#SBATCH -J yl.bwa.index
#SBATCH -D /home/a-m/deewan2/2021_yl_dnaseq_jin/out_slurm

cd ../genome/

module purge
module load BWA/0.7.17-IGB-gcc-4.9.4
module load GATK/4.1.4.0-Java-1.8.0_152
module load SAMtools/1.11-IGB-gcc-8.2.0

echo "start align of genome"

bwa index -p PO1f -a is GCA_009372015.1_YarliW29_genomic.fna
bwa index -p PO1f_plasmid -a is PO1f_plasmid.fna

gatk CreateSequenceDictionary -R GCA_009372015.1_YarliW29_genomic.fna
gatk CreateSequenceDictionary -R PO1f_plasmid.fna

samtools faidx PO1f_plasmid.fna
samtools faidx GCA_009372015.1_YarliW29_genomic.fna

echo "\n\n end align of genome \n\n "
