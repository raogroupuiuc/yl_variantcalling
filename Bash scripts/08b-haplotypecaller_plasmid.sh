#!/bin/bash
#SBATCH -N 1
#SBATCH -n 3
#SBATCH --ntasks=2
#SBATCH --mem=40g
#SBATCH -p normal
#SBATCH --mail-user=deewan2@illinois.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -J vcf_3
#SBATCH --array=1-3%2
#SBATCH -D /home/a-m/deewan2/2021_yl_dnaseq_jin/out_slurm

cd ../

line=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" raw_data/filename_plasmid.txt)

module purge 
module load GATK/4.1.4.0-Java-1.8.0_152
cd results/gatk

echo "gatk - remove dups start  "

gatk BuildBamIndex -I ${line}_aligned.marked_dup.sort.RG.bam

gatk HaplotypeCaller -I ${line}_aligned.marked_dup.sort.RG.bam -O ${line}.GVCF.vcf -R ../../genome/PO1f_plasmid.fna -ERC GVCF -ploidy 1

echo "end"

