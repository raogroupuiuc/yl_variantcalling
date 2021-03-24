#!/bin/bash
#SBATCH -N 1
#SBATCH -n 3
#SBATCH --ntasks=2
#SBATCH --mem=40g
#SBATCH -p normal
#SBATCH --mail-user=deewan2@illinois.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -J vcf_5
#SBATCH --array=1-5%2
#SBATCH -D /home/a-m/deewan2/2021_yl_dnaseq_jin/out_slurm

cd ../

line=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" raw_data/filename_PO1f.txt)

module purge 
module load GATK/4.1.4.0-Java-1.8.0_152

echo "gatk - remove dups start  "
cd results/bwa_align/ 

# gatk BuildBamIndex -I ${line}_aligned.marked_dup.sort.RG.bam #
gatk HaplotypeCaller -I ${line}_aligned.marked_dup.sort.RG.bam -O ${line}.GVCF.vcf -R ../../genome/GCA_009372015.1_YarliW29_genomic.fna -ERC GVCF -ploidy 1

echo "end"

