#!/bin/bash
#SBATCH -N 1
#SBATCH -n 3
#SBATCH --ntasks=2
#SBATCH --mem=30g
#SBATCH -p normal
#SBATCH --mail-user=deewan2@illinois.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -J YL_SAM_rem_dup
#SBATCH --array=1-8%2
#SBATCH -D /home/a-m/deewan2/2021_yl_dnaseq_jin/out_slurm

cd ../

line=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" raw_data/shortnames.txt)

module purge 
module load GATK/4.1.4.0-Java-1.8.0_152

cd results/bwa_align/ 

echo "gatk - Sort"


gatk AddOrReplaceReadGroups -I ${line}_aligned.marked_dup.sort.bam -O ${line}_aligned.marked_dup.sort.RG.bam -ID 1 -LB lib1 -PL ILLUMINA -PU A00835:204:HWJVHDRXX:1 -SM ${line}

echo "end"


	
