#!/bin/bash
#SBATCH -N 1
#SBATCH -n 3
#SBATCH --ntasks=2
#SBATCH --mem=30g
#SBATCH -p normal
#SBATCH --mail-user=deewan2@illinois.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -J YL_SAM_coverage
#SBATCH --array=1-3%2
#SBATCH -D /home/a-m/deewan2/2021_yl_dnaseq_jin/out_slurm

cd ../

line=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" raw_data/filename_plasmid.txt)

module purge
module load SAMtools/1.11-IGB-gcc-8.2.0

echo " samtools - get coverage "

cd results/gatk/

samtools coverage -q 30 \
        --reference ../../genome/PO1f_plasmid.fna \
        -o ${line}_aligned.sorted.coverage.allchr.txt \
         ${line}_aligned.sorted.bam

echo "end"


	
