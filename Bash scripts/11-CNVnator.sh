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
module load CNVnator/0.3.3-IGB-gcc-4.9.4

echo " CNVnator - read depth  "

cd results/CNVnator/

echo "1-done"

cnvnator -root ${line}.root -chrom pINT03-XYL123.gb -tree ../gatk/${line}_aligned.sorted.bam

echo "2-done"

cnvnator -root ${line}.root -chrom pINT03-XYL123.gb -his 1000

echo "3-done"

cnvnator -root ${line}.root -stat 1000 

echo "4-done"

cnvnator -root ${line}.root -partition 1000

echo "5-done"

cnvnator -root ${line}.root -call 1000

echo "6-done"

# plotcircular.py ${line}.root

echo "end"


	
