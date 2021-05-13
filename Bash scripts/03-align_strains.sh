#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks=2
#SBATCH --mem=10g
#SBATCH -p normal
#SBATCH --mail-user=deewan2@illinois.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -J align
#SBATCH -D /home/a-m/deewan2/2021_yl_dnaseq_jin/out_slurm
#SBATCH --array=1-5%2

cd ../ 

line=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" raw_data/filename_PO1f.txt)

module purge
module load BWA/0.7.17-IGB-gcc-4.9.4
module load SAMtools/1.7-IGB-gcc-4.9.4

echo "Alignment - start"

bwa mem -v 1 genome/PO1f \
        raw_data/${line}_L001_R1_001.fastq \
        raw_data/${line}_L001_R2_001.fastq \
         > results/gatk/${line}_aligned.sam

echo "Alignment - end"
echo " "

cd /home/a-m/deewan2/2021_yl_dnaseq_jin/results/gatk

echo "start bam"

samtools view -bS ${line}_aligned.sam > ${line}_aligned.bam
echo "bam file done"
samtools sort ${line}_aligned.bam > ${line}_aligned.sorted.bam
echo "sorting bam file done"
samtools index ${line}_aligned.sorted.bam

echo "end bam indexing done "


	
