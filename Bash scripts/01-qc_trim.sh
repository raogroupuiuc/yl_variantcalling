#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks=4
#SBATCH --mem=4g
#SBATCH -p normal
#SBATCH --mail-user=deewan2@illinois.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -J YL_jin_QC
#SBATCH --array=1-16%4
#SBATCH -D /home/a-m/deewan2/2021_yl_dnaseq_jin/out_slurm


module load FastQC/0.11.8-Java-1.8.0_152

line=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" ../raw_data/filenames.txt)

echo " "
echo "Start FastQC Job"
fastqc -o ../results/fastqc_raw ../raw_data/${line}*.fastq
echo "Finish FastQC Job"
echo " "

module load Trimmomatic/0.38-Java-1.8.0_152

echo "start trim"

trimmomatic SE -phred33 -threads 1 ../raw_data/${line}*.fastq ../results/trim/${line}_trimmed.fastq \
	ILLUMINACLIP:../raw_data/TruSeq3-SE.fa:2:30:10 LEADING:28 TRAILING:28 MINLEN:30 

echo "end trim"
echo "Start fastqc"
fastqc -o ../results/fastqc_trim ../results/trim/${line}_trimmed.fastq
echo "Finish FastQC"
echo " "
