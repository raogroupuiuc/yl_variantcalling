#!/bin/bash
#SBATCH -N 1
#SBATCH -n 3
#SBATCH --ntasks=1
#SBATCH --mem=30g
#SBATCH -p normal
#SBATCH --mail-user=deewan2@illinois.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -J YL_CNVkit
#SBATCH -D /home/a-m/deewan2/2021_yl_dnaseq_jin/out_slurm



cd ../

module purge
module load CNVkit/0.9.8-IGB-gcc-8.2.0-Python-3.7.2

echo " CNVkit - read depth  "

cd results/CNVkit/

#cnvkit.py batch ../gatk/XEV_GTATGTTC-TTCCTGTT_aligned.sorted.bam \
#	-n ../gatk/X123_CGCTATGT-GTGTCGGA_aligned.sorted.bam \
#	-m wgs -f ../../genome/PO1f_plasmid.fna \
#	--annotate ../../genome/PO1f_plasmid.gff3 \
#	--output-reference XEV_relative \
#	--scatter --diagram 

echo "1-done"


cnvkit.py batch ../gatk/GEV_ACGCACCT-CCTTCACC_aligned.sorted.bam \
       -n ../gatk/X123_CGCTATGT-GTGTCGGA_aligned.sorted.bam \
       -m wgs -f ../../genome/PO1f_plasmid.fna \
       --annotate ../../genome/PO1f_plasmid.gff3 \
       --output-reference GEV_relative \
       --scatter --diagram

echo "end"


	
