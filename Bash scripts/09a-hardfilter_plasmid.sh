#!/bin/bash
#SBATCH -N 1
#SBATCH -n 3
#SBATCH --ntasks=2
#SBATCH --mem=40g
#SBATCH -p normal
#SBATCH --mail-user=deewan2@illinois.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -J fltr_3
#SBATCH --array=1-3%2
#SBATCH -D /home/a-m/deewan2/2021_yl_dnaseq_jin/out_slurm

cd ../

line=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" raw_data/filename_plasmid.txt)

-R ../../genome/PO1f_plasmid.fna

module purge 
module load GATK/4.1.4.0-Java-1.8.0_152

echo "gatk - hard filter"
cd results/bwa_align/ 

gatk VariantFiltration \
	-V ${line}.vcf \
	-O ${line}.filtered.vcf \
	-R ../../genome/PO1f_plasmid.fna \
	--cluster-size 3 \
	--cluster-window-size 10 \
	--filter-expression "QD < 2.0" \
	--filter-name "QDFilter" \
	--filter-expression "MQ < 40.0" \
	--filter-name "MQFilter" \
	--filter-expression "FS > 60.0" \
	--filter-name "FSFilter" \
        --filter-expression "SOR > 3.0" \
        --filter-name "SORfilter" \
	--filter-expression "HaplotypeScore > 13.0" \
	--filter-name "HaplotypeScoreFilter" \
	--filter-expression "MQRankSum < -12.5" \
	--filter-name "MQRankSumFilter" \
	--filter-expression "ReadPosRankSum < -8.0" \
	--filter-name "ReadPosRankSumFilter" 

echo "end"

