#!/bin/bash
#SBATCH -N 1
#SBATCH -n 3
#SBATCH --ntasks=1
#SBATCH --mem=10g
#SBATCH -p normal
#SBATCH --mail-user=deewan2@illinois.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -J gcvf
#SBATCH -D /home/a-m/deewan2/2021_yl_dnaseq_jin/out_slurm

cd ../

module purge 
module load GATK/4.1.4.0-Java-1.8.0_152

cd results/bwa_align/

gatk CombineGVCFs \
	-R ../../genome/PO1f_plasmid.fna \
        -V X123_CGCTATGT-GTGTCGGA.GVCF.vcf \
        -V XEV_GTATGTTC-TTCCTGTT.GVCF.vcf \
        -V GEV_ACGCACCT-CCTTCACC.GVCF.vcf \
	-O plasmid.g.vcf 

gatk CombineGVCFs \
        -R ../../genome/GCA_009372015.1_YarliW29_genomic.fna \
        -V 16EV_TCGATATC-ACTCGTGT.GVCF.vcf \
        -V 2EV_CGTCTGCG-ATTGTGAA.GVCF.vcf \
        -V 1EV_TACTCATA-GCCACAGG.GVCF.vcf \
        -V 41EV_CTAGCGCT-GTCTACAC.GVCF.vcf \
        -V PO1f_TATCGCAC-ACACTAAG.GVCF.vcf \
        -O strains.g.vcf

gatk GenotypeGVCFs \
        -R ../../genome/PO1f_plasmid.fna \
        -V plasmid.g.vcf \
	-O plasmid.vcf

gatk GenotypeGVCFs \
        -R ../../genome/GCA_009372015.1_YarliW29_genomic.fna \
        -V strains.g.vcf \
        -O strains.vcf
