#!/bin/bash

# This shell script is a simple pipline thats aligns pair end reads with the reference genome provided by tiny-test-data online.

cd ./tiny-test-data/wgs ||

#Runs bwa mem to align reads.
bwa mem ../genomes/Hsapiens/hg19/bwa/hg19.fa mt_1.fq.gz mt_2.fq.gz > mt.sam

#Sorts the SAM file to BAM file.
samtools sort mt.sam > sorted_mt.bam


#Generate Pile Ups and Variant Calling
bcftools mpileup -Ou -f ../genomes/Hsapiens/hg19/seq/hg19.fa sorted_mt.bam | bcftools call -vmO v -o mt.vcf

# Compress and index VCF file
bgzip mt.vcf
bcftools index -f mt.vcf.gz

# Filter variants with quality > 20
bcftools filter -i'QUAL>20' mt.vcf.gz | bcftools stats > variant_stats.txt



