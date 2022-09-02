#!/bin/bash

samples_folder="raw_fastq"
output_folder="outputs"

#------------------------------------------------
# Extract read group information from each folder
#------------------------------------------------
# usage: extract_group_info.sh arg1 [arg2]
# arg1: name of the folder that contains all samples
# arg2: folders that have specified prefix are processed in this batch

../codes/extract_group_info.sh $samples_folder  > outputs/read_group_info.csv
#../codes/extract_group_info.sh $samples_folder CDCF > outputs/read_group_info_CDCF.csv
#../codes/extract_group_info.sh $samples_folder Sample > outputs/read_group_info_Sample.csv

#------------------------------------------------------------------------------
# Elign reads to reference genome CDC317, sort sam files, and markup duplicates
#------------------------------------------------------------------------------
# usage: run_bwa.sh arg1 arg2 arg3 [arg4]
# arg1: name of the folder that contains all samples
# arg2: output folder name
# arg3: number of cpus
# arg4: folders that have specified prefix are processed in this batch

../codes/run_bwa.sh $samples_folder $output_folder 36
#../codes/run_bwa.sh $samples_folder $output_folder 36 CDCF
#../codes/run_bwa.sh $samples_folder $output_folder 36 Sample

#--------------------------------------------------------------------------------
# Calling SNPs and indels simultaneously via local de-novo assembly of haplotypes
#--------------------------------------------------------------------------------
# usage: run_gatk_haplotypecaller.sh arg1 arg2 [arg3]
# arg1: name of the folder that contains all samples
# arg2: output folder name
# arg3: folders that have specified prefix are processed in this batch

#../codes/run_gatk_haplotypecaller.sh $samples_folder $output_folder
../codes/run_gatk_haplotypecaller.sh $samples_folder $output_folder CDCF
../codes/run_gatk_haplotypecaller.sh $samples_folder $output_folder Sample

#--------------------------------
# Combine GVCF and call genotypes
#--------------------------------
# usage: run_gatk_genotype.sh arg1 arg2 [arg3]
# arg1: name of the folder that contains all samples
# arg2: output folder name
# arg3: folders that have specified prefix are processed in this batch

../codes/run_gatk_genotype.sh $samples_folder $output_folder
#../codes/run_gatk_genotype.sh $samples_folder $output_folder CDCF
#../codes/run_gatk_genotype.sh $samples_folder $output_folder Sample

#-----------------------------------
# Filter variants and run annotation
#-----------------------------------
# usage: run_gatk_variant_filter.sh arg1
# usage: run_variant_annotation.sh arg1
# arg1: output folder name

../codes/run_gatk_variant_filter.sh $output_folder
../codes/run_variant_annotation.sh $output_folder
mv $PWD/snpEff_genes.txt $output_folder
mv $PWD/snpEff_summary.html $output_folder

#-----------------------------------
# Compute genome coverage statistics
#-----------------------------------
# usage: run_bedtools_coverage.sh arg1 arg2 [arg3]
# arg1: name of the folder that contains all samples
# arg2: output folder name
# arg3: folders that have specified prefix are processed in this batch

../codes/run_bedtools_coverage.sh $samples_folder $output_folder
#../codes/run_bedtools_coverage.sh $samples_folder $output_folder CDCF
#../codes/run_bedtools_coverage.sh $samples_folder $output_folder Sample

#-------------------------------
# Compute copy number variations
#-------------------------------
# usage: python compute_cnv.py arg1 arg2 [arg3]
# arg1: name of the folder that contains all samples
# arg2: output folder name
# arg3: folders that have specified prefix are processed in this batch

python ../codes/compute_cnv.py $samples_folder $output_folder
#python ../codes/compute_cnv.py $samples_folder $output_folder CDCF
#python ../codes/compute_cnv.py $samples_folder $output_folder Sample
