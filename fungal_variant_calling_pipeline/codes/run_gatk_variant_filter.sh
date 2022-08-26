#---------------------------------------------------
# Exclude variants marked by flags other than "PASS"
#---------------------------------------------------
$PWD/../resources/GATK/gatk-4.1.9.0/gatk SelectVariants \
    -V $PWD/$1/cohort.unfiltered.gt.vcf.gz \
    --exclude-filtered true\
    -O $PWD/$1/cohort.filtered.raw.gt.vcf.gz
num_variants_step0=$(bcftools view -H $PWD/$1/cohort.unfiltered.gt.vcf.gz | wc -l)
num_variants_step1=$(bcftools view -H $PWD/$1/cohort.filtered.raw.gt.vcf.gz | wc -l)

#------------------------------------------
# select variants: SNPs, INDELs, Mixed-type
#------------------------------------------
$PWD/../resources/GATK/gatk-4.1.9.0/gatk SelectVariants \
    -V $PWD/$1/cohort.filtered.raw.gt.vcf.gz \
    -select-type SNP \
    --exclude-filtered true\
    -O $PWD/$1/cohort.filtered.raw.gt.snp.vcf.gz
num_snps_step1=$(bcftools view -H $PWD/$1/cohort.filtered.raw.gt.snp.vcf.gz | wc -l)

$PWD/../resources/GATK/gatk-4.1.9.0/gatk SelectVariants \
    -V $PWD/$1/cohort.filtered.raw.gt.vcf.gz \
    -select-type INDEL \
    --exclude-filtered true\
    -O $PWD/$1/cohort.filtered.raw.gt.indel.vcf.gz
num_indels_step1=$(bcftools view -H $PWD/$1/cohort.filtered.raw.gt.indel.vcf.gz | wc -l)

$PWD/../resources/GATK/gatk-4.1.9.0/gatk SelectVariants \
    -V $PWD/$1/cohort.filtered.raw.gt.vcf.gz \
    -select-type MIXED \
    --exclude-filtered true\
    -O $PWD/$1/cohort.filtered.raw.gt.mixed.vcf.gz
num_mixed_step1=$(bcftools view -H $PWD/$1/cohort.filtered.raw.gt.mixed.vcf.gz | wc -l)

#------------------------------------------
# Hard-filtering using GATK recommendations
# For SNP, label SNP clusters and filter
#------------------------------------------
$PWD/../resources/GATK/gatk-4.1.9.0/gatk VariantFiltration \
  -V $PWD/$1/cohort.filtered.raw.gt.snp.vcf.gz \
  -filter "QD < 2.00" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "SOR > 3.00" --filter-name "SOR3" \
  -filter "FS > 60.00" --filter-name "FS60" \
  -filter "MQ < 40.0" --filter-name "MQ40" \
  -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  -O $PWD/$1/cohort.filtered.hf.gt.snp.vcf.gz

$PWD/../resources/GATK/gatk-4.1.9.0/gatk VariantFiltration \
  -V $PWD/$1/cohort.filtered.raw.gt.indel.vcf.gz \
  -filter "QD < 2.00" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "FS > 200.0" --filter-name "FS200" \
  -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
  -O $PWD/$1/cohort.filtered.hf.gt.indel.vcf.gz

# mixed variants are evaluated with the indel model
$PWD/../resources/GATK/gatk-4.1.9.0/gatk VariantFiltration \
  -V $PWD/$1/cohort.filtered.raw.gt.mixed.vcf.gz \
  -filter "QD < 2.00" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "FS > 200.0" --filter-name "FS200" \
  -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
  -O $PWD/$1/cohort.filtered.hf.gt.mixed.vcf.gz

#-------------------------
# Merge hard filtered VCFs
#-------------------------
java -jar $PWD/../resources/GATK/picard.jar MergeVcfs \
  I=$PWD/$1/cohort.filtered.hf.gt.snp.vcf.gz \
  I=$PWD/$1/cohort.filtered.hf.gt.indel.vcf.gz \
  I=$PWD/$1/cohort.filtered.hf.gt.mixed.vcf.gz \
  O=$PWD/$1/cohort.filtered.hf.tmp.gt.vcf.gz
$PWD/../resources/GATK/gatk-4.1.9.0/gatk SelectVariants \
     -V $PWD/$1/cohort.filtered.hf.tmp.gt.vcf.gz \
     --exclude-filtered true\
     -O $PWD/$1/cohort.filtered.hf.gt.vcf.gz
num_variants_step2=$(bcftools view -H $PWD/$1/cohort.filtered.hf.gt.vcf.gz | wc -l)

#-------------------------------------------------------------------------------------------------------------------------------------------
# Further filtering 1: minimum read depth (DP) and genotype quality (GQ), and ref-to-alt AD ratio
# Note that the -S option changes genotype of samples that do not fulfill the requirement to ./. but does not remove the entire variant site
# & and | apply multiple filters to the same sample simultaneously, while && and || apply to different samples independently
#-------------------------------------------------------------------------------------------------------------------------------------------
bcftools filter -S . -e 'FMT/DP<10 | FMT/GQ<20' -O z -o $PWD/$1/cohort.filtered.hf.DP10.GQ20.gt.vcf.gz $PWD/$1/cohort.filtered.hf.gt.vcf.gz
num_variants_step3=$(bcftools view -H $PWD/$1/cohort.filtered.hf.DP10.GQ20.gt.vcf.gz | wc -l)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------
# Further filtering 2: remove monomorphic SNPs/INDELs, multiallelic SNPs and indels, SNPs in the close proximity of INDELS, and clusters of INDELs with a window
# & and | apply multiple filters to the same sample simultaneously, while && and || apply to different samples independently
#---------------------------------------------------------------------------------------------------------------------------------------------------------------
# AC==0: no variants (all the same as the reference)
# AC==AN: only alternative alleles are called. Here we keep a site as a variant if only alternative alleles are called
#bcftools filter -e 'AC==0 || AC==AN' --IndelGap 5 --SnpGap 10 $PWD/$1/cohort.filtered.hf.DP10.GQ20.gt.vcf.gz | bcftools view -m2 -M2 -O z -o $PWD/$1/cohort.filtered.hf.DP10.GQ20.allele.gt.vcf.gz
bcftools filter -e 'AC==0' --IndelGap 5 --SnpGap 10 $PWD/$1/cohort.filtered.hf.DP10.GQ20.gt.vcf.gz | bcftools view -m2 -M2 -O z -o $PWD/$1/cohort.filtered.hf.DP10.GQ20.allele.gt.vcf.gz
num_variants_step4=$(bcftools view -H $PWD/$1/cohort.filtered.hf.DP10.GQ20.allele.gt.vcf.gz | wc -l)

#--------------------------------------------------------------
# Further filtering 3: remove variants in repetitive regions
# Note that tabix would fail if gzip, instead of bgzip, is used
#--------------------------------------------------------------
bedtools intersect -v -a $PWD/$1/cohort.filtered.hf.DP10.GQ20.allele.gt.vcf.gz \
    -b $PWD/../resources/TRF/cpar_v3_inc_mito.fasta.2.5.7.80.10.50.2000.bed -wa -header \
    > $PWD/$1/cohort.filtered.hf.DP10.GQ20.allele.trf.gt.vcf
bgzip $PWD/$1/cohort.filtered.hf.DP10.GQ20.allele.trf.gt.vcf
tabix -fp vcf $PWD/$1/cohort.filtered.hf.DP10.GQ20.allele.trf.gt.vcf.gz
num_variants_step5=$(bcftools view -H $PWD/$1/cohort.filtered.hf.DP10.GQ20.allele.trf.gt.vcf.gz | wc -l)

#----------------------------------------
# Further filtering 4: remove SNPclusters
#----------------------------------------
$PWD/../resources/GATK/gatk-4.1.9.0/gatk VariantFiltration \
    -V $PWD/$1/cohort.filtered.hf.DP10.GQ20.allele.trf.gt.vcf.gz\
    -cluster 3 -window 10\
    -O $PWD/$1/cohort.filtered.hf.DP10.GQ20.allele.trf.sc.tmp.gt.vcf.gz
$PWD/../resources/GATK/gatk-4.1.9.0/gatk SelectVariants \
    -V $PWD/$1/cohort.filtered.hf.DP10.GQ20.allele.trf.sc.tmp.gt.vcf.gz \
    --exclude-filtered true\
    -select "FILTER == SnpCluster" --invertSelect \
    -O $PWD/$1/cohort.filtered.hf.DP10.GQ20.allele.trf.sc.gt.vcf.gz
num_variants_step6=$(bcftools view -H $PWD/$1/cohort.filtered.hf.DP10.GQ20.allele.trf.sc.gt.vcf.gz | wc -l)

#-----------------------------------------------------------
# Create vcf that includes only variants with filter PASS
#-----------------------------------------------------------
bcftools view -f PASS $PWD/$1/cohort.filtered.hf.DP10.GQ20.allele.trf.sc.gt.vcf.gz -O z -o $PWD/$1/cohort.filtered.gt.vcf.gz
tabix -fp vcf $PWD/$1/cohort.filtered.gt.vcf.gz
num_variants_step7=$(bcftools view -H $PWD/$1/cohort.filtered.gt.vcf.gz | wc -l)
num_snps_step7=$(bcftools view -H --types snps $PWD/$1/cohort.filtered.gt.vcf.gz | wc -l)
num_indels_step7=$(bcftools view -H --types indels $PWD/$1/cohort.filtered.gt.vcf.gz | wc -l)
num_mixed_step7=$(bcftools view -H --types other $PWD/$1/cohort.filtered.gt.vcf.gz | wc -l)

# print number of variants at each filtering step
echo "number of variants [step0: raw variants] = $num_variants_step0"
echo "number of variants [step1: exclude filtered variants] = $num_variants_step1($num_snps_step1,$num_indels_step1,$num_mixed_step1)"
echo "number of variants [step2: hard filtering] = $num_variants_step2"
echo "number of variants [step3: filter by DP/GQ] = $num_variants_step3"
echo "number of variants [step4: filter by alleles 1] = $num_variants_step4"
echo "number of variants [step5: filter by repetitive regions] = $num_variants_step5"
echo "number of variants [step6: remove sap cluster] = $num_variants_step6"
echo "number of variants [step7: exclude variants failed to pass our filters] = $num_variants_step7($num_snps_step7,$num_indels_step7,$num_mixed_step7)"

# convert to table
$PWD/../resources/GATK/gatk-4.1.9.0/gatk VariantsToTable \
     -V $PWD/$1/cohort.filtered.gt.vcf.gz \
     -F CHROM -F POS -F ID -F TYPE -F REF -F ALT -F HET -F HOM-REF -F HOM-VAR -F NO-CALL -F VAR -F NSAMPLES -F NCALLED -F MULTI-ALLELIC -GF AD -GF GT \
     --error-if-missing-data \
     -O $PWD/$1/cohort.filtered.gt.txt

# remove temporary files
rm $PWD/$1/cohort.filtered.raw.gt.vcf.gz
rm $PWD/$1/cohort.filtered.raw.gt.snp.vcf.gz
rm $PWD/$1/cohort.filtered.raw.gt.indel.vcf.gz
rm $PWD/$1/cohort.filtered.raw.gt.mixed.vcf.gz
rm $PWD/$1/cohort.filtered.hf.gt.snp.vcf.gz
rm $PWD/$1/cohort.filtered.hf.gt.indel.vcf.gz
rm $PWD/$1/cohort.filtered.hf.gt.mixed.vcf.gz
rm $PWD/$1/cohort.filtered.hf.tmp.gt.vcf.gz
rm $PWD/$1/cohort.filtered.hf.gt.vcf.gz
rm $PWD/$1/cohort.filtered.hf.DP10.GQ20.gt.vcf.gz
rm $PWD/$1/cohort.filtered.hf.DP10.GQ20.allele.gt.vcf.gz
rm $PWD/$1/cohort.filtered.hf.DP10.GQ20.allele.trf.gt.vcf.gz
rm $PWD/$1/cohort.filtered.hf.DP10.GQ20.allele.trf.sc.tmp.gt.vcf.gz
rm $PWD/$1/cohort.filtered.hf.DP10.GQ20.allele.trf.sc.gt.vcf.gz
rm $PWD/$1/cohort.filtered.raw.gt.vcf.gz.tbi
rm $PWD/$1/cohort.filtered.raw.gt.snp.vcf.gz.tbi
rm $PWD/$1/cohort.filtered.raw.gt.indel.vcf.gz.tbi
rm $PWD/$1/cohort.filtered.raw.gt.mixed.vcf.gz.tbi
rm $PWD/$1/cohort.filtered.hf.gt.snp.vcf.gz.tbi
rm $PWD/$1/cohort.filtered.hf.gt.indel.vcf.gz.tbi
rm $PWD/$1/cohort.filtered.hf.gt.mixed.vcf.gz.tbi
rm $PWD/$1/cohort.filtered.hf.tmp.gt.vcf.gz.tbi
rm $PWD/$1/cohort.filtered.hf.gt.vcf.gz.tbi
rm $PWD/$1/cohort.filtered.hf.DP10.GQ20.allele.trf.gt.vcf.gz.tbi
rm $PWD/$1/cohort.filtered.hf.DP10.GQ20.allele.trf.sc.tmp.gt.vcf.gz.tbi
rm $PWD/$1/cohort.filtered.hf.DP10.GQ20.allele.trf.sc.gt.vcf.gz.tbi

#-------------------
# Variant annotation
#-------------------
java -jar $PWD/../resources/snpEff/snpEff.jar ann -fi $PWD/../resources/CDC317/cpar_v3_no_mito.genes.positions.bed -v CDC317 $PWD/$1/cohort.filtered.gt.vcf.gz -c ../snpEff_database/snpEff.config > cdc.genes.variants.         annotation.vcf
bgzip < cdc.genes.variants.annotation.vcf > cdc.genes.variants.annotation.vcf.gz
tabix -fp vcf cdc.genes.variants.annotation.vcf.gz
../gatk-4.1.9.0/gatk VariantsToTable \
    -V cdc.genes.variants.annotation.vcf.gz \
    -F CHROM -F POS -F ID -F TYPE -F REF -F ALT -F HET -F HOM-REF -F HOM-VAR -F NO-CALL -F VAR -F NSAMPLES -F NCALLED -F MULTI-ALLELIC -F ANN -GF AD -GF GT \
    --error-if-missing-data \
    -O cdc.genes.variants.annotation.txt

