awk '$3 == "gene"' ../CDC317/C_parapsilosis_CDC317_current_features.gff | convert2bed -i gff > CDC317.genes.positions.bed
java -jar ../snpEff_database/snpEff.jar ann -fi CDC317.genes.positions.bed -v CDC317 cohort.filtered.gt.vcf.gz -c ../snpEff_database/snpEff.config > cdc.genes.variants.annotation.vcf
bgzip < cdc.genes.variants.annotation.vcf > cdc.genes.variants.annotation.vcf.gz
tabix -fp vcf cdc.genes.variants.annotation.vcf.gz
../gatk-4.1.9.0/gatk VariantsToTable \
    -V cdc.genes.variants.annotation.vcf.gz \
    -F CHROM -F POS -F ID -F TYPE -F REF -F ALT -F HET -F HOM-REF -F HOM-VAR -F NO-CALL -F VAR -F NSAMPLES -F NCALLED -F MULTI-ALLELIC -F ANN -GF AD -GF GT \
    --error-if-missing-data \
    -O cdc.genes.variants.annotation.txt

