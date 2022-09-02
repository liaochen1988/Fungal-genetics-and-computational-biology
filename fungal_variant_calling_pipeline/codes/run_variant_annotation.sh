awk '$3 == "gene"' $PWD/../resources/CDC317/cpar_v3.gff3 | convert2bed -i gff > $PWD/../resources/CDC317/cpar_v3.genes.positions.bed
java -jar $PWD/../resources/snpEff/snpEff.jar ann -fi $PWD/../resources/CDC317/cpar_v3.genes.positions.bed -v CDC317 $PWD/$1/cohort.filtered.gt.vcf.gz -c $PWD/../resources/snpEff/snpEff.config > $PWD/$1/cdc.genes.variants.annotation.vcf
bgzip < $PWD/$1/cdc.genes.variants.annotation.vcf > $PWD/$1/cdc.genes.variants.annotation.vcf.gz
tabix -fp vcf $PWD/$1/cdc.genes.variants.annotation.vcf.gz
$PWD/../resources/GATK/gatk-4.1.9.0/gatk VariantsToTable \
    -V $PWD/$1/cdc.genes.variants.annotation.vcf.gz \
    -F CHROM -F POS -F ID -F TYPE -F REF -F ALT -F HET -F HOM-REF -F HOM-VAR -F NO-CALL -F VAR -F NSAMPLES -F NCALLED -F MULTI-ALLELIC -F ANN -GF AD -GF GT \
    --error-if-missing-data \
    -O $PWD/$1/cdc.genes.variants.annotation.txt

