curr_folder=$PWD
command="$curr_folder/../resources/GATK/gatk-4.1.9.0/gatk CombineGVCFs -R $curr_folder/../resources/CDC317/cpar_v3_inc_mito.fasta"

# iterate all files in the folder specified by the first argument ($1)
for path in $curr_folder/$1/*; do  # note that path is the absolute directory
    # if not a directory, skip
    [ -d "${path}" ] || continue
    base_folder=`basename $path`

    # the second argment allows you to run the pipeline for specific cohort with prefix defined by the second argument ($3)
    if [ ! -z "$3" ]
    then
        # if the folder name does not start with $3, skip
        [[ $base_folder == $3* ]] || continue
    fi

    # enter the folder
    cd $curr_folder/$2/$base_folder/

    # get bam file with marked duplicates
    bamfile=(*.marked_duplicates.bam)
    base_fastq=${bamfile%.marked_duplicates.bam}
    command="$command --variant $curr_folder/$2/$base_folder/$base_fastq.unfiltered.vcf.gz"
done

command="$command -O $curr_folder/$2/cohort.unfiltered.vcf.gz"
#echo $command
eval $command

$curr_folder/../resources/GATK/gatk-4.1.9.0/gatk --java-options "-Xms4g -Xmx196g" GenotypeGVCFs \
    -R $curr_folder/../resources/CDC317/cpar_v3_inc_mito.fasta \
    -V $curr_folder/$2/cohort.unfiltered.vcf.gz \
    -O $curr_folder/$2/cohort.unfiltered.gt.vcf.gz
gunzip -k $curr_folder/$2/cohort.unfiltered.gt.vcf.gz
