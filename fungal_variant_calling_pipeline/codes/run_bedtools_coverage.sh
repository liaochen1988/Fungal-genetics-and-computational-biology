# iterate all files in the folder specified by the first argument ($1)
curr_folder=$PWD
command="wait"
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

    # continue if genome_coverage.txt exists
    if test -f "$curr_folder/$2/$base_folder/$base_fastq.genome_coverage.txt"
    then
        echo "$base_fastq.genome_coverage.txt exists."
        continue
    fi

    # .read_density.hist.txt and .genome_coverage.txt are needed for computing copy numbers
    # coverage analysis per gene (gene positions are listed in coverage_ORFs_v3.txt)
    bedtools coverage -a $curr_folder/../resources/CDC317/coverage_ORFs_v3.txt -b $bamfile -hist > $curr_folder/$2/$base_folder/$base_fastq.read_density.hist.txt &
    # coverage at the chromosome level
    bedtools genomecov -ibam $bamfile > $curr_folder/$2/$base_folder/$base_fastq.genome_coverage.txt &
    command="$command $!"

    # .genome_coverage.per_base.txt can be used to plot genome copy number at nucleotide level
    # coverage at the per-site level
    bedtools genomecov -ibam $bamfile -d > $curr_folder/$2/$base_folder/$base_fastq.genome_coverage.per_base.txt &
    command="$command $!"

done

echo $command
eval $command
