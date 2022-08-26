# iterate all files in the folder specified by the first argument ($1)
curr_folder=$PWD
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

    bedtools coverage -a $curr_folder/../resources/CDC317/coverage_ORFs_v3.txt -b $bamfile > $curr_folder/$2/$base_folder/$base_fastq.read_density.txt &
    bedtools coverage -a $curr_folder/../resources/CDC317/coverage_ORFs_v3.txt -b $bamfile -hist > $curr_folder/$2/$base_folder/$base_fastq.read_density.hist.txt &
    bedtools genomecov -ibam $bamfile > $curr_folder/$2/$base_folder/$base_fastq.genome_coverage.txt &
    bedtools genomecov -ibam $bamfile -d > $curr_folder/$2/$base_folder/$base_fastq.genome_coverage.per_base.txt &

done
