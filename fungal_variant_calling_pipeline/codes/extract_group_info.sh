# iterate all files in the folder specified by the first argument ($1)
for path in $PWD/$1/*; do  # note that path is the absolute directory
    # if not a directory, skip
    [ -d "${path}" ] || continue
    base_folder=`basename $path`

    # the second argment allows you to run the pipeline for specific cohort with prefix defined by the second argument ($2)
    if [ ! -z "$2" ]
    then
        # if the folder name does not start with $2, skip
        [[ $base_folder == $2* ]] || continue
    fi

    # enter the folder
    cd $path

    # get sample names
    Reads1=(*R1.fastq.gz)
    Reads2=(*R2.fastq.gz)
    base_fastq=${Reads1%_R*}

    # get read group information
    header=$(zcat < $Reads1 | head -n 1)
    id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g')
    sm=$(echo $header | head -n 1 | grep -Eo "[ATGCN]+\+[ATGCN]+$" | sed 's/+/_/g')
    echo "$base_fastq,$id"_"$sm"
done
