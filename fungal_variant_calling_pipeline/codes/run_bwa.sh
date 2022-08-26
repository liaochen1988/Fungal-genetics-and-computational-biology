# iterate all files in the folder specified by the first argument ($1)
curr_folder=$PWD
for path in $curr_folder/$1/*; do  # note that path is the absolute directory
    # if not a directory, skip
    [ -d "${path}" ] || continue
    base_folder=`basename $path`

    # the second argment allows you to run the pipeline for specific cohort with prefix defined by the second argument ($2)
    if [ ! -z "$4" ]
    then
        # if the folder name does not start with $2, skip
        [[ $base_folder == $4* ]] || continue
    fi

    # enter the folder
    cd $path

    # get pair-end read samples
    Reads1=(*R1.fastq.gz)
    Reads2=(*R2.fastq.gz)
    base_fastq=${Reads1%_R*} # remove the shortest path in Reads1 after _R*

    # continue if $base_fastq.sam.gz exists
    if test -f "$curr_folder/$2/$base_folder/$base_fastq.sam.gz"
    then
        echo "$base_fastq.sam.gz exists."
        continue
    else
        mkdir $curr_folder/$2/$base_folder
    fi

    # get read group information
    header=$(zcat < $Reads1 | head -n 1)
    id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g')
    sm=$(echo $header | head -n 1 | grep -Eo "[ATGCN]+\+[ATGCN]+$" | sed 's/+/_/g')
    #echo "Read Group @RG\tID:$id\tSM:$id"_"$sm\tLB:$id"_"$sm\tPL:ILLUMINA"

    # run bwa-mem
    bwa mem \
     -M \
     -t $3 \
     -R $(echo "@RG\tID:$id\tSM:$id"_"$sm\tLB:$id"_"$sm\tPL:ILLUMINA") \
     $curr_folder/../resources/CDC317/cpar_v3_inc_mito.fasta \
     $Reads1 $Reads2 | gzip > $curr_folder/$2/$base_folder/$base_fastq.sam.gz

    # sort sam file
    samtools sort -O bam -T $base_fastq.sort -o $curr_folder/$2/$base_folder/$base_fastq.sort.bam -@ $3 $curr_folder/$2/$base_folder/$base_fastq.sam.gz

    # mark duplicates
    java -jar $curr_folder/../resources/GATK/picard.jar MarkDuplicates I=$curr_folder/$2/$base_folder/$base_fastq.sort.bam O=$curr_folder/$2/$base_folder/$base_fastq.marked_duplicates.bam M=$curr_folder/$2/$base_folder/$base_fastq.marked_dup_metrics.txt

done
