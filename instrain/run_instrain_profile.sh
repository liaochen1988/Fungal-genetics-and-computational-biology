for path in */; do
    # iterate each sample folder
    [ -d "${path}" ] || continue # if not a directory, skip
    [[ $path == Sample* ]] || continue # if not start with "Sample", skip
    echo $path
    cd $path

    # get bam file with marked duplicates
    bamfile=(*.marked_duplicates.bam)
    base=${bamfile%.marked_duplicates.bam} # remove the shortest path in Reads1 after _R*
    echo $bamfile
    echo $base

    # remove $base.ID in case the task fails
    if [[ -d "$base.IS" ]]; then
        rm -r $base.IS
     fi

    # continue if $base.ID (a folder) exists
    if [[ -d "$base.IS" ]]; then
      echo "$base.IS exists."
      cd ..
      continue
    fi

    # run instrain
    inStrain profile $bamfile ../C_parapsilosis_CDC317_current_chromosomes.fasta -o $base.IS -p 36 --skip_mm_profiling &

    cd ..
done
