 command="inStrain compare -i"

 for path in */; do
     # iterate each sample folder
     [ -d "${path}" ] || continue # if not a directory, skip
     [[ $path == Sample* ]] || continue # if not start with "Sample", skip
     echo $path
     cd $path

     # get base name from bamfile
     bamfile=(*.marked_duplicates.bam)
     base=${bamfile%.marked_duplicates.bam} # remove the shortest path in Reads1 after _R*

     # continue if $base.ID (a folder) exists
     if [[ -d "$base.IS" ]]; then
         echo "$base.IS exists."
         command="$command $path$base.IS/"
     fi

     cd ..
 done

 command="$command -p 36"
 echo $command
 eval $command
