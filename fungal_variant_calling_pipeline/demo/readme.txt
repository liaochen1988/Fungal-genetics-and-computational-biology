Please follow the steps below:
1. For all fastq files, name them as *R1.fastq.gz and *R2.fastq.gz.
2. Put the forward and reverse reads of each sample in separate folders. Then put all different folders (samples) under the directory "raw_fastq".
3. Create an output folder with name "outputs"
4. Run call_variants.sh

Key output fies:
1. Filtered genetic variant profiles: cohort.filtered.gt.vcf/.vcf.gz/.txt
2. Annotated genetic variant profiles: cdc.genes.variants.annotation.vcf/.vcf.gz/.txt
3. Copy number of ORFs: ORF_v3_copy_number.csv
4. Read group information: read_group_info.csv

Other notes:
1. Read group is identified by a number of tags that allow us to differentiate not only samples, but also various technical features that are associated with artifacts. With this information, we can mitigate the effects of those artifacts during the duplicate marking and base recalibration steps.
