awk '$3 == "gene"' cpar_v3_no_mito.gff3 | convert2bed -i gff > cpar_v3_no_mito.genes.positions.bed
