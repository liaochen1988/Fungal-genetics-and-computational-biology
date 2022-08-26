trf ../CDC317/cpar_v3_inc_mito.fasta 2 5 7 80 10 50 2000 -d -h
python3 TRFdat_to_bed.py --dat cpar_v3_inc_mito.fasta.2.5.7.80.10.50.2000.dat --bed  cpar_v3_inc_mito.fasta.2.5.7.80.10.50.2000.bed
python3 TRFdat_to_txt.py --dat cpar_v3_inc_mito.fasta.2.5.7.80.10.50.2000.dat --txt  cpar_v3_inc_mito.fasta.2.5.7.80.10.50.2000.txt
