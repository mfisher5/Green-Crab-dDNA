#!/usr/bin/env bash
cd /mnt/nfs/home/KellyCEG/blastdb/
PATH=$PATH:/mnt/nfs/home/KellyCEG/ncbi-blast-2.13.0+/bin
BLAST_DB='/mnt/nfs/home/KellyCEG/blastdb/nt'
blast_input="/mnt/nfs/home/KellyCEG/mcf/dada2/run3_Hash_Key_clean.fasta"
blast_output='/mnt/nfs/home/KellyCEG/mcf/blast/run3_hash_key_clean_blast_2023-05-10.fasta'
# BLAST PARAMETERS
PERCENT_IDENTITY="94"   # ran at 85% for 5/10/2023
WORD_SIZE="15"
EVALUE="1e-30"
# number of matches recorded in the alignment:
MAXIMUM_MATCHES="50"
CULLING="5"
################################################################################
# BLAST CLUSTERS
################################################################################
echo $(date +%H:%M) "BLASTing..."
/mnt/nfs/home/KellyCEG/ncbi-blast-2.13.0+/bin/blastn -query "${blast_input}" -db "${BLAST_DB}" -num_threads 16 -perc_identity "${PERCENT_IDENTITY}" -word_size "${WORD_SIZE}" -evalue "${EVALUE}" -max_target_seqs "${MAXIMUM_MATCHES}" -culling_limit="${CULLING}" -outfmt "6 qseqid sseqid sacc pident length mismatch gapopen qcovus qstart qend sstart send evalue bitscore staxids qlen sscinames sseq" -out "${blast_output}"