#!/bin/bash

# phage databases
GPD_db=/home/team/phage_databases/GutPhageDatabase/GPD_sequences.fa
MVP_db=/home/team/phage_databases/MicrobeVersusPhage/mvp_viral_cluster_representative_seqs.fasta
RVDB_db=/home/team/phage_databases/RVDB/C-RVDBv22.0.fasta

# plasmid databases
COMPASS_db=/home/team/plasmid_databases/COMPASS/COMPASS.fasta
PLSDB_db=/home/team/plasmid_databases/PLSDB/plsdb.fna

INDIR=/home/team/Bfragilis/spacer-graph-new/
OUTDIR=/home/team/Bfragilis_spacer_analysis/blast_outputs/

fasta=/home/team/Bfragilis/spacer-graph-new/spacer_nr.fa

echo "Blasting $header against phage"

# blast pphage
blastn -query $fasta -db "${GPD_db} ${MVP_db} ${RVDB_db} ${mMGE_phage}" -evalue 0.001 -outfmt "6 qseqid sseqid pident nident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs"  -num_threads 8 -out ${OUTDIR}/spacer_nr_phages.m8

echo "Blasting $header against plasmids"

# blast plasmids
blastn -query $fasta -db "${COMPASS_db} ${PLSDB_db} ${mMGE_plasmid}" -evalue 0.001 -outfmt "6 qseqid sseqid pident nident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs"  -num_threads 8 -out ${OUTDIR}/spacer_nr_plasmids.m8
