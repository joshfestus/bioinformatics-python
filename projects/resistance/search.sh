#!/usr/bin/env bash
# Search every protein in the genome against the reference catalog of
# resistance proteins. The catalog is small, so it becomes the database
# and the genome's 5,779 proteins are the queries.
set -euo pipefail
makeblastdb -in amr_proteins.fasta -dbtype prot -out amr_db > /dev/null
blastp -query proteins.fasta -db amr_db \
    -outfmt "6 qseqid sseqid pident length qlen slen qstart qend \
              sstart send evalue bitscore" \
    -evalue 1e-10 -max_target_seqs 5 -num_threads 2 > hits.tsv
wc -l < hits.tsv
