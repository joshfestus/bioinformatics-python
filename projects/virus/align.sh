#!/usr/bin/env bash
# Align every clean genome to the reference, keeping reference numbering.
set -euo pipefail
mafft --6merpair --thread 2 --keeplength --addfragments clean.fasta \
    reference.fasta > aligned.fasta
grep -c '>' aligned.fasta
