#!/bin/bash
#SBATCH --job-name=concat_fasta
#SBATCH --output=concat_fasta.out
#SBATCH --time=12:00:00
#SBATCH --mem=50G

#SBATCH --mail-user=jkwon@caltech.edu   # your email
#SBATCH --mail-type=ALL        # when to get emails (BEGIN,END,FAIL,ALL)

indir=/oak/stanford/groups/quake/jkwon/isoseq
outfile=/oak/stanford/groups/quake/jkwon/isoseq/all_tissues_slurm.fasta

> "$outfile"
for f in "$indir"/*.fasta; do
    sed -e '$a\' "$f" >> "$outfile"
done
