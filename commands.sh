# Concatenate fastas for all the tissues
indir=/oak/stanford/groups/quake/jkwon/isoseq
outfile=/oak/stanford/groups/quake/jkwon/isoseq/all_tissues.fasta

> "$outfile"
for f in "$indir"/*.fasta; do
    cat "$f"
    echo
done >> "$outfile"

# pbmm2 alignment with pigeon environment
pbmm2 align /oak/stanford/groups/quake/mmantri/group.quake/references/gencode_GRCh38.p13_release41 /oak/stanford/groups/quake/jkwon/isoseq/all_tissues.fasta /oak/stanford/groups/quake/jkwon/isoseq/all_tissues_mapped.bam --sort --preset ISOSEQ -j 24 --report-json mapping_stats.report.json

# isoseq collapse with pigeon environment
isoseq collapse /oak/stanford/groups/quake/jkwon/isoseq/all_tissues_mapped.bam /oak/stanford/groups/quake/jkwon/isoseq/all_tissues_mapped.bam/concatenated_collapsed.gff
