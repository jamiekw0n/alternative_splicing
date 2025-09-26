#!/usr/bin/env python3
import argparse, re
from pathlib import Path
from collections import defaultdict
import pysam

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--bams", required=True)
    ap.add_argument("--bp", type=int, default=3)
    ap.add_argument("--out", required=True)
    return ap.parse_args()

def load_model_junctions(gtf):
    by_tx = defaultdict(list)
    with open(gtf) as f:
        for line in f:
            if not line.strip() or line.startswith("#"): continue
            chrom, src, feat, start, end, score, strand, frame, attrs = line.rstrip().split("\t")
            if feat != "exon": continue
            mtx = re.search(r'transcript_id "([^"]+)"', attrs)
            if not mtx: continue
            by_tx[(chrom, strand, mtx.group(1))].append((int(start), int(end)))
    jset = defaultdict(set)  # chrom -> set((donor, acceptor))
    for (chrom, strand, tx), exons in by_tx.items():
        exons = sorted(exons)
        for i in range(len(exons)-1):
            d = exons[i][1]; a = exons[i+1][0]
            jset[chrom].add((d, a))
    return jset

def bam_junctions(bam):
    obs = defaultdict(set)
    with pysam.AlignmentFile(bam, "rb") as b:
        for r in b.fetch(until_eof=True):
            if r.is_unmapped: continue
            chrom = b.get_reference_name(r.reference_id)
            refpos = r.reference_start
            for op, ln in (r.cigartuples or []):
                if op in (0,7,8,2):  # M, =, X, D advance ref
                    refpos += ln
                elif op == 3:  # N
                    donor = refpos
                    acceptor = refpos + ln
                    obs[chrom].add((donor, acceptor))
                    refpos += ln
    return obs

def main():
    a = parse_args()
    model = load_model_junctions(a.gtf)
    observed = defaultdict(set)
    for bam in [x for x in a.bams.split(",") if x.strip()]:
        bj = bam_junctions(bam)
        for c, ss in bj.items():
            observed[c].update(ss)
    outp = Path(a.out); outp.parent.mkdir(parents=True, exist_ok=True)
    with open(outp, "w") as w:
        w.write("chrom	donor	acceptor	matched_model")
        for c, js in observed.items():
            for (d, ac) in sorted(js):
                is_match = any(abs(d - md) <= a.bp and abs(ac - ma) <= a.bp for (md, ma) in model.get(c, set()))
                w.write(f"{c}	{d}	{ac}	{int(is_match)}\n")

if __name__ == "__main__":
    main()
