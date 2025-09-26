#!/usr/bin/env python3
"""
Concatenate lane-split FASTQ.GZ files per tissue using a *simple cat*
strategy that guarantees exactly one newline between files (only if needed).

Folder layout assumed:
  <run_dir>/<tissue_dir>/*_L00?_R1_*.fastq.gz  (and analogous R2)

Usage:
  python concat_fastqs_simple.py -o /out/parent [-p 8] <run1> [<run2> …]
"""

import argparse
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm


# --------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------- #
def build_concat_cmd(sources, dest_path):
    """
    Return a shell one-liner that, for each fastq.gz:
      - decompress with gzip -cd
      - check if the last byte is already '\n'
      - echo a newline only if it's missing
    and then recompress the whole stream to dest_path.
    """
    if not sources:
        return None

    srcs_str = " ".join(str(p) for p in sources)
    return (
        f"(for f in {srcs_str}; do "
        f"  gzip -cd \"$f\"; "
        f"  gzip -cd \"$f\" | tail -c1 | grep -q $'\\n' || echo; "
        f"done) | gzip -c > {dest_path}"
    )


def concat_one_tissue(tissue_dir: Path, out_parent: Path):
    tissue = tissue_dir.name
    out_dir = out_parent / tissue 
    out_dir.mkdir(parents=True, exist_ok=True)
    tqdm.write(f"Processing tissue: {tissue}")
    r1_files = sorted(tissue_dir.glob("*_L00?_R1_*.fastq.gz"))
    r2_files = sorted(tissue_dir.glob("*_L00?_R2_*.fastq.gz"))

    if not (r1_files or r2_files):
        tqdm.write(f"[SKIP] {tissue}: no FASTQ.GZ found")
        return

    # ----- R1 -----
    if r1_files:
        r1_out = out_dir / f"{tissue}_R1_all.fastq.gz"
        cmd = build_concat_cmd(r1_files, r1_out)
        subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)
        tqdm.write(f"[OK] {tissue}: {len(r1_files)} R1 file(s) → {r1_out.name}")

    # ----- R2 -----
    if r2_files:
        r2_out = out_dir / f"{tissue}_R2_all.fastq.gz"
        cmd = build_concat_cmd(r2_files, r2_out)
        subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)
        tqdm.write(f"[OK] {tissue}: {len(r2_files)} R2 file(s) → {r2_out.name}")


# --------------------------------------------------------------------- #
# CLI / main
# --------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description="FASTQ.GZ concatenation per tissue")
    ap.add_argument(
        "run_dirs", nargs="+", type=Path,
        help="Run directories that each contain tissue sub-folders"
    )
    ap.add_argument(
        "-o", "--output", required=True, type=Path,
        help="Parent directory for concatenated FASTQ.GZ files"
    )
    ap.add_argument(
        "-p", "--parallel", type=int, default=8,
        help="Thread pool size (default 8)"
    )
    args = ap.parse_args()

    # Discover tissue folders (one level below each run dir)
    tissues = []
    for run in args.run_dirs:
        if not run.is_dir():
            tqdm.write(f"[WARN] {run} is not a directory — skipping")
            continue
        for sub in run.iterdir():
            if sub.is_dir() and any(sub.glob("*.fastq.gz")):
                tissues.append(sub)

    if not tissues:
        print("No tissue directories containing *.fastq.gz were found. Nothing to do.")
        return

    tqdm.write(f"→ Found {len(tissues)} tissue directories. Starting …")

    with ThreadPoolExecutor(max_workers=args.parallel) as tp:
        futures = {tp.submit(concat_one_tissue, t, args.output): t for t in tissues}
        for _ in tqdm(as_completed(futures), total=len(futures), desc="Concat"):
            pass


if __name__ == "__main__":
    main()
