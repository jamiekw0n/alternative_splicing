#!/usr/bin/env python3
import os
import sys
import csv
import shutil
import argparse
import tempfile
import subprocess
import pandas as pd
from pathlib import Path
from collections import defaultdict
from typing import Optional, Dict, List, Set

# Where the tissue-aligned BAMs live
BAM_ROOT = "/hpc/projects/genomics/TSP_kinnex_access/no_bc_rank_filter/aligned_bams"

def run(cmd: List[str], check: bool = True) -> subprocess.CompletedProcess:
    """Run a shell command and capture output."""
    print(">>", " ".join(cmd), flush=True)
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if check and proc.returncode != 0:
        print(proc.stdout)
        print(proc.stderr, file=sys.stderr)
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")
    return proc

def find_tissue_bam(tissue_dir: str) -> Optional[str]:
    """
    Given a tissue_dir like ".../TSP33_ascending-colon_23_kinnexsc",
    return the absolute path to its BAM:
        BAM_ROOT/<tissue_dir_name>/0_perc.bam
    """
    if not isinstance(tissue_dir, str) or not tissue_dir:
        return None
    tissue_name = os.path.basename(tissue_dir.rstrip("/"))
    bam_path = os.path.join(BAM_ROOT, tissue_name, "0_perc.bam")
    return bam_path

def write_readname_list(readnames: Set[str], out_path: str) -> None:
    with open(out_path, "w") as f:
        for rn in sorted(r.strip() for r in readnames if r and r.strip()):
            f.write(f"{rn}\n")

def extract_from_bam_by_names(bam_path: str, names_file: str, out_bam: str, threads: int = 4) -> bool:
    """
    Use samtools to extract reads whose QNAME is listed in names_file.
    Returns True if output BAM is non-empty, else False.
    """
    if not os.path.exists(bam_path):
        print(f"[WARN] BAM not found: {bam_path}")
        return False

    # Extract to BAM: -N name list, -b BAM output, -@ threads
    cmd = ["samtools", "view", "-@", str(threads), "-b", "-N", names_file, "-o", out_bam, bam_path]
    run(cmd, check=False)

    # Check if non-empty (has alignments)
    if not os.path.exists(out_bam) or os.path.getsize(out_bam) == 0:
        # samtools may create 0-byte file if none matched
        try:
            os.remove(out_bam)
        except FileNotFoundError:
            pass
        return False
    return True

def merge_bams(out_bam: str, bam_list: List[str], threads: int = 4) -> bool:
    """
    Merge multiple BAMs into out_bam. If there is only one, just move it.
    """
    if len(bam_list) == 0:
        return False
    if len(bam_list) == 1:
        shutil.move(bam_list[0], out_bam)
        return True

    # Try samtools merge first
    cmd = ["samtools", "merge", "-@", str(threads), "-f", out_bam] + bam_list
    proc = run(cmd, check=False)
    if proc.returncode == 0 and os.path.exists(out_bam) and os.path.getsize(out_bam) > 0:
        return True

    # Fallback: samtools cat (concatenate)
    cmd = ["samtools", "cat", "-o", out_bam] + bam_list
    run(cmd)
    return os.path.exists(out_bam) and os.path.getsize(out_bam) > 0

def index_bam(bam_path: str, threads: int = 4) -> None:
    run(["samtools", "index", "-@", str(threads), bam_path])

def main():
    ap = argparse.ArgumentParser(description="Extract molecule reads per new_id into one BAM per new_id.")
    ap.add_argument("--map_csv", required=True,
                    help="CSV from previous step with columns: new_id,old_id,tissue_dir,pb_id,molecules")
    ap.add_argument("--out_dir", required=True,
                    help="Output directory for per-new_id BAMs")
    ap.add_argument("--threads", type=int, default=4, help="Threads for samtools")
    ap.add_argument("--dry_run", action="store_true", help="Plan only; don’t run samtools")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.map_csv, dtype=str)
    required_cols = {"new_id", "tissue_dir", "molecules"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Input CSV missing columns: {missing}")

    # Build: new_id -> tissue_name -> set(readnames)
    # molecules column is comma-separated "molecule/####"
    per_newid_reads: Dict[str, Dict[str, Set[str]]] = defaultdict(lambda: defaultdict(set))

    for _, row in df.iterrows():
        new_id = row["new_id"]
        tissue_dir = row["tissue_dir"]
        mols = row.get("molecules", "") or ""

        if not isinstance(new_id, str) or not isinstance(tissue_dir, str) or not mols:
            continue

        tissue_name = os.path.basename(tissue_dir.rstrip("/"))
        # split molecule IDs
        for m in mols.split(","):
            m = m.strip()
            if not m:
                continue
            per_newid_reads[new_id][tissue_name].add(m)

    print(f"Loaded new_id → tissue → readnames for {len(per_newid_reads)} new_ids.")

    # For each new_id, extract reads from each tissue BAM and merge
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)

        for new_id, tissue_map in per_newid_reads.items():
            out_bam = out_dir / f"{new_id}.bam"
            if out_bam.exists() and out_bam.stat().st_size > 0:
                print(f"[SKIP] Output exists: {out_bam}")
                continue

            temp_bams: List[str] = []
            print(f"\n=== Processing {new_id} ===")
            for tissue_name, readnames in tissue_map.items():
                if not readnames:
                    continue

                bam_path = os.path.join(BAM_ROOT, tissue_name, "0_perc.bam")
                if not os.path.exists(bam_path):
                    print(f"[WARN] BAM missing for tissue {tissue_name}: {bam_path}")
                    continue

                # Write readnames file
                names_file = tmpdir_path / f"{new_id}.{tissue_name}.names.txt"
                write_readname_list(readnames, str(names_file))

                # Extract into temp BAM
                tissue_tmp_bam = tmpdir_path / f"{new_id}.{tissue_name}.bam"
                if args.dry_run:
                    print(f"(dry-run) would extract from {bam_path} using {names_file} -> {tissue_tmp_bam}")
                else:
                    ok = extract_from_bam_by_names(bam_path, str(names_file), str(tissue_tmp_bam), threads=args.threads)
                    if ok:
                        temp_bams.append(str(tissue_tmp_bam))
                    else:
                        print(f"[INFO] no reads matched in {tissue_name} for {new_id}")

            if args.dry_run:
                print(f"(dry-run) would merge {len(temp_bams)} temp BAMs -> {out_bam}")
                continue

            if not temp_bams:
                print(f"[INFO] No reads found for {new_id}; no output written.")
                continue

            # Merge (or move if only one)
            ok = merge_bams(str(out_bam), temp_bams, threads=args.threads)
            if not ok:
                print(f"[WARN] Failed to produce BAM for {new_id}")
                continue

            # Index final BAM
            index_bam(str(out_bam), threads=args.threads)
            print(f"[DONE] Wrote {out_bam} (+ .bai)")

    print("\nAll done.")

if __name__ == "__main__":
    main()
