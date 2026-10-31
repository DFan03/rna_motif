#!/usr/bin/env python3
"""Redundancy clustering with MMseqs2 (off-the-shelf). Runs several identity thresholds for sensitivity.
Output: clusters_id{X}.tsv (acc, cluster_rep, cluster_size, weight=1/size)."""
import argparse, subprocess, shutil, tempfile
from pathlib import Path
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--fasta", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--ids", type=float, nargs="+", default=[0.8, 0.9, 0.95])
ap.add_argument("--cov", type=float, default=0.8, help="min alignment coverage of the shorter... see --cov-mode")
ap.add_argument("--mmseqs", default="mmseqs"); a = ap.parse_args()
out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
for pid in a.ids:
    tmp = Path(tempfile.mkdtemp())
    pref = tmp / "c"
    # --dbtype 2 = nucleotide; --cov-mode 1 = coverage of the target (representative) -> allows partial genomes to join
    subprocess.run([a.mmseqs, "easy-cluster", a.fasta, str(pref), str(tmp/"t"), "--min-seq-id", str(pid),
                    "-c", str(a.cov), "--cov-mode", "1", "--dbtype", "2", "--threads", "4", "-v", "1"],
                   check=True, capture_output=True)
    cl = pd.read_csv(f"{pref}_cluster.tsv", sep="\t", header=None, names=["rep", "acc"])
    cl["cluster_size"] = cl.groupby("rep").acc.transform("size"); cl["weight"] = 1 / cl.cluster_size
    cl = cl[["acc", "rep", "cluster_size", "weight"]].sort_values(["rep", "acc"])
    cl.to_csv(out / f"clusters_id{int(pid*100)}.tsv", sep="\t", index=False)
    print(f"id>={pid}: {cl.rep.nunique()} clusters from {len(cl)} seqs; largest={cl.cluster_size.max()}")
    shutil.rmtree(tmp)
