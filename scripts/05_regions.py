#!/usr/bin/env python3
"""Cut 5'UTR / CDS / 3'UTR from the longest annotated CDS and write FASTAs per split, for all sequences and for
one representative per 90% MMseqs2 cluster (the cluster representative chosen by MMseqs2)."""
import argparse
from pathlib import Path
import pandas as pd
from Bio import SeqIO
ap = argparse.ArgumentParser(); ap.add_argument("--dir", required=True); a = ap.parse_args()
d = Path(a.dir); out = d / "regions"; out.mkdir(exist_ok=True)
meta = pd.read_csv(d/"metadata_clean.tsv", sep="\t").merge(pd.read_csv(d/"split.tsv", sep="\t")[["acc","split","cluster90"]], on="acc")
meta = meta[meta.cds_start.notna() & meta.cds_strand.eq(1)].set_index("acc")   # needs a + strand CDS annotation
seqs = {r.id: str(r.seq) for r in SeqIO.parse(d/"viral_clean.fasta", "fasta")}
reps = set(pd.read_csv(d/"clusters/clusters_id90.tsv", sep="\t").rep)
handles = {}
def w(region, split, scheme, acc, s):
    if len(s) < 30: return
    k = (region, split, scheme)
    if k not in handles: handles[k] = open(out / f"{region}_{split}_{scheme}.fa", "w")
    handles[k].write(f">{acc}\n{s}\n")
for acc, r in meta.iterrows():
    s = seqs[acc]; a0, a1 = int(r.cds_start), int(r.cds_end)
    parts = {"utr5": s[:a0], "cds": s[a0:a1], "utr3": s[a1:]}
    for reg, sub in parts.items():
        w(reg, r.split, "all", acc, sub)
        if acc in reps: w(reg, r.split, "reps90", acc, sub)
for h in handles.values(): h.close()
for k in sorted(handles): 
    n = sum(1 for l in open(out / f"{k[0]}_{k[1]}_{k[2]}.fa") if l[0] == ">")
    print(k, n)
