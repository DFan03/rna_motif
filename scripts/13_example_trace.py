#!/usr/bin/env python3
"""Trace ONE record through the whole pipeline and print the real values (used verbatim in the README worked example)."""
import re, sys
from pathlib import Path
import pandas as pd
from Bio import SeqIO
acc = sys.argv[1] if len(sys.argv) > 1 else "MN566104.1"
D = Path("data/processed/orthoflavi_v1")
rec = next(r for f in Path("data/raw/orthoflavi_v1").glob("*.gb") for r in SeqIO.parse(f, "genbank") if r.id == acc)
src = next(f for f in rec.features if f.type == "source").qualifiers
print("== 1 PARSE (GenBank record -> fields)")
print(f"LOCUS {rec.name} {len(rec.seq)} bp {rec.annotations['molecule_type']} | {rec.description}")
print({k: src[k][0] for k in ("organism", "host", "collection_date", "isolate", "strain") if k in src}, "| taxon:", [x for x in src["db_xref"] if x.startswith("taxon")])
cds = max([f for f in rec.features if f.type == "CDS"], key=lambda f: len(f.location))
print(f"longest CDS (polyprotein): {int(cds.location.start)+1}..{int(cds.location.end)}  -> 5'UTR {int(cds.location.start)} nt, 3'UTR {len(rec.seq)-int(cds.location.end)} nt")
seq = str(rec.seq).upper(); print("non-ACGT bases:", sum(c not in "ACGT" for c in seq), "| T->U applied for RNA")
m = pd.read_csv(D/"metadata_clean.tsv", sep="\t"); sp = m[m.species_query == m[m.acc == acc].species_query.iloc[0]]
import numpy as np
med = sp.length.median(); scale = max((sp.length - med).abs().median() / 0.6745, 50); z = (m[m.acc == acc].length.iloc[0] - med) / scale
print("\n== 2 QC checks")
print(f"length {len(seq)} nt vs species median {med:.0f}: robust z = {z:.2f} (limit 3.5) -> pass | ambiguity 0 % (limit 1 %) -> pass | unique sequence hash -> pass | has CDS and UTR flanks >= 10 nt -> pass")
log = pd.read_csv("docs/qc_log_orthoflavi_v1.tsv", sep="\t"); print("entries for this accession in qc_log:", log[log.accession == acc].reason.tolist() or "none (kept, no flags)")
s = pd.read_csv(D/"split.tsv", sep="\t"); r = s[s.acc == acc].iloc[0]; cl = pd.read_csv(D/"clusters/clusters_id90.tsv", sep="\t")
c = cl[cl.acc == acc].iloc[0]; print("\n== 3 CLUSTER (MMseqs2, >=90 % identity)")
print(f"cluster representative {c.rep}, cluster size {c.cluster_size} -> this record's weight = 1/{c.cluster_size} = {c.weight:.3f}; members: {', '.join(cl[cl.rep == c.rep].acc.head(6))} ...")
print("\n== 4 SPLIT"); print(f"serogroup = {r.group}; fixed split (seed 42): {r.split}; in leave-one-group-out fold 'hold out {r.group}' this record is TEST data, in the 4 other folds it is TRAINING data")
print("\n== 5 REGIONS"); print(f"5'UTR {int(cds.location.start)} nt | CDS {len(cds.location)} nt | 3'UTR {len(rec.seq)-int(cds.location.end)} nt (file data/processed/orthoflavi_v1/regions/utr3_*_reps90.fa)")
F = pd.read_csv("results/families_v1/families_ranked.tsv", sep="\t"); print("\n== 6 MOTIF SITES FOUND IN THIS 3'UTR (family representatives, FIMO p<1e-4)")
rows = []
for _, f in F[F.region == "utr3"].iterrows():
    p = Path(f"results/families_v1/utr3/fimo_fam{int(f.fam_idx)}/fimo.tsv")
    try: h = pd.read_csv(p, sep="\t", comment="#").dropna(subset=["sequence_name"])
    except Exception: continue
    for _, x in h[h.sequence_name == acc].iterrows(): rows.append((int(f["rank"]), f.rep_consensus, int(x.start), int(x.stop), x.matched_sequence, x["p-value"]))
for t in sorted(rows, key=lambda t: t[2]): print(f"family #{t[0]} {t[1]:14s} 3'UTR nt {t[2]}-{t[3]}  matched {t[4]}  p={t[5]:.1e}")
v = pd.read_csv("results/logo_v2/logo_validation.tsv", sep="\t"); top = F[(F.region == "utr3")].iloc[0]; fold, mid = top.rep_uid.split("__")
z = v[(v.region == "utr3") & (v.heldout == fold.replace("_", " ")) & (v.motif_id == mid)].iloc[0]
print("\n== 7 EVIDENCE FOR THE #1 FAMILY (held out:", z.heldout + ")")
print(f"discovered on the other 4 groups with STREME E={z.discovery_E:g}; in held-out {z.heldout}: {int(z.clusters_with_hit)}/{int(z.n_heldout_clusters)} clusters have a site vs {z.null_mean:.2f} expected in dinucleotide shuffles; empirical p={z.emp_p:.4f}, BH FDR={z.fdr_bh:.3f}, replicated={z.replicated}")
