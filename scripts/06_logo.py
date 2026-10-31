#!/usr/bin/env python3
"""Leave-one-serogroup-out motif discovery + replication (decisions C-15..C-18).
For each serogroup g: STREME on all sequences of the other groups (discovery); FIMO scan of the discovered motifs on
the held-out group's cluster representatives (one per 90% cluster); empirical p-value vs 200 dinucleotide-shuffled
copies (MEME Suite fasta-shuffle-letters -kmer 2); Benjamini-Hochberg across all tests (statsmodels)."""
import argparse, collections, subprocess, sys, shutil
from pathlib import Path
import pandas as pd, numpy as np
from Bio import SeqIO
from statsmodels.stats.multitest import multipletests
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rnamotif.streme import run_streme
BIN = "/opt/anaconda3/envs/rnamotif/bin"

ap = argparse.ArgumentParser()
ap.add_argument("--dir", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--regions", nargs="+", default=["utr5", "utr3"])
ap.add_argument("--n-shuffle", type=int, default=200); ap.add_argument("--fimo-p", type=float, default=1e-4)
ap.add_argument("--seed", type=int, default=1); a = ap.parse_args()
d, out = Path(a.dir), Path(a.out); out.mkdir(parents=True, exist_ok=True)

meta = pd.read_csv(d/"metadata_clean.tsv", sep="\t").merge(pd.read_csv(d/"split.tsv", sep="\t")[["acc","group","cluster90"]], on="acc")
meta = meta[meta.cds_start.notna() & meta.cds_strand.eq(1)].set_index("acc")
reps = set(pd.read_csv(d/"clusters/clusters_id90.tsv", sep="\t").rep)
seqs = {r.id: str(r.seq) for r in SeqIO.parse(d/"viral_clean.fasta", "fasta")}
def region(acc, reg):
    s, a0, a1 = seqs[acc], int(meta.cds_start[acc]), int(meta.cds_end[acc])
    return {"utr5": s[:a0], "utr3": s[a1:]}[reg]
def write(path, items):
    with open(path, "w") as f:
        for k, s in items:
            if len(s) >= 30: f.write(f">{k}\n{s}\n")
def run(cmd): 
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode: raise RuntimeError(" ".join(map(str, cmd)) + "\n" + r.stderr[-600:])

rows_m, rows_v = [], []
for reg in a.regions:
    for g in sorted(meta.group.unique()):
        w = out / "work" / f"{reg}_{g}"; w.mkdir(parents=True, exist_ok=True)
        tr = meta[meta.group != g].index; te = [x for x in meta[meta.group == g].index if x in reps]
        write(w/"train.fa", [(x, region(x, reg)) for x in tr])
        write(w/"test_reps.fa", [(x, region(x, reg)) for x in te])
        motifs = run_streme(w/"train.fa", w/"streme", seed=a.seed)
        for m in motifs: rows_m.append(dict(region=reg, heldout=g, n_train_seq=len(tr), **m))
        if not motifs: continue
        n_te = sum(1 for _ in SeqIO.parse(w/"test_reps.fa", "fasta"))
        # shuffled copies of the held-out representatives (dinucleotide-preserving)
        run([f"{BIN}/fasta-shuffle-letters", "-kmer", "2", "-rna", "-copies", str(a.n_shuffle), "-seed", str(a.seed), str(w/"test_reps.fa"), str(w/"shuf.fa")])
        def scan(fa, tag):
            run([f"{BIN}/fimo", "--oc", str(w/f"fimo_{tag}"), "--thresh", str(a.fimo_p), "--norc", "--verbosity", "1",
                 str(w/"streme"/"streme.txt"), str(fa)])
            try: t = pd.read_csv(w/f"fimo_{tag}"/"fimo.tsv", sep="\t", comment="#")
            except pd.errors.EmptyDataError: return pd.DataFrame(columns=["motif_alt_id", "sequence_name"])
            return t.dropna(subset=["sequence_name"]) if len(t) else t
        obs = scan(w/"test_reps.fa", "obs"); shf = scan(w/"shuf.fa", "shuf")
        # fimo names shuffled copies like "<acc>_shuf_<i>" -> copy index i
        if len(shf): shf["copy"] = shf.sequence_name.str.extract(r"_shuf_(\d+)$")[0]
        for m in motifs:
            o = obs[obs.motif_alt_id == m["id"]].sequence_name.nunique() if len(obs) else 0
            per_copy = (shf[shf.motif_alt_id == m["id"]].groupby("copy").sequence_name.nunique().reindex(
                        [str(i) for i in range(1, a.n_shuffle + 1)], fill_value=0) if len(shf) else pd.Series([0]*a.n_shuffle))
            p = (1 + (per_copy >= o).sum()) / (1 + a.n_shuffle)
            rows_v.append(dict(region=reg, heldout=g, motif=m["consensus"], motif_id=m["id"], discovery_E=m["evalue"],
                               n_heldout_clusters=n_te, clusters_with_hit=o, null_mean=float(per_copy.mean()), emp_p=p))
M = pd.DataFrame(rows_m); V = pd.DataFrame(rows_v)
if len(V):
    V["fdr_bh"] = multipletests(V.emp_p, method="fdr_bh")[1]
    V["replicated"] = (V.fdr_bh < 0.05) & (V.clusters_with_hit >= 2)
M.to_csv(out/"logo_discovered_motifs.tsv", sep="\t", index=False); V.to_csv(out/"logo_validation.tsv", sep="\t", index=False)
print(f"discovered {len(M)} motifs; validated {len(V)}; replicated {int(V.replicated.sum()) if len(V) else 0}")
if len(V): print(V.sort_values("emp_p").head(15).to_string(index=False))
