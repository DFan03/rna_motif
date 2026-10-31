#!/usr/bin/env python3
"""Positive control (C-19): compare replicated motifs with known flavivirus elements from Rfam (Infernal cmscan, --cut_ga).
(a) presence of each Rfam family per serogroup in the 3'UTRs; (b) overlap of replicated-motif FIMO sites with Rfam hits
vs the chance expectation (fraction of 3'UTR nt covered by Rfam hits); (c) recall: share of Rfam-covered nt hit by a replicated motif."""
import argparse, subprocess, collections
from pathlib import Path
import pandas as pd
BIN = "/opt/anaconda3/envs/rnamotif/bin"
ap = argparse.ArgumentParser(); ap.add_argument("--dir", required=True); ap.add_argument("--logo", required=True)
ap.add_argument("--cms", default="reference/known_elements"); a = ap.parse_args()
cm = Path(a.cms); d = Path(a.dir); logo = Path(a.logo)
db = cm / "flavi_rfam.cm"
with open(db, "w") as o:
    for f in sorted(cm.glob("RF*.cm")): o.write(f.read_text())
subprocess.run([f"{BIN}/cmpress", "-F", str(db)], check=True, capture_output=True)

def cmscan(fa, tag):
    tbl = logo / f"cmscan_{tag}.tbl"
    r = subprocess.run([f"{BIN}/cmscan", "--cut_ga", "--rfam", "--nohmmonly", "--fmt", "2", "--tblout", str(tbl), "-o", "/dev/null",
                        "--cpu", "4", str(db), str(fa)], capture_output=True, text=True)
    if r.returncode: raise RuntimeError(r.stderr[-500:])
    rows = []
    for l in tbl.read_text().splitlines():
        if l.startswith("#") or not l.strip(): continue
        p = l.split()   # fmt2: idx, target name, target acc, clan, query name, query acc, mdl, mdl from, mdl to, seq from, seq to, strand,...
        rows.append(dict(family=p[1], seq=p[3], s=min(int(p[9]), int(p[10])), e=max(int(p[9]), int(p[10])), score=float(p[16]), evalue=float(p[17])))
    return pd.DataFrame(rows, columns=["family", "seq", "s", "e", "score", "evalue"])

# (a) presence per group over ALL 3'UTRs
meta = pd.read_csv(d/"split.tsv", sep="\t")
allfa = logo / "utr3_all_seqs.fa"
with open(allfa, "w") as o:
    for f in [d/"regions"/"utr3_train_all.fa", d/"regions"/"utr3_heldout_all.fa"]: o.write(f.read_text())
hits = cmscan(allfa, "utr3_all").merge(meta[["acc", "group"]], left_on="seq", right_on="acc")
n_by_group = meta.groupby("group").acc.size()
pres = hits.groupby(["group", "family"]).seq.nunique().unstack(fill_value=0)
print("== (a) 3'UTRs with an Rfam hit (GA cutoff), by group (n sequences in parentheses)")
pres.index = [f"{g} ({n_by_group[g]})" for g in pres.index]; print(pres.to_string())
pres.to_csv(logo / "rfam_presence_by_group.tsv", sep="\t")

# (b,c) overlap with replicated motif sites in held-out representatives
val = pd.read_csv(logo / "logo_validation.tsv", sep="\t"); rep = val[(val.region == "utr3") & val.replicated]
tot_nt = cov_nt = site_total = site_overlap = rec_nt = rfam_nt = 0; per_motif = []
for g, sub in rep.groupby("heldout"):
    w = logo / "work" / f"utr3_{g}"
    fa = w / "test_reps.fa"; seqlen = {}; cur = None
    for l in fa.read_text().splitlines():
        if l[0] == ">": cur = l[1:].split()[0]; seqlen[cur] = 0
        else: seqlen[cur] += len(l.strip())
    rf = cmscan(fa, f"fold_{g.replace(' ', '_')}")
    fim = pd.read_csv(w/"fimo_obs"/"fimo.tsv", sep="\t", comment="#").dropna(subset=["sequence_name"])
    fim = fim[fim.motif_alt_id.isin(sub.motif_id)]
    covered = {s: set() for s in seqlen}; mcov = {s: set() for s in seqlen}
    for _, r in rf.iterrows(): covered[r.seq].update(range(r.s, r.e + 1))
    for _, r in fim.iterrows(): mcov[r.sequence_name].update(range(int(r.start), int(r.stop) + 1))
    tot_nt += sum(seqlen.values()); cov_nt += sum(len(v) for v in covered.values())
    rfam_nt += sum(len(v) for v in covered.values()); rec_nt += sum(len(covered[s] & mcov[s]) for s in seqlen)
    for _, r in fim.iterrows():
        site = set(range(int(r.start), int(r.stop) + 1)); site_total += 1
        site_overlap += bool(site & covered[r.sequence_name])
print("\n== (b) replicated 3'UTR motif sites (held-out reps) overlapping an Rfam-defined known element")
print(f"sites: {site_total}; overlapping a known element: {site_overlap} ({site_overlap/max(site_total,1):.1%}); "
      f"chance baseline = fraction of 3'UTR nt covered by Rfam hits: {cov_nt/max(tot_nt,1):.1%}")
print(f"== (c) recall: of Rfam-covered nt in held-out reps, {rec_nt/max(rfam_nt,1):.1%} lie under a replicated motif site")
pd.DataFrame([dict(sites=site_total, sites_overlapping_known=site_overlap, frac_overlap=site_overlap/max(site_total,1),
                   chance_frac=cov_nt/max(tot_nt,1), recall_nt=rec_nt/max(rfam_nt,1))]).to_csv(logo/"known_element_overlap.tsv", sep="\t", index=False)
