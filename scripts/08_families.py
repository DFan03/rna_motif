#!/usr/bin/env python3
"""Merge replicated motifs into families (Tomtom, MEME Suite) and rank families on separate, visible criteria (C-20, C-21)."""
import argparse, re, subprocess, collections
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
BIN = "/opt/anaconda3/envs/rnamotif/bin"
ap = argparse.ArgumentParser(); ap.add_argument("--dir", required=True); ap.add_argument("--logo", required=True)
ap.add_argument("--out", required=True); ap.add_argument("--q", type=float, default=0.05); ap.add_argument("--min-overlap", type=int, default=6)
a = ap.parse_args()
d, logo, out = Path(a.dir), Path(a.logo), Path(a.out); out.mkdir(parents=True, exist_ok=True)
val = pd.read_csv(logo/"logo_validation.tsv", sep="\t")
split = pd.read_csv(d/"split.tsv", sep="\t")
HDR = "MEME version 5\n\nALPHABET= ACGU\n\nstrands: +\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 U 0.25\n\n"

def blocks(path):
    L = Path(path).read_text().splitlines(); res = {}
    for i, l in enumerate(L):
        m = re.match(r"MOTIF (\S+) (\S+)", l)
        if m:
            w = int(re.search(r"w= (\d+)", L[i+1]).group(1)); res[m.group(2)] = (L[i+1], L[i+2:i+2+w])
    return res

def consensus(rows):
    iu = {"A":"A","C":"C","G":"G","U":"U","AG":"R","CU":"Y","CG":"S","AU":"W","GU":"K","AC":"M"}
    s = ""
    for r in rows:
        p = np.array(r.split(), float); order = np.argsort(-p)
        top = p[order[0]]
        if top >= 0.6: s += "ACGU"[order[0]]
        else:
            two = "".join(sorted("ACGU"[j] for j in order[:2])); s += iu.get(two, "N") if p[order[0]] + p[order[1]] >= 0.8 else "N"
    return s

fam_rows, site_rows = [], []
for reg in ["utr3", "utr5"]:
    rep = val[(val.region == reg) & val.replicated].copy()
    if rep.empty: continue
    rep["uid"] = rep.heldout.str.replace(" ", "_") + "__" + rep.motif_id
    wdir = out / reg; wdir.mkdir(exist_ok=True); pooled = wdir / "pooled.meme"; mats = {}
    with open(pooled, "w") as f:
        f.write(HDR)
        for _, r in rep.iterrows():
            hdr, rows = blocks(logo/"work"/f"{reg}_{r.heldout}"/"streme"/"streme.txt")[r.motif_id]
            mats[r.uid] = rows; f.write(f"MOTIF {r.uid}\n{hdr}\n" + "\n".join(rows) + "\n\n")
    subprocess.run([f"{BIN}/tomtom", "-oc", str(wdir/"tomtom"), "-thresh", str(a.q), "-min-overlap", str(a.min_overlap), "-norc",
                    "-verbosity", "1", str(pooled), str(pooled)], check=True, capture_output=True)
    t = pd.read_csv(wdir/"tomtom"/"tomtom.tsv", sep="\t", comment="#").dropna(subset=["Target_ID"])
    t = t[t["q-value"] < a.q]; ok = set(zip(t.Query_ID, t.Target_ID))
    # edge only if significant in BOTH directions (conservative; limits chaining)
    parent = {u: u for u in rep.uid}
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for q_, t_ in ok:
        if q_ != t_ and (t_, q_) in ok: parent[find(q_)] = find(t_)
    rep["fam"] = rep.uid.map(find)
    # family-level table
    for fam, g in rep.groupby("fam"):
        best = g.sort_values(["discovery_E", "clusters_with_hit"], ascending=[True, False]).iloc[0]
        fam_rows.append(dict(region=reg, family_key=fam, n_motifs=len(g), n_groups_replicated=g.heldout.nunique(),
            groups_replicated=",".join(sorted(g.heldout.unique())), best_discovery_E=g.discovery_E.min(), best_fdr=g.fdr_bh.min(),
            rep_uid=best.uid, rep_consensus=consensus(mats[best.uid]), member_consensus=";".join(sorted(set(g.motif))[:6]),
            rep_fold=best.heldout, rep_motif_id=best.motif_id, rep_region=reg))
F = pd.DataFrame(fam_rows); F["fam_idx"] = F.index   # index used for the per-family fimo_fam{idx} folders

# --- descriptive coverage of each family's representative on cluster representatives of ALL groups (includes discovery data)
reps90 = {}
for reg in ["utr3", "utr5"]:
    p = out / f"{reg}_all_reps.fa"
    with open(p, "w") as o:
        for s_ in ("train", "heldout"): o.write((d/"regions"/f"{reg}_{s_}_reps90.fa").read_text())
    reps90[reg] = p
grp = split.set_index("acc").group
def read_rfam(path):
    rows = []
    for l in Path(path).read_text().splitlines():
        if l.startswith("#") or not l.strip(): continue
        p = l.split(); rows.append((p[1], p[3], min(int(p[9]), int(p[10])), max(int(p[9]), int(p[10]))))
    return pd.DataFrame(rows, columns=["family", "seq", "s", "e"])
RF = read_rfam(logo/"cmscan_utr3_all.tbl")

cov_cols = []
for i, r in F.iterrows():
    reg = r.region; wdir = out / reg; sel = wdir / f"rep_{i}.meme"
    hdr, rows = blocks(logo/"work"/f"{reg}_{r.rep_fold}"/"streme"/"streme.txt")[r.rep_motif_id]
    sel.write_text(HDR + f"MOTIF fam{i}\n{hdr}\n" + "\n".join(rows) + "\n")
    subprocess.run([f"{BIN}/fimo", "--oc", str(wdir/f"fimo_fam{i}"), "--thresh", "1e-4", "--norc", "--verbosity", "1", str(sel), str(reps90[reg])],
                   check=True, capture_output=True)
    try: h = pd.read_csv(wdir/f"fimo_fam{i}"/"fimo.tsv", sep="\t", comment="#").dropna(subset=["sequence_name"])
    except pd.errors.EmptyDataError: h = pd.DataFrame(columns=["sequence_name", "start", "stop"])
    seqs = [l[1:].strip() for l in reps90[reg].read_text().splitlines() if l.startswith(">")]
    slen = {}; cur = None
    for l in reps90[reg].read_text().splitlines():
        if l[0] == ">": cur = l[1:].strip(); slen[cur] = 0
        else: slen[cur] += len(l.strip())
    hit = set(h.sequence_name); byg = collections.defaultdict(lambda: [0, 0])
    for s in seqs: byg[grp[s]][1] += 1; byg[grp[s]][0] += s in hit
    for gname, (k, n) in byg.items(): F.loc[i, f"cov_{gname}"] = k / n
    F.loc[i, "clusters_with_hit"] = len(hit); F.loc[i, "clusters_total"] = len(seqs)
    F.loc[i, "groups_with_hit"] = sum(1 for k, n in byg.values() if k > 0)
    if len(h):
        rel = ((h.start + h.stop) / 2) / h.sequence_name.map(slen); F.loc[i, "rel_pos_median"] = rel.median(); F.loc[i, "rel_pos_iqr"] = rel.quantile(.75) - rel.quantile(.25)
    if reg == "utr3" and len(h):
        ov = collections.Counter(); tot = 0
        for _, s in h.iterrows():
            tot += 1; rf = RF[(RF.seq == s.sequence_name) & (RF.s <= s.stop) & (RF.e >= s.start)]
            for fam_ in set(rf.family): ov[fam_] += 1
        F.loc[i, "frac_sites_in_rfam"] = sum(1 for _, s in h.iterrows() if len(RF[(RF.seq == s.sequence_name) & (RF.s <= s.stop) & (RF.e >= s.start)])) / tot
        F.loc[i, "rfam_families"] = ",".join(f"{k}:{v}" for k, v in ov.most_common(3))
F = F.sort_values(["n_groups_replicated", "groups_with_hit", "clusters_with_hit", "best_discovery_E"], ascending=[False, False, False, True]).reset_index(drop=True)
F.insert(0, "rank", F.index + 1)
keep = ["rank", "fam_idx", "rep_uid", "region", "rep_consensus", "n_motifs", "n_groups_replicated", "groups_replicated", "clusters_with_hit", "clusters_total", "groups_with_hit",
        "best_discovery_E", "best_fdr", "rel_pos_median", "rel_pos_iqr", "frac_sites_in_rfam", "rfam_families", "member_consensus"] + [c for c in F if c.startswith("cov_")]
F[[c for c in keep if c in F]].to_csv(out/"families_ranked.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); print(F[["rank","region","rep_consensus","n_motifs","n_groups_replicated","clusters_with_hit","clusters_total","groups_with_hit","best_discovery_E","rel_pos_median","frac_sites_in_rfam","rfam_families"]].to_string(index=False))
print(f"\n{int(val.replicated.sum())} replicated motifs -> {len(F)} families")
# heatmap of coverage by group
cc = [c for c in F if c.startswith("cov_")]
fig, ax = plt.subplots(figsize=(7, 0.35*len(F) + 1.5)); im = ax.imshow(F[cc].fillna(0).values, vmin=0, vmax=1, cmap="viridis", aspect="auto")
ax.set_xticks(range(len(cc))); ax.set_xticklabels([c[4:] for c in cc], rotation=30, ha="right")
ax.set_yticks(range(len(F))); ax.set_yticklabels([f"{r['rank']}. {r.region} {r.rep_consensus}" for _, r in F.iterrows()], fontsize=8)
plt.colorbar(im, label="fraction of clusters with a site"); ax.set_title("Family coverage by serogroup (descriptive)")
plt.tight_layout(); plt.savefig(out/"family_coverage.png", dpi=150)
