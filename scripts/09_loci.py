#!/usr/bin/env python3
"""Merge 3'UTR motif families into loci by the Rfam element their sites fall in (C-22).
A family is assigned to the Rfam family that overlaps >= --min-frac of its sites; families sharing an Rfam family form one locus.
Families without a majority Rfam overlap stay separate ('unannotated'): no principled basis to merge them."""
import argparse, collections
from pathlib import Path
import pandas as pd
ap = argparse.ArgumentParser(); ap.add_argument("--fam", required=True); ap.add_argument("--logo", required=True)
ap.add_argument("--min-frac", type=float, default=0.5); a = ap.parse_args()
fam_dir, logo = Path(a.fam), Path(a.logo)
F = pd.read_csv(fam_dir/"families_ranked.tsv", sep="\t"); F = F[F.region == "utr3"].copy()
rows = []
for l in (logo/"cmscan_utr3_all.tbl").read_text().splitlines():
    if l.startswith("#") or not l.strip(): continue
    p = l.split(); rows.append((p[1], p[3], min(int(p[9]), int(p[10])), max(int(p[9]), int(p[10]))))
RF = pd.DataFrame(rows, columns=["rfam", "seq", "s", "e"])
reps = {l[1:].strip() for l in (fam_dir/"utr3_all_reps.fa").read_text().splitlines() if l.startswith(">")}
RF = RF[RF.seq.isin(reps)]

sites = {}
for _, r in F.iterrows():
    try: h = pd.read_csv(fam_dir/"utr3"/f"fimo_fam{int(r.fam_idx)}"/"fimo.tsv", sep="\t", comment="#").dropna(subset=["sequence_name"])
    except Exception: h = pd.DataFrame(columns=["sequence_name", "start", "stop"])
    sites[r.fam_idx] = h
def overlaps(h):
    c = collections.Counter()
    for _, s in h.iterrows():
        for rf in set(RF[(RF.seq == s.sequence_name) & (RF.s <= s.stop) & (RF.e >= s.start)].rfam): c[rf] += 1
    return c
F["rfam_assigned"] = None; F["rfam_frac"] = float("nan")
for i, r in F.iterrows():
    c = overlaps(sites[r.fam_idx]); n = len(sites[r.fam_idx])
    if n and c:
        top, k = c.most_common(1)[0]
        F.loc[i, "rfam_frac"] = k / n
        if k / n >= a.min_frac: F.loc[i, "rfam_assigned"] = top
F["locus"] = F.rfam_assigned.where(F.rfam_assigned.notna(), "unannotated_" + F.rep_consensus)
out = []
for loc, g in F.groupby("locus"):
    allsites = pd.concat([sites[x].assign(fam=x) for x in g.fam_idx])
    seqs = set(allsites.sequence_name); grp_rep = set(",".join(g.groups_replicated).split(","))
    row = dict(locus=loc, n_families=len(g), families=";".join(g.rep_consensus), n_motifs=int(g.n_motifs.sum()),
               n_groups_replicated=len(grp_rep), groups_replicated=",".join(sorted(grp_rep)), clusters_with_hit=len(seqs),
               groups_with_hit=int(g.groups_with_hit.max()), best_discovery_E=g.best_discovery_E.min(),
               mean_rfam_frac=g.rfam_frac.mean(), median_site_pos_nt=allsites.assign(m=(allsites.start+allsites.stop)/2).m.median() if len(allsites) else float("nan"))
    if not loc.startswith("unannotated"):
        el = RF[RF.rfam == loc]; tot = cov = 0
        mc = {s: set() for s in el.seq}
        for _, s in allsites.iterrows():
            if s.sequence_name in mc: mc[s.sequence_name].update(range(int(s.start), int(s.stop) + 1))
        for _, e in el.iterrows():
            span = set(range(e.s, e.e + 1)); tot += len(span); cov += len(span & mc[e.seq])
        row.update(rfam_hits_on_reps=len(el), rfam_nt_covered_by_locus=cov / tot if tot else float("nan"),
                   rfam_hits_with_motif=sum(1 for _, e in el.iterrows() if set(range(e.s, e.e+1)) & mc[e.seq]) / max(len(el), 1))
    out.append(row)
L = pd.DataFrame(out).sort_values(["n_groups_replicated", "groups_with_hit", "clusters_with_hit", "best_discovery_E"], ascending=[False, False, False, True]).reset_index(drop=True)
L.insert(0, "rank", L.index + 1); L.to_csv(fam_dir/"loci_ranked.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
print(L[["rank", "locus", "n_families", "n_motifs", "n_groups_replicated", "clusters_with_hit", "groups_with_hit", "best_discovery_E", "rfam_nt_covered_by_locus", "rfam_hits_with_motif", "families"]].to_string(index=False))
print(f"\n{len(F)} 3'UTR families -> {len(L)} loci ({(~L.locus.str.startswith('unannotated')).sum()} Rfam-defined)")
print("\nfamily -> assigned element (fraction of sites):"); print(F[["rep_consensus", "rfam_assigned", "rfam_frac"]].to_string(index=False))
