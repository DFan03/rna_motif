#!/usr/bin/env python3
"""Phase 3 (descriptive, option D4-1): ribosome vs RNA density along ZIKV / DENV-2 genomes, with motif loci and Rfam elements marked.
NO statistical motif->TE model (n = 2 donors per virus; see docs/TE_DATA_SURVEY.md)."""
import argparse, subprocess, collections, json, re
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from Bio import SeqIO
RN = "/opt/anaconda3/envs/rnamotif/bin"; RB = "/opt/anaconda3/envs/riboseq/bin"
ap = argparse.ArgumentParser(); ap.add_argument("--bams", required=True); ap.add_argument("--ref", required=True)
ap.add_argument("--fam", required=True); ap.add_argument("--out", required=True); ap.add_argument("--psite", type=int, default=12)
a = ap.parse_args(); bams, ref, fam, out = Path(a.bams), Path(a.ref), Path(a.fam), Path(a.out); out.mkdir(parents=True, exist_ok=True)

G = {}   # genome -> dict(seq, cds_start0, cds_end)
for name in ["ZIKV_IbH30656", "DENV2_16681"]:
    r = SeqIO.read(ref/f"{name}.gb", "genbank"); c = max([f for f in r.features if f.type == "CDS"], key=lambda f: len(f.location))
    G[name] = dict(seq=str(r.seq).replace("T", "U"), L=len(r.seq), c0=int(c.location.start), c1=int(c.location.end))
virus_of = {"1": None, "2": "ZIKV_IbH30656", "3": "DENV2_16681"}

# ---------- read 5' ends / lengths from BAMs
def load(label):
    cnt = {g: {"+": collections.defaultdict(lambda n=G[g]["L"] + 60: np.zeros(n)), "-": collections.defaultdict(lambda n=G[g]["L"] + 60: np.zeros(n))} for g in G}   # n bound per genome (no late binding)
    tot = 0
    for l in subprocess.run([f"{RB}/samtools", "view", "-F", "4", str(bams/f"{label}.bam")], capture_output=True, text=True).stdout.splitlines():
        p = l.split("\t"); flag, g, pos, seq = int(p[1]), p[2], int(p[3]) - 1, p[9]; n = len(seq); tot += 1
        if flag & 16: cnt[g]["-"][n][pos + n - 1] += 1          # 5' end of a reverse read = rightmost
        else: cnt[g]["+"][n][pos] += 1
    return cnt
def total_reads(label):
    t = (bams/f"{label}.cutadapt.txt").read_text(); return int(re.search(r"Total reads processed:\s+([\d,]+)", t).group(1).replace(",", ""))
labels = sorted(p.name[:-4] for p in bams.glob("*.bam"))
D = {l: load(l) for l in labels}; N = {l: total_reads(l) for l in labels}

# ---------- 1. QC: counts, read-length profile, strand, periodicity
rows = []
for l in labels:
    donor, cond, kind = l[0], l[1], l.split("_")[1]
    for g in G:
        pl = sum(sum(v.sum() for v in D[l][g]["+"].values()) for _ in [0]); mi = sum(v.sum() for v in D[l][g]["-"].values())
        rows.append(dict(label=l, donor=donor, cond={"1": "uninfected", "2": "ZIKV", "3": "DENV2"}[cond], kind=kind, genome=g, total_reads=N[l],
                         plus=int(pl), minus=int(mi), plus_per_M=pl / N[l] * 1e6))
Q = pd.DataFrame(rows); Q.to_csv(out/"viral_read_counts.tsv", sep="\t", index=False)
print(Q.pivot_table(index=["cond", "donor", "kind"], columns="genome", values="plus", aggfunc="sum").to_string())
# periodicity of RPF P-sites in CDS (infected ribo samples, matching genome)
per = []
for l in labels:
    if not l.endswith("ribo") or l[1] == "1": continue
    g = virus_of[l[1]]; c0, c1 = G[g]["c0"], G[g]["c1"]; fr = np.zeros(3)
    for n, arr in D[l][g]["+"].items():
        if 27 <= n <= 31:
            for i in range(c0 + 30, c1 - 30): fr[(i + a.psite - c0) % 3] += arr[i]
    per.append(dict(label=l, f0=fr[0], f1=fr[1], f2=fr[2], frac_best=fr.max() / fr.sum() if fr.sum() else np.nan))
P = pd.DataFrame(per); P.to_csv(out/"periodicity.tsv", sep="\t", index=False); print("\nperiodicity (P-site offset %d, 27-31 nt):" % a.psite); print(P.to_string(index=False))
lenprof = {l: {n: int(v.sum()) for n, v in D[l][virus_of[l[1]]]["+"].items()} for l in labels if l[1] != "1" and l.endswith("ribo")}

# ---------- 2. region density ratios (RPF vs RNA), relative to CDS
def psite(l, g):
    v = np.zeros(G[g]["L"] + 60)
    for n, arr in D[l][g]["+"].items():
        if l.endswith("ribo"):
            if 27 <= n <= 31: v[a.psite:] += arr[:len(v) - a.psite]
        else: v += arr                              # RNA: 5' ends of all reads
    return v[:G[g]["L"]]
reg = []
for g in G:
    c0, c1, L = G[g]["c0"], G[g]["c1"], G[g]["L"]; bounds = {"5UTR": (0, c0), "CDS": (c0, c1), "3UTR": (c1, L)}
    for donor in "AB":
        cond = "2" if g.startswith("ZIKV") else "3"
        rb, rn = psite(f"{donor}{cond}_ribo", g), psite(f"{donor}{cond}_RNA", g)
        dens = {k: (rb[s:e].sum() / (e - s), rn[s:e].sum() / (e - s)) for k, (s, e) in bounds.items()}
        for k, (s, e) in bounds.items():
            reg.append(dict(genome=g, donor=donor, region=k, length=e - s, ribo_reads=int(rb[s:e].sum()), rna_reads=int(rn[s:e].sum()),
                            ribo_per_nt=dens[k][0], rna_per_nt=dens[k][1], ribo_over_rna=dens[k][0] / dens[k][1] if dens[k][1] else np.nan,
                            rel_to_CDS=(dens[k][0] / dens[k][1]) / (dens["CDS"][0] / dens["CDS"][1]) if dens[k][1] and dens["CDS"][1] else np.nan))
R = pd.DataFrame(reg); R.to_csv(out/"region_density_ratios.tsv", sep="\t", index=False)
print("\nribo/RNA density ratio by region, relative to CDS (=1):"); print(R.pivot_table(index=["genome", "donor"], columns="region", values="rel_to_CDS").round(3).to_string())

# ---------- 3. motif loci and Rfam elements on the two genomes
fa = out/"genomes.fa"; fa.write_text("".join(f">{g}\n{d['seq']}\n" for g, d in G.items()))
F = pd.read_csv(fam/"families_ranked.tsv", sep="\t"); Lc = pd.read_csv(fam/"loci_ranked.tsv", sep="\t")
sites = []
for _, r in F[F.region == "utr3"].iterrows():
    m = fam/"utr3"/f"rep_{int(r.fam_idx)}.meme"
    subprocess.run([f"{RN}/fimo", "--oc", str(out/f"fimo_{int(r.fam_idx)}"), "--thresh", "1e-4", "--norc", "--verbosity", "1", str(m), str(fa)], check=True, capture_output=True)
    try: h = pd.read_csv(out/f"fimo_{int(r.fam_idx)}"/"fimo.tsv", sep="\t", comment="#").dropna(subset=["sequence_name"])
    except pd.errors.EmptyDataError: continue
    for _, s in h.iterrows(): sites.append(dict(genome=s.sequence_name, start=int(s.start), stop=int(s.stop), family=r.rep_consensus, rank=int(r["rank"])))
S = pd.DataFrame(sites)
db = Path("reference/known_elements/flavi_rfam.cm")
subprocess.run([f"{RN}/cmscan", "--cut_ga", "--rfam", "--nohmmonly", "--fmt", "2", "--tblout", str(out/"genomes_cmscan.tbl"), "-o", "/dev/null", str(db), str(fa)], check=True, capture_output=True)
els = []
for l in (out/"genomes_cmscan.tbl").read_text().splitlines():
    if l.startswith("#") or not l.strip(): continue
    p = l.split(); els.append(dict(genome=p[3], element=p[1], start=min(int(p[9]), int(p[10])), stop=max(int(p[9]), int(p[10]))))
E = pd.DataFrame(els); E.to_csv(out/"rfam_elements_on_genomes.tsv", sep="\t", index=False); print("\nRfam elements on genomes:"); print(E.to_string(index=False))
# locus -> density
loc_map = {r["rep_consensus"]: next((x.locus for x in Lc.itertuples() if r["rep_consensus"] in x.families.split(";")), None) for _, r in F[F.region == "utr3"].iterrows()}
S["locus"] = S.family.map(loc_map); S.to_csv(out/"motif_sites_on_genomes.tsv", sep="\t", index=False)
dens_rows = []
for _, s in S.iterrows():
    g = s.genome; cond = "2" if g.startswith("ZIKV") else "3"
    for donor in "AB":
        rb, rn = psite(f"{donor}{cond}_ribo", g), psite(f"{donor}{cond}_RNA", g)
        u = (s.start - 1, s.stop); c1, L = G[g]["c1"], G[g]["L"]
        dens_rows.append(dict(genome=g, donor=donor, locus=s.locus, family=s.family, start=s.start, stop=s.stop, in_region="3UTR" if s.start - 1 >= c1 else ("CDS" if s.start - 1 >= G[g]["c0"] else "5UTR"),
                              ribo_in_site=rb[u[0]:u[1]].sum(), rna_in_site=rn[u[0]:u[1]].sum()))
Dn = pd.DataFrame(dens_rows); Dn.to_csv(out/"motif_site_densities.tsv", sep="\t", index=False)
print("\nmotif sites by genome/region:"); print(S.assign(reg=lambda x: [("3UTR" if r.start-1 >= G[r.genome]["c1"] else "CDS/5UTR") for r in x.itertuples()]).groupby(["genome", "reg", "locus"]).size().to_string())

# ---------- 4. tracks
fig, axes = plt.subplots(len(G) * 2, 1, figsize=(13, 9), sharex=False)
for gi, g in enumerate(G):
    cond = "2" if g.startswith("ZIKV") else "3"; c0, c1, L = G[g]["c0"], G[g]["c1"], G[g]["L"]
    for ki, kind in enumerate(["ribo", "RNA"]):
        ax = axes[gi * 2 + ki]
        for donor, col in (("A", "#1f77b4"), ("B", "#ff7f0e")):
            v = psite(f"{donor}{cond}_{kind}", g); sm = np.convolve(v, np.ones(25) / 25, mode="same"); ax.plot(np.arange(L), sm / max(sm.sum(), 1) * 1e4, lw=0.8, color=col, label=f"donor {donor}")
        ax.axvspan(c0, c1, color="0.93", zorder=0); ax.set_ylabel(f"{kind}\n(per 1e4 viral reads)"); ax.set_xlim(0, L)
        if ki == 0: ax.set_title(f"{g}  (grey = CDS; red ticks = motif-locus sites; blue boxes = Rfam elements)", fontsize=9)
        if kind == "ribo": ax.legend(fontsize=7, loc="upper right")
        for _, s in S[S.genome == g].iterrows(): ax.axvline((s.start + s.stop) / 2, ymin=0.9, ymax=1, color="red", lw=0.8)
        for _, e in E[E.genome == g].iterrows(): ax.axvspan(e.start, e.stop, ymin=0, ymax=0.08, color="tab:blue", alpha=0.7)
plt.tight_layout(); plt.savefig(out/"tracks.png", dpi=140)
json.dump({"psite_offset": a.psite, "read_length_profile_ribo": lenprof}, open(out/"qc.json", "w"))
