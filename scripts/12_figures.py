#!/usr/bin/env python3
"""All README figures, regenerated from the result tables (docs/figures/*.png). Palette: validated default (dataviz skill),
categorical slots 1-3 (blue/orange/aqua) + neutral gray; single-hue sequential for heatmaps."""
import re, sys, collections, warnings
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib.transforms import Affine2D
from Bio import SeqIO
warnings.filterwarnings("ignore")
ROOT = Path("."); D = ROOT/"data/processed/orthoflavi_v1"; OUT = ROOT/"docs/figures"; OUT.mkdir(parents=True, exist_ok=True)
BLUE, ORANGE, AQUA, GRAY, VIOLET, RED = "#2a78d6", "#eb6834", "#1baf7a", "#9a9993", "#6250d6", "#e34948"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
SEQ = LinearSegmentedColormap.from_list("seq", ["#eaf2fc", "#9cc2ee", "#2a78d6", "#143f80"])
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF,
    "savefig.facecolor": SURF, "axes.titleweight": "bold", "axes.titlesize": 10, "axes.titlelocation": "left"})
def save(fig, name): fig.savefig(OUT/f"{name}.png", dpi=150, bbox_inches="tight"); plt.close(fig); print("wrote", name)
def grid(ax, axis="x"): ax.grid(axis=axis, color=GRID, lw=0.8); ax.set_axisbelow(True)
meta = pd.read_csv(D/"metadata_clean.tsv", sep="\t"); split = pd.read_csv(D/"split.tsv", sep="\t")
qlog = pd.read_csv(ROOT/"docs/qc_log_orthoflavi_v1.tsv", sep="\t")
val = pd.read_csv(ROOT/"results/logo_v2/logo_validation.tsv", sep="\t")
fam = pd.read_csv(ROOT/"results/families_v1/families_ranked.tsv", sep="\t"); loci = pd.read_csv(ROOT/"results/families_v1/loci_ranked.tsv", sep="\t")
short = lambda s: s.replace(" virus", "").replace("encephalitis", "enceph.").replace("Japanese", "Jap.").replace("Tick-borne", "Tick-b.").replace("Murray Valley", "Murray V.").replace("Saint Louis", "St. Louis")

# ---------------- F1 workflow diagram
fig, ax = plt.subplots(figsize=(11, 5.2)); ax.axis("off"); ax.set_xlim(0, 11); ax.set_ylim(0, 5.2)
def box(x, y, w, h, title, sub, col):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08", fc="white", ec=col, lw=1.6))
    ax.text(x + w/2, y + h - 0.2, title, ha="center", va="top", fontsize=8.6, fontweight="bold", color=INK)
    ax.text(x + w/2, y + h - 0.52, sub, ha="center", va="top", fontsize=7.4, color=INK2, linespacing=1.25)
def arrow(x0, y0, x1, y1): ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=9, color=INK2, lw=1.1))
ax.text(0.05, 5.05, "Phase 1 - data and QC (task 1)", fontsize=9, fontweight="bold", color=BLUE)
B = [(0.05, "1 Fetch", "NCBI GenBank\n(Biopython Entrez)\n12 species x <=20"), (2.2, "2 Parse + QC", "organism, host, CDS\nambiguity, length,\nduplicates, flags"),
     (4.35, "3 Cluster", "MMseqs2\n80 / 90 / 95 % id\n207 -> 25 clusters"), (6.5, "4 Split", "5 serogroups\nseed 42\nclusters never cross"), (8.65, "5 Regions", "5'UTR / CDS / 3'UTR\nfrom longest CDS")]
for x, t, s in B: box(x, 3.55, 2.05, 1.3, t, s, BLUE)
for i in range(4): arrow(B[i][0] + 2.05, 4.2, B[i+1][0], 4.2)
ax.text(0.05, 3.3, "Phase 2 - motif discovery, validation, ranking (task 2)", fontsize=9, fontweight="bold", color=ORANGE)
C = [(0.05, "6 Discover", "STREME on 4 serogroups\n(leave-one-group-out)\n246 motifs"), (2.2, "7 Test held-out", "FIMO on 5th group;\nnull = 200 dinuc. shuffles\nBH-FDR -> 41 replicated"),
     (4.35, "8 Families", "Tomtom all-vs-all\n41 -> 16 families"), (6.5, "9 Loci", "merge by Rfam element\n(Infernal, 12 -> 10 loci)"), (8.65, "10 Rank", "groups replicated,\nclusters, discovery E")]
for x, t, s in C: box(x, 1.75, 2.05, 1.3, t, s, ORANGE)
for i in range(4): arrow(C[i][0] + 2.05, 2.4, C[i+1][0], 2.4)
ax.plot([9.7, 9.7, 1.05], [3.55, 3.13, 3.13], color=INK2, lw=1.1); arrow(1.05, 3.13, 1.05, 3.06)
ax.text(0.05, 1.5, "Phase 3 - translation context (descriptive)   |   Phase 4 - reliability checks", fontsize=9, fontweight="bold", color=AQUA)
Dd = [(0.05, "11 Ribo-seq", "GSE207347: 12 libraries\ncutadapt, bowtie2 ->\nZIKV, DENV-2 BAMs"), (2.2, "12 Density", "ribosome / RNA per\nregion; motif sites\nmarked on genomes"),
      (4.35, "Controls", "planted motif, null RNA,\nRfam positive control,\nperiodicity"), (6.5, "Leakage rules", "group split, cluster\ncounting, no joint\npreprocessing"), (8.65, "Reports", "README, DECISIONS,\nRESULTS_LOG")]
for x, t, s in Dd: box(x, 0.1, 2.05, 1.3, t, s, AQUA)
save(fig, "f01_workflow")

# ---------------- F2 QC funnel
reasons = qlog.reason.str.replace(r"[:(].*", "", regex=True).value_counts()
n_raw, n_keep = 238, len(meta)
items = [("Downloaded", n_raw, INK2), ("Non-ACGU bases > 1 %", int(reasons.get("non_ACGU_frac>0.01", 13)), ORANGE), ("Length outlier (robust z)", int(reasons.get("length_outlier_MADz>3.5", 16)), ORANGE),
         ("No UTR flanks", int(reasons.get("UTR_incomplete", 1)), ORANGE), ("Exact duplicate", int(reasons.get("exact_duplicate_of", 1)), ORANGE), ("Kept", n_keep, BLUE),
         ("  of which flagged: no CDS annotation", int(reasons.get("FLAG_no_CDS_annotation", 12)), GRAY), ("  of which flagged: GenBank 'DNA' (cDNA) type", int(reasons.get("FLAG_moltype_DNA_cDNA", 10)), GRAY)]
fig, ax = plt.subplots(figsize=(8, 3.6)); ys = np.arange(len(items))[::-1]
ax.barh(ys, [i[1] for i in items], color=[i[2] for i in items], height=0.62)
for y, (n, v, c) in zip(ys, items): ax.text(v + 3, y, str(v), va="center", fontsize=8.5, color=INK)
ax.set_yticks(ys); ax.set_yticklabels([i[0] for i in items]); ax.set_xlabel("records"); grid(ax); ax.set_title("QC: 238 downloaded -> 207 kept (orange removed, gray flagged but kept)")
save(fig, "f02_qc_funnel")

# ---------------- F3 lengths
raw = {}
for f in sorted(Path("data/raw/orthoflavi_v1").glob("*.gb")):
    for r in SeqIO.parse(f, "genbank"): raw[r.id] = (f.stem.replace("_", " "), len(r.seq))
rawdf = pd.DataFrame([(k, *v) for k, v in raw.items()], columns=["acc", "sp", "len"]); rawdf["status"] = np.where(rawdf.acc.isin(meta.acc), "kept", "removed")
order = rawdf.groupby("sp").len.median().sort_values().index.tolist()
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1.25, 1]})
rng = np.random.RandomState(0)
for i, sp in enumerate(order):
    for st, col, off in (("kept", BLUE, 0), ("removed", ORANGE, 0)):
        v = rawdf[(rawdf.sp == sp) & (rawdf.status == st)].len.values; a1.scatter(v, i + rng.uniform(-0.28, 0.28, len(v)), s=14, color=col, alpha=0.85, lw=0, label=st if i == 0 else None)
a1.set_yticks(range(len(order))); a1.set_yticklabels([short(s) for s in order]); a1.set_xlabel("genome length (nt)"); grid(a1); a1.legend(frameon=False, fontsize=8, loc="lower right"); a1.set_title("Genome length per species (all 238 downloaded)")
ut = meta.dropna(subset=["utr3_len"]); o2 = ut.groupby("species_query").utr3_len.median().sort_values().index.tolist()
for i, sp in enumerate(o2):
    v = ut[ut.species_query == sp].utr3_len.values; a2.scatter(v, i + rng.uniform(-0.28, 0.28, len(v)), s=14, color=BLUE, alpha=0.85, lw=0)
a2.set_yticks(range(len(o2))); a2.set_yticklabels([short(s) for s in o2]); a2.set_xlabel("3'UTR length (nt)"); grid(a2); a2.set_title("3'UTR length per species (kept records)")
fig.subplots_adjust(wspace=0.5); save(fig, "f03_lengths")

# ---------------- F4 redundancy
cl = {i: pd.read_csv(D/f"clusters/clusters_id{i}.tsv", sep="\t") for i in (80, 90, 95)}
sp = meta[["acc", "species_query"]]
per = sp.groupby("species_query").size().rename("sequences").to_frame().join(cl[90].merge(sp, on="acc").groupby("species_query").rep.nunique().rename("clusters90")).sort_values("sequences")
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1.5, 1]})
y = np.arange(len(per)); a1.barh(y + 0.19, per.sequences, 0.36, color=GRAY, label="sequences kept"); a1.barh(y - 0.19, per.clusters90, 0.36, color=BLUE, label="clusters at 90 % identity")
for yy, s, c in zip(y, per.sequences, per.clusters90): a1.text(s + 0.3, yy + 0.19, str(s), va="center", fontsize=7.5, color=INK2); a1.text(c + 0.3, yy - 0.19, str(c), va="center", fontsize=7.5, color=INK2)
a1.set_yticks(y); a1.set_yticklabels([short(s) for s in per.index]); a1.set_xlabel("count"); grid(a1); a1.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2); a1.set_title("Redundancy per species: many genomes, few clusters")
n = [207] + [cl[i].rep.nunique() for i in (80, 90, 95)]; labs = ["sequences", "clusters\n>=80 %", "clusters\n>=90 %", "clusters\n>=95 %"]
a2.bar(range(4), n, color=[GRAY, BLUE, BLUE, BLUE], width=0.6); [a2.text(i, v + 3, str(v), ha="center", fontsize=9) for i, v in enumerate(n)]
a2.set_xticks(range(4)); a2.set_xticklabels(labs); a2.set_ylabel("count"); grid(a2, "y"); a2.set_title("Effective sample size depends on the threshold")
save(fig, "f04_redundancy")

# ---------------- F5 split design
groups = ["Dengue", "Zika", "JE-serogroup", "Yellow fever", "TBEV"]; gi = split.groupby("group").agg(seq=("acc", "size"), clu=("cluster90", "nunique"))
rows = [f"fold {i+1}: hold out {g}" for i, g in enumerate(groups)] + ["fixed split (seed 42)"]
fig, ax = plt.subplots(figsize=(9, 3.6))
for r in range(6):
    for c, g in enumerate(groups):
        held = (g == groups[r]) if r < 5 else g in ("JE-serogroup", "Zika")
        ax.add_patch(Rectangle((c, 5 - r), 0.94, 0.86, fc=ORANGE if held else BLUE, ec="none", alpha=0.95))
        ax.text(c + 0.47, 5 - r + 0.43, f"{'TEST' if held else 'train'}\n{gi.clu[g]} cl.", ha="center", va="center", color="white", fontsize=7.8, fontweight="bold")
ax.set_xlim(0, 5); ax.set_ylim(0, 6); ax.set_xticks(np.arange(5) + 0.47); ax.set_xticklabels([f"{g}\n({gi.seq[g]} seq)" for g in groups]); ax.xaxis.tick_top()
ax.set_yticks(np.arange(6) + 0.43); ax.set_yticklabels(rows[::-1]); [s.set_visible(False) for s in ax.spines.values()]; ax.tick_params(length=0)
ax.set_title("Leave-one-serogroup-out design: whole groups (never single sequences) are held out", pad=34)
save(fig, "f05_split_design")

# ---------------- F6 calibration (planted vs null)
sys.path.insert(0, "src"); sys.path.insert(0, "tests")
from rnamotif.streme import run_streme
import test_planted_motif as T
tmp = Path("data/interim/fig_calib"); tmp.mkdir(parents=True, exist_ok=True)
T.make(str(tmp/"null.fa"), 150, 0.0, 2); T.make(str(tmp/"planted.fa"), 150, 0.7, 1)
res = {k: run_streme(str(tmp/f"{k}.fa"), str(tmp/f"o_{k}")) for k in ("null", "planted")}
from matplotlib.patches import Patch
fig, ax = plt.subplots(figsize=(8, 3.6)); x = 0
for k, col in (("planted", BLUE), ("null", GRAY)):
    for m in res[k]:
        ax.bar(x, m["evalue"], color=col, width=0.7); ax.text(x, m["evalue"] * 1.2, m["consensus"], rotation=60, fontsize=7, ha="left", va="bottom", color=INK2); x += 1
    x += 1
ax.set_yscale("log"); ax.set_ylim(0.005, 400); ax.axhline(0.05, color=INK, lw=1, ls="--"); ax.text(x - 0.4, 0.056, "E = 0.05 cut-off", ha="right", fontsize=8); ax.set_xticks([]); ax.set_ylabel("STREME E-value (log)")
ax.legend(handles=[Patch(color=BLUE, label="synthetic RNA, 70 % of sequences carry GCUAGUCAGUA"), Patch(color=GRAY, label="synthetic RNA, nothing planted")], frameon=False, fontsize=8, loc="upper left"); grid(ax, "y")
ax.set_title("Calibration on synthetic RNA: STREME lists chance motifs too, so I keep only E < 0.05")
save(fig, "f06_calibration")

# ---------------- F7 discovery -> replication
val["frac"] = val.clusters_with_hit / val.n_heldout_clusters; val["nullfrac"] = val.null_mean / val.n_heldout_clusters
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1.3, 1]})
nr = val[~val.replicated]; r = val[val.replicated]
a1.scatter(nr.discovery_E, nr.frac + rng.uniform(-0.015, 0.015, len(nr)), s=16, color=GRAY, alpha=0.7, lw=0, label=f"not replicated (n={len(nr)})")
a1.scatter(r.discovery_E, r.frac + rng.uniform(-0.015, 0.015, len(r)), s=22, color=BLUE, lw=0, label=f"replicated, FDR<0.05 (n={len(r)})")
a1.set_xscale("log"); a1.axvline(0.05, color=INK2, lw=0.8, ls="--"); a1.text(0.052, 1.03, "discovery E = 0.05", fontsize=7.5, color=INK2); a1.set_xlabel("STREME discovery E-value (log; training groups)"); a1.set_ylabel("fraction of held-out clusters with a site")
a1.legend(frameon=False, fontsize=8, loc="center right"); grid(a1, "y"); a1.set_title("Stronger discovery evidence -> more likely to replicate")
t = val.groupby(["region", "heldout"]).agg(tested=("motif", "size"), rep=("replicated", "sum")).reset_index(); t["lab"] = t.region.map({"utr3": "3'", "utr5": "5'"}) + "UTR " + t.heldout.str.replace("-serogroup", "").str.replace("Yellow fever", "YF")
t = t.sort_values(["region", "heldout"], ascending=[False, True]); yy = np.arange(len(t))[::-1]
a2.barh(yy, t.tested, color=GRAY, height=0.6, label="tested"); a2.barh(yy, t["rep"], color=BLUE, height=0.6, label="replicated"); [a2.text(tt + 0.5, y_, f"{int(rr)}/{tt}", va="center", fontsize=7.5) for y_, tt, rr in zip(yy, t.tested, t["rep"])]
a2.set_yticks(yy); a2.set_yticklabels(t.lab); a2.set_xlabel("motifs"); a2.legend(frameon=False, fontsize=8); grid(a2); a2.set_title("Replicated / tested per held-out group")
save(fig, "f07_replication")

# ---------------- F8 observed vs null (top 3'UTR families)
top = fam[fam.region == "utr3"].head(4); fig, axs = plt.subplots(1, 4, figsize=(12, 3.0), sharey=False)
for ax, (_, f) in zip(axs, top.iterrows()):
    fold, mid = f.rep_uid.split("__"); fold_n = fold.replace("_", " ")
    row = val[(val.region == "utr3") & (val.heldout == fold_n) & (val.motif_id == mid)].iloc[0]
    try: h = pd.read_csv(f"results/logo_v2/work/utr3_{fold_n}/fimo_shuf/fimo.tsv", sep="\t", comment="#").dropna(subset=["sequence_name"])
    except Exception: h = pd.DataFrame(columns=["motif_alt_id", "sequence_name"])
    h = h[h.motif_alt_id == mid]; h["copy"] = h.sequence_name.str.extract(r"_shuf_(\d+)$")[0]
    nul = h.groupby("copy").sequence_name.nunique().reindex([str(i) for i in range(1, 201)], fill_value=0)
    bins = np.arange(0, row.n_heldout_clusters + 2) - 0.5; ax.hist(nul, bins=bins, color=GRAY, rwidth=0.8)
    ax.axvline(row.clusters_with_hit, color=BLUE, lw=2); ax.text(row.clusters_with_hit - 0.1, ax.get_ylim()[1] * 0.92, f"observed\n{int(row.clusters_with_hit)}", color=BLUE, ha="right", fontsize=8, fontweight="bold")
    ax.set_title(f"{f.rep_consensus}\nheld out: {fold_n}", fontsize=8.5); ax.set_xlabel(f"clusters with a site (of {int(row.n_heldout_clusters)})"); ax.set_ylabel("shuffles (of 200)" if ax is axs[0] else ""); grid(ax, "y")
fig.suptitle("Observed cluster count in the held-out group vs 200 dinucleotide-shuffled copies (empirical p = 1/201 when no shuffle reaches it)", x=0.01, y=1.1, ha="left", fontsize=9.5, fontweight="bold")
save(fig, "f08_observed_vs_null")

# ---------------- F9 logos
def blocks(path):
    L = Path(path).read_text().splitlines(); res = {}
    for i, l in enumerate(L):
        m = re.match(r"MOTIF (\S+) (\S+)", l)
        if m: w = int(re.search(r"w= (\d+)", L[i+1]).group(1)); res[m.group(2)] = np.array([r.split() for r in L[i+2:i+2+w]], float)
    return res
LC = {"A": AQUA, "C": BLUE, "G": ORANGE, "U": VIOLET}
def logo(ax, P):
    for j, p in enumerate(P):
        ic = 2 + np.sum(p[p > 0] * np.log2(p[p > 0])); y = 0
        for k in np.argsort(p):
            h = p[k] * ic
            if h < 0.02: continue
            tp = TextPath((0, 0), "ACGU"[k], size=1, prop=matplotlib.font_manager.FontProperties(family="DejaVu Sans", weight="bold")); bb = tp.get_extents()
            tr = Affine2D().translate(-bb.x0, -bb.y0).scale(0.92 / bb.width, h / bb.height).translate(j + 0.04, y)
            ax.add_patch(PathPatch(tr.transform_path(tp), fc=LC["ACGU"[k]], ec="none")); y += h
    ax.set_xlim(0, len(P)); ax.set_ylim(0, 2); ax.set_xticks([]); ax.set_yticks([0, 1, 2]); ax.set_ylabel("bits")
sel = fam.head(8); fig, axs = plt.subplots(len(sel), 1, figsize=(6.5, 1.25 * len(sel)))
for ax, (_, f) in zip(axs, sel.iterrows()):
    fold, mid = f.rep_uid.split("__"); P = blocks(f"results/logo_v2/work/{f.region}_{fold.replace('_', ' ')}/streme/streme.txt")[mid]; logo(ax, P)
    ax.set_title(f"#{int(f['rank'])} {f.region} {f.rep_consensus} - replicated in {int(f.n_groups_replicated)} held-out groups", fontsize=8, pad=2)
fig.tight_layout(); save(fig, "f09_logos")

# ---------------- F10 coverage heatmap
cc = [c for c in fam if c.startswith("cov_")]; ordc = ["cov_Dengue", "cov_Zika", "cov_JE-serogroup", "cov_Yellow fever", "cov_TBEV"]; cc = [c for c in ordc if c in cc]
fig, ax = plt.subplots(figsize=(7.2, 0.34 * len(fam) + 1.6)); M = fam[cc].fillna(0).values; im = ax.imshow(M, cmap=SEQ, vmin=0, vmax=1, aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]): ax.text(j, i, f"{M[i, j]:.0%}", ha="center", va="center", fontsize=7.2, color="white" if M[i, j] > 0.55 else INK)
ax.set_xticks(range(len(cc))); ax.set_xticklabels([c[4:] for c in cc]); ax.xaxis.tick_top(); ax.set_yticks(range(len(fam))); ax.set_yticklabels([f"{int(r['rank'])}. {r.region} {r.rep_consensus}" for _, r in fam.iterrows()], fontsize=8)
[s.set_visible(False) for s in ax.spines.values()]; ax.tick_params(length=0); ax.set_title("Share of clusters per serogroup with a site of each family (descriptive; rank = replication order)", pad=26)
fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="share of clusters"); save(fig, "f10_family_coverage")

# ---------------- F11 loci ranking
L3 = loci.copy(); L3["lab"] = [(l + "  (" + f.split(";")[0] + ")") if not l.startswith("unannotated") else (l.replace("unannotated_", "") + "  (unannotated)") for l, f in zip(L3.locus, L3.families)]; L3 = L3.iloc[::-1]
fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True); yy = np.arange(len(L3))
col = [BLUE if not l.startswith("unannotated") else GRAY for l in L3.locus]
a1.barh(yy, L3.n_groups_replicated, color=col, height=0.6); a1.set_xlim(0, 5); a1.set_xlabel("held-out groups where it replicated (of 5)")
a2.barh(yy, L3.clusters_with_hit, color=col, height=0.6); a2.set_xlim(0, 25); a2.set_xlabel("clusters with a site (of 25; descriptive)")
a3.scatter(L3.best_discovery_E, yy, color=col, s=40, lw=0); a3.set_xscale("log"); a3.axvline(0.05, color=INK2, ls="--", lw=0.8); a3.set_xlabel("best discovery E-value (log; dashed = 0.05)"); 
a1.set_yticks(yy); a1.set_yticklabels(L3.lab, fontsize=8); [grid(a) for a in (a1, a2, a3)]; a1.set_title("Ranked 3'UTR loci (blue = defined by an Rfam element, gray = unannotated)", fontsize=9.5)
save(fig, "f11_loci_ranking")

# ---------------- F12 3'UTR map with Rfam elements and motif sites
S = pd.read_csv(ROOT/"results/te_descriptive_v1/motif_sites_on_genomes.tsv", sep="\t"); E = pd.read_csv(ROOT/"results/te_descriptive_v1/rfam_elements_on_genomes.tsv", sep="\t")
c1 = {}
for g in ("ZIKV_IbH30656", "DENV2_16681"):
    r = SeqIO.read(f"data/raw/te_ribo/ref/{g}.gb", "genbank"); c = max([f for f in r.features if f.type == "CDS"], key=lambda f: len(f.location)); c1[g] = (int(c.location.end), len(r.seq))
loc_col = {"Flavivirus_DB": BLUE, "Flavi_CRE": ORANGE, "Flavi_TBFV_xrRNA": AQUA}
fig, axs = plt.subplots(2, 1, figsize=(11, 4.2))
for ax, g in zip(axs, c1):
    s0, L = c1[g]; ax.set_xlim(0, L - s0); ax.set_ylim(-1.4, 0.3 + 0.45 * S[(S.genome == g) & (S.start - 1 >= s0)].family.nunique() + 0.3)
    for _, e in E[E.genome == g].iterrows():
        ax.add_patch(Rectangle((e.start - s0, -1.2), e.stop - e.start, 0.9, fc=loc_col.get(e.element, GRAY), alpha=0.25, ec=loc_col.get(e.element, GRAY))); ax.text((e.start + e.stop) / 2 - s0, -0.75, e.element, ha="center", va="center", fontsize=7.5)
    s3 = S[(S.genome == g) & (S.start - 1 >= s0)]
    for k, (fam_, d) in enumerate(s3.groupby("family")):
        y = 0.3 + k * 0.45; lc = d.locus.iloc[0]; ax.scatter((d.start + d.stop) / 2 - s0, [y] * len(d), s=26, marker="s", color=loc_col.get(lc, GRAY), lw=0)
        ax.text(L - s0 + 4, y, fam_, va="center", fontsize=7, color=INK2)
    ax.set_yticks([]); ax.set_xlabel("position in 3'UTR (nt)" if g.startswith("DENV") else ""); ax.set_title(g.replace("_", " ") + " 3'UTR: Rfam elements (boxes) and motif-family sites (squares)", fontsize=9); [ax.spines[s].set_visible(False) for s in ("left",)]
fig.legend(handles=[Patch(color=BLUE, label="Flavivirus_DB locus"), Patch(color=ORANGE, label="Flavi_CRE locus"), Patch(color=AQUA, label="Flavi_TBFV_xrRNA locus"), Patch(color=GRAY, label="unannotated families")], frameon=False, fontsize=8, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.04))
fig.tight_layout(rect=(0, 0.04, 1, 1)); save(fig, "f12_utr_map")

# ---------------- F13 Rfam control
pres = pd.read_csv(ROOT/"results/logo_v2/rfam_presence_by_group.tsv", sep="\t", index_col=0); ov = pd.read_csv(ROOT/"results/logo_v2/known_element_overlap.tsv", sep="\t").iloc[0]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 3.6), gridspec_kw={"width_ratios": [1.6, 1]})
nn = pres.index.map(lambda s: int(re.search(r"\((\d+)\)", s).group(1))); P = pres.values / np.array(nn)[:, None]
im = a1.imshow(P, cmap=SEQ, vmin=0, vmax=1, aspect="auto")
for i in range(P.shape[0]):
    for j in range(P.shape[1]): a1.text(j, i, f"{P[i, j]:.0%}", ha="center", va="center", fontsize=8, color="white" if P[i, j] > 0.55 else INK)
a1.set_xticks(range(P.shape[1])); a1.set_xticklabels(pres.columns, rotation=20, ha="right", fontsize=8); a1.set_yticks(range(P.shape[0])); a1.set_yticklabels([s.split(" (")[0] for s in pres.index]); [s.set_visible(False) for s in a1.spines.values()]; a1.tick_params(length=0)
a1.set_title("Rfam models: share of 3'UTRs with a hit, by serogroup")
a2.bar([0, 1, 2], [ov.frac_overlap, ov.chance_frac, ov.recall_nt], color=[BLUE, GRAY, ORANGE], width=0.6); [a2.text(i, v + 0.015, f"{v:.0%}", ha="center") for i, v in enumerate([ov.frac_overlap, ov.chance_frac, ov.recall_nt])]
a2.set_xticks([0, 1, 2]); a2.set_xticklabels(["motif sites inside\na known element", "chance level\n(Rfam coverage)", "Rfam nt covered\nby a motif (recall)"], fontsize=8); a2.set_ylim(0, 0.8); a2.set_yticks([]); a2.spines["left"].set_visible(False); a2.set_title("Replicated 3'UTR motifs vs known elements")
save(fig, "f13_rfam_control")

# ---------------- F14 translation context
Q = pd.read_csv(ROOT/"results/te_descriptive_v1/viral_read_counts.tsv", sep="\t"); Pm = pd.read_csv(ROOT/"results/te_descriptive_v1/periodicity.tsv", sep="\t"); Rg = pd.read_csv(ROOT/"results/te_descriptive_v1/region_density_ratios.tsv", sep="\t")
fig, axs = plt.subplots(1, 3, figsize=(12.5, 3.8), gridspec_kw={"width_ratios": [1.3, 0.8, 1]})
a = axs[0]; Q["lib"] = Q.donor + "-" + Q.cond.str.replace("uninfected", "uninf.").str.replace("DENV2", "DENV-2") + " " + Q.kind
tot = Q.groupby(["lib", "kind", "cond"]).plus.sum().reset_index(); tot["o"] = tot.cond.map({"uninfected": 0, "ZIKV": 1, "DENV2": 2}); tot = tot.sort_values(["kind", "o", "lib"], ascending=[False, True, True]); yy = np.arange(len(tot))[::-1]
a.barh(yy, tot.plus.clip(lower=1), color=[GRAY if c == "uninfected" else (BLUE if k == "ribo" else ORANGE) for c, k in zip(tot.cond, tot.kind)], height=0.65); a.set_xscale("log")
[a.text(max(v, 1) * 1.15, y_, f"{int(v):,}", va="center", fontsize=7.5) for y_, v in zip(yy, tot.plus)]; a.set_yticks(yy); a.set_yticklabels(tot.lib, fontsize=7.5); a.set_xlabel("reads mapped to the two viral genomes (log)"); grid(a); a.set_title("Viral reads per library\n(gray = uninfected = noise floor)")
b = axs[1]; Pm["f1"] = Pm.f1 / (Pm.f0 + Pm.f1 + Pm.f2); Pm["f2"] = Pm.f2 / (Pm.f0 + Pm.f1 + 0 + Pm.f2); Pm["f0"] = 1 - Pm.f1 - Pm.f2; xb = np.arange(len(Pm))
b.bar(xb - 0.27, Pm.f0, 0.27, color=BLUE, label="frame 0"); b.bar(xb, Pm.f1, 0.27, color=GRAY, label="frame 1"); b.bar(xb + 0.27, Pm.f2, 0.27, color=GRAY, alpha=0.55, label="frame 2")
b.set_xticks(xb); b.set_xticklabels(Pm.label.str.replace("_ribo", "").map({"A2": "A\nZIKV", "A3": "A\nDENV-2", "B2": "B\nZIKV", "B3": "B\nDENV-2"}), fontsize=8); b.set_ylabel("share of P-sites in CDS"); b.set_ylim(0, 1); b.legend(frameon=False, fontsize=7.5); grid(b, "y"); b.set_title("3-nt periodicity\n(consistent with Ribo-seq)")
c = axs[2]; Rr = Rg[Rg.region != "CDS"].copy(); Rr["lab"] = Rr.genome.str.split("_").str[0] + " " + Rr.donor; xs = {"5UTR": 0, "3UTR": 1}; w = 0.2
for i, (lab, d) in enumerate(Rr.groupby("lab")):
    c.bar([xs[r] + (i - 1.5) * w for r in d.region], d.rel_to_CDS, w * 0.9, color=[BLUE, AQUA, ORANGE, VIOLET][i], label=lab)
c.axhline(1, color=INK, lw=1, ls="--"); c.text(1.45, 1.03, "CDS = 1", ha="right", fontsize=7.5); c.set_xticks([0, 1]); c.set_xticklabels(["5'UTR", "3'UTR"]); c.set_ylabel("ribosome/RNA density, relative to CDS"); c.legend(frameon=False, fontsize=7.5); grid(c, "y"); c.set_title("The 3'UTR is mostly not translated")
fig.tight_layout(); save(fig, "f14_translation_context")

# ---------------- F15 motif funnel
n_loci = len(loci); n3 = int((fam.region == "utr3").sum()); n5 = int((fam.region == "utr5").sum())
steps = [("STREME motifs reported (5 folds x 2 regions)", len(val), GRAY), ("with discovery E < 0.05", int((val.discovery_E < 0.05).sum()), GRAY), ("replicated in a held-out serogroup (FDR < 0.05)", int(val.replicated.sum()), BLUE),
         (f"motif families (Tomtom): {n3} in 3'UTR + {n5} in 5'UTR", len(fam), BLUE), (f"units to rank: {n_loci} 3'UTR loci (Rfam merge) + {n5} 5'UTR families", n_loci + n5, BLUE)]
fig, ax = plt.subplots(figsize=(8, 3.3)); yy = np.arange(len(steps))[::-1]
ax.barh(yy, [t[1] for t in steps], color=[t[2] for t in steps], height=0.62); [ax.text(v + 3, y_, str(v), va="center") for y_, (_, v, _) in zip(yy, steps)]
ax.set_yticks(yy); ax.set_yticklabels([t[0] for t in steps], fontsize=8.5); grid(ax); ax.set_xlabel("number of motifs / families / loci"); ax.set_title("How many motifs? 246 candidates -> 41 replicated -> 16 families")
save(fig, "f15_motif_funnel")
