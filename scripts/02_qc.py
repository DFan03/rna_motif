#!/usr/bin/env python3
"""QC of GenBank records -> clean RNA FASTA + metadata + qc log. Raw input is never modified."""
import argparse, hashlib, re, statistics as st, sys, platform
from collections import Counter
from pathlib import Path
import pandas as pd, Bio
from Bio import SeqIO

ap = argparse.ArgumentParser()
ap.add_argument("--raw", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--max-amb-frac", type=float, default=0.01, help="drop record if non-ACGU fraction above this")
ap.add_argument("--min-scale", type=float, default=50.0, help="floor (nt) on robust length scale")
ap.add_argument("--len-mad-z", type=float, default=3.5, help="robust z (per species) for length outlier")
a = ap.parse_args()
raw, out = Path(a.raw), Path(a.out); out.mkdir(parents=True, exist_ok=True)

recs, log = [], []
def drop(acc, reason): log.append((acc, reason))

for gb in sorted(raw.glob("*.gb")):
    for r in SeqIO.parse(gb, "genbank"):
        src = next((f for f in r.features if f.type == "source"), None)
        q = src.qualifiers if src else {}
        taxid = next((x.split(":")[1] for x in q.get("db_xref", []) if x.startswith("taxon:")), None)
        seq = str(r.seq).upper().replace("T", "U")
        cds = [f for f in r.features if f.type == "CDS"]
        main = max(cds, key=lambda f: len(f.location)) if cds else None
        recs.append(dict(acc=r.id, species_query=gb.stem.replace("_", " "), organism=r.annotations.get("organism"),
            taxid=taxid, host=(q.get("host") or [None])[0], country=(q.get("country") or [None])[0],
            collection_date=(q.get("collection_date") or [None])[0],
            moltype=r.annotations.get("molecule_type"), length=len(seq), seq=seq,
            n_cds=len(cds), cds_start=int(main.location.start) if main else None,
            cds_end=int(main.location.end) if main else None,
            cds_strand=main.location.strand if main else None,
            has_5utr_feat=any(f.type == "5'UTR" for f in r.features), has_3utr_feat=any(f.type == "3'UTR" for f in r.features),
            definition=r.description))
n0 = len(recs); df = pd.DataFrame(recs)

# 1 molecule type: GenBank stores some RNA-virus genomes as cDNA ("DNA"); these are valid, so flag only.
# Drop only if moltype is missing (cannot tell what the record is).
bad = df[df.moltype.isna()]
for acc in bad.acc: drop(acc, "moltype_missing")
df = df.drop(bad.index)
for acc in df[~df.moltype.str.contains("RNA")].acc: log.append((acc, "FLAG_moltype_DNA_cDNA(kept)"))
# 2 alphabet
df["n_nonacgu"] = df.seq.map(lambda s: sum(c not in "ACGU" for c in s))
df["amb_frac"] = df.n_nonacgu / df.length
bad = df[df.amb_frac > a.max_amb_frac]
for acc in bad.acc: drop(acc, f"non_ACGU_frac>{a.max_amb_frac}")
df = df.drop(bad.index)
# 3 length outliers (robust z within species query)
def mz(g):
    med = g.median(); scale = max((g - med).abs().median() / 0.6745, a.min_scale)
    return (g - med) / scale
df["len_z"] = df.groupby("species_query").length.transform(mz)
bad = df[df.len_z.abs() > a.len_mad_z]
for acc in bad.acc: drop(acc, f"length_outlier_MADz>{a.len_mad_z}")
df = df.drop(bad.index)
# 4 missing annotation (no CDS => cannot define UTRs); keep flagged, not dropped
df["no_cds"] = df.cds_start.isna()
for acc in df[df.no_cds].acc: log.append((acc, "FLAG_no_CDS_annotation(kept)"))
# 4b CDS spanning (almost) the whole record => UTRs unusable; drop (UTR-based analysis needs both flanks)
df["utr5_len"] = df.cds_start; df["utr3_len"] = df.length - df.cds_end
bad = df[(~df.no_cds) & ((df.utr5_len < 10) | (df.utr3_len < 10))]
for acc in bad.acc: drop(acc, "UTR_incomplete(<10nt flank)")
df = df.drop(bad.index)
# 5 exact duplicates (sequence hash); keep first accession, record members
df["sha"] = df.seq.map(lambda s: hashlib.sha256(s.encode()).hexdigest())
dup = df[df.duplicated("sha", keep="first")]
members = df.groupby("sha").acc.apply(list)
for _, r in dup.iterrows(): drop(r.acc, f"exact_duplicate_of:{members[r.sha][0]}")
df = df.drop(dup.index)
# 6 missing metadata flags (kept)
for col in ["host", "collection_date", "country", "taxid"]:
    df[f"missing_{col}"] = df[col].isna()

# regions from the longest CDS (polyprotein): 5'UTR = [0,cds_start), 3'UTR = [cds_end,L)
df = df.reset_index(drop=True)
with open(out/"viral_clean.fasta", "w") as f:
    for _, r in df.iterrows(): f.write(f">{r.acc} {r.organism}\n{r.seq}\n")
df.drop(columns=["seq"]).to_csv(out/"metadata_clean.tsv", sep="\t", index=False)
pd.DataFrame(log, columns=["accession", "reason"]).to_csv(out/"qc_log.tsv", sep="\t", index=False)

c = Counter(r.split(":")[0].split("(")[0] for _, r in log)
rep = [f"# QC report", f"- raw dir: {raw}", f"- records read: {n0}", f"- records kept: {len(df)}",
       f"- removed/flagged by reason: {dict(c)}",
       f"- params: max_amb_frac={a.max_amb_frac}, len_mad_z={a.len_mad_z}, min_scale={a.min_scale}",
       f"- versions: python {platform.python_version()}, biopython {Bio.__version__}, pandas {pd.__version__}",
       "", "## Kept per species", df.species_query.value_counts().to_string(),
       "", "## Length (nt) per species", df.groupby("species_query").length.describe()[["count","min","50%","max"]].to_string(),
       "", "## UTR lengths (nt)", df[["utr5_len","utr3_len"]].describe().to_string(),
       "", "## Missing metadata (kept, flagged)", df[[c for c in df if c.startswith("missing_")]].sum().to_string(),
       "", "## Hosts", df.host.fillna("NA").value_counts().head(8).to_string()]
(out/"qc_report.md").write_text("\n".join(rep) + "\n"); print("\n".join(rep))
