#!/usr/bin/env python3
"""Fetch complete Orthoflavivirus genomes (GenBank flat files) from NCBI, N per species.
Raw output is write-once: refuses to overwrite an existing outdir. Writes MANIFEST.tsv."""
import argparse, hashlib, sys, time, datetime
from pathlib import Path
from Bio import Entrez

SPECIES = ["Dengue virus 1", "Dengue virus 2", "Dengue virus 3", "Dengue virus 4",
           "Zika virus", "West Nile virus", "Japanese encephalitis virus",
           "Yellow fever virus", "Tick-borne encephalitis virus", "Usutu virus",
           "Saint Louis encephalitis virus", "Murray Valley encephalitis virus"]

ap = argparse.ArgumentParser()
ap.add_argument("--outdir", required=True)
ap.add_argument("--per-species", type=int, default=20)
ap.add_argument("--n-species", type=int, default=len(SPECIES))
ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()
out = Path(a.outdir)
if out.exists(): sys.exit(f"{out} exists; raw data is read-only, choose a new dir")
out.mkdir(parents=True)
Entrez.tool = "rna_motif_project"
rows = []
for sp in SPECIES[:a.n_species]:
    term = (f'"{sp}"[Organism] AND viruses[filter] '
            f'AND 9000:12500[SLEN] AND "complete genome"[Title]')
    h = Entrez.esearch(db="nuccore", term=term, retmax=10000, idtype="acc"); r = Entrez.read(h); h.close()
    ids = r["IdList"]
    # deterministic subsample (sorted accessions, evenly spaced) -> reproducible, spans submission time
    ids = sorted(ids)
    step = max(1, len(ids)//a.per_species) if ids else 1
    pick = ids[::step][:a.per_species]
    print(f"{sp}: {len(ids)} hits, fetching {len(pick)}", flush=True)
    if not pick: continue
    h = Entrez.efetch(db="nuccore", id=",".join(pick), rettype="gb", retmode="text")
    txt = h.read(); h.close()
    f = out / (sp.replace(" ", "_") + ".gb"); f.write_text(txt)
    rows.append((sp, term, len(ids), len(pick), f.name, hashlib.sha256(txt.encode()).hexdigest()))
    time.sleep(0.4)
with open(out/"MANIFEST.tsv", "w") as m:
    m.write("# source=NCBI nuccore (GenBank) via Entrez; downloaded=%s; per_species=%d\n" %
            (datetime.date.today(), a.per_species))
    try:
        import urllib.request; rel = urllib.request.urlopen("https://ftp.ncbi.nlm.nih.gov/genbank/GB_Release_Number", timeout=20).read().decode().strip()
    except Exception: rel = "unavailable"
    m.write("# GenBank release number at download: %s\n" % rel)
    m.write("species\tquery\tn_hits\tn_fetched\tfile\tsha256\n")
    for r in rows: m.write("\t".join(map(str, r)) + "\n")
print("done ->", out)
