#!/usr/bin/env python3
"""Group-level train/held-out split, fixed BEFORE any motif discovery (leakage standard LK-3, LK-4, LK-8).
Unit of splitting = flavivirus serogroup (cluster ids at 90% identity never cross species, see docs)."""
import argparse, numpy as np, pandas as pd
GROUP = {"Dengue virus 1": "Dengue", "Dengue virus 2": "Dengue", "Dengue virus 3": "Dengue", "Dengue virus 4": "Dengue",
         "Zika virus": "Zika",
         "West Nile virus": "JE-serogroup", "Japanese encephalitis virus": "JE-serogroup", "Usutu virus": "JE-serogroup",
         "Saint Louis encephalitis virus": "JE-serogroup", "Murray Valley encephalitis virus": "JE-serogroup",
         "Yellow fever virus": "Yellow fever", "Tick-borne encephalitis virus": "TBEV"}
ap = argparse.ArgumentParser(); ap.add_argument("--dir", required=True)
ap.add_argument("--n-heldout", type=int, default=2); ap.add_argument("--seed", type=int, default=42); a = ap.parse_args()
m = pd.read_csv(f"{a.dir}/metadata_clean.tsv", sep="\t")
cl = pd.read_csv(f"{a.dir}/clusters/clusters_id90.tsv", sep="\t").rename(columns={"rep": "cluster90", "weight": "weight90"})
m = m.merge(cl[["acc", "cluster90", "weight90"]], on="acc"); m["group"] = m.species_query.map(GROUP)
assert m.group.notna().all()
groups = sorted(m.group.unique()); rng = np.random.RandomState(a.seed)
held = list(rng.choice(groups, a.n_heldout, replace=False))
m["split"] = np.where(m.group.isin(held), "heldout", "train")
assert m.groupby("cluster90").split.nunique().max() == 1, "a cluster spans both splits"
m[["acc", "species_query", "group", "cluster90", "weight90", "split"]].to_csv(f"{a.dir}/split.tsv", sep="\t", index=False)
print("seed", a.seed, "heldout groups:", held)
print(m.groupby(["split", "group"]).agg(n_seq=("acc", "size"), n_clusters=("cluster90", "nunique")))
