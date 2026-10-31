"""Sanity (GP-8, C-13): STREME must recover a planted motif from synthetic RNA and report ~nothing on pure background."""
import random, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rnamotif.streme import run_streme, significant
MOTIF = "GCUAGUCAGUA"

def make(path, n, plant_frac, seed, L=400):
    rnd = random.Random(seed)
    with open(path, "w") as f:
        for i in range(n):
            s = [rnd.choice("ACGU") for _ in range(L)]
            if rnd.random() < plant_frac:
                m = list(MOTIF)
                for j in range(len(m)):
                    if rnd.random() < 0.1: m[j] = rnd.choice("ACGU")  # 10% per-base mutation
                p = rnd.randrange(0, L - len(m)); s[p:p+len(m)] = m
            f.write(f">s{i}\n{''.join(s)}\n")

def run(fa, out):
    return [m["consensus"] for m in significant(run_streme(fa, out))]

def test_recovers_planted(tmp_path):
    fa = tmp_path / "planted.fa"; make(fa, 150, 0.7, 1)
    motifs = run(str(fa), str(tmp_path / "o1"))
    assert motifs, "no motif found"
    # best motif's consensus should overlap the planted one
    best = motifs[0]
    assert any(MOTIF[i:i+6] in best for i in range(len(MOTIF)-5)), (motifs, MOTIF)

def test_null_gives_nothing(tmp_path):
    fa = tmp_path / "null.fa"; make(fa, 150, 0.0, 2)
    motifs = run(str(fa), str(tmp_path / "o2"))
    assert len(motifs) == 0, motifs   # E<0.05 threshold must give no false positives on random RNA
