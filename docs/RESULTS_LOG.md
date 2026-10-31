# Results log (chronological; negative results included)

## 2026-10-31 — Phase 2 dry run: STREME on 3'UTR, training groups (Dengue, TBEV, YFV)
- Input: `utr3_train_reps90.fa` (13 cluster representatives, one per 90% cluster) vs `utr3_train_all.fa` (103 sequences).
- STREME 5.5.9 `--rna --minw 6 --maxw 12 --order 2` (control = 3-mer-preserving shuffle), filter E < 0.05 (C-09b).
- reps90: 0 motifs reported. all: 11 motifs reported, **0 with E < 0.05** (best: GGGAGGC, E = 0.11; nsites 82/103 seqs).
- Interpretation: with only 13 independent sequences of ~450 nt (~6 kb), power is too low to call anything; the 103-sequence
  run is pseudo-replicated (sequences share ancestry) yet still shows nothing significant, so there is no sign of an inflated
  false positive here. Nothing was tested on the held-out groups (not yet justified: no discovered motif to test).
- Sanity tests (planted motif recovered; null gives 0 motifs at E<0.05) passed — see tests/.
- Status: dry run passed technically; next-step options are in the chat log / PLAN revision (leave-one-group-out discovery).
