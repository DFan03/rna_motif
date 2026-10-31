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

## 2026-10-31 — Leave-one-serogroup-out (LOGO) run `results/logo_v2` (scripts 06, 07)
Setup: C-15..C-19. 5 folds x {utr5, utr3}; 246 motifs reported by STREME (45 with E<0.05); each tested on the held-out group.
- **Bug found and fixed before trusting any result:** first run (`results/logo_v1`, invalid, kept but git-ignored) joined FIMO output to
  STREME motifs on the wrong ID column (`motif_id` vs `motif_alt_id`), so every held-out count was 0 ("0 replicated"). Caught because
  all counts, including the null mean, were exactly 0; confirmed by scanning motifs on their own training sequences (many hits).
  v2 uses the correct column. (No automated regression test yet for this join.)
- **Replication:** 41 of 246 tested motifs replicate (BH FDR<0.05, hits in >=2 held-out clusters): 35 in 3'UTR, 6 in 5'UTR.
  Per fold (replicated/tested): 3'UTR Dengue 9/33, JE 9/34, TBEV 2/39, YFV 6/42, Zika 9/42; 5'UTR Dengue 1/9, JE 4/12, TBEV 0/11, YFV 0/11, Zika 1/13.
- Discovery strength predicts replication: 14/45 (31%) of E<0.05 motifs replicate vs 27/201 (13%) of the rest.
- Many replicated motifs are overlapping variants of the same locus (e.g. GGACUAGAGGUU / GACUAGAGGU / UAGAGGAGAC), so 41 is not 41 independent elements; they still need merging into families.
- **Known elements (Rfam/Infernal):** Flavi_CRE hits in all 5 groups; Flavivirus_DB in Dengue 64/67, JE 70/83, Zika 17/18, not YFV/TBEV; Flavi_xrRNA only in JE-serogroup (31/83);
  TBFV-specific models only in TBEV (sanity check passes). Dengue xrRNA not detected by the general CM at GA cutoff (CM sensitivity limit).
  Replicated 3'UTR motif sites: 120/221 (54.3%) overlap a known element vs 41.1% of 3'UTR nt covered by Rfam hits (descriptive; sites overlap each other, no significance test);
  recall: 19.4% of Rfam-covered nt are under a replicated motif site.
- Interpretation: the pipeline finds sequence that is conserved across serogroups and is modestly enriched at known cis-elements. This is evidence of conservation,
  not of regulatory function or of effect on TE.
