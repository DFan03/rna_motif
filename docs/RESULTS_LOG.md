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

## 2026-10-31 — Motif families and ranking `results/families_v1` (script 08)
- 41 replicated motifs -> **16 families** (12 in 3'UTR, 4 in 5'UTR; corrected from an earlier miscount of 13/3). Table: `results/families_v1/families_ranked.tsv`; figure: `family_coverage.png`.
- Top 3'UTR families (representative consensus; groups replicated of 5; clusters with a site of 25; main Rfam overlap):
  1. GUUAGAGGAGAC (4; 20/25; 86% of sites in Flavivirus_DB)  2. GUCAGGCCA (4; 20/25; early 3'UTR, rel. pos 0.24; 19% in Flavi_xrRNA)
  3. GGACUAGAGGUU (4; 20/25; 81% in Flavivirus_DB)  4. GUGCAGYCUGY (4; 20/25)  5. CCAGGCACAGA (3; 23/25; 96% in Flavi_CRE; all 5 groups; end of 3'UTR, rel. pos 0.91).
- Families 1, 3 and 7 all fall in Flavivirus_DB: Tomtom keeps them apart because they cover different parts of the element -> effectively one *locus*, not three.
- Strong pattern: the top families are present in Dengue, Yellow fever, JE-serogroup and Zika but absent from TBEV (tick-borne vs mosquito-borne divergence); only family 5 (CRE) and a few weak ones occur in TBEV.
- Caveats: coverage numbers are descriptive (include discovery sequences); rank order depends on C-21; the 5'UTR families are weaker (best discovery E 0.002-0.08 but low-complexity-looking, e.g. UAACAGUUUUUU); no TE data yet, so nothing here says anything about translation.

## 2026-10-31 — 3'UTR loci by Rfam element `results/families_v1/loci_ranked.tsv` (script 09)
- 12 3'UTR families -> **10 loci**, 3 Rfam-defined. Only one merge happened: **Flavivirus_DB** = 3 families (GUUAGAGGAGAC, GGACUAGAGGUU, SGAAGCUGUANG; 13 motifs;
  replicated in 4 groups; 20/25 clusters). Every DB hit on the cluster representatives contains a site of this locus; the locus covers 39% of DB-hit nucleotides.
- Flavi_CRE locus (CCAGGCACAGA): 3 groups replicated, all 5 groups with a site, 96% of hits contain a site but only 11% of CRE nt are covered (a short motif inside a ~96-nt element).
- Flavi_TBFV_xrRNA locus (GUGGCGACG): weak (discovery E 1.3; 1 group replicated; 8/25 clusters).
- Seven families remain unannotated (e.g. GUCAGGCCA, GUGCAGYCUGY, GCAUAUUGAC), all replicated in 1-4 groups. GUCAGGCCA / GUGCAGYCUGY sit early in the 3'UTR where xrRNA is expected, but the Rfam xrRNA model detects xrRNA only in JE-serogroup, so this is NOT confirmed.
- Correction: the earlier family count split was 12 (3'UTR) / 4 (5'UTR), not 13/3.
- Ranking of loci uses the C-21 criteria applied to loci (union of groups/clusters).
