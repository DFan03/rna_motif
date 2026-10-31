# PLAN — Viral RNA motif discovery (small-scale test)

Status: DRAFT v0.1 (2026-10-31) — needs user sign-off on the "Decisions to confirm" section.
Constraint: everything must run on a local laptop CPU (M2 Pro) in minutes, not hours.

## 0. Scope and honest framing

Goal: (1) get a small, clean set of viral RNA sequences, (2) discover candidate motifs, rank them with
statistics against a sensible background, (3) later, test whether motifs relate to TE.

- Viral sequences from public databases have NO functional or TE labels. A discovered motif is a
  *statistically over-represented pattern*, not a validated element. Over-representation can come from
  codon usage, dinucleotide bias, shared ancestry (homology), or conserved protein-coding sequence, as
  well as real regulation. The report must say so explicitly (task 2 requirement).
- Validated functional elements (e.g. IRES, 5'/3' UTR structures, frameshift signals, packaging signals,
  poly(A)/AU-rich elements) have experimental evidence (mutagenesis, reporter assay, SHAPE/structure
  probing, conservation + covariation). Our candidates only get that label if they match such evidence.
  We will compare top candidates to known elements (Rfam, literature) as a plausibility check only.
- TE readout: TE (ribosome-profiling RPF / mRNA) is not in sequence databases. Phase 3 needs an external
  dataset (see Decision D4). Phases 1-2 do not depend on it.

## 1. Phases and deliverables

### Phase 0 — Setup (done in this session)
Git repo, folder layout, .gitignore, requirements, this plan.

### Phase 1 — Data acquisition and QC  (task 1)
1. Source: NCBI Virus / GenBank-RefSeq via Entrez (Biopython) OR the offline FASTA + metadata if provided
   (Decision D1). Target: ~100-500 records, one clear virus group, so k-mer stats are meaningful.
2. Raw download goes in `data/raw/` (never edited). Write `data/raw/MANIFEST.tsv`: accession, source URL/query,
   download date, database + release (e.g. RefSeq release number), sha256 of file. Small, committed.
3. Parse FASTA + metadata (organism, taxid, host, molecule type, collection date, length, annotations/CDS/UTR
   features if GenBank format is fetched).
4. QC filters, each logged with counts removed and reason (`results/tables/qc_log.tsv`):
   - molecule type: keep RNA-genome / mRNA-derived records; note ssRNA(+), ssRNA(-), dsRNA, DNA-virus transcripts.
   - alphabet: convert T→U consistently; non-ACGU bases (N, IUPAC): drop records above a threshold
     (e.g. >1% ambiguous) and mask the rest; thresholds recorded, not hidden.
   - length outliers: flag by robust z-score / IQR within group; partial/incomplete genomes handled by metadata flags.
   - exact duplicates: collapse by sequence hash; keep one representative, record cluster members.
   - missing metadata: keep but flag; exclude from any metadata-stratified analysis.
5. Redundancy / sampling bias (required explanation): many records are near-identical isolates of the same
   virus (outbreak sequencing), which would inflate k-mer counts and make "recurrent" motifs just
   sequence copies. Mitigation:
   - cluster at nucleotide identity (e.g. 90% / 95%, MMseqs2 or CD-HIT-EST if installed, otherwise a
     k-mer/Mash-style Jaccard clustering in pure Python for this small set), keep one representative per cluster;
   - and/or weight sequences by 1/cluster size;
   - count motif support as number of *independent clusters / taxa* containing it, not number of sequences;
   - split for validation at the cluster/species level so homologs never sit on both sides (leakage rule).
6. Output: `data/processed/viral_clean.fasta`, `data/processed/metadata_clean.tsv`, QC report in `docs/qc_report.md`.
- Dry run: first run on ~20 records; check the QC log and counts by eye before pulling the full set.

### Phase 2 — Motif discovery  (task 2)
Candidate method (simple, interpretable, CPU-cheap): **k-mer enrichment** (k = 4-8, optionally with
degenerate/IUPAC extension), then **PWM refinement** by merging similar enriched k-mers.
1. Count k-mers per sequence (presence/absence and counts), weighted by cluster weights.
2. Background models (compare several; report sensitivity of the ranking to the choice):
   - B1: shuffled sequences preserving mononucleotide composition;
   - B2: shuffled preserving dinucleotide composition (Altschul-Erikson style) — default, controls for dinucleotide bias;
   - B3: order-k Markov model (k=2-3) fit on the data;
   - B4 (optional): codon-preserving shuffle for coding regions, to separate motifs from codon usage.
3. Statistics per k-mer/motif: observed vs expected count, fold enrichment (log2), empirical p-value
   from N shuffles (N ~ 1000), and/or binomial/Poisson test; Benjamini-Hochberg FDR across all tested motifs;
   effect size + bootstrap CI over clusters (not over sequences).
4. Coverage metrics: number of sequences/clusters/taxa containing the motif, fraction of total sequence length,
   positional distribution (5'UTR / CDS / 3'UTR, relative position), strand/orientation.
5. Merge overlapping significant k-mers into PWMs (information content reported); re-scan sequences with the
   PWM to get final counts (scan threshold chosen on training half only).
6. Ranking: composite of FDR, log2 enrichment, cluster-level support, and robustness across backgrounds;
   components are shown separately (no opaque single score).
7. Optional secondary layer (later): local RNA structure (ViennaRNA RNAfold, MFE / base-pair probability,
   stem-loop context) as an extra feature; skipped if ViennaRNA install is a hurdle.
8. Output: ranked table `results/tables/motifs_ranked.tsv`, logos/positional plots in `results/figures/`.
- Validation within the data: discover on half of the clusters, test enrichment on the held-out half
  (split by cluster/species, set before any fitting). Report how many motifs replicate.
- Negative control: run the identical pipeline on dinucleotide-shuffled sequences; the number of "significant"
  motifs must be ~0 (calibration of false positives).
- Dry run: tiny subset (k=6, 50 shuffles) end-to-end before the full run.

### Phase 3 — Linking motifs to function / TE  (later; after Phase 2 review)
- Needs a TE dataset (D4). Options: public Ribo-seq of virus-infected cells (TE per viral/host transcript),
  or MPRA-style reporter data. Model: motif presence/count + covariates (length, GC, UTR length, uORFs)
  → TE; start with regression with cross-validation split by cluster/gene family; report effect sizes with CIs.
- Compare top motifs to Rfam / literature on known viral elements; classify each as
  "matches known element" vs "novel/unsupported".
- Honest limits: correlation with TE is not causation; wet-lab (mutagenesis/reporter) would be needed.

### Phase 4 — Reliability evaluation and reporting
- Reliability checklist per motif: FDR, replication in held-out clusters, robustness to background choice,
  robustness to clustering threshold, not explained by codon/dinucleotide bias, positional/structural
  plausibility, overlap with known elements.
- Run the `stat-check` skill to audit the finished run against the golden standards.
- Write `docs/report.md`: methods, parameters, seeds, software versions, results, limitations.

## 2. Repository layout

```
rna_motif/
  PLAN.md  README.md  requirements.txt  .gitignore
  data/raw/        read-only downloads + MANIFEST.tsv
  data/interim/    intermediate files (not committed)
  data/processed/  cleaned FASTA/metadata (not committed; regenerable)
  reference/golden_standard/   project methodological references (papers, SUMMARY.md, checklist.md)
  src/rnamotif/    library code (io, qc, cluster, kmer, background, stats, pwm)
  scripts/         numbered steps: 01_fetch.py 02_qc.py 03_cluster.py 04_kmer_enrich.py 05_pwm.py 06_report.py
  tests/           small unit tests (toy sequences with a planted motif)
  results/{tables,figures}/    outputs, each with params + seed + versions
  docs/            qc_report.md, report.md
  logs/            run logs (not committed)
```

## 3. Cross-cutting rules (from user conventions)

- Raw data read-only; every processing step scripted; fixed random seeds; parameters and software versions saved with every result.
- Dry run → report "dry run passed" + what to inspect → ask before any full run.
- Unit test with a planted motif in synthetic sequences must recover it before touching real data.
- At most 2-3 debug attempts, then stop and summarize.
- Commit small and often; no data/secrets in git.
- Compute budget: all steps designed for a few thousand sequences max, minutes on M2 Pro; no GPU.

## 4. Decisions to confirm

Golden standards now exist in `reference/golden_standard/` (SUMMARY.md = DRAFT pending your approval; rule IDs LK-n = leakage
paper, MV-n = JBC mRNA review). They cover splitting/leakage and known-element context but are silent on the items below,
so I will NOT improvise on these; please decide or point me to a reference:

- D1. Data source: NCBI Virus/RefSeq download via Entrez, or do you have the offline FASTA + metadata from the task? Which virus group (e.g. coronaviruses, flaviviruses, picornaviruses, influenza)? Suggest one family with RNA genomes and annotated UTRs (e.g. Flaviviridae or Coronaviridae).
- D2. Redundancy handling: identity threshold (suggest 90%, sensitivity at 80/95%) and clustering tool (MMseqs2/CD-HIT install vs pure-Python).
- D3. Background model: default dinucleotide-preserving shuffle (B2) with B1/B3 as sensitivity checks — OK?
- D4. TE data: do you already have a TE dataset (Ribo-seq/MPRA) in mind, or should Phase 3 start with a literature/data search?
- D5. Region scope: whole genome, or split into 5'UTR / CDS / 3'UTR using annotations (needs GenBank features)?
- D6. If you have a review paper on RNA motif discovery / viral RNA data processing, put it in `reference/golden_standard/` so I can summarize it and follow it. (Neither current standard covers motif-enrichment statistics.)

## 5. Immediate next steps (after your sign-off)
1. Answer D1-D6 (even partially).
2. Write `01_fetch.py` + `02_qc.py`, dry run on ~20 records, show you the QC log.
3. Planted-motif unit test, then Phase 2 dry run.
