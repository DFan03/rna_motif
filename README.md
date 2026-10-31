# rna_motif — candidate regulatory motifs in viral RNA

**Status: work in progress (Phases 1-2 first pass done).** This README is rewritten at the end for review.

## What this project does
Viral RNA genomes contain local sequence/structure patterns that may regulate stability, translation and protein
binding, but public databases hold almost no functional labels for them. This project builds a small, laptop-scale
pipeline that (1) collects and cleans public viral RNA sequences, (2) discovers *candidate* motifs, ranks them with
statistics against a sensible background, and (3) later asks whether they relate to translation efficiency (TE).
A discovered motif is a statistically over-represented pattern, **not** a validated regulatory element.

## Where things stand
| Phase | What | State |
|---|---|---|
| 0 | Project set-up, standards, plan | done |
| 1 | Data download, QC, redundancy clustering, train/held-out split | **done** (207 genomes, 12 species) |
| 2 | Motif discovery (STREME), leave-one-serogroup-out validation (FIMO + shuffle null), Rfam positive control | **first pass done** (41 replicated motifs, mostly 3'UTR) |
| 3 | Link motifs to TE: descriptive Ribo-seq/RNA-seq density on ZIKV + DENV-2 (GSE207347) | **done (descriptive only)** - no statistical link possible (n=2 donors/virus) |
| 4 | Reliability assessment, final report | not started |

## Data (Phase 1)
- Source: NCBI GenBank, genus *Orthoflavivirus*, complete genomes, <=20 per species (12 species); downloaded 2026-10-31.
- 238 downloaded -> 207 after QC (removed: 13 high ambiguity, 16 length outliers, 1 UTR-incomplete, 1 exact duplicate).
- Highly redundant: the 207 genomes form only 14 / 25 / 47 clusters at 80 / 90 / 95 % identity. Evidence is therefore
  counted per cluster, never per sequence.
- Split (fixed before any motif work, seed 42): **train** = Dengue, TBEV, Yellow fever; **held-out** = JE-serogroup, Zika.
- Details: `docs/qc_report_orthoflavi_v1.md`, `docs/qc_log_orthoflavi_v1.tsv`, `docs/MANIFEST_orthoflavi_v1.tsv`.

## First results (Phase 2, details in `docs/RESULTS_LOG.md`)
- Discovery on 4 serogroups, test on the 5th (all 5 rotations); 246 candidate motifs, **41 replicate** (FDR<0.05, found in >=2 held-out clusters), 35 of them in the 3'UTR.
- Replicated motifs are partly redundant (variants of the same locus); they sit on known Rfam flavivirus elements 54% of the time vs 41% expected by chance (descriptive).
- **Families:** the 41 motifs merge (Tomtom) into **16 families** (12 in 3'UTR, 4 in 5'UTR). Top families lie on the Rfam Flavivirus_DB, xrRNA and CRE elements and are conserved in all mosquito-borne groups but not TBEV. Ranked table: `results/families_v1/families_ranked.tsv`; coverage figure: `results/families_v1/family_coverage.png`.
- Sanity checks: planted motif recovered; random RNA gives 0 motifs at E<0.05; a bug in the validation join was caught and fixed (see RESULTS_LOG).
- This shows conservation across serogroups, **not** regulatory function or effect on TE.

## Phase 3 result (descriptive; `docs/RESULTS_LOG.md`, `docs/TE_DATA_SURVEY.md`)
- Public infection Ribo-seq gives one TE-like value per virus, so a motif-to-TE model is not possible. Instead I mapped Ribo-seq and RNA-seq from ZIKV/DENV-2-infected neural progenitor cells (GSE207347; 12 libraries) onto the viral genomes.
- Libraries pass QC (3-nt periodicity 73-77%, background <=200 reads). The 3'UTR has only 6-22% of the CDS ribosome/RNA density, so the top motif loci (Flavivirus_DB, CRE, ...) sit in a region that is mostly **not translated**; with only ~350 3'UTR footprints per virus nothing finer can be concluded.

## Why each method was chosen
See **`docs/DECISIONS.md`** (every choice, its rationale, alternatives, status). Principle: use established tools
(MMseqs2, MEME Suite, ViennaRNA, Biopython) rather than home-made code. Methodological references are in
`reference/golden_standard/` (data-leakage guidelines, Bernett et al. 2024; mRNA design review, Jin et al. 2025).

## Reproduce
```bash
conda env create -f environment.yml -n rnamotif && conda activate rnamotif
python scripts/01_fetch.py --outdir data/raw/orthoflavi_v1 --per-species 20
python scripts/02_qc.py    --raw data/raw/orthoflavi_v1 --out data/processed/orthoflavi_v1
python scripts/03_cluster.py --fasta data/processed/orthoflavi_v1/viral_clean.fasta --out data/processed/orthoflavi_v1/clusters
python scripts/04_split.py --dir data/processed/orthoflavi_v1
python scripts/05_regions.py --dir data/processed/orthoflavi_v1
python scripts/06_logo.py --dir data/processed/orthoflavi_v1 --out results/logo_v2
python scripts/07_known_elements.py --dir data/processed/orthoflavi_v1 --logo results/logo_v2
python scripts/08_families.py --dir data/processed/orthoflavi_v1 --logo results/logo_v2 --out results/families_v1
pytest tests   # planted-motif / null sanity tests
# Phase 3 (needs conda env `riboseq`: sra-tools cutadapt bowtie2 samtools biopython):
scripts/10c_run_parallel.sh   # 12 SRA runs of GSE207347 -> viral BAMs (about 6 GB download, streamed)
python scripts/11_te_descriptive.py --bams data/processed/te_ribo/bams --ref data/raw/te_ribo/ref --fam results/families_v1 --out results/te_descriptive_v1
```
Data files are not stored in git (see `.gitignore`); they are regenerated by the commands above.

## Known limitations (so far)
- One virus genus; conclusions do not extend to other viruses.
- Related species still share ancestry; the serogroup split reduces but cannot eliminate this.
- No functional/TE labels in the data; any claim of function needs external evidence.

## Layout
`PLAN.md` (workflow) · `docs/DECISIONS.md` (rationale log) · `scripts/` (numbered steps) · `reference/golden_standard/` ·
`data/` (untracked) · `results/` · `environment.yml`
