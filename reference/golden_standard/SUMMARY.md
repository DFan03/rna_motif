# Golden standards — summaries  (STATUS: DRAFT, awaiting user approval; do not rely on until approved)

Two sources. Rule IDs (LK-n, MV-n) are cited in PLAN.md / results to justify methodological choices.

---
## S1. Bernett, Blumenthal, Grimm, Haselbeck, Joeres, Kalinina & List. "Guiding questions to avoid data
##     leakage in biological machine learning applications." Nature Methods 21:1444-1453 (Aug 2024).
##     doi:10.1038/s41592-024-02362-y   File: Guiding-questions-...pdf  (10 PDF pages = journal pp. 1444-1453)
Read in full from the local PDF. Page numbers below = journal page numbers.

Definition (p.1445): leakage = performance inflated because the model uses information that generalizes
train->test but NOT train->inference-time (the real use case). Taxonomy (Kapoor & Narayanan): L1 no clean
train/test separation (joint preprocessing/normalization/feature selection/clustering on train+test; duplicates
in both); L2 illegitimate features; L3 test set not from the distribution of interest (L3.1 temporal, L3.2
dependencies absent at inference, L3.3 sampling bias).
Why biology is hard (p.1446): complex dependencies (homology between sequences) and data biased to well-studied molecules.

Seven questions and the rule each implies for THIS project:
- LK-1 Study bias (p.1446). Overstudied elements inflate performance. -> Viral databases are dominated by
  a few well-sequenced viruses (SARS-CoV-2, HIV, influenza, HCV, dengue). Report taxon composition; cap
  per-species contribution; state which viruses the conclusions do/do not cover.
- LK-2 Label distribution (p.1446). Class imbalance / naive resampling creates shortcuts. -> Applies to
  Phase 3 (TE labels): inspect label distribution per split; no majority-class shortcuts; stratify.
- LK-3 Data shift (p.1446, p.1449). Test set must represent inference-time data; keep whole groups
  (e.g. same protein/species) in ONE split. -> Split by cluster/species, never by individual sequence.
- LK-4 Inter-sample similarities (p.1446, 1448-1450). Hidden homology/dependencies let models memorize.
  Use homology-aware splitting tools (GraphPart, DataSAIL named at p.1448-1449, p.1450); also treat
  similarity in embedding/feature space. -> Cluster by sequence identity BEFORE any split; the cluster is
  the unit of splitting and of counting evidence. (Directly supports the "avoid duplicate sampling bias" answer.)
- LK-5 Unavailable features (p.1446). Features absent at inference time are illegitimate. -> Any feature used
  to relate motifs to TE must be computable from sequence alone (e.g. predicted, not experimentally probed, structure) if the
  intended use is prediction for new viruses.
- LK-6 Features as illegitimate surrogates (p.1446, 1448). Engineered dataset-level features (e.g. counts
  per protein/taxon) can proxy the label. -> No per-taxon/per-cluster summary features (e.g. "number of
  isolates of this virus") in TE models; check that motif features are not proxies for virus identity/length/GC.
- LK-7 Pretraining (p.1448). Embeddings from pretrained models (RNA-FM, UTR-LM etc.) may have seen test
  sequences. -> If a pretrained RNA language model is ever used, check overlap of its pretraining data with
  our test sequences. (Currently not planned.)
Other: general leakage items from the taxonomy (p.1445): split BEFORE preprocessing/normalization/feature
selection/clustering of the fit data (LK-8, no joint preprocessing); remove duplicates across splits (LK-9);
keep a real hold-out test set looked at once (LK-10). Box 1 (p.1446): report so that method, result and
inferential reproducibility are possible (LK-11).
Note: the paper is about supervised ML. Phase 2 (unsupervised enrichment) is not "trained", but the same logic
applies when a threshold/PWM is tuned or a motif is "discovered then tested" -> discover on one set of
clusters, test on held-out clusters (our reading, not stated in the paper).

---
## S2. Jin, Zhou, Zhang & Chen. "mRNA vaccine sequence and structure design and optimization: Advances and
##     challenges." J Biol Chem 301(1):108015 (2025; online 26 Nov 2024). doi:10.1016/j.jbc.2024.108015
##     PMC11728972 (open access). No PDF stored locally.
CAUTION on provenance: summarized from the PMC web text via a summarizing fetch of roughly the first 100k of 232k
characters, not from a PDF; no page numbers, so citations are by section. Numbers should be re-checked against the
paper before being quoted in a report. If you drop the PDF here I will re-verify.
This is a design review for therapeutic mRNA, NOT a motif-discovery protocol. Its value for this project is
(a) a catalogue of known regulatory elements/features to compare candidates against, (b) which features relate to
TE and stability, (c) structure-prediction tooling.

- MV-1 (Intro / four principles). Mechanistic frame: UTRs -> ribosome loading & TE; codon usage -> TE;
  global structure -> in-solution stability; local motifs -> in-cell lifetime & translation fidelity.
  Use these four categories to annotate candidate motifs.
- MV-2 (5' UTR). Known inhibitory elements: uORFs, iron-responsive elements; structured 5'UTR impairs loading.
  Human 5'UTRs average ~200 nt (~70% >100 nt). ML predictors (Optimus 5-Prime CNN, FramePool, MTtrans,
  UTR-LM, RNA-FM) predict mean ribosome load from sequence; FramePool reports r=0.11-0.25 against measured
  TE in human data -> sequence-only TE prediction is weak in native contexts; expect small effects in Phase 3.
- MV-3 (3' UTR / poly(A)). AU-rich elements (50-150 nt, multiple AUUUA) promote decay; miRNA response elements
  repress translation; poly(A) tail/PABP protects and circularizes. These are the first "known element"
  references for Phase 4. Review notes no comprehensive 3'UTR database.
- MV-4 (Codon optimization). CAI (geometric mean of relative codon usage; vaccine ORFs >0.9) and CSC
  (codon stabilization coefficient from mRNA half-life) are TE/stability covariates. -> k-mer enrichment
  inside CDS is confounded by codon usage: use codon-preserving background (our B4) and report CAI/GC
  as covariates.
- MV-5 (Global structure). Lower MFE / higher GC / fewer unpaired nts -> more in-solution stable; but
  over-stable ORF structure can stall ribosomes, and MFE correlates weakly with in-cell expression (LinearDesign
  benchmarks). Do not treat MFE as a TE proxy. Best designs balance CAI and MFE (MV-7).
- MV-6 (Local structure / fidelity). Stem-loops/pseudoknots can induce frameshifting; single-stranded loops
  are endonuclease targets; consensus across tools or base-pair probabilities (partition function) recommended
  to call local structure. -> If we add a structure layer, use RNAfold partition-function bpp, not MFE only.
- MV-7 (Multi-objective). No single metric (CAI, MFE, AUP) suffices; the best-performing sequences were not
  the extremes of either (Fig. 3). -> Rank motifs with several shown components, not one score (matches PLAN).
- MV-8 (Tools). Thermodynamic folders (RNAfold, RNAstructure, mfold, Vfold2D; Turner parameters) struggle with
  multibranch loops and pseudoknots; ML structure predictors have not surpassed them on long mRNAs
  (OpenVaccine, CASP15). LinearFold/LinearPartition for long sequences. -> Prefer RNAfold/LinearFold; avoid ML 3D models.
- MV-9 (Modifications / data). m1Psi, Psi alter energetics (~0.25/0.18 kcal/mol per pair); irrelevant for viral
  genomic RNA analysis. PERSIST-seq (233 vaccine mRNAs) is the data resource mentioned for stability; the review
  states TE/fidelity/lifetime databases are lacking (supports Decision D4: TE data must come from elsewhere).
- Limit: vaccine-design context (modified, IVT, human 5'/3'UTR reuse) differs from natural viral genomes
  (cap-independent translation, IRES, frameshift signals, genome-wide structure). Do not transfer numeric
  thresholds directly.

---
## What these standards do NOT settle (still need your decision)
- Identity threshold & clustering tool (LK-4 says homology-aware, gives no number).
- Choice of background model for k-mer enrichment (neither source addresses it).
- Multiple-testing method: only general_principles.md #3 (FDR/Bonferroni) applies.
