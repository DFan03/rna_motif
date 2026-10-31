# TE data survey for Phase 3 (decision D4)  — 2026-11-01

Question: can any public dataset give a translation-efficiency (TE) readout that lets us relate the discovered
flavivirus UTR motif families/loci to TE? Search was by web search + reading the primary papers; nothing downloaded yet.

## Candidate datasets (flavivirus infection Ribo-seq)
| # | Accession | Study | Design | TE readout | Notes |
|---|-----------|-------|--------|-----------|-------|
| 1 | GEO GSE207347 | Singh et al., Viruses 2022 ([PMC9316442](https://pmc.ncbi.nlm.nih.gov/articles/PMC9316442)) | Human neural progenitor cells (2 donor clones); ZIKV IbH 30656 and DENV-2 16681, MOI 1, 72 hpi, mock controls; paired Ribo-seq + RNA-seq (HiSeq 2000) | Host genes: RiboDiff TE (19,821 genes). Viruses: descriptive ribosome/RNA ratio only, no formal per-ORF TE or test | Replicate number inconsistent in the paper (n=2 vs "three biological replicates") — check GEO sample list. DENV 3'UTR: high RNA, low ribosome coverage; 5'UTR high coverage; candidate uORFs near capsid start |
| 2 | GEO GSE133111 | DENV-2 NGC in A549, RPLP1/2 knockdown ([bioRxiv 713016](https://www.biorxiv.org/content/10.1101/713016v1.full)) | A549, MOI 10, 2.5 hpi; Ribo-seq + RNA-seq in triplicate; non-targeting vs RPLP1/RPLP2 siRNA; no cycloheximide pre-treatment | Ribosome density (RPF/RNA), CDS-focused; pause site at nt 1025-1035 (start of E) | UTRs not analysed separately |
| 3 | ArrayExpress E-MTAB-5418 (accession reported by a search summary, NOT verified by me) | ZIKV Ribo-seq, mammalian + insect cells ([bioRxiv 112904](https://www.biorxiv.org/content/10.1101/112904v1.full); follow-up [Nat Commun 2024](https://www.nature.com/articles/s41467-024-53085-9)) | ZIKV in Vero and C6/36; Asian/American vs African lineage strains in the 2024 paper | Ribosome density along genome, uORFs | Verify accession before use |

Host-gene TE resources (not viral): [RPFdb v3.0](https://academic.oup.com/nar/article/53/D1/D293/7774864) (5018 Ribo-seq datasets, 34 species, gene-level TE) and
[TranslatomeDB](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5753366/). Neither covers viral genomes.
Not found: any massively parallel reporter assay (MPRA) of viral UTRs; the only MPRA hit was human disease-UTR polysome profiling ([eLife reviewed preprint 98814](https://elifesciences.org/reviewed-preprints/98814)).
One search (for flavivirus UTR mutant/reporter libraries) was blocked by the search tool, so that gap is **not** a negative result.
Only part of papers 1 and 2 was read by me (7.6k and 37.9k characters at the end of each page not read).

## Feasibility (the important part)
- Infection Ribo-seq gives **one TE-like number per virus per condition** (the viral genome is one polyprotein ORF): at most n = 2 viruses
  (ZIKV, DENV-2) in dataset 1 and one in dataset 2. A regression "motif content -> TE" needs many sequences with varying motif content; n<=3
  cannot support it, and any virus-to-virus difference is confounded by strain, cell type, time post infection, and host response.
- Leakage standard LK-5/LK-6 would also bite: virus identity itself is a perfect surrogate for the label.
- Therefore a **statistical** motif->TE model is not supportable with these data. What they can support:
  (a) descriptive: where do replicated motif loci sit relative to ribosome density along the 5'UTR (uORFs/IRES region) and 3'UTR (sfRNA region: expected ~no ribosomes);
  (b) checking that top loci (Flavivirus_DB, Flavi_CRE, xrRNA candidates) lie in the non-translated 3'UTR, so any TE effect would be indirect (e.g. via sfRNA, replication, circularisation), not a direct translation signal;
  (c) host-side: infection-induced TE changes of host genes (dataset 1, RiboDiff already computed in the paper) — but our motifs are viral and rare in host mRNA, so linking them is a new hypothesis, not recommended without your decision.
- For a real motif->TE model the data would need many UTR sequence variants with a TE/ribosome-load readout (MPRA / polysome profiling of a viral-UTR library, or many viruses with Ribo-seq). None was found.

## Options for you (D4)
1. **Descriptive Phase 3 (recommended):** download datasets 1 and 2 (processed/raw after checking sizes), map ribosome density onto the viral genomes, and annotate each motif locus as translated (5'UTR/CDS) or non-translated (3'UTR/sfRNA). No regression; n stated plainly.
2. **Targeted further search** for a viral-UTR library/reporter dataset (different wording / specific labs), then decide.
3. **Redefine Phase 3** with host-cell data (human UTR MPRA TE) to learn generic TE motifs and check whether any resemble our viral motifs. Different question, needs your approval.
Nothing has been downloaded; GEO raw FASTQs are typically several GB per series, so I will ask before any download.
