# rna_motif — viral RNA regulatory motif discovery (small-scale test project)

Goal: discover candidate sequence/structure motifs in viral mRNA/RNA, rank them with statistical
evidence, and (later phase) relate them to translation efficiency (TE). Scale: laptop CPU (M2 Pro).

- Plan / workflow: `PLAN.md`  (start here)
- Methodological references: `reference/golden_standard/`, `~/.claude/standards/`
- Layout: `data/{raw,interim,processed}` (raw is read-only), `src/rnamotif` (library),
  `scripts/` (numbered pipeline steps), `results/`, `docs/`, `tests/`, `logs/`
