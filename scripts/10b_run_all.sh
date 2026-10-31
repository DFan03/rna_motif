#!/bin/bash
# GSE207347: label -> SRR (A/B = donor clone; 1 uninfected, 2 ZIKV, 3 DENV-2)  [GEO sample table, see docs/TE_DATA_SURVEY.md]
cd /Users/a000000/rna_motif
S=/private/tmp/claude-501/-Users-a000000-rna-motif/528baf19-8494-4844-b24e-800973afb3f7/scratchpad/ribo_full
PY=/opt/anaconda3/envs/riboseq/bin/python
while read label run; do
  [ -f data/processed/te_ribo/bams/$label.bam.bai ] && { echo "skip $label"; continue; }
  $PY scripts/10_ribo_map.py --run $run --label $label --out data/processed/te_ribo/bams --scratch $S || { echo "FAILED $label $run"; exit 1; }
done <<'LIST'
A1_ribo SRR19918553
B1_ribo SRR19918545
A2_ribo SRR19918537
B2_ribo SRR19918527
A3_ribo SRR19918487
B3_ribo SRR19918519
A1_RNA SRR19918511
B1_RNA SRR19918503
A2_RNA SRR19918493
B2_RNA SRR19918481
A3_RNA SRR19918473
B3_RNA SRR19918465
LIST
echo ALL_DONE
