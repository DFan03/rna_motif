#!/bin/bash
# Same work as 10b_run_all.sh but 4 runs at a time (user: use all CPU). Skips finished runs (a .bai file marks completion).
cd /Users/a000000/rna_motif
S=/private/tmp/claude-501/-Users-a000000-rna-motif/528baf19-8494-4844-b24e-800973afb3f7/scratchpad/ribo_full
PY=/opt/anaconda3/envs/riboseq/bin/python
export S PY
cat <<'LIST' | xargs -P 4 -L 1 bash -c '
label=$0; run=$1
[ -f data/processed/te_ribo/bams/$label.bam.bai ] && { echo "skip $label"; exit 0; }
$PY scripts/10_ribo_map.py --run $run --label $label --out data/processed/te_ribo/bams --scratch $S/$label --threads 3 || echo "FAILED $label $run"
'
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
