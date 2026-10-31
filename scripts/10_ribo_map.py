#!/usr/bin/env python3
"""Stream SRA runs -> trim (cutadapt) -> map to the two viral genomes (bowtie2) -> sorted BAM + stats (C-24).
Disk is tight, so each run is processed alone and its .sra deleted afterwards. --max-spots for dry runs (no prefetch)."""
import argparse, json, shutil, subprocess, sys, collections
from pathlib import Path
R = "/opt/anaconda3/envs/riboseq/bin"
ap = argparse.ArgumentParser()
ap.add_argument("--run", required=True); ap.add_argument("--label", required=True)
ap.add_argument("--index", default="data/interim/te_ribo/viruses_idx"); ap.add_argument("--out", required=True)
ap.add_argument("--max-spots", type=int, default=0); ap.add_argument("--scratch", required=True)
ap.add_argument("--min-len", type=int, default=18); ap.add_argument("--threads", type=int, default=4); a = ap.parse_args()
out = Path(a.out); out.mkdir(parents=True, exist_ok=True); scratch = Path(a.scratch); scratch.mkdir(parents=True, exist_ok=True)
sra = None
if a.max_spots:
    src = [f"{R}/fastq-dump", "-X", str(a.max_spots), "-Z", "--split-files", a.run]
else:
    subprocess.run([f"{R}/prefetch", "-O", str(scratch), "--max-size", "5G", a.run], check=True, capture_output=True)
    sra = scratch / a.run / f"{a.run}.sra"; src = [f"{R}/fastq-dump", "-Z", "--split-files", str(sra)]
bam = out / f"{a.label}.bam"; logf = out / f"{a.label}.cutadapt.txt"
# library: SMARTer smRNA (5' template-switch nts + polyA tail + Illumina adapter); 3 nt 5' trim per kit guidance
p1 = subprocess.Popen(src, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
p2 = subprocess.Popen([f"{R}/cutadapt", "-u", "3", "-a", "AGATCGGAAGAGC", "-a", "A{10}", "--times", "2", "-m", str(a.min_len), "-q", "10", "-j", "1", "-"],
                      stdin=p1.stdout, stdout=subprocess.PIPE, stderr=open(logf, "w"))
p3 = subprocess.Popen([f"{R}/bowtie2", "-x", a.index, "-U", "-", "--no-unal", "-p", str(a.threads), "--very-sensitive", "-L", "16", "--no-1mm-upfront"],
                      stdin=p2.stdout, stdout=subprocess.PIPE, stderr=open(out / f"{a.label}.bowtie2.txt", "w"))
p4 = subprocess.Popen([f"{R}/samtools", "sort", "-o", str(bam), "-"], stdin=p3.stdout)
for p in (p1, p2, p3, p4): p.wait()
subprocess.run([f"{R}/samtools", "index", str(bam)], check=True)
st = subprocess.run([f"{R}/samtools", "idxstats", str(bam)], capture_output=True, text=True).stdout
if sra and not a.max_spots: shutil.rmtree(scratch / a.run, ignore_errors=True)
ca = Path(logf).read_text(); tot = [l for l in ca.splitlines() if "Total reads processed" in l]
bt = Path(out / f"{a.label}.bowtie2.txt").read_text().splitlines()
print(a.label, a.run, "|", tot[0].strip() if tot else "", "|", bt[-1] if bt else "", "|", st.replace("\n", " ; "))
