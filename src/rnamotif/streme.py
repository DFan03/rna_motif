"""Run STREME (MEME Suite, off-the-shelf) and parse its motif list. Significance filter is applied here on STREME's
reported E-value, because STREME's own --thresh does not suppress chance motifs (see docs/DECISIONS.md C-09b)."""
import re, shutil, subprocess
from pathlib import Path

ENV_BIN = "/opt/anaconda3/envs/rnamotif/bin"

def streme_bin():
    return shutil.which("streme") or f"{ENV_BIN}/streme"

def run_streme(primary_fa, outdir, control_fa=None, minw=6, maxw=12, seed=1, order=2):
    cmd = [streme_bin(), "--rna", "--p", str(primary_fa), "--oc", str(outdir), "--minw", str(minw),
           "--maxw", str(maxw), "--seed", str(seed), "--order", str(order)]
    if control_fa: cmd += ["--n", str(control_fa)]
    subprocess.run(cmd, check=True, capture_output=True)
    return parse_streme(Path(outdir) / "streme.txt")

def parse_streme(txt_path):
    """-> list of dict(id, consensus, width, nsites, evalue) in STREME's order."""
    out, cur = [], None
    for line in Path(txt_path).read_text().splitlines():
        m = re.match(r"MOTIF (\S+) (\S+)", line)
        if m: cur = dict(id=m.group(2), consensus=re.sub(r"^\d+-", "", m.group(1))); continue
        m = re.match(r"letter-probability matrix: alength= 4 w= (\d+) nsites= (\d+) E= (\S+)", line)
        if m and cur:
            cur.update(width=int(m.group(1)), nsites=int(m.group(2)), evalue=float(m.group(3))); out.append(cur); cur = None
    return out

def significant(motifs, max_e=0.05):
    return [m for m in motifs if m["evalue"] < max_e]
