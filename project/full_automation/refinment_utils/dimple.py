import os
import subprocess
import re
from pathlib import Path

def generate_density_image(pdb_path, map_path, output_image):
    script = f"""
load {pdb_path}, model
load {map_path}, map
isomesh mesh, map, level=1.0, carve=1.6, model
orient
ray 800, 600
png {output_image}
quit
"""
    with open("pymol_script.pml", "w") as f:
        f.write(script)

    subprocess.run(["pymol", "-cq", "pymol_script.pml"])
    os.remove("pymol_script.pml")


def dimple_execution(mtz, pdb, highres_cutoff=None):
    """
    Run DIMPLE rough refinement on input MTZ + PDB files.

    :param hkl_input_file: path to .mtz file
    :param highres_cutoff: optional resolution cutoff (Å), e.g. 2.0
    :returns: (rwork, rfree, resolution_high, resolution_low)
    """

    if mtz.endswith('.hkl'):
        # Convert to MTZ if UC file exists
        mtz_file = mtz.replace("hkl", "mtz")
        mtz_command = f'get_hkl -i {mtz} -p {pdb} -o {mtz_file} --output-format=mtz | tee conversion_log.txt'
        os.system(mtz_command)
        mtz = mtz_file

    
    if not os.path.exists(pdb):
        raise FileNotFoundError(f"PDB file not found: {pdb}")

    base = Path(mtz).stem
    out_dir = Path(mtz).parent / f"{base}_dimple_out"
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = ['dimple', '--mapconvert=1', '--no-cleanup', '--loglevel=debug']
    if highres_cutoff:
        cmd.append(f'--highres={highres_cutoff}')
    cmd += [mtz, pdb, str(out_dir)]


    print("Running:", " ".join(cmd))
    proc = subprocess.run(cmd, cwd=out_dir, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"DIMPLE failed:\n{proc.stderr}")

    log = out_dir / 'dimple.log'
    if not log.exists():
        log = out_dir / f"{base}.log"  # fallback

    rwork = rfree = res_hi = res_lo = None
    pattern_rf = re.compile(r"R_work\s*=\s*([\d\.]+)")
    pattern_rfree = re.compile(r"R_free\s*=\s*([\d\.]+)")
    pattern_res = re.compile(r"Resolution used\s*([\d\.]+)\s*-\s*([\d\.]+)\s\AA")

    with open(log) as f:
        for line in f:
            if not rwork:
                m = pattern_rf.search(line)
                if m:
                    rwork = float(m.group(1))
            if not rfree:
                m = pattern_rfree.search(line)
                if m:
                    rfree = float(m.group(1))
            if not (res_hi and res_lo):
                m = pattern_res.search(line)
                if m:
                    res_lo = float(m.group(1))
                    res_hi = float(m.group(2))
    generate_density_image(pdb, out_dir / f"{base}_map_coefficients.mtz", out_dir / f"{base}_density_image.png")
    
    return rwork, rfree, res_hi, res_lo
