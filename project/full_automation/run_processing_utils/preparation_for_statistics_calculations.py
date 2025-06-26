import os
import re
import glob
import shlex
import subprocess
from partialator_utils.partialator_execution import run_partialator
from unit_cell_utils.parsing_UC_files import parse_UC_file
from partialator_utils.resolution_cutoff_determination import get_d_at_snr_one, get_d_at_cc_threshold
from refinment_utils.parsing_phenix_pdb_file import parsing_phenix_pdb_file

def extract_first_match(pattern, filepath):
    try:
        cmd = f"grep -e {pattern} {filepath}"
        result = subprocess.check_output(shlex.split(cmd)).decode().strip().split('\n')[0]
        return result
    except subprocess.CalledProcessError:
        return None


def get_UC(hkl_input_file):
    result = extract_first_match("indexamajig", hkl_input_file)
    if result:
        cell_or_pdb = re.findall(r"\b\S+\.(?:cell|pdb)", result)
        if cell_or_pdb:
            return os.path.join("/", cell_or_pdb[0])
    return None


def get_pg(hkl_input_file):
    result = extract_first_match("Symmetry", hkl_input_file)
    return result.split(":")[-1].strip() if result else None


def find_latest_file(pattern):
    files = glob.glob(pattern, recursive=True)
    return max(files, key=os.path.getctime) if files else None


def get_resolution_cutoff_from_logs(hkl_name, Rfree_Rwork_path, run_name):
    pdb_path = f'{Rfree_Rwork_path}/*{run_name.strip()}*/**/*.pdb'
    log_path = f'{Rfree_Rwork_path}/*{run_name.strip()}*/**/*.log'

    phenix_file = find_latest_file(pdb_path) or find_latest_file(log_path)
    if phenix_file:
        return (*parsing_phenix_pdb_file(phenix_file), phenix_file)
    return (None, None, 1.5, None, None)


def fallback_resolution_cutoff(hkl_name):
    SNR_file = hkl_name.replace('.hkl', '_SNR.dat')
    CC_file = hkl_name.replace('.hkl', '_CC.dat')

    snr_cutoff = get_d_at_snr_one(SNR_file) if os.path.exists(SNR_file) else 10.0
    cc_cutoff = get_d_at_cc_threshold(CC_file) if os.path.exists(CC_file) else 10.0

    return max(1.5, snr_cutoff, cc_cutoff)


def resolve_pdb_path(pdb, hkl_input_file, cell_path):
    if os.path.exists(pdb):
        return pdb

    local_dir = os.path.join(os.path.dirname(hkl_input_file), os.path.basename(pdb))
    if os.path.exists(local_dir):
        return local_dir

    if cell_path:
        base_name = os.path.basename(hkl_input_file)
        for ext in ['cell', 'pdb']:
            alt_path = os.path.join(cell_path, base_name.replace('hkl', ext))
            if os.path.exists(alt_path):
                return alt_path
    return None


def prep_for_calculating_overall_statistics(
    hkl_input_file, offset, cell_path, Rfree_Rwork_path=None, nsh=10):
    if not os.path.exists(hkl_input_file):
        print(f"{os.path.basename(hkl_input_file)} does not exist.")
        return (None, ) * 6

    hkl_name = os.path.basename(hkl_input_file)
    run_name = hkl_name.replace('.hkl', '')

    Rwork = Rfree = resolution_low = None

    # Phenix results
    if Rfree_Rwork_path:
        Rwork, Rfree, resolution_cut_off_high, resolution_low, _ = get_resolution_cutoff_from_logs(
            hkl_name, Rfree_Rwork_path, run_name
        )
    else:
        resolution_cut_off_high = fallback_resolution_cutoff(hkl_name)

    # Get unit cell or PDB file
    pdb = get_UC(hkl_input_file)
    pdb = resolve_pdb_path(pdb, hkl_input_file, cell_path)

    if not pdb or not os.path.exists(pdb):
        print(f"No cell/pdb file exists for {hkl_input_file}")
        return (None, ) * 6

    pg = get_pg(hkl_input_file)
    resolution_cut_off_new = resolution_cut_off_high + offset

    CCstar_dat_file, error_file = run_partialator(
        hkl_input_file, resolution_cut_off_new, pg, pdb, nsh, str(offset)
    )

    return CCstar_dat_file, error_file, Rwork, Rfree, resolution_cut_off_new, resolution_low
