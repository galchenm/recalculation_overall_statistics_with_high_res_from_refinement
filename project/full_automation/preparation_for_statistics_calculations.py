import subprocess
import shlex
import re
import os
import glob
from partialator_execution import run_partialator


def get_UC(hkl_input_file):
    UC_filename = None
    try:
        command = f'grep -e indexamajig {hkl_input_file}'
        result = subprocess.check_output(shlex.split(command)).decode('utf-8').strip().split('\n')[0]
        UC_filename = '/'+re.findall(r"\b\S+\.cell", result)[0] if len(re.findall(r"\b\S+\.cell", result)) > 0 else re.findall(r"\b\S+\.pdb", result)[0]
    except subprocess.CalledProcessError:
        pass
    return UC_filename

def get_pg(hkl_input_file):
    point_group = None
    try:
        command = f'grep -e Symmetry {hkl_input_file}'
        result = subprocess.check_output(shlex.split(command)).decode('utf-8').strip().split('\n')[0]
        point_group = result.split(': ')[-1].strip()
    except subprocess.CalledProcessError:
        pass
    return point_group



def prep_for_calculating_overall_statistics(hkl_input_file, phenix_pdb_file, resolution_cut_off_high, offset, cell_path, nsh):
    """Prepare for calculating overall statistics by checking the existence of necessary files and running the partialator.
    This function checks if the provided hkl_input_file exists, retrieves the unit cell file, and runs the partialator
    with the specified parameters. If the unit cell file is not found, it attempts to locate
    it in the same directory as the hkl_input_file or in a specified cell_path. If successful, it returns the path to the
    CCstar.dat file and any error filename to parse. If any required file is missing,
    it returns (None, None).
    Parameters:
        hkl_input_file (str): Path to the input hkl file.
        phenix_pdb_file (str): Path to the Phenix PDB file.
        resolution_cut_off_high (float): High resolution cutoff for the calculations.
        offset (float): Offset to be applied to the resolution cutoff.
        cell_path (str): Path to the directory where the unit cell file might be located.
        nsh (int): Number of shells for the calculations.
    Returns:
        tuple: A tuple containing the path to the CCstar.dat file and any error filename to
        parse, or (None, None) if any required file is missing.
    """
    
    # Check if hkl_input_file exists
    if not hkl_input_file:
        print(f"No hkl file provided for {phenix_pdb_file}")
        return (None, None)
    
    hkl_name = os.path.basename(hkl_input_file)

    # Get the unit cell file path
    pdb = get_UC(hkl_input_file)

    # Check if the file exists in the original location
    if not os.path.exists(pdb):
        # Check if the file exists in the same directory as `hkl_input_file`
        potential_path = os.path.join(os.path.dirname(hkl_input_file), os.path.basename(pdb))
        pdb = potential_path if os.path.exists(potential_path) else None

    # Adjust the path if `cell_path` is provided
    if cell_path is not None:
        # Attempt to replace the extension with '.cell'
        cell_path_file = os.path.join(cell_path, os.path.basename(hkl_input_file).replace('hkl', 'cell'))

        if os.path.exists(cell_path_file):
            pdb = cell_path_file
        else:
            # Fallback to replacing the extension with '.pdb'
            pdb_path_file = os.path.join(cell_path, hkl_input_file.replace('hkl', 'pdb'))
            pdb = pdb_path_file if os.path.exists(pdb_path_file) else None

    pg = get_pg(hkl_input_file)
    resolution_cut_off_new = resolution_cut_off_high + offset
    if pdb and os.path.exists(pdb):
        CCstar_dat_file, error_filename_to_parse = run_partialator(
            hkl_input_file, resolution_cut_off_new, pg, pdb, nsh, str(offset)
        )
        return (CCstar_dat_file, error_filename_to_parse)
    else:
        print(f'No cell/pdb file exists for {hkl_input_file}')
        return (None, None)
    