import os
import re
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


def parsing_phenix_pdb_file(phenix_pdb_file):
    """
    Parse the Phenix PDB file to extract resolution cut-off values and R-factors from REMARK 3 lines.

    Parameters:
        phenix_pdb_file (str): Path to the Phenix PDB file.

    Returns:
        tuple: (Rwork, Rfree, resolution_cut_off_high, resolution_cut_off_low)

    Raises:
        FileNotFoundError: If the specified file does not exist.

    Example:
        >>> parsing_phenix_pdb_file("refine_001.pdb")
        (0.2103, 0.2512, 1.67, 50.0)
    """
    values = {
        "RESOLUTION RANGE HIGH (ANGSTROMS)": None,
        "RESOLUTION RANGE LOW  (ANGSTROMS)": None,
        "R VALUE            (WORKING SET)": None,
        "FREE R VALUE": None,
    }

    with open(phenix_pdb_file, 'r') as file:
        for line in file:
            if "REMARK   3" not in line:
                continue
            for key in values:
                if key in line:
                    try:
                        values[key] = float(line.split(":")[1].strip())
                    except (IndexError, ValueError):
                        continue

    return (
        values["R VALUE            (WORKING SET)"],
        values["FREE R VALUE"],
        values["RESOLUTION RANGE HIGH (ANGSTROMS)"],
        values["RESOLUTION RANGE LOW  (ANGSTROMS)"],
    )

def parsing_phenix_log_file(phenix_log_file):
    """
    Parse the Phenix log file to extract R-factors and resolution cut-off.

    Parameters:
        phenix_log_file (str): Path to the Phenix log file to be parsed.

    Returns:
        tuple: A tuple containing Rwork (float), Rfree (float), and resolution_cut_off (float or None).

    Raises:
        FileNotFoundError: If the Phenix log file does not exist.
    """
    with open(phenix_log_file, 'r') as file:
        content = file.read()

    # Extract final R-work and R-free values
    r_match = re.search(r"Final R-work = ([\d.+-]+), R-free = ([\d.+-]+)", content)
    if not r_match:
        return None, None, None

    Rwork = float(r_match.group(1))
    Rfree = float(r_match.group(2))

    # Extract resolution range
    res_match = re.search(r"Resolution range: ([\d.+-]+) ([\d.+-]+)", content)
    if not res_match:
        return Rwork, Rfree, None

    resolution_cut_off = round(min(float(res_match.group(1)), float(res_match.group(2))), 3)
    return Rwork, Rfree, resolution_cut_off


def parse_phenix_results(pdb_file=None, log_file=None):
    """
    Determine which Phenix output format to use and parse accordingly.
    
    Parameters:
        pdb_file (str): Path to the Phenix PDB file.
        log_file (str): Path to the Phenix log file.
    
    Returns:
        tuple: (Rwork, Rfree, resolution_high, resolution_low)
               If using old format, resolution_low will be None.
    Raises:
        FileNotFoundError: If neither file is provided or found.
    """

    # Try parsing from the PDB if provided
    if pdb_file and os.path.exists(pdb_file):
        with open(pdb_file, 'r') as f:
            content = f.read()
            if "REMARK   3   RESOLUTION RANGE HIGH" in content and "REMARK   3   R VALUE" in content:
                print("Detected new Phenix format. Parsing PDB...")
                return parsing_phenix_pdb_file(pdb_file)
    
    # Try parsing from the LOG if provided
    if log_file and os.path.exists(log_file):
        print("Falling back to old Phenix format. Parsing log...")
        Rwork, Rfree, resolution = parsing_phenix_log_file(log_file)
        return float(Rwork), float(Rfree), resolution, None  # resolution_low is None

    # If nothing found
    raise FileNotFoundError("Could not determine valid Phenix output format. Check input files.")
