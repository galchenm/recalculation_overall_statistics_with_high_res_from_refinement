import os
import re
import time
from partialator_utils.resolution_cutoff_determination import get_d_at_snr_one, get_d_at_cc_threshold
from stream_utils.parsing_stream import parse_stream
from partialator_utils.wait_for_file import wait_for_line
from collections import defaultdict
from partialator_utils.resolution_cutoff_determination import calculating_max_res_from_Rsplit_CCstar_dat
from partialator_utils.resolution_cutoff_determination import get_UC
from partialator_utils.wait_for_file import wait_for_file
from unit_cell_utils.parsing_UC_files import parse_UC_file

def outer_shell(CCstar_dat_file, is_extended=False):
    """
    Extracts outer shell statistics from the CCstar.dat file and related data files.

    Args:
        CCstar_dat_file (str): Path to the CCstar.dat file.
        is_extended (bool): Currently unused but preserved for interface compatibility.

    Returns:
        tuple: (
            shell, CCstar_shell, Rsplit_shell, CC_shell,
            max_shell, min_shell, SNR_shell, Completeness_shell,
            unique_refs_shell, multiplicity_shell
        )
    """

    def wait_for_file_ready(filepath, check_line=''):
        """Waits for file existence and optional specific content."""
        if check_line:
            return wait_for_line(filepath, check_line, max_attempts=20, delay=2.0)
        while not os.path.exists(filepath) or os.stat(filepath).st_size == 0:
            time.sleep(2)

    def read_single_value(filepath, value_index, shell_index=None):
        with open(filepath, 'r') as f:
            for line in f:
                parts = re.sub(r' +', ' ', line.strip()).split()
                try:
                    value = round(float(parts[value_index]), 3)
                    shell = parts[shell_index] if shell_index is not None else ''
                    return value, shell
                except (IndexError, ValueError):
                    continue
        return '', ''

    # Wait and parse CCstar
    wait_for_file_ready(CCstar_dat_file)
    CCstar_shell, shell = read_single_value(CCstar_dat_file, 1, shell_index=3)

    # Wait and parse Rsplit
    rsplit_file = CCstar_dat_file.replace("CCstar", "Rsplit")
    wait_for_file_ready(rsplit_file)
    Rsplit_shell, _ = read_single_value(rsplit_file, 1)

    # Wait and parse CC
    cc_file = CCstar_dat_file.replace("CCstar", "CC")
    wait_for_file_ready(cc_file)
    CC_shell, _ = read_single_value(cc_file, 1)

    # Wait and parse SNR and related values
    snr_file = CCstar_dat_file.replace("CCstar", "SNR")
    wait_for_file_ready(snr_file)

    max_shell = min_shell = SNR_shell = Completeness_shell = unique_refs_shell = multiplicity_shell = ''
    with open(snr_file, 'r') as f:
        for line in f:
            parts = re.sub(r' +', ' ', line.strip()).split()
            try:
                unique_refs_shell = parts[1]
                Completeness_shell = round(float(parts[3]), 3)
                multiplicity_shell = round(float(parts[5]), 3)
                SNR_shell = round(float(parts[6]), 3)
                max_shell = round(10 / float(parts[-2]), 2)
                min_shell = round(10 / float(parts[-1]), 2)
                break
            except (IndexError, ValueError, ZeroDivisionError):
                continue

    return (
        shell,
        CCstar_shell,
        Rsplit_shell,
        CC_shell,
        max_shell,
        min_shell,
        SNR_shell,
        Completeness_shell,
        unique_refs_shell,
        multiplicity_shell
    )

def parse_err(data_info, name_of_run, filename, CCstar_dat_file, is_extended=False):
    """
    Parses the error file to extract overall statistics and updates the data_info dictionary.

    Args:
        data_info (dict): Dictionary to store parsed data.
        name_of_run (str): Name of the run used as a key in data_info.
        filename (str): Path to the error file.
        CCstar_dat_file (str): Path to the CCstar.dat file.
        is_extended (bool): If True, includes additional statistics.

    Returns:
        dict: Updated data_info dictionary with parsed statistics.
    """

    def extract_value(line, pattern):
        match = re.match(pattern, line)
        return str(round(float(match.group(1)), 4)) if match else ''

    if not wait_for_line(filename, "B =", max_attempts=10, delay=2.0):
        print("The file is not updating, and the required line did not appear. Exiting function.")
        data_info[name_of_run]['Comment'] = 'Something odd happened with calculation overall statistics. Did not finish calculating B-factor'
        return data_info

    # Initialize variables
    metrics = {
        'CCstar': '', 'Rsplit': '', 'CC': '', 'CCano': '', 'snr': '',
        'completeness': '', 'multiplicity': '', 'total_measurements': '',
        'unique_reflections': '', 'Wilson_B_factor': '', 'resolution': ''
    }

    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('Overall CC* ='):
                metrics['CCstar'] = extract_value(line, r'Overall CC\* =\s*(\d+\.\d+)')
            elif line.startswith('Overall Rsplit ='):
                metrics['Rsplit'] = extract_value(line, r'Overall Rsplit =\s*(\d+\.\d+)')
            elif line.startswith('Overall CC ='):
                metrics['CC'] = extract_value(line, r'Overall CC =\s*(\d+\.\d+)')
            elif line.startswith('Overall CCano ='):
                metrics['CCano'] = extract_value(line, r'Overall CCano =\s*(\d+\.\d+)')
            elif line.startswith('Fixed resolution range:'):
                metrics['resolution'] = line.split('(')[1].split(')')[0].replace('to', '-').replace('Angstroms', '').strip()
            elif ' measurements in total.' in line:
                m = re.match(r'(\d+) measurements in total\.', line)
                if m:
                    metrics['total_measurements'] = m.group(1)
            elif ' reflections in total.' in line:
                m = re.match(r'(\d+) reflections in total\.', line)
                if m:
                    metrics['unique_reflections'] = m.group(1)
            elif line.startswith('Overall <snr> ='):
                metrics['snr'] = extract_value(line, r'Overall <snr> =\s*(\d+\.\d+)')
            elif line.startswith('Overall redundancy ='):
                metrics['multiplicity'] = extract_value(line, r'Overall redundancy =\s*(\d+\.\d+)')
            elif line.startswith('Overall completeness ='):
                metrics['completeness'] = extract_value(line, r'Overall completeness =\s*(\d+\.\d+)')
            elif line.startswith('B ='):
                metrics['Wilson_B_factor'] = extract_value(line, r'B =\s*(\d+\.\d+)')

    # File paths for extra data
    SNR_dat_file = CCstar_dat_file.replace("CCstar", "SNR")
    CC_dat_file = CCstar_dat_file.replace("CCstar", "CC")

    # Shell statistics
    shell, CCstar_shell, Rsplit_shell, CC_shell, max_shell, min_shell, \
    SNR_shell, Completeness_shell, unique_refs_shell, multiplicity_shell = outer_shell(CCstar_dat_file, is_extended)

    # Populate primary data
    run_data = {
        'Resolution': f"{metrics['resolution']} ({max_shell} - {min_shell})",
        'Rsplit(%)': f"{metrics['Rsplit']} ({Rsplit_shell})",
        'CC1/2': f"{metrics['CC']} ({CC_shell})",
        'CC*': f"{metrics['CCstar']} ({CCstar_shell})",
        'CCano': metrics['CCano'],
        'SNR': f"{metrics['snr']} ({SNR_shell})",
        'Completeness(%)': f"{metrics['completeness']} ({Completeness_shell})",
        'Multiplicity': f"{metrics['multiplicity']} ({multiplicity_shell})",
        'Total Measurements': metrics['total_measurements'],
        'Unique Reflections': f"{metrics['unique_reflections']} ({unique_refs_shell})",
        'Wilson B-factor': metrics['Wilson_B_factor'],
        'Resolution CC>=0.3': get_d_at_cc_threshold(CC_dat_file),
        'Resolution SNR=1': get_d_at_snr_one(SNR_dat_file)
    }
    data_info[name_of_run].update(run_data)

    # Extended stats if needed
    if is_extended:
        extended = {
            'Resolution_overall': metrics['resolution'],
            'Rsplit(%)_overall': metrics['Rsplit'],
            'CC1/2_overall': metrics['CC'],
            'CC*_overall': metrics['CCstar'],
            'SNR_overall': metrics['snr'],
            'Completeness(%)_overall': metrics['completeness'],
            'Multiplicity_overall': metrics['multiplicity'],
            'Unique Reflections_overall': metrics['unique_reflections'],
            'Resolution_outer_shell': f'{max_shell} - {min_shell}',
            'Rsplit(%)_outer_shell': Rsplit_shell,
            'CC1/2_outer_shell': CC_shell,
            'CC*_outer_shell': CCstar_shell,
            'SNR_outer_shell': SNR_shell,
            'Completeness(%)_outer_shell': Completeness_shell,
            'Multiplicity_outer_shell': multiplicity_shell,
            'Unique_Reflections_outer_shell': unique_refs_shell
        }
        data_info[name_of_run].update(extended)

    return data_info
