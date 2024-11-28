#!/usr/bin/env python3
# coding: utf8
# Written by Galchenkova M.

import os
import sys
import glob
import re
import pandas as pd
import numpy as np
import argparse
from collections import defaultdict
import subprocess
import shlex
import time
import concurrent.futures


os.nice(0)
indexes = ['Num. patterns/hits', 'Indexed patterns/crystals', 'Resolution', 'Rsplit(%)', 'CC1/2', 'CC*', 'CCano', 'SNR', 'Completeness(%)', 'Multiplicity','Total Measurements' ,'Unique Reflections', 'Wilson B-factor', 'Resolution SNR=1', 'Resolution CC>=0.3', 'a,b,c,alpha,betta,gamma']

x_arg_name = 'd'
y_arg_name = 'CC*'
y_arg_name2 = 'Rsplit/%'

USER='galchenm'

data_info_all = defaultdict(dict)

class CustomFormatter(argparse.RawDescriptionHelpFormatter,
                      argparse.ArgumentDefaultsHelpFormatter):
    pass

def parse_cmdline_args():
    parser = argparse.ArgumentParser(
        description=sys.modules[__name__].__doc__,
        formatter_class=CustomFormatter)
    parser.add_argument('path_from', type=str, help="The path of folder/s that contain/s files")
    parser.add_argument('output', type=str, help="Path to lists of files")
    parser.add_argument('-f','--f', type=str, help='File with blocks')
    parser.add_argument('-p','--p', type=str, help='Pattern in filename or path')
    parser.add_argument('-c', '--cell', type=str, help='Path to the folder with cell/pdb files. Check that naming is the same as hkl files!')
    parser.add_argument('-a','--a', type=str, help='Additional path with Rfree/Rwork')
    parser.add_argument('-n', '--nshells', default=10, type=int,  help="Number of shells")
    parser.add_argument('-offset', '--offset', default=None, nargs='+', type=float, help="Value/s for offset parameter for resolution [Angstrom]")
    parser.add_argument('--e', action="store_true", help='Use this flag if you want to get an extended version of results')
    return parser.parse_args()


def parsing_stream(stream):
    try:
        res_hits = subprocess.check_output(['grep', '-rc', 'hit = 1', stream]).decode('utf-8').strip().split('\n')
        hits = int(res_hits[0])
        chunks = int(subprocess.check_output(['grep', '-c', 'Image filename',stream]).decode('utf-8').strip().split('\n')[0]) #len(res_hits)
    except subprocess.CalledProcessError:
        hits = 0
        chunks = 0

    try:
        res_indexed = subprocess.check_output(['grep', '-rc', 'Begin crystal',stream]).decode('utf-8').strip().split('\n')
        indexed = int(res_indexed[0])
    except subprocess.CalledProcessError:
        indexed = 0

    try:
        res_none_indexed_patterns = subprocess.check_output(['grep', '-rc', 'indexed_by = none',stream]).decode('utf-8').strip().split('\n')
        none_indexed_patterns = int(res_none_indexed_patterns[0])
    except subprocess.CalledProcessError:
        none_indexed_patterns = 0

    indexed_patterns = chunks - none_indexed_patterns

    return chunks, hits, indexed_patterns, indexed 


def run_partialator(hkl_input_file, highres, pg, pdb, nsh=10, suffix=''):
    path = os.path.dirname(os.path.abspath(hkl_input_file))
    os.chdir(path)
    print(f'We are in {os.getcwd()}')
    data = os.path.basename(hkl_input_file).split('.')[0]
    data_output_name = data if len(suffix) == 0 else f"{data}_offset_{suffix.replace('.', '_')}"

    if os.path.exists(f'{data}.hkl1') and os.path.exists(f'{data}.hkl2'):
        
        job_file = os.path.join(path,"%s.sh" % data_output_name)
        
        with open(job_file, 'w+') as fh:
            fh.writelines("#!/bin/sh\n")
            fh.writelines("#SBATCH --job=%s\n" % data_output_name)
            fh.writelines("#SBATCH --partition=short,upex,allcpu\n")
            fh.writelines("#SBATCH --time=12:00:00\n")
            fh.writelines("#SBATCH --nodes=1\n")
            fh.writelines("#SBATCH --nice=100\n")
            fh.writelines("#SBATCH --mem=500000\n")
            fh.writelines("#SBATCH --output=%s.out\n" % data_output_name)
            fh.writelines("#SBATCH --error=%s.err\n" % data_output_name)
            fh.writelines("source /etc/profile.d/modules.sh\n")
            fh.writelines("module load xray\n")

            fh.writelines("module load hdf5/1.10.5\n")
            fh.writelines("module load anaconda3/5.2\n")
            fh.writelines("module load maxwell crystfel\n")
            fh.writelines("export QT_QPA_PLATFORM=offscreen\n") 

            command = f"compare_hkl -p {pdb} -y {pg} --highres={highres} --nshells={nsh} --fom=CCstar --shell-file={data_output_name}_CCstar.dat {data}.hkl1 {data}.hkl2\n"
            fh.writelines(command)

            command = f"compare_hkl -p {pdb} -y {pg} --highres={highres} --nshells={nsh} --fom=Rsplit --shell-file={data_output_name}_Rsplit.dat {data}.hkl1 {data}.hkl2\n"
            fh.writelines(command)

            command = f"compare_hkl -p {pdb} -y {pg} --highres={highres} --nshells={nsh} --fom=CC --shell-file={data_output_name}_CC.dat {data}.hkl1 {data}.hkl2\n"
            fh.writelines(command)

            command = f"compare_hkl -p {pdb} -y {pg} --highres={highres} --nshells={nsh} --fom=CCano --shell-file={data_output_name}_CCano.dat {data}.hkl1 {data}.hkl2\n"
            fh.writelines(command)

            command = f"check_hkl -p {pdb} -y {pg} --highres={highres} --nshells={nsh} --shell-file={data_output_name}_SNR.dat {data}.hkl\n"
            fh.writelines(command)

            command = f"check_hkl -p {pdb} -y {pg} --highres={highres} --nshells={nsh} --wilson --shell-file={data_output_name}_Wilson.dat {data}.hkl\n"
            fh.writelines(command)
            
            max_dd = round(10./highres,3)

            command = f"python3 /gpfs/cfel/group/cxi/scratch/2020/EXFEL-2019-Schmidt-Mar-p002450/scratch/galchenm/scripts_for_work/plot_func/many_plots-upt-v2.py -i {data_output_name}_CCstar.dat -x '1/d' -y 'CC*' -o {data_output_name}.png -add_nargs {data_output_name}_Rsplit.dat -yad 'Rsplit/%' -x_lim_dw 1. -x_lim_up {max_dd} -t {data_output_name} -legend {data_output_name} >> output.err\n"
            fh.writelines(command)
        print(f'The {job_file} is going to be submitted')    
        os.system("sbatch %s" % job_file)
        
        return "%s_CCstar.dat" % os.path.join(path, data_output_name), "%s.err" % os.path.join(path, data_output_name)
    else:
        print(f'You do not have hkl1 and/or hkl2 files for {hkl_input_file}')
        return None, None

def get_d_at_snr_one(file_path):
    """Get resolution at SNR=1 using linear interpolation, handling non-monotonic SNR data."""
    data = []
    
    # Read data from file
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines[1:]:
            parts = line.split()
            if len(parts) < 9:
                continue
            try:
                snr = float(parts[6])
                d_value = float(parts[8])
                data.append((snr, d_value))
            except ValueError:
                continue

    # Filter only monotonically decreasing SNR
    filtered_data = []
    for i, (snr, d_value) in enumerate(data):
        if i == 0 or snr < filtered_data[-1][0]:  # SNR must decrease
            filtered_data.append((snr, d_value))

    # Find d(A) at SNR=1 using linear interpolation
    for i in range(len(filtered_data) - 1):
        snr1, d1 = filtered_data[i]
        snr2, d2 = filtered_data[i + 1]
        if snr1 == 1.0:
            return d1
        if snr1 > 1.0 > snr2:  # Interpolate
            return round(d1 + (d2 - d1) * (1.0 - snr1) / (snr2 - snr1), 3)

    return None  # Return None if SNR=1 not found or data insufficient


def get_d_at_cc_threshold(file_path, target_cc=0.3):
    """
    Parse the file and return the d(A) value corresponding to a specific CC value using linear interpolation if needed.

    :param file_path: Path to the CC.dat file.
    :param target_cc: The CC value for which d(A) is required. Default is 0.3.
    :return: The d(A) value corresponding to the target CC.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Skip the header and parse the rows
    data = []
    for line in lines[1:]:  # Skip the header
        parts = line.split()
        if len(parts) < 5:  # Ensure there are enough columns
            continue
        try:
            cc = float(parts[1])  # CC column
            d_value = float(parts[3])  # d(A) column
            data.append((cc, d_value))
        except ValueError:
            continue

    # Filter only monotonically decreasing CC values
    filtered_data = []
    for i, (cc, d_value) in enumerate(data):
        if i == 0 or cc < filtered_data[-1][0]:  # CC must decrease
            filtered_data.append((cc, d_value))

    # Find the exact or interpolate for the target CC
    for i in range(len(filtered_data) - 1):
        cc1, d1 = filtered_data[i]
        cc2, d2 = filtered_data[i + 1]
        
        if cc1 == target_cc:  # Exact match
            return d1
        if cc1 > target_cc > cc2:  # Interpolation
            # Linear interpolation formula: d = d1 + (d2 - d1) * (target_cc - cc1) / (cc2 - cc1)
            return round(d1 + (d2 - d1) * (target_cc - cc1) / (cc2 - cc1), 3)
    return None  # Return None if target CC is out of range


def wait_for_line(filename, line_for_checking, max_attempts=10, delay=2.0):
    """
    Waits for a specific line to appear in the file. Stops if max_attempts are exceeded.
    
    Args:
        filename (str): Path to the file being monitored.
        line_for_checking (str): The line to look for in the file.
        max_attempts (int): Maximum number of attempts to check the file.
        delay (float): Time in seconds to wait between each attempt.

    Returns:
        bool: True if the line is found, False if the line does not appear within max_attempts.
    """
    attempts = 0

    while attempts < max_attempts:
        with open(filename) as f:
            if any(line_for_checking in line for line in f):
                print(f"Line '{line_for_checking}' found.")
                return True  # Exit early if the line is found
        
        attempts += 1
        print(f"Attempt {attempts}/{max_attempts}: Line '{line_for_checking}' not found. Retrying in {delay} seconds...")
        time.sleep(delay)
    
    print(f"Line '{line_for_checking}' not found after {max_attempts} attempts. Exiting.")
    return False  # Exit if the line isn't found after max_attempts


def parse_err(data_info, name_of_run, filename, CCstar_dat_file):
    global is_extended
    print(220)
    resolution = ''
    Rsplit = ''
    CC = ''
    CCano = ''
    CCstar = ''
    snr = ''
    completeness = ''
    multiplicity = '' # it is the same as redanduncy
    total_measuremenets = ''
    unique_reflections = ''
    Wilson_B_factor = ''

    is_checked = False
    line_for_checking = 'B ='

    SNR_dat_file = CCstar_dat_file.replace("CCstar","SNR")
    CC_dat_file = CCstar_dat_file.replace("CCstar","CC")
    shell = CCstar_shell = Rsplit_shell = CC_shell = max_shell = min_shell = SNR_shell = Completeness_shell = unique_refs_shell = multiplicity_shell = ''

    if not wait_for_line(filename, line_for_checking, max_attempts=10, delay=2.0):
        print("The file is not updating, and the required line did not appear. Exiting function.")
        data_info[name_of_run]['Comment'] = 'Something odd happened with calculation overall statistics. Did not finish calculating B-factor'

    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('Overall CC* = '):
                CCstar = re.search(r'\d+\.\d+',line).group(0)
                CCstar = str(round(float(CCstar),4))
            if line.startswith('Overall Rsplit = '):
                Rsplit = re.search(r'\d+\.\d+',line).group(0)
                Rsplit = str(round(float(Rsplit),4))
            if line.startswith('Overall CC = '):
                CC = re.search(r'\d+\.\d+',line).group(0)
                CC = str(round(float(CC),4))
            if line.startswith('Overall CCano = '):
                CCano = re.search(r'\d+\.\d+',line).group(0)
                CCano = str(round(float(CCano),4)) if len(CCano) > 0 else ''           
            if line.startswith('Fixed resolution range: '):
                resolution = line[line.find("(")+1:line.find(")")].replace('to','-').replace('Angstroms','').strip()
                print(261)
            if ' measurements in total.' in line:
                total_measuremenets = re.search(r'\d+', line).group(0)
            if ' reflections in total.' in line:
                unique_reflections = re.search(r'\d+', line).group(0)
            if line.startswith('Overall <snr> ='):
                snr = re.search(r'\d+\.\d+',line).group(0)
                snr = str(round(float(snr),4))
            if line.startswith('Overall redundancy ='):
                multiplicity = re.search(r'\d+\.\d+',line).group(0)
                multiplicity = str(round(float(multiplicity),4))
            if line.startswith('Overall completeness ='):
                completeness = re.search(r'\d+\.\d+',line).group(0)
                completeness = str(round(float(completeness),4))
            if line.startswith('B ='):
                Wilson_B_factor = re.search(r'\d+\.\d+',line).group(0) if re.search(r'\d+\.\d+',line) is not None else ''
                Wilson_B_factor = str(round(float(Wilson_B_factor),4)) if len(Wilson_B_factor) > 0 else ''
    
    shell, CCstar_shell, Rsplit_shell, CC_shell, max_shell, min_shell, SNR_shell, Completeness_shell, unique_refs_shell, multiplicity_shell = outer_shell(CCstar_dat_file)

    data_info[name_of_run]['Resolution'] = resolution + f' ({max_shell} - {min_shell})'
    data_info[name_of_run]['Rsplit(%)'] = Rsplit + f' ({Rsplit_shell})'
    data_info[name_of_run]['CC1/2'] =  CC + f' ({CC_shell})'
    data_info[name_of_run]['CC*'] = CCstar  + f' ({CCstar_shell})'
    data_info[name_of_run]['CCano'] = CCano
    data_info[name_of_run]['SNR'] =  snr  + f' ({SNR_shell})'
    data_info[name_of_run]['Completeness(%)'] = completeness  + f' ({Completeness_shell})'
    data_info[name_of_run]['Multiplicity'] =  multiplicity  + f' ({multiplicity_shell})'
    data_info[name_of_run]['Total Measurements'] =  total_measuremenets
    data_info[name_of_run]['Unique Reflections'] =  unique_reflections  + f' ({unique_refs_shell})'
    data_info[name_of_run]['Wilson B-factor'] = Wilson_B_factor
    data_info[name_of_run]['Resolution CC>=0.3'] = str(get_d_at_cc_threshold(CC_dat_file))
    data_info[name_of_run]['Resolution SNR=1'] = str(get_d_at_snr_one(SNR_dat_file))

    if is_extended:
        data_info[name_of_run]['Resolution_overall'] = resolution
        data_info[name_of_run]['Rsplit(%)_overall'] = Rsplit 
        data_info[name_of_run]['CC1/2_overall'] =  CC 
        data_info[name_of_run]['CC*_overall'] = CCstar  
        data_info[name_of_run]['SNR_overall'] =  snr  
        data_info[name_of_run]['Completeness(%)_overall'] = completeness 
        data_info[name_of_run]['Multiplicity_overall'] =  multiplicity  
        data_info[name_of_run]['Unique Reflections_overall'] =  unique_reflections

        data_info[name_of_run]['Resolution_outer_shell'] =  f'{max_shell} - {min_shell}'
        data_info[name_of_run]['Rsplit(%)_outer_shesll'] =  f'{Rsplit_shell}'
        data_info[name_of_run]['CC1/2_outer_shell'] =  f'{CC_shell}'
        data_info[name_of_run]['CC*_outer_shell'] = f'{CCstar_shell}'
        data_info[name_of_run]['SNR_outer_shell'] =  f'{SNR_shell}'
        data_info[name_of_run]['Completeness(%)_outer_shell'] = f'{Completeness_shell}'
        data_info[name_of_run]['Multiplicity_outer_shell'] =  f'{multiplicity_shell}'
        data_info[name_of_run]['Unique_Reflections_outer_shell'] =  f'{unique_refs_shell}'

    return data_info    

def outer_shell(CCstar_dat_file):
    shell, CCstar_shell = '',''
    while not(os.path.exists(CCstar_dat_file)):
        time.sleep(5)
    while os.stat(CCstar_dat_file).st_size == 0:
        time.sleep(5)
           
    with open(CCstar_dat_file, 'r') as file:
        for line in file:
            line = re.sub(' +',' ', line.strip()).split(' ')
            CCstar_shell = line[1]
            shell = line[3]
    
    Rsplit_dat_file = CCstar_dat_file.replace("CCstar","Rsplit")
    Rsplit_shell = ''
    while not(os.path.exists(Rsplit_dat_file)):
        time.sleep(5)
    while os.stat(Rsplit_dat_file).st_size == 0:
        time.sleep(5)
  
    with open(Rsplit_dat_file, 'r') as file:
        for line in file:
            line = re.sub(' +',' ', line.strip()).split(' ')
            Rsplit_shell = line[1]
            
    CC_dat_file = CCstar_dat_file.replace("CCstar","CC")
    CC_shell = ''
    while not(os.path.exists(CC_dat_file)): 
        time.sleep(5)   
    while os.stat(CC_dat_file).st_size == 0:
        time.sleep(5)
    with open(CC_dat_file, 'r') as file:
        for line in file:
            line = re.sub(' +',' ', line.strip()).split(' ')
            CC_shell = line[1]
    
    max_shell, min_shell, SNR_shell, Completeness_shell, unique_refs_shell, multiplicity_shell = '','', '', '', '', ''
    SNR_dat_file = CCstar_dat_file.replace('CCstar', 'SNR')
    
    while not(os.path.exists(SNR_dat_file)):
        time.sleep(5)
    while os.stat(SNR_dat_file).st_size == 0:
        time.sleep(5)
    
    with open(SNR_dat_file, 'r') as file:
        for line in file:
            line = re.sub(' +',' ', line.strip()).split(' ')
            unique_refs_shell, Completeness_shell, multiplicity_shell, SNR_shell, max_shell, min_shell = line[1], line[3], line[5], line[6], line[-2], line[-1]
                
    return shell, round(float(CCstar_shell),3), round(float(Rsplit_shell),3), round(float(CC_shell),3), round(10/float(max_shell),2), round(10/float(min_shell),2), round(float(SNR_shell),3), round(float(Completeness_shell),3), unique_refs_shell, multiplicity_shell  


def get_xy(file_name, x_arg_name, y_arg_name):
    x = []
    y = []

    with open(file_name, 'r') as stream:
        for line in stream:
            if y_arg_name in line:
                tmp = line.replace('1/nm', '').replace('# ', '').replace('centre', '').replace('/ A', '').replace(' dev','').replace('(A)','')
                tmp = tmp.split()
                y_index = tmp.index(y_arg_name)
                x_index = tmp.index(x_arg_name)

            else:
                tmp = line.split()
                
                x.append(float(tmp[x_index]) if not np.isnan(float(tmp[x_index])) else 0. )
                y.append(float(tmp[y_index]) if not np.isnan(float(tmp[y_index])) else 0. )

    x = np.array(x)
    y = np.array(y)

    list_of_tuples = list(zip(x, y))
    
    df = pd.DataFrame(list_of_tuples, 
                  columns = [x_arg_name, y_arg_name])
    
    df = df[df[y_arg_name].notna()]
    df = df[df[y_arg_name] >= 0.]
    return df[x_arg_name], df[y_arg_name]

def calculating_max_res_from_Rsplit_CCstar_dat(CCstar_dat_file, Rsplit_dat_file):
    d_CCstar, CCstar = get_xy(CCstar_dat_file, x_arg_name, y_arg_name)
    CCstar *= 100
    
    d_Rsplit, Rsplit = get_xy(Rsplit_dat_file, x_arg_name, y_arg_name2)
    
    i = 0

    CC2, d2 = CCstar[0], d_CCstar[0]
    CC1, d1 = 0., 0.

    Rsplit2 = Rsplit[0]
    Rsplit1 = 0.

    while Rsplit[i]<=CCstar[i] and i < len(d_CCstar):
        CC1, d1, Rsplit1 = CC2, d2, Rsplit2
        i+=1
        try:
            CC2, d2, Rsplit2 = CCstar[i], d_CCstar[i], Rsplit[i]
        except:
            return -1000            
        if Rsplit[i]==CCstar[i]:
            resolution = d_CCstar[i]
            return resolution
            
    k1 = round((CC2-CC1)/(d2-d1),3)
    b1 = round((CC1*d2-CC2*d1)/(d2-d1),3)     

    k2 = round((Rsplit2-Rsplit1)/(d2-d1),3)
    b2 = round((Rsplit1*d2-Rsplit2*d1)/(d2-d1),3)

    resolution = round(0.98*(b2-b1)/(k1-k2),3)
    return resolution

def parse_cryst1_from_pdb(file_path):
    # Regular expression to match the CRYST1 line and capture the parameters
    pattern = r"CRYST1\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)"
    
    # Open and read the PDB file
    with open(file_path, 'r') as file:
        for line in file:
            # Search for the CRYST1 line
            if line.startswith("CRYST1"):
                # Match and capture the unit cell parameters
                match = re.match(pattern, line.strip())
                if match:
                    # Extract the values from the match
                    a = float(match.group(1))
                    b = float(match.group(2))
                    c = float(match.group(3))
                    alpha = float(match.group(4))
                    beta = float(match.group(5))
                    gamma = float(match.group(6))
                    
                    return a, b, c, alpha, beta, gamma
    # If no CRYST1 line is found
    raise ValueError("No CRYST1 line found in the provided PDB file.")

def parse_UC_file(UC_file):
    if UC_file.endswith('pdb'):
        return parse_cryst1_from_pdb(UC_file)
    else:
        print(461)
        # Regular expressions to capture the unit cell parameters
        pattern = r"a = ([\d.]+) A.*?b = ([\d.]+) A.*?c = ([\d.]+) A.*?al = ([\d.]+) deg.*?be = ([\d.]+) deg.*?ga = ([\d.]+) deg"
        
        # Read the content of the file
        with open(UC_file, 'r') as file:
            file_content = file.read()

        # Search for the pattern in the file content
        match = re.search(pattern, file_content, re.DOTALL)

        if match:
            # Extract values from the match
            a = float(match.group(1))
            b = float(match.group(2))
            c = float(match.group(3))
            al = float(match.group(4))
            be = float(match.group(5))
            ga = float(match.group(6))
            return a, b, c, al, be, ga
        else:
            raise ValueError("Unit cell parameters not found in the provided file.")
            return None, None, None, None, None, None

def processing(input_tuple):
    global is_extended
    global cell_path

    data_info = defaultdict(dict)
    CCstar_dat_file, error_filename_to_parse, Rwork, Rfree, resolution_cut_off_high, resolution_cut_off_low = input_tuple
    
    name_of_run = os.path.basename(CCstar_dat_file).split(".")[0].replace("_CCstar","")
    print(f'Processing {name_of_run}')
    
    data_info[name_of_run] = {i:'' for i in indexes}
    stream = CCstar_dat_file.replace("_CCstar.dat",".stream") if "_offset_" not in CCstar_dat_file else CCstar_dat_file.split("_offset_")[0]+'.stream'
    chunks, hits, indexed_patterns, indexed = parsing_stream(stream)
    

    Rsplit_dat_file = CCstar_dat_file.replace("CCstar","Rsplit")
    hkl_name = CCstar_dat_file.replace("_CCstar.dat",".hkl") if "_offset_" not in CCstar_dat_file else CCstar_dat_file.split("_offset_")[0]+'.hkl'

    # Get the unit cell file path
    UC_file = get_UC(hkl_name)
    print('UC_file = ', UC_file)
    # Check if the file exists in the original location
    if not os.path.exists(UC_file):
        # Check if the file exists in the same directory as `hkl_name`
        potential_path = os.path.join(os.path.dirname(hkl_name), os.path.basename(UC_file))
        print('potential_path = ', potential_path)
        if os.path.exists(potential_path):
            UC_file = potential_path
        else:
            UC_file = None

    # Adjust the path if `cell_path` is provided
    if cell_path is not None:
        # Attempt to replace the extension with '.cell'
        cell_path_file = os.path.join(cell_path, os.path.basename(hkl_name).replace('hkl', 'cell'))
        print('cell_path_file = ', cell_path_file)
        if os.path.exists(cell_path_file):
            
            UC_file = cell_path_file
        else:
            # Fallback to replacing the extension with '.pdb'
            pdb_path_file = os.path.join(cell_path, hkl_name.replace('hkl', 'pdb'))
            UC_file = pdb_path_file if os.path.exists(pdb_path_file) else None

    while not(os.path.exists(Rsplit_dat_file)):
        time.sleep(5)
    while os.stat(Rsplit_dat_file).st_size == 0:
        time.sleep(5)
    data_info[name_of_run]['CC* intersects with Rsplit at'] = f'{calculating_max_res_from_Rsplit_CCstar_dat(CCstar_dat_file, Rsplit_dat_file)}'
    
    data_info[name_of_run]['Num. patterns/hits'] = str(chunks)+"/"+str(hits)
    data_info[name_of_run]['Indexed patterns/crystals'] = str(indexed_patterns)+"/"+str(indexed)
    data_info[name_of_run]['Rwork/Rfree'] = str(Rwork)  +"/"+ str(Rfree)
    data_info[name_of_run]['Refinement resolution cut-off high'] = str(resolution_cut_off_high)
    data_info[name_of_run]['Refinement resolution cut-off low'] = str(resolution_cut_off_low)   
    data_info = parse_err(data_info, name_of_run, error_filename_to_parse, CCstar_dat_file)
    a, b, c, al, be, ga = (parse_UC_file(UC_file) if UC_file is not None else (None, None, None, None, None, None))
    data_info[name_of_run]['a,b,c,alpha,betta,gamma'] = a, b, c, al, be, ga

    if is_extended:
        data_info[name_of_run]['UC_file'] = UC_file
        data_info[name_of_run]['N_patterns'] = str(chunks)
        data_info[name_of_run]['Indexed_patterns'] = str(indexed_patterns)
        data_info[name_of_run]['N_hits'] = str(hits)
        data_info[name_of_run]['Indexed_crystals'] = str(indexed)
        data_info[name_of_run]['Rwork'] = str(Rwork) 
        data_info[name_of_run]['Rfree'] = str(Rfree)

        data_info[name_of_run]['a'] = a
        data_info[name_of_run]['b'] = b
        data_info[name_of_run]['c'] = c
        data_info[name_of_run]['alpha'] = al
        data_info[name_of_run]['betta'] = be 
        data_info[name_of_run]['gamma'] = ga
        data_info[name_of_run]['CC* intersects with Rsplit at'] = f'{calculating_max_res_from_Rsplit_CCstar_dat(CCstar_dat_file, Rsplit_dat_file)}'

    return data_info


def parsing_phenix_pdb_file(phenix_pdb_file):
    resolution_cut_off_high = None
    resolution_cut_off_low = None
    Rfree = None
    Rwork = None

    # Open and read the file line by line
    with open(phenix_pdb_file, 'r') as file:
        for line in file:
            # Check for each specific line and extract the number using split
            if "REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS)" in line:
                resolution_cut_off_high = float(line.split(":")[1].strip())
            elif "REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS)" in line:
                resolution_cut_off_low = float(line.split(":")[1].strip())
            elif "REMARK   3   R VALUE            (WORKING SET)" in line:
                Rwork = float(line.split(":")[1].strip())
            elif "REMARK   3   FREE R VALUE                     :" in line:
                Rfree = float(line.split(":")[1].strip())

    return (Rwork, Rfree, resolution_cut_off_high, resolution_cut_off_low)

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

def prep_for_calculating_overall_statistics(input_prior_information):
    global cell_path

    phenix_pdb_file, resolution_cut_off_high, resolution_cut_off_low, offset, Rwork, Rfree, cell_path = input_prior_information
    pattern_in_hkl_input_filename = os.path.basename(os.path.dirname(phenix_pdb_file))
    hkl_input_file = glob.glob(f'{main_path}/{pattern_in_hkl_input_filename}.hkl', recursive=True)
    
    if len(hkl_input_file) == 0:
        print(f'In {main_path} there is no hkl file associated with {phenix_pdb_file}')
        return (None, None) 
    hkl_input_file = hkl_input_file[0]
    hkl_name = os.path.basename(hkl_input_file)

    # Get the unit cell file path
    pdb = get_UC(hkl_input_file)

    # Check if the file exists in the original location
    if not os.path.exists(pdb):
        # Check if the file exists in the same directory as `hkl_input_file`
        potential_path = os.path.join(os.path.dirname(hkl_input_file), os.path.basename(pdb))
        print('potential_path = ', potential_path)
        if os.path.exists(potential_path):
            pdb = potential_path
        else:
            pdb = None

    # Adjust the path if `cell_path` is provided
    if cell_path is not None:
        # Attempt to replace the extension with '.cell'
        cell_path_file = os.path.join(cell_path, os.path.basename(hkl_input_file).replace('hkl', 'cell'))
        print('cell_path_file = ', cell_path_file)
        if os.path.exists(cell_path_file):
            
            pdb = cell_path_file
        else:
            # Fallback to replacing the extension with '.pdb'
            pdb_path_file = os.path.join(cell_path, hkl_input_file.replace('hkl', 'pdb'))
            pdb = pdb_path_file if os.path.exists(pdb_path_file) else None
    
    
    pg = get_pg(hkl_input_file)
    resolution_cut_off_new = resolution_cut_off_high + offset
    if os.path.exists(pdb):
        CCstar_dat_file, error_filename_to_parse = run_partialator(hkl_input_file, resolution_cut_off_new, pg, pdb, nsh, str(offset))
        return (CCstar_dat_file, error_filename_to_parse)
    else:
        print(f'No cell/pdb file exists for {hkl_input_file}')  
        return (None, None)   

if __name__ == "__main__":

    args = parse_cmdline_args()
    main_path = args.path_from
    output = args.output
    Rfree_Rwork_path = args.a
    nsh = args.nshells
    offsets = args.offset
    is_extended = args.e
    cell_path = args.cell
    
    if offsets is None:
        offsets = [0.]
    else:
        offsets.append(0.)
    
    Rwork, Rfree, resolution_cut_off_high, resolution_cut_off_low = '', '', '', ''
    phenix_pdb_files = []
    
    if args.f is not None:
        with open(args.f, 'r') as file:
            runs = file.read().split('\n')
        
        
        for run in runs:
            phenix_pdb_file = glob.glob(f'{Rfree_Rwork_path}/*{run.strip()}*/**/*.pdb', recursive=True)
            phenix_pdb_file = max(phenix_pdb_file, key=os.path.getctime)
            phenix_pdb_file = list(filter(lambda x: (args.p in x), phenix_pdb_file)) if args.p is not None else phenix_pdb_file
            phenix_pdb_files.append(phenix_pdb_file)
       
        CC_dat_with_err_files_to_parse = []

        with concurrent.futures.ProcessPoolExecutor() as executor:
            for phenix_pdb_file, tuple_information in zip(phenix_pdb_files, executor.map(parsing_phenix_pdb_file, phenix_pdb_files)):
                Rwork, Rfree, resolution_cut_off_high, resolution_cut_off_low = tuple_information
                
                if resolution_cut_off_high is not None:
                    input_prior_information = list(map(lambda offset: (phenix_pdb_file, resolution_cut_off_high, resolution_cut_off_low, offset, Rwork, Rfree, cell_path), offsets))
                    with concurrent.futures.ProcessPoolExecutor() as executor2:
                        for tuple_information_CCstar_err_file in executor2.map(prep_for_calculating_overall_statistics, input_prior_information):
                            CCstar_dat_file, error_filename_to_parse = tuple_information_CCstar_err_file
                            
                            if CCstar_dat_file is None and error_filename_to_parse is None:
                                continue
                            #CCstar_dat_file, error_filename_to_parse, Rwork, Rfree, resolution_cut_off
                            CCstar_err_file_Rflags_res_cut_off = (*tuple_information_CCstar_err_file, Rwork, Rfree, resolution_cut_off_high, resolution_cut_off_low)
                            CC_dat_with_err_files_to_parse.append(CCstar_err_file_Rflags_res_cut_off)
        
        #check how many processes are pending in order not to submit
        pending_command = f'squeue -u {USER} -t pending'
        number_of_pending_processes = subprocess.check_output(shlex.split(pending_command)).decode('utf-8').strip().split('\n')
        while len(number_of_pending_processes) > 1:
            time.sleep(5)
        
        existed_files_to_check = [ CC_dat_with_err_file for CC_dat_with_err_file in CC_dat_with_err_files_to_parse if os.path.exists(CC_dat_with_err_file[0]) and os.path.exists(CC_dat_with_err_file[1])]
        for CC_dat_with_err_file in existed_files_to_check:
            processing(*CC_dat_with_err_file)
    else:

        for directory in os.listdir(Rfree_Rwork_path):
            path = os.path.join(Rfree_Rwork_path, directory)
            if not os.path.isdir(path):
                continue
            phenix_pdb_file = glob.glob(os.path.join(path, '*.pdb')) #, recursive=True)
            if len(phenix_pdb_file) == 0:
                continue
            phenix_pdb_file = [i for i in phenix_pdb_file if 'dimple' not in os.path.basename(i)]    
            phenix_pdb_file = max(phenix_pdb_file, key=os.path.getctime)

            #phenix_pdb_file = list(filter(lambda x: (args.p in x), phenix_pdb_file)) if args.p is not None else phenix_pdb_file
            phenix_pdb_files.append(phenix_pdb_file)
         
        
        CC_dat_with_err_files_to_parse = []

        input_prior_information = []
        with concurrent.futures.ProcessPoolExecutor() as executor:
            for phenix_pdb_file, tuple_information in zip(phenix_pdb_files, executor.map(parsing_phenix_pdb_file, phenix_pdb_files)):
                
                Rwork, Rfree, resolution_cut_off_high, resolution_cut_off_low = tuple_information
                
                if resolution_cut_off_high is not None:
                    input_prior_information+=list(map(lambda offset: (phenix_pdb_file, resolution_cut_off_high, resolution_cut_off_low, offset, Rwork, Rfree, cell_path), offsets))
        
        
        with concurrent.futures.ProcessPoolExecutor() as executor2:
            for info, tuple_information_CCstar_err_file in zip(input_prior_information, executor2.map(prep_for_calculating_overall_statistics, input_prior_information)):
                
                phenix_pdb_file, resolution_cut_off_high, resolution_cut_off_low, offset, Rwork, Rfree, cell_path = info
                CCstar_dat_file, error_filename_to_parse = tuple_information_CCstar_err_file
                if CCstar_dat_file is None and error_filename_to_parse is None:
                    continue
                CCstar_err_file_Rflags_res_cut_off = (*tuple_information_CCstar_err_file, Rwork, Rfree, resolution_cut_off_high, resolution_cut_off_low)
                CC_dat_with_err_files_to_parse.append(CCstar_err_file_Rflags_res_cut_off)   
        
        
        pending_command = f'squeue -u {USER} -t pending'
        number_of_pending_processes = subprocess.check_output(shlex.split(pending_command)).decode('utf-8').strip().split('\n')
        
        while len(number_of_pending_processes) > 2:
            time.sleep(15)
            pending_command = f'squeue -u {USER} -t pending'
            number_of_pending_processes = subprocess.check_output(shlex.split(pending_command)).decode('utf-8').strip().split('\n')
            
        time.sleep(5)
        existed_files_to_check = [ CC_dat_with_err_file for CC_dat_with_err_file in CC_dat_with_err_files_to_parse if os.path.exists(CC_dat_with_err_file[0]) and os.path.exists(CC_dat_with_err_file[1])]
        
        #for CC_dat_with_err_file in existed_files_to_check:
        #    tmp_data_info = processing(CC_dat_with_err_file)
        #    data_info_all.update(tmp_data_info)

        
        with concurrent.futures.ProcessPoolExecutor() as executor:
            for result in executor.map(processing, existed_files_to_check):
                for name_of_run in result:
                    data_info_all[name_of_run] = result[name_of_run]         

    print(len(data_info_all))           
    df = pd.DataFrame.from_dict(data_info_all)
    df.to_csv(output, sep=';')
