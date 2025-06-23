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

os.nice(0)

#indexes = ['Num. patterns/Num. hits', 'Indexed patterns/Indexed crystals', 'Resolution', 'Rsplit (%)', 'CC1/2', 'CC*', 'CCano', 'SNR', 'Completeness (%)', 'Multiplicity','Total Measurements' ,'Unique Reflections', 'Wilson B-factor', 'Resolution SNR=1', 'Resolution CC>=0.3', 'a, b, c, alpha, betta, gamma', 'a', 'b', 'c', 'alpha', 'betta', 'gamma', 'Num. patterns', 'Indexed patterns', 'Num. hits', 'Indexed crystals', 'Resolution overall', 'Rsplit (%) overall', 'CC1/2 overall', 'CC* overall', 'SNR overall', 'Completeness (%) overall', 'Multiplicity overall','Total Measurements overall' ,'Unique Reflections overall' ,'Resolution outer shell', 'Rsplit (%) outer shell', 'CC1/2 outer shell', 'CC* outer shell', 'SNR outer shell', 'Completeness (%) outer shell', 'Multiplicity outer shell','Total Measurements outer shell' ,'Unique Reflections outer shell', 'UC file']
indexes = ['Num. patterns/Num. hits', 'Indexed patterns/Indexed crystals', 'Resolution', 'Rsplit(%)', 'CC1/2', 'CC*', 'CCano', 'SNR', 'Completeness(%)', 'Multiplicity','Total Measurements' ,'Unique Reflections', 'Wilson B-factor', 'Resolution SNR=1', 'Resolution CC>=0.3', 'a,b,c,alpha,betta,gamma']

x_arg_name = 'd'
y_arg_name = 'CC*'
y_arg_name2 = 'Rsplit/%'


data_info = defaultdict(dict)

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
    parser.add_argument('--e', action="store_true", help='Use this flag if you want to get an extended version of results')
    return parser.parse_args()

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



def parse_err(name_of_run, filename, CCstar_dat_file, is_extended):
    resolution = ''
    Rsplit = ''
    CC = ''
    CCstar = ''
    snr = ''
    completeness = ''
    multiplicity = '' # it is the same as redanduncy
    total_measuremenets = ''
    unique_reflections = ''
    Wilson_B_factor = ''
    CCano = ''

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
            if line.startswith('Overall Rsplit = '):
                Rsplit = re.search(r'\d+\.\d+',line).group(0)
            if line.startswith('Overall CC = '):
                CC = re.search(r'\d+\.\d+',line).group(0)
            if line.startswith('Fixed resolution range: '):
                resolution = line[line.find("(")+1:line.find(")")].replace('to','-').replace('Angstroms','').strip()
            if ' measurements in total.' in line:
                total_measuremenets = re.search(r'\d+', line).group(0)
            if ' reflections in total.' in line:
                unique_reflections = re.search(r'\d+', line).group(0)
            if line.startswith('Overall <snr> ='):
                snr = re.search(r'\d+\.\d+',line).group(0)
            if line.startswith('Overall redundancy ='):
                multiplicity = re.search(r'\d+\.\d+',line).group(0)                
            if line.startswith('Overall completeness ='):
                completeness = re.search(r'\d+\.\d+',line).group(0)
            if line.startswith('B ='):
                
                Wilson_B_factor = re.search(r'\d+\.\d+',line).group(0) if re.search(r'\d+\.\d+',line) is not None else ''
            if line.startswith('Overall CCano ='):
                CCano = line.split('=')[-1].strip()
    
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
        data_info[name_of_run]['Rsplit(%)_outer_shell'] =  f'{Rsplit_shell}'
        data_info[name_of_run]['CC1/2_outer_shell'] =  f'{CC_shell}'
        data_info[name_of_run]['CC*_outer_shell'] = f'{CCstar_shell}'
        data_info[name_of_run]['SNR_outer_shell'] =  f'{SNR_shell}'
        data_info[name_of_run]['Completeness(%)_outer_shell'] = f'{Completeness_shell}'
        data_info[name_of_run]['Multiplicity_outer_shell'] =  f'{multiplicity_shell}'
        data_info[name_of_run]['Unique_Reflections_outer_shell'] =  f'{unique_refs_shell}'

    return min_shell   

def outer_shell(CCstar_dat_file):
    shell, CCstar_shell = '',''
    with open(CCstar_dat_file, 'r') as file:
        for line in file:
            line = re.sub(' +',' ', line.strip()).split(' ')
            CCstar_shell = line[1]
            shell = line[3]
    
    max_shell, min_shell, SNR_shell, Completeness_shell, unique_refs_shell, multiplicity_shell = '','', '', '', '', ''
    SNR_dat_file = CCstar_dat_file.replace('CCstar', 'SNR')
    with open(SNR_dat_file, 'r') as file:
        for line in file:
            line = re.sub(' +',' ', line.strip()).split(' ')
            unique_refs_shell, Completeness_shell, multiplicity_shell, SNR_shell, max_shell, min_shell = line[1], line[3], line[5], line[6], line[-2], line[-1]

    
    Rsplit_dat_file = CCstar_dat_file.replace("CCstar","Rsplit")
    Rsplit_shell = ''
    with open(Rsplit_dat_file, 'r') as file:
        for line in file:
            line = re.sub(' +',' ', line.strip()).split(' ')
            Rsplit_shell = line[1]
            
    CC_dat_file = CCstar_dat_file.replace("CCstar","CC")
    CC_shell = ''
    with open(CC_dat_file, 'r') as file:
        for line in file:
            line = re.sub(' +',' ', line.strip()).split(' ')
            CC_shell = line[1]
    
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

    resolution = round(0.98*(b2-b1)/(k1-k2), 3)
    return resolution

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

def get_UC(hkl_input_file):
    UC_filename = None
    try:
        command = f'grep -e indexamajig {hkl_input_file}'
        result = subprocess.check_output(shlex.split(command)).decode('utf-8').strip().split('\n')[0]
        UC_filename = '/'+re.findall(r"\b\S+\.cell", result)[0] if len(re.findall(r"\b\S+\.cell", result)) > 0 else re.findall(r"\b\S+\.pdb", result)[0]
    except subprocess.CalledProcessError:
        pass
    return UC_filename


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

def processing(CCstar_dat_files, is_extended):
    global main_path
    
    for CCstar_dat_file in CCstar_dat_files:    
        name_of_run = os.path.basename(CCstar_dat_file).split(".")[0].replace("_CCstar","")
        data_info[name_of_run] = {i:'' for i in indexes}
        
        Rsplit_dat_file = CCstar_dat_file.replace("CCstar","Rsplit")
        SNR_dat_file = CCstar_dat_file.replace("CCstar","SNR")
        CC_dat_file = CCstar_dat_file.replace("CCstar","CC")
        hkl_file = CCstar_dat_file.replace("_CCstar.dat",".hkl")
        UC_file = get_UC(hkl_file)
        UC_file = UC_file if os.path.exists(UC_file) else (os.path.join(os.path.dirname(hkl_file), os.path.basename(UC_file)) if os.path.exists(os.path.join(os.path.dirname(hkl_file), os.path.basename(UC_file))) else None)
        
        
        stream = CCstar_dat_file.replace("_CCstar.dat",".stream")
        chunks, hits, indexed_patterns, indexed = parsing_stream(stream)
        
        data_info[name_of_run]['Num. patterns/Num. hits'] = str(chunks)+"/"+str(hits)
        data_info[name_of_run]['Indexed patterns/Indexed crystals'] = str(indexed_patterns)+"/"+str(indexed)
        
        d_at_snr_one = get_d_at_snr_one(SNR_dat_file)
        data_info[name_of_run]['Resolution SNR=1'] = str(d_at_snr_one)

        d_at_cc_threshold = get_d_at_cc_threshold(CC_dat_file)
        data_info[name_of_run]['Resolution CC>=0.3'] = str(d_at_cc_threshold) 

        a, b, c, al, be, ga = (parse_UC_file(UC_file) if UC_file is not None else (None, None, None, None, None, None))
        data_info[name_of_run]['a,b,c,alpha,betta,gamma'] = a, b, c, al, be, ga

        filename_to_parse = glob.glob(f'{main_path}/**/{os.path.basename(CCstar_dat_file).replace("_CCstar.dat","*.err")}', recursive=True)
        min_res = None
        if len(filename_to_parse)>0:
            filename_to_parse = max(filename_to_parse, key=os.path.getctime)
            
            min_res = parse_err(name_of_run, filename_to_parse, CCstar_dat_file, is_extended)


        if is_extended:
            data_info[name_of_run]['UC_file'] = UC_file
            data_info[name_of_run]['N_patterns'] = str(chunks)
            data_info[name_of_run]['Indexed_patterns'] = str(indexed_patterns)
            data_info[name_of_run]['N_hits'] = str(hits)
            data_info[name_of_run]['Indexed_crystals'] = str(indexed)
            

            data_info[name_of_run]['a'] = a
            data_info[name_of_run]['b'] = b
            data_info[name_of_run]['c'] = c
            data_info[name_of_run]['alpha'] = al
            data_info[name_of_run]['betta'] = be 
            data_info[name_of_run]['gamma'] = ga
            data_info[name_of_run]['CC* intersects with Rsplit at'] = f'{calculating_max_res_from_Rsplit_CCstar_dat(CCstar_dat_file, Rsplit_dat_file)}'
            get_hkl_command = f'get_hkl -i {hkl_file} -p {UC_file}  -o {hkl_file.replace("hkl", "mtz")} --output-format=mtz | tee conversion_log.txt'
            if UC_file:
                os.system(get_hkl_command)
                data_info[name_of_run]['mtz_command'] = get_hkl_command
                data_info[name_of_run]['mtz_with_the_defined_resolution'] = f'{hkl_file.replace("hkl", "mtz")} {d_at_snr_one}' if d_at_snr_one else (f'{hkl_file.replace("hkl", "mtz")} {min_res}' if min_res else '')
            else:
                data_info[name_of_run]['mtz_command'] = ''
                data_info[name_of_run]['mtz_with_the_defined_resolution'] = ''
        
        
if __name__ == "__main__":

    args = parse_cmdline_args()
    main_path = args.path_from
    output = args.output
    is_extended = args.e
    
    if args.f is not None:
        with open(args.f, 'r') as file:
            for line in file:

                
                CCstar_dat_files = glob.glob(f'{main_path}/*{line.strip()}*/**/*CCstar.dat', recursive=True)
                CCstar_dat_files = list(filter(lambda x: (args.p in x), CCstar_dat_files)) if args.p is not None else CCstar_dat_files
                
                processing(CCstar_dat_files, is_extended)
    else:
        CCstar_dat_files = glob.glob(os.path.join(main_path, '*CCstar.dat'))
        if len(CCstar_dat_files) > 0:
            
            processing(CCstar_dat_files, is_extended)
        else:    
            for path, dirs, all_files in os.walk(main_path):
                
                CCstar_dat_files = glob.glob(os.path.join(path, '*CCstar.dat'), recursive=True)
                CCstar_dat_files = list(filter(lambda x: (args.p in x), CCstar_dat_files)) if args.p is not None else CCstar_dat_files
                
                if len(CCstar_dat_files)>0:
                    processing(CCstar_dat_files, is_extended)
    
    df = pd.DataFrame.from_dict(data_info)
    print(df)
    df.to_csv(output, sep=';')
