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

indexes = ['Resolution SNR=1', 'Resolution CC>=0.3', 'CC* intersects Rsplit at this resolution']
x_arg_name = 'd'
y_arg_name = 'CC*'
y_arg_name2 = 'Rsplit/%'

data_info_all = defaultdict(dict)

class CustomFormatter(argparse.RawDescriptionHelpFormatter,
                      argparse.ArgumentDefaultsHelpFormatter):
    pass

def parse_cmdline_args():
    """Parse the command-line arguments."""
    parser = argparse.ArgumentParser(
        description=sys.modules[__name__].__doc__,
        formatter_class=CustomFormatter)
    parser.add_argument('path_from', type=str, help="The path of folder/s that contain/s files")
    parser.add_argument('output', type=str, help="Path to lists of files")

    return parser.parse_args()

def get_xy(file_name, x_arg_name, y_arg_name):
    """Extract the x and y values from the given file."""
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
                x.append(float(tmp[x_index]) if not np.isnan(float(tmp[x_index])) else np.nan)
                y.append(float(tmp[y_index]) if not np.isnan(float(tmp[y_index])) else np.nan)

    x = np.array(x)
    y = np.array(y)

    list_of_tuples = list(zip(x, y))
    df = pd.DataFrame(list_of_tuples, columns=[x_arg_name, y_arg_name])
    df = df[df[y_arg_name].notna()]
    df = df[df[y_arg_name] >= 0.]
    return df[x_arg_name], df[y_arg_name]

def calculating_max_res_from_Rsplit_CCstar_dat(CCstar_dat_file, Rsplit_dat_file):
    """Calculate the maximum resolution where CC* intersects with Rsplit."""
    d_CCstar, CCstar = get_xy(CCstar_dat_file, x_arg_name, y_arg_name)
    CCstar *= 100

    d_Rsplit, Rsplit = get_xy(Rsplit_dat_file, x_arg_name, y_arg_name2)
    
    i = 0
    CC2, d2 = CCstar[0], d_CCstar[0]
    CC1, d1 = 0., 0.
    Rsplit2 = Rsplit[0]
    Rsplit1 = 0.

    while Rsplit[i] <= CCstar[i] and i < len(d_CCstar):
        CC1, d1, Rsplit1 = CC2, d2, Rsplit2
        i += 1
        try:
            CC2, d2, Rsplit2 = CCstar[i], d_CCstar[i], Rsplit[i]
        except IndexError:
            return -1000
            
        if Rsplit[i] == CCstar[i]:
            resolution = d_CCstar[i]
            return resolution
            
    k1 = round((CC2 - CC1) / (d2 - d1), 3)
    b1 = round((CC1 * d2 - CC2 * d1) / (d2 - d1), 3)     
    k2 = round((Rsplit2 - Rsplit1) / (d2 - d1), 3)
    b2 = round((Rsplit1 * d2 - Rsplit2 * d1) / (d2 - d1), 3)
    resolution = round(0.98 * (b2 - b1) / (k1 - k2), 3)
    return resolution

def get_d_at_snr_one(file_path):
    """Get resolution at SNR=1 using linear interpolation, handling non-monotonic SNR data."""
    data = []
    
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
    """Parse the file and return the d(A) value corresponding to a specific CC value."""
    with open(file_path, 'r') as file:
        lines = file.readlines()

    data = []
    for line in lines[1:]:
        parts = line.split()
        if len(parts) < 5:
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
            return round(d1 + (d2 - d1) * (target_cc - cc1) / (cc2 - cc1), 3)
    return None  # Return None if target CC is out of range

def processing(CCstar_dat_file):
    """Process a single file and extract relevant data."""
    data_info = defaultdict(dict)
    
    name_of_run = os.path.basename(CCstar_dat_file).split(".")[0].replace("_CCstar","")
    data_info[name_of_run] = {i: '' for i in indexes}
    
    print(f'Processing {name_of_run}')
    
    Rsplit_dat_file = CCstar_dat_file.replace("CCstar", "Rsplit")
    
    # Wait for Rsplit_dat_file to exist and be non-empty
    start_time = time.time()
    timeout = 60  # Timeout after 60 seconds
    while not os.path.exists(Rsplit_dat_file) and time.time() - start_time < timeout:
        time.sleep(5)
    while os.stat(Rsplit_dat_file).st_size == 0 and time.time() - start_time < timeout:
        time.sleep(5)
    
    # Similarly, wait for other necessary files
    CC_dat_file = CCstar_dat_file.replace("CCstar", "CC")
    while not os.path.exists(CC_dat_file) and time.time() - start_time < timeout:
        time.sleep(5)
    while os.stat(CC_dat_file).st_size == 0 and time.time() - start_time < timeout:
        time.sleep(5)
    
    SNR_dat_file = CCstar_dat_file.replace("CCstar", "SNR")
    while not os.path.exists(SNR_dat_file) and time.time() - start_time < timeout:
        time.sleep(5)
    while os.stat(SNR_dat_file).st_size == 0 and time.time() - start_time < timeout:
        time.sleep(5)
    
    # Process the files and collect results
    data_info[name_of_run]['CC* intersects with Rsplit at'] = f'{calculating_max_res_from_Rsplit_CCstar_dat(CCstar_dat_file, Rsplit_dat_file)}'
    data_info[name_of_run]['Resolution SNR=1'] = str(get_d_at_snr_one(SNR_dat_file))
    data_info[name_of_run]['Resolution CC>=0.3'] = str(get_d_at_cc_threshold(CC_dat_file))
    
    return data_info

if __name__ == "__main__":
    args = parse_cmdline_args()
    main_path = args.path_from
    output = args.output

    CC_dat_to_parse = []

    # Collect all CCstar.dat files to process
    for path, subdirs, files in os.walk(main_path):
        for name in files:
            if name.endswith("CCstar.dat"):
                CC_dat_to_parse.append(os.path.join(path, name))

    # Use multiprocessing to process files in parallel
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for result in executor.map(processing, CC_dat_to_parse):
            for name_of_run in result:
                data_info_all[name_of_run] = result[name_of_run]         
    
    print(len(data_info_all))           
    df = pd.DataFrame.from_dict(data_info_all)
    df.to_csv(output, sep=';')
    print(f"Results saved to {output}")
    print("Done!")