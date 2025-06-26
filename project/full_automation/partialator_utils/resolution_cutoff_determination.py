import os
import sys
import numpy as np
import pandas as pd

x_arg_name = 'd'
y_arg_name = 'CC*'
y_arg_name2 = 'Rsplit/%'

def get_d_at_snr_one(file_path):
    """
    Get resolution at SNR=1 using 
    linear interpolation, handling 
    non-monotonic SNR data.
    
    :param file_path: Path to the SNR.dat file.
    :return: The d(A) value at SNR=1, or None if not
            found or data is insufficient."""
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

def get_xy(file_name, x_arg_name, y_arg_name):
    """
    Extract the x and y values from the given file.
    :param file_name: Path to the file containing x and y data.
    :param x_arg_name: Name of the x argument (e.g., 'd(A)').
    :param y_arg_name: Name of the y argument (e.g., 'CC').
    :return: Two numpy arrays containing the x and y values.
    This function reads a file, extracts the x and y values based on the provided argument names,
    and returns them as numpy arrays. It handles cases where the values might be NaN or zero.
    It also constructs a DataFrame from the extracted values, ensuring that only valid entries are retained.
    
    """
    
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
    
    """Calculate the maximum resolution from the provided CCstar and Rsplit data files.
    This function reads the CCstar and Rsplit data files, extracts the relevant x and y values,
    and computes the maximum resolution based on the relationship between CCstar and Rsplit.
    :param CCstar_dat_file: Path to the CCstar data file.
    :param Rsplit_dat_file: Path to the Rsplit data file.
    :return: The maximum resolution as a float value.
    This function assumes that the CCstar data file contains columns for 'd(A)' and 'CC*',
    and the Rsplit data file contains columns for 'd(A)' and 'Rsplit'.
    It uses linear interpolation to find the resolution at which the CCstar value is equal to the Rsplit value.
    If the CCstar value is not found or the data is insufficient, it returns -1000.
    """
    
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
