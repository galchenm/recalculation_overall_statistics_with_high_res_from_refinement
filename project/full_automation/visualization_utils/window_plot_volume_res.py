import math
import re
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd
import shlex

import subprocess
import glob
import os
import sys


#

scan_types = [
                "left-to-right top-down",  # number_of_pattern = index_line * NUM_PORES_PER_LINE + index_pore
                "right-to-left top-down",  # number_of_pattern = (index_line + 1) * NUM_PORES_PER_LINE - index_pore - 1
                "snake start top-left",    # number_of_pattern = (index_line * NUM_PORES_PER_LINE + index_pore) if index_line%2 == 0 else ((index_line + 1) * NUM_PORES_PER_LINE - index_pore - 1)
                "snake start top-right",   # number_of_pattern = ((index_line + 1) * NUM_PORES_PER_LINE - index_pore - 1) if index_line%2 == 0 else (index_line * NUM_PORES_PER_LINE + index_pore)
                "left-to-right bottom-up", # number_of_pattern = TOTAL_NUMBER_OF_PORES - ((index_line + 1) * NUM_PORES_PER_LINE - index_pore - 1)
                "right-to-left bottom-up", # number_of_pattern = TOTAL_NUMBER_OF_PORES - (index_line * NUM_PORES_PER_LINE + index_pore) 
                "snake start bottom-left", # number_of_pattern = TOTAL_NUMBER_OF_PORES - ((index_line + 1) * NUM_PORES_PER_LINE - index_pore - 1) if index_line%2 == 0 else (TOTAL_NUMBER_OF_PORES - (index_line * NUM_PORES_PER_LINE + index_pore))
                "snake start bottom-right" # number_of_pattern = (TOTAL_NUMBER_OF_PORES - (index_line * NUM_PORES_PER_LINE + index_pore)) if index_line%2 == 0 else (TOTAL_NUMBER_OF_PORES - ((index_line + 1) * NUM_PORES_PER_LINE - index_pore - 1))
             ]



def reading_point_info_file(points_info_filename):
    if not(os.path.exists(points_info_filename)):
        return None, None, None, None, None
        
    WINDOW_NUM = 0
    NUM_PORES_PER_LINE = 0
    NUM_LINES_PER_WINDOW = 0
    with open(points_info_filename) as points_info:
        points_info_stream = points_info.read()
    
    scan_type = None
    with open(os.path.join(os.path.dirname(points_info_filename),'info.txt')) as txt_info:
        for line in txt_info:
            if line.startswith("scan type:"):
                scan_type = line.split(':')[-1].strip()
                break
    
    p = re.compile("WINDOWNUM=([\+\-\d\.]*)")
    WINDOW_NUM = p.findall(points_info_stream)
    WINDOW_NUM = list(map(int, WINDOW_NUM))[0] if len(WINDOW_NUM) > 0 else None
    
    
    p = re.compile("NUM_PORES_PER_LINE=([\+\-\d\.]*)")
    NUM_PORES_PER_LINE = p.findall(points_info_stream)
    NUM_PORES_PER_LINE = list(map(int, NUM_PORES_PER_LINE))[0] if len(NUM_PORES_PER_LINE) > 0 else None

    p = re.compile("NUM_LINES_PER_WINDOW=([\+\-\d\.]*)")
    NUM_LINES_PER_WINDOW = p.findall(points_info_stream)
    NUM_LINES_PER_WINDOW = list(map(int, NUM_LINES_PER_WINDOW))[0] if len(NUM_LINES_PER_WINDOW) > 0 else None

    p = re.compile("LINE_IDX=([\+\-\d\.]*)")
    LINES_IDX = p.findall(points_info_stream)
    LINES_IDX = list(map(int, LINES_IDX)) if len(LINES_IDX) > 0 else None
    
    return scan_type, WINDOW_NUM, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW, LINES_IDX


def calculate_indices_left_to_right_bottom_up(number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW):
    TOTAL_NUMBER = NUM_LINES_PER_WINDOW * NUM_PORES_PER_LINE - 1
    for index_line in range(NUM_LINES_PER_WINDOW):
        for index_pore in range(NUM_PORES_PER_LINE):
            calculated_number_of_pattern = TOTAL_NUMBER - ((index_line + 1) * NUM_PORES_PER_LINE - index_pore - 1)
            if calculated_number_of_pattern == number_of_pattern:
                return index_line, index_pore
    return None, None

def calculate_indices_right_to_left_bottom_up(number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW):
    TOTAL_NUMBER = NUM_LINES_PER_WINDOW * NUM_PORES_PER_LINE - 1
    return (TOTAL_NUMBER_OF_PORES - number_of_pattern)//NUM_PORES_PER_LINE, (TOTAL_NUMBER_OF_PORES - number_of_pattern)%NUM_PORES_PER_LINE

def calculate_indices_left_to_right_top_down(number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW):
    return number_of_pattern//NUM_PORES_PER_LINE, number_of_pattern%NUM_PORES_PER_LINE

def calculate_indices_right_to_left_top_down(number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW):
    for index_line in range(NUM_LINES_PER_WINDOW):
        for index_pore in range(NUM_PORES_PER_LINE):
            calculated_number_of_pattern = (index_line + 1) * NUM_PORES_PER_LINE - index_pore - 1
            if calculated_number_of_pattern == number_of_pattern:
                return index_line, index_pore
    return None, None

def calculate_indices_right_to_left_top_down(number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW):
    for index_line in range(NUM_LINES_PER_WINDOW):
        for index_pore in range(NUM_PORES_PER_LINE):
            calculated_number_of_pattern = (index_line + 1) * NUM_PORES_PER_LINE - index_pore - 1
            if calculated_number_of_pattern == number_of_pattern:
                return index_line, index_pore
    return None, None

def calculate_indices_snake_start_top_left(number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW):
    for index_line in range(NUM_LINES_PER_WINDOW):
        for index_pore in range(NUM_PORES_PER_LINE):
            calculated_number_of_pattern = (index_line * NUM_PORES_PER_LINE + index_pore) if index_line%2 == 0 else ((index_line + 1) * NUM_PORES_PER_LINE - index_pore - 1)
            if calculated_number_of_pattern == number_of_pattern:
                return index_line, index_pore
    return None, None

def calculate_indices_snake_start_bottom_left(number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW):
    TOTAL_NUMBER = NUM_LINES_PER_WINDOW * NUM_PORES_PER_LINE - 1
    for index_line in range(NUM_LINES_PER_WINDOW):
        for index_pore in range(NUM_PORES_PER_LINE):
            calculated_number_of_pattern = TOTAL_NUMBER - ((index_line + 1) * NUM_PORES_PER_LINE - index_pore - 1) if index_line%2 == 0 else (TOTAL_NUMBER - (index_line * NUM_PORES_PER_LINE + index_pore))
            if calculated_number_of_pattern == number_of_pattern:
                return index_line, index_pore
    return None, None

def calculate_indices_snake_start_top_right(number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW):
    for index_line in range(NUM_LINES_PER_WINDOW):
        for index_pore in range(NUM_PORES_PER_LINE):
            calculated_number_of_pattern = ((index_line + 1) * NUM_PORES_PER_LINE - index_pore - 1) if index_line%2 == 0 else (index_line * NUM_PORES_PER_LINE + index_pore)
            if calculated_number_of_pattern == number_of_pattern:
                return index_line, index_pore
    return None, None

def calculate_indices_snake_start_bottom_right(number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW):
    for index_line in range(NUM_LINES_PER_WINDOW):
        for index_pore in range(NUM_PORES_PER_LINE):
            calculated_number_of_pattern = (TOTAL_NUMBER - (index_line * NUM_PORES_PER_LINE + index_pore)) if index_line%2 == 0 else (TOTAL_NUMBER - ((index_line + 1) * NUM_PORES_PER_LINE - index_pore - 1))
            if calculated_number_of_pattern == number_of_pattern:
                return index_line, index_pore
    return None, None

def get_id_line_id_pores_from_filename(raw_filename, extension='cbf'):
    path_dir = os.path.dirname(raw_filename)
    points_info_filename = [os.path.join(path_dir, i) for i in os.listdir(path_dir) if i == 'pointsinfo.txt']
    
    if len(points_info_filename)==0:
        return None
    
    points_info_filename = points_info_filename[0]

    scan_type, WINDOW_NUM, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW, LINES_IDX = reading_point_info_file(points_info_filename)
    if scan_type is None:
        return None
    if WINDOW_NUM is None:
        return None
    d_scan_types_for_id_line_id_pore = {
            "left-to-right top-down": (lambda number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW:(calculate_indices_left_to_right_top_down(number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW))), 
            "right-to-left top-down": (lambda number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW:(calculate_indices_right_to_left_top_down(number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW))),
            "snake start top-left": (lambda number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW: (calculate_indices_snake_start_top_left(number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW))), 
            "snake start top-right": (lambda number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW: (calculate_indices_snake_start_top_right(number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW))),
            "left-to-right bottom-up": (lambda number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW: (calculate_indices_left_to_right_bottom_up(number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW))), 
            "right-to-left bottom-up": (lambda number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW: (calculate_indices_right_to_left_bottom_up(number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW))),
            "snake start bottom-left": (lambda number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW: (calculate_indices_snake_start_bottom_left(number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW))), 
            "snake start bottom-right": (lambda number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW: (calculate_indices_snake_start_bottom_right(number_of_pattern, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW))),
            }
    
    pattern_number = int(re.search(fr'\d+.{extension}', os.path.basename(raw_filename)).group().split('.')[0])
    
    return  d_scan_types_for_id_line_id_pore[scan_type](pattern_number, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW)   

def calculate_number_of_pattern_left_to_right_bottom_up(index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW):
    TOTAL_NUMBER = NUM_LINES_PER_WINDOW * NUM_PORES_PER_LINE - 1
    return TOTAL_NUMBER - ((index_line + 1) * NUM_PORES_PER_LINE - index_pore - 1)


def calculate_number_of_pattern_right_to_left_bottom_up(index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW):
    TOTAL_NUMBER = NUM_LINES_PER_WINDOW * NUM_PORES_PER_LINE - 1
    return TOTAL_NUMBER - (index_line * NUM_PORES_PER_LINE + index_pore)

def calculate_number_of_pattern_left_to_right_top_down(index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW):
    return index_line * NUM_PORES_PER_LINE + index_pore

def calculate_number_of_pattern_right_to_left_top_down(index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW):
    return (index_line + 1) * NUM_PORES_PER_LINE - index_pore - 1


def calculate_number_of_pattern_right_to_left_top_down(index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW):
    return (index_line + 1) * NUM_PORES_PER_LINE - index_pore - 1


def calculate_number_of_pattern_snake_start_top_left(index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW):
    return (index_line * NUM_PORES_PER_LINE + index_pore) if index_line%2 == 0 else ((index_line + 1) * NUM_PORES_PER_LINE - index_pore - 1)

def calculate_number_of_pattern_snake_start_bottom_left(index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW):
    TOTAL_NUMBER = NUM_LINES_PER_WINDOW * NUM_PORES_PER_LINE - 1
    return TOTAL_NUMBER - ((index_line + 1) * NUM_PORES_PER_LINE - index_pore - 1) if index_line%2 == 0 else (TOTAL_NUMBER - (index_line * NUM_PORES_PER_LINE + index_pore))


def calculate_number_of_pattern_snake_start_top_right(index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW):
    return ((index_line + 1) * NUM_PORES_PER_LINE - index_pore - 1) if index_line%2 == 0 else (index_line * NUM_PORES_PER_LINE + index_pore)

def calculate_number_of_pattern_snake_start_bottom_right(index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW):
    TOTAL_NUMBER = NUM_LINES_PER_WINDOW * NUM_PORES_PER_LINE - 1
    return TOTAL_NUMBER - (index_line * NUM_PORES_PER_LINE + index_pore) if index_line%2 == 0 else (TOTAL_NUMBER - ((index_line + 1) * NUM_PORES_PER_LINE - index_pore - 1))

def obtain_coordinates_for_current_folder(path_dir, extension='cbf'):
    points_info_filename = [os.path.join(path_dir, i) for i in os.listdir(path_dir) if i == 'pointsinfo.txt']
    
    if len(points_info_filename)==0:
        return None
    
    points_info_filename = points_info_filename[0]

    scan_type, WINDOW_NUM, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW, LINES_IDX = reading_point_info_file(points_info_filename)
    if scan_type is None:
        return None
    if WINDOW_NUM is None:
        return None
    d_scan_types_for_calculate_number_of_pattern = {
            "left-to-right top-down": (lambda index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW:calculate_number_of_pattern_left_to_right_top_down(index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW)), 
            "right-to-left top-down": (lambda index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW:calculate_number_of_pattern_right_to_left_top_down(index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW)),
            "snake start top-left": (lambda index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW: calculate_number_of_pattern_snake_start_top_left(index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW)), 
            "snake start top-right": (lambda index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW: calculate_number_of_pattern_snake_start_top_right(index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW)),
            "left-to-right bottom-up": (lambda index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW: calculate_number_of_pattern_left_to_right_bottom_up(index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW)), 
            "right-to-left bottom-up": (lambda index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW: calculate_number_of_pattern_right_to_left_bottom_up(index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW)),
            "snake start bottom-left": (lambda index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW: calculate_number_of_pattern_snake_start_bottom_left(index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW)), 
            "snake start bottom-right": (lambda index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW: calculate_number_of_pattern_snake_start_bottom_right(index_line, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW)),
            }
      
    raw_filename = glob.glob(f'{path_dir}/*{extension}')[0]
    
    template_file_name = re.sub(r'\d+\.', lambda m: "?" *(len(m.group())-1) +'.', raw_filename)
    
    
    dict_line_info = defaultdict(list)
    
    min_LINES_IDX = min(LINES_IDX)
    
    
    for index_line in LINES_IDX:
        dict_line_info[index_line] = {}
        for index_pore in range(NUM_PORES_PER_LINE):
            
            number_of_pattern = d_scan_types_for_calculate_number_of_pattern[scan_type](index_line - min_LINES_IDX, index_pore, NUM_PORES_PER_LINE, NUM_LINES_PER_WINDOW)
            
            filename = re.sub('\?\?\?\?\?', "%05d" % number_of_pattern, template_file_name)
            
            if os.path.exists(filename):
                dict_line_info[index_line][index_pore] = filename
    return dict_line_info

def generate_value_pairs(dictionary):
    value_pairs = {}
    for key1, inner_dict in dictionary.items():
        for key2, value in inner_dict.items():
            value_pairs[value] = (key1, key2)         
    return value_pairs
    
def reading_streamfile(stream_filename):
    path_dir = None
    try:
        command = f'grep -e "Image filename" {stream_filename}'
        result = subprocess.check_output(shlex.split(command)).decode('utf-8').strip().split('\n')[0]
        path_dir = os.path.dirname(result.split(":")[-1].strip())
    except subprocess.CalledProcessError:
        pass
    
    if path_dir is None:
        return None, None, None

    dict_line_info = obtain_coordinates_for_current_folder(path_dir, extension='cbf')
    
    if dict_line_info is None:
        return None, None, None
    pairs = generate_value_pairs(dict_line_info)
    
    coordinates = list(pairs.values())
    max_x = max(map(lambda x: x[0], coordinates))
    max_y = max(map(lambda x: x[1], coordinates))
    
    volume = np.zeros((max_x+1, max_y+1))
    res = np.zeros((max_x+1, max_y+1))
    num_crystals = np.zeros((max_x+1, max_y+1), dtype=int)
    cell = {}  # Use a dictionary to store crystal information for each image filename
    resolution_limit = defaultdict(list)
    
    with open(stream_filename, 'r') as stream:
        for line in stream:
            if line.startswith('----- Begin chunk -----'):
                cell = {}
                indexed = False
            elif line.startswith('Image filename:'):
                image_filename = line.split(': ')[-1].strip()
                cell[image_filename] = []
            elif line.startswith('--- End crystal'):
                if indexed:
                    volume[pairs[image_filename]] += 1 / np.linalg.det(cell[image_filename])
                    res[pairs[image_filename]] += resolution_limit[image_filename]
                    num_crystals[pairs[image_filename]] += 1  # Count the number of crystals for this file
                    #print(f'The number of crystals for {os.path.basename(image_filename)} is {num_crystals[pairs[image_filename]]}', end="\r")
                cell[image_filename] = []
                indexed = False
            elif line[:5] in ('astar', 'bstar', 'cstar'):
                cell[image_filename].append([float(i) for i in line.split()[2:5]])
            elif line.startswith('diffraction_resolution_limit'):
                resolution_limit[image_filename] = float(line.split()[-2])
                indexed = True

    # Average the volume and resolution by the number of crystals found for each file
    num_crystals_total = num_crystals.sum()
    num_crystals[num_crystals == 0] = 1  # Avoid division by zero
    volume /= num_crystals
    res /= num_crystals
    
    return num_crystals_total, volume, res

    
def plot_over_resolution(data, picture_filename='foo.png', block_size = 4, label='UC volume, nm^3', title='Volume'):
    # Determine the block dimensions (block_size x block_size)
    
    block_rows = data.shape[0] // block_size
    block_cols = data.shape[1] // block_size

    # Reshape the 'data' array into (block_size x block_size) blocks and calculate the mean of each block
    block_data_means = data[:block_rows * block_size, :block_cols * block_size].reshape(block_rows, block_size, block_cols, block_size).mean(axis=(1, 3))

    # Determine the min and max values of 'block_data_means' for setting color levels
    vmin = np.nanmin(block_data_means)
    vmax = np.nanmax(block_data_means)

    # Create the plot with a size of 16x16
    plt.figure(figsize=(16, 16))

    # Plot the spatial distribution of volume
    plt.imshow(block_data_means, cmap='viridis', aspect='auto', origin='lower', vmin=vmin, vmax=vmax)
    plt.colorbar(label=label)

    plt.title(f'Spatial Distribution of {title} (Averaged in {block_size}x{block_size} Bins)')
    
    plt.savefig(picture_filename)
    return picture_filename