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
    parser.add_argument('output', type=str, help="Path to csv file with final results")
    parser.add_argument('-f','--f', type=str, help='File with blocks')
    parser.add_argument('-p','--p', type=str, help='Pattern in filename or path')
    parser.add_argument('-c', '--cell', type=str, help='Path to the folder with cell/pdb files. Check that naming is the same as in/or exactly as hkl files!')
    parser.add_argument('-a','--a', type=str, help='Additional path with Rfree/Rwork')
    parser.add_argument('-n', '--nshells', default=10, type=int,  help="Number of shells")
    parser.add_argument('-offset', '--offset', default=None, nargs='+', type=float, help="Value/s for offset parameter for resolution [Angstrom]")
    parser.add_argument('--e', action="store_true", help='Use this flag if you want to get an extended version of results')
    parser.add_argument('--r', action="store_true", help='Use this flag if you want to run the robust refinement')
    return parser.parse_args()

if __name__ == "__main__":

    args = parse_cmdline_args()
    
    if not os.path.exists(args.path_from):
        print(f"Path {args.path_from} does not exist. Please check the path and try again.")
        sys.exit(1)
        
    main_path = args.path_from
    output = args.output
    
    Rfree_Rwork_path = args.a
    
    nsh = args.nshells
    
    offsets = args.offset
    if offsets is None:
        offsets = [0.]
    else:
        offsets.append(0.)
    
    is_extended = args.e
    cell_path = args.cell
    
    is_refining = args.r
    
    pattern = args.p
    
    block_file = args.f
    if block_file:
        if not os.path.exists(block_file):
            print(f"Block file {block_file} does not exist. Please check the path and try again.")
            sys.exit(1)
        with open(block_file, 'r') as f:
            hkl_files = [os.path.join(main_path, line.strip()) for line in f if line.strip() and os.path.exists(os.path.join(main_path, line.strip()))] 
    else:
        hkl_files = glob.glob(f'{main_path}/*{line.strip()}*/**/*.hkl', recursive=True)
        
        if not hkl_files:
            print(f"No valid blocks found in {main_path}. Please check the directory.")
            sys.exit(1)
        
    
    
    hkl_files = list(filter(lambda x: (pattern in x), hkl_files)) if pattern is not None else hkl_files