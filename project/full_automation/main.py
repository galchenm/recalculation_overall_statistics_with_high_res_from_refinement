#!/usr/bin/env python3
# coding: utf8
# Written by Galchenkova M.

import os
import sys
import csv
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
from run_processing_utils.preparation_for_statistics_calculation import prep_for_calculating_overall_statistics
from partialator_utils.partialator_execution import run_partialator


os.nice(0)
indexes = ['Num. patterns/hits', 'Indexed patterns/crystals', 'Resolution', 'Rsplit(%)', 'CC1/2', 'CC*', 'CCano', 'SNR', 'Completeness(%)', 'Multiplicity','Total Measurements' ,'Unique Reflections', 'Wilson B-factor', 'Resolution SNR=1', 'Resolution CC>=0.3', 'a,b,c,alpha,betta,gamma']

USER='galchenm'
SLEEP_TIME = 10
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
    parser.add_argument('-s', '--single', action='store_true', help='Process a single .hkl file (use with -p or specify exact path)')
    parser.add_argument('--online', action='store_true', help='Monitor folder and process each new .hkl file as it appears')
    parser.add_argument('--offline', action='store_true', help='Process all .hkl files at once (default behavior)')
    return parser.parse_args()


def discover_hkl_files():
    hkl_files = []
    if block_file:
        with open(block_file, 'r') as f:
            hkl_files = [os.path.join(main_path, line.strip()) for line in f if line.strip()]
    else:
        hkl_files = glob.glob(f'{main_path}/**/*.hkl', recursive=True)
    if pattern:
        hkl_files = list(filter(lambda x: pattern in x, hkl_files))
    return [f for f in hkl_files if os.path.exists(f)]

def wait_for_jobs_to_finish():
    while True:
        pending_command = f'squeue -u {USER} -t pending'
        number_of_pending = subprocess.check_output(shlex.split(pending_command)).decode().strip().split('\n')
        if len(number_of_pending) <= 1:
            break
        time.sleep(SLEEP_TIME)


def append_to_csv(data_dict, output_file):
    file_exists = os.path.isfile(output_file)
    
    with open(output_file, 'a', newline='') as csvfile:
        fieldnames = ['Run'] + list(next(iter(data_dict.values())).keys())
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        if not file_exists:
            writer.writeheader()
        
        for run_name, data in data_dict.items():
            row = {'Run': run_name}
            row.update(data)
            writer.writerow(row)

def write_to_csv(all_data_info, output_file):
    all_flattened = []
    for run_name, data in all_data_info.items():
        row = {'Run': run_name}
        row.update(data)
        all_flattened.append(row)
    df = pd.DataFrame(all_flattened)
    df.to_csv(output_file, index=False)
    print(f"Results written to {output_file}")


if __name__ == "__main__":
    args = parse_cmdline_args()

    if not os.path.exists(args.path_from):
        print(f"Path {args.path_from} does not exist. Please check the path and try again.")
        sys.exit(1)

    main_path = args.path_from
    output = args.output
    pattern = args.p
    Rfree_Rwork_path = args.a
    nsh = args.nshells
    offsets = args.offset or [0.0]
    if 0.0 not in offsets:
        offsets.append(0.0)
    is_extended = args.e
    cell_path = args.cell
    is_refining = args.r
    block_file = args.f
    is_single = args.single
    is_online = args.online
    is_offline = args.offline or not is_online  # default to offline

    data_info_all = defaultdict(dict)

    if is_single:
        hkl_files = discover_hkl_files()
        if not hkl_files:
            print("No .hkl files found.")
            sys.exit(1)
        hkl_file = hkl_files[0]
        submission_data = {}

        def prep_wrapper(offset):
            run_name = hkl_file.split('.')[0] + f'_{str(offset).replace(".", "p")}'
            try:
                results = prep_for_calculating_overall_statistics(
                    hkl_file, offset, cell_path, Rfree_Rwork_path, nsh)
                return (run_name, {
                    'hkl_file': hkl_file,
                    'data': {
                        'CCstar_dat_file': results[0],
                        'error_file': results[1],
                        'Rwork': results[2],
                        'Rfree': results[3],
                        'resolution_cut_off_high': results[4],
                        'resolution_low': results[5]
                    }
                })
            except Exception as e:
                print(f"[ERROR] prep failed for offset {offset}: {e}")
                return None

        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = [executor.submit(prep_wrapper, offset) for offset in offsets]
            for future in concurrent.futures.as_completed(futures):
                result = future.result()
                if result:
                    run_name, run_info = result
                    submission_data[run_name] = run_info

        wait_for_jobs_to_finish()

        for run_name, item in submission_data.items():
            run_data = processing_statistics_for_run(
                run_name, {run_name: item['data']},
                item['hkl_file'], main_path, is_extended, cell_path, is_refining)
            data_info_all.update(run_data)

        write_to_csv(data_info_all, output)

    elif is_online:
        seen = set()
        print("Watching for new .hkl files...")
        try:
            while True:
                hkl_files = discover_hkl_files()
                new_files = [f for f in hkl_files if f not in seen]
                for hkl_file in new_files:
                    seen.add(hkl_file)
                    for offset in offsets:
                        run_name = hkl_file.split('.')[0] + f'_{str(offset).replace(".", "p")}'
                        CCstar_dat_file, error_file, Rwork, Rfree, resolution_cut_off_high, resolution_low = prep_for_calculating_overall_statistics(
                            hkl_file, offset, cell_path, Rfree_Rwork_path, nsh)
                        data_info = {
                            'CCstar_dat_file': CCstar_dat_file,
                            'error_file': error_file,
                            'Rwork': Rwork,
                            'Rfree': Rfree,
                            'resolution_cut_off_high': resolution_cut_off_high,
                            'resolution_low': resolution_low
                        }
                        wait_for_jobs_to_finish()
                        run_data = processing_statistics_for_run(run_name, {run_name: data_info},
                                                                 hkl_file, main_path,
                                                                 is_extended, cell_path, is_refining)
                        append_to_csv(run_data, output)
                        data_info_all.update(run_data)

                time.sleep(SLEEP_TIME)
        except KeyboardInterrupt:
            print("Exiting online mode.")

    elif is_offline:
        hkl_files = discover_hkl_files()
        if not hkl_files:
            print("No .hkl files found.")
            sys.exit(1)

        submission_data = {}

        def prep_wrapper_offline(args):
            hkl_file, offset = args
            run_name = hkl_file.split('.')[0] + f'_{str(offset).replace(".", "p")}'
            try:
                results = prep_for_calculating_overall_statistics(
                    hkl_file, offset, cell_path, Rfree_Rwork_path, nsh)
                return (run_name, {
                    'hkl_file': hkl_file,
                    'data': {
                        'CCstar_dat_file': results[0],
                        'error_file': results[1],
                        'Rwork': results[2],
                        'Rfree': results[3],
                        'resolution_cut_off_high': results[4],
                        'resolution_low': results[5]
                    }
                })
            except Exception as e:
                print(f"[ERROR] prep failed for {hkl_file} offset {offset}: {e}")
                return None

        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = [executor.submit(prep_wrapper_offline, (hkl_file, offset)) for hkl_file in hkl_files for offset in offsets]
            for future in concurrent.futures.as_completed(futures):
                result = future.result()
                if result:
                    run_name, run_info = result
                    submission_data[run_name] = run_info

        wait_for_jobs_to_finish()

        for run_name, item in submission_data.items():
            run_data = processing_statistics_for_run(
                run_name, {run_name: item['data']},
                item['hkl_file'], main_path, is_extended, cell_path, is_refining)
            data_info_all.update(run_data)
        data_info_all = data_info_all.T
        write_to_csv(data_info_all, output)
