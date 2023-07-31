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
indexes = ['Num. patterns/hits', 'Indexed patterns/crystals', 'CC* intersects with Rsplit at', 'Resolution', 'Rsplit (%)', 'CC1/2', 'CC*', 'CCano', 'SNR', 'Completeness (%)', 'Multiplicity','Total Measurements' ,'Unique Reflections', 'Wilson B-factor']

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
    parser.add_argument('-a','--a', type=str, help='Additional path with Rfree/Rwork')
    parser.add_argument('-n', '--nshells', default=10, type=int,  help="Number of shells")
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


def run_partialator(hkl_input_file, highres, pg, pdb, nsh):
    path = os.path.dirname(os.path.abspath(hkl_input_file))
    os.chdir(path)
    print(f'We are in {os.getcwd()}')
    data = os.path.basename(hkl_input_file).split('.')[0]
    
    if os.path.exists(f'{data}.hkl1') and os.path.exists(f'{data}.hkl2'):
        print(77)
        
        print('SBATCH PROCESS\n')
        job_file = os.path.join(path,"%s.sh" % data)
        with open(job_file, 'w+') as fh:
            fh.writelines("#!/bin/sh\n")
            fh.writelines("#SBATCH --job=%s\n" % data)
            fh.writelines("#SBATCH --partition=upex\n")
            fh.writelines("#SBATCH --time=12:00:00\n")
            fh.writelines("#SBATCH --nodes=1\n")
            fh.writelines("#SBATCH --nice=100\n")
            fh.writelines("#SBATCH --mem=500000\n")
            fh.writelines("#SBATCH --output=%s.out\n" % data)
            fh.writelines("#SBATCH --error=%s.err\n" % data)
            fh.writelines("source /etc/profile.d/modules.sh\n")
            fh.writelines("module load xray\n")


            fh.writelines("module load hdf5/1.10.5\n")
            fh.writelines("module load anaconda3/5.2\n")
            fh.writelines("module load maxwell crystfel\n")
            fh.writelines("export QT_QPA_PLATFORM=offscreen\n") 


            command = f"compare_hkl -p {pdb} -y {pg} --highres={highres} --nshells={nsh} --fom=CCstar --shell-file={data}_CCstar.dat {data}.hkl1 {data}.hkl2\n"
            fh.writelines(command)

            command = f"compare_hkl -p {pdb} -y {pg} --highres={highres} --nshells={nsh} --fom=Rsplit --shell-file={data}_Rsplit.dat {data}.hkl1 {data}.hkl2\n"
            fh.writelines(command)

            command = f"compare_hkl -p {pdb} -y {pg} --highres={highres} --nshells={nsh} --fom=CC --shell-file={data}_CC.dat {data}.hkl1 {data}.hkl2\n"
            fh.writelines(command)

            command = f"compare_hkl -p {pdb} -y {pg} --highres={highres} --nshells={nsh} --fom=CCano --shell-file={data}_CCano.dat {data}.hkl1 {data}.hkl2\n"
            fh.writelines(command)

            command = f"check_hkl -p {pdb} -y {pg} --highres={highres} --nshells={nsh} --shell-file={data}_SNR.dat {data}.hkl\n"
            fh.writelines(command)

            command = f"check_hkl -p {pdb} -y {pg} --highres={highres} --nshells={nsh} --wilson --shell-file={data}_Wilson.dat {data}.hkl\n"
            fh.writelines(command)
            
            max_dd = round(10./highres,3)

            command = f"python3 /gpfs/cfel/group/cxi/scratch/data/2020/EXFEL-2019-Schmidt-Mar-p002450/scratch/galchenm/scripts_for_work/plot_func/many_plots-upt-v2.py -i {data}_CCstar.dat -x '1/d' -y 'CC*' -o {data}.png -add_nargs {data}_Rsplit.dat -yad 'Rsplit/%' -x_lim_dw 1. -x_lim_up {max_dd} -t {data} -legend {data} >> output.err\n"
            fh.writelines(command)
            
        os.system("sbatch %s" % job_file)
        
        return "%s_CCstar.dat" % data, "%s.err" % data
    else:
        print(f'You do not have hkl1 and/or hkl2 files for {hkl_input_file}')
        return None, None

def parse_err(name_of_run, filename, CCstar_dat_file):
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
    
    shell, CCstar_shell, Rsplit_shell, CC_shell, max_shell, min_shell, SNR_shell, Completeness_shell, unique_refs_shell, multiplicity_shell = outer_shell(CCstar_dat_file)
    
    data_info[name_of_run]['Resolution'] = resolution + f' ({max_shell} - {min_shell})'
    data_info[name_of_run]['Rsplit (%)'] = Rsplit + f' ({Rsplit_shell})'
    data_info[name_of_run]['CC1/2'] =  CC + f' ({CC_shell})'
    data_info[name_of_run]['CC*'] = CCstar  + f' ({CCstar_shell})'
    data_info[name_of_run]['CCano'] = CCano
    data_info[name_of_run]['SNR'] =  snr  + f' ({SNR_shell})'
    data_info[name_of_run]['Completeness (%)'] = completeness  + f' ({Completeness_shell})'
    data_info[name_of_run]['Multiplicity'] =  multiplicity  + f' ({multiplicity_shell})'
    data_info[name_of_run]['Total Measurements'] =  total_measuremenets
    data_info[name_of_run]['Unique Reflections'] =  unique_reflections  + f' ({unique_refs_shell})'
    data_info[name_of_run]['Wilson B-factor'] = Wilson_B_factor
    print(data_info)    

def outer_shell(CCstar_dat_file):
    shell, CCstar_shell = '',''
    while not(os.path.exists(CCstar_dat_file)):
        time.sleep(25)
    while os.stat(CCstar_dat_file).st_size == 0:
        time.sleep(25)
           
    with open(CCstar_dat_file, 'r') as file:
        for line in file:
            line = re.sub(' +',' ', line.strip()).split(' ')
            CCstar_shell = line[1]
            shell = line[3]
    
    Rsplit_dat_file = CCstar_dat_file.replace("CCstar","Rsplit")
    Rsplit_shell = ''
    while not(os.path.exists(Rsplit_dat_file)):
        time.sleep(25)
    while os.stat(Rsplit_dat_file).st_size == 0:
        time.sleep(25)
  
    with open(Rsplit_dat_file, 'r') as file:
        for line in file:
            line = re.sub(' +',' ', line.strip()).split(' ')
            Rsplit_shell = line[1]
            
    CC_dat_file = CCstar_dat_file.replace("CCstar","CC")
    CC_shell = ''
    while not(os.path.exists(CC_dat_file)): 
        time.sleep(25)   
    while os.stat(CC_dat_file).st_size == 0:
        time.sleep(25)
    with open(CC_dat_file, 'r') as file:
        for line in file:
            line = re.sub(' +',' ', line.strip()).split(' ')
            CC_shell = line[1]
    
    max_shell, min_shell, SNR_shell, Completeness_shell, unique_refs_shell, multiplicity_shell = '','', '', '', '', ''
    SNR_dat_file = CCstar_dat_file.replace('CCstar', 'SNR')
    
    while not(os.path.exists(SNR_dat_file)):
        time.sleep(25)
    while os.stat(SNR_dat_file).st_size == 0:
        time.sleep(25)
    
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

def processing(CCstar_dat_file, error_filename_to_parse):

    name_of_run = os.path.basename(CCstar_dat_file).split(".")[0].replace("_CCstar","")
    data_info[name_of_run] = {i:'' for i in indexes}
    
    Rsplit_dat_file = CCstar_dat_file.replace("CCstar","Rsplit")
    
    data_info[name_of_run]['CC* intersects with Rsplit at'] = f'{calculating_max_res_from_Rsplit_CCstar_dat(CCstar_dat_file, Rsplit_dat_file)}'
    
    stream = CCstar_dat_file.replace("_CCstar.dat",".stream")
    chunks, hits, indexed_patterns, indexed = parsing_stream(stream)
        
    data_info[name_of_run]['Num. patterns/hits'] = str(chunks)+"/"+str(hits)
    data_info[name_of_run]['Indexed patterns/crystals'] = str(indexed_patterns)+"/"+str(indexed)
        
    parse_err(name_of_run, error_filename_to_parse, CCstar_dat_file)


def parsing_phenix_log_file(phenix_log_file):
    file = open(phenix_log_file, 'r')
    content = file.read()
    file.close()
    p = re.compile("Final R-work = ([\+\-\d\.]*), R-free = ([\+\-\d\.]*)")
    result = p.findall(content)
    
    if len(result) == 0:
        return None, None, None
    Rwork, Rfree = result[0]
    
    #Resolution range: 56.3397 1.70002
    p = re.compile("Resolution range: ([\+\-\d\.]*) ([\+\-\d\.]*)")
    result = p.findall(content)
    if len(result) == 0:
        return float(Rwork), float(Rfree), None
    resolution_cut_off = min(list(map(lambda x: float(x), result[0])))
    
    return Rwork, Rfree, resolution_cut_off

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


if __name__ == "__main__":

    args = parse_cmdline_args()
    main_path = args.path_from
    output = args.output
    Rfree_Rwork_path = args.a
    nsh = args.nshells
    if args.f is not None:
        with open(args.f, 'r') as file:
            for line in file:
                phenix_log_file = glob.glob(f'{Rfree_Rwork_path}/*{line.strip()}*/**/*.log', recursive=True)
                phenix_log_file = max(phenix_log_file, key=os.path.getctime)
                
                phenix_log_file = list(filter(lambda x: (args.p in x), phenix_log_file)) if args.p is not None else phenix_log_file
                Rwork, Rfree, resolution_cut_off = parsing_phenix_log_file(phenix_log_file)
                if resolution_cut_off is not None:
                    pattern_in_hkl_input_filename = os.path.basename(os.path.dirname(phenix_log_file))
                    
                    hkl_input_file = glob.glob(f'{main_path}/*{pattern_in_hkl_input_filename}*.hkl', recursive=True)
                    if len(hkl_input_file) == 0:
                        print(f'In {main_path} there is no hkl file associated with {phenix_log_file}')
                        continue
                    hkl_input_file = hkl_input_file[0]
                    pdb = get_UC(hkl_input_file)
                    pg = get_pg(hkl_input_file)
                    
                    CCstar_dat_file, error_filename_to_parse = run_partialator(hkl_input_file, resolution_cut_off, pg, pdb, nsh)
                    
                    while CCstar_dat_file is not None and error_filename_to_parse is not None and not(os.path.exists(CCstar_dat_file)) and not(os.path.exists(error_filename_to_parse)):
                        time.sleep(5.0)
                    
                    if CCstar_dat_file is not None and error_filename_to_parse is not None:
                        processing(CCstar_dat_file, error_filename_to_parse)                
    else:
        for path, dirs, all_files in os.walk(Rfree_Rwork_path):
            
            phenix_log_file = glob.glob(os.path.join(path, '*.log'), recursive=True)
            if len(phenix_log_file) == 0:
                continue
            phenix_log_file = max(phenix_log_file, key=os.path.getctime)
            
            phenix_log_file = list(filter(lambda x: (args.p in x), phenix_log_file)) if args.p is not None else phenix_log_file
            Rwork, Rfree, resolution_cut_off = parsing_phenix_log_file(phenix_log_file)
            if resolution_cut_off is not None:
                pattern_in_hkl_input_filename = os.path.basename(os.path.dirname(phenix_log_file))
                
                hkl_input_file = glob.glob(f'{main_path}/*{pattern_in_hkl_input_filename}*.hkl', recursive=True)
                if len(hkl_input_file) == 0:
                    print(f'In {main_path} there is no hkl file associated with {phenix_log_file}')
                    continue
                hkl_input_file = hkl_input_file[0]
                pdb = get_UC(hkl_input_file)
                pg = get_pg(hkl_input_file)
                
                CCstar_dat_file, error_filename_to_parse = run_partialator(hkl_input_file, resolution_cut_off, pg, pdb, nsh)
                
                while CCstar_dat_file is not None and error_filename_to_parse is not None and not(os.path.exists(CCstar_dat_file)) and not(os.path.exists(error_filename_to_parse)):
                    time.sleep(5.0)
                
                if CCstar_dat_file is not None and error_filename_to_parse is not None:
                    processing(CCstar_dat_file, error_filename_to_parse)
                
                
    df = pd.DataFrame.from_dict(data_info)
    df.to_csv(output, sep=';')