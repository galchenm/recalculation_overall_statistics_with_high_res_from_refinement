import os
import sys

USER='galchenm'

def find_script(script_name, search_path):
    """Search for a script in the given directory and its subdirectories.
    This function looks for a script with the specified name in the provided search path.
    Parameters:
        script_name (str): The name of the script to search for.
        search_path (str): The directory path where the script should be searched.
    Returns:
        str: The full path to the script if found, otherwise None.
    Raises:
        FileNotFoundError: If the script is not found in the specified path.
    Usage:
    >>> script_path = find_script('my_script.py', '/path/to/search')
    This will search for 'my_script.py' in the specified directory and return its full path.
    """
    for root, dirs, files in os.walk(search_path):
        if script_name in files:
            return os.path.join(root, script_name)
    return None

def run_partialator(hkl_input_file, highres, pg, pdb, nsh=10, suffix=''):
    """Run the partialator to compare two hkl files and generate statistics.
    This function prepares a job script to run the `compare_hkl` and `check_hkl` commands
    for the provided hkl input file, high resolution cutoff, point group, and pdb file
    with the specified number of shells. It generates output files for various statistics
    such as CCstar, Rsplit, CC, CCano, SNR, and Wilson
    statistics. The job script is submitted to the SLURM scheduler.
    It also generates a plot combining the CCstar and Rsplit statistics.
    If the hkl1 and hkl2 files already exist, it uses them directly.
    If the files do not exist, it will not run the job and will print a message.
    This function assumes that the necessary modules for `compare_hkl` and `check_hkl`
    are available in the environment, and it sets up the job script accordingly.
    It returns the path to the CCstar.dat file and the error filename to parse.
    Parameters:
        hkl_input_file (str): Path to the input hkl file.
        highres (float): High resolution cutoff for the calculations.
        pg (str): Point group for the calculations.
        pdb (str): Path to the pdb file for unit cell information.
        nsh (int): Number of shells for the calculations.
        suffix (str): Suffix to be added to the output file names.
    Returns:
        tuple: A tuple containing the path to the CCstar.dat file and the error filename to
        parse, or (None, None) if the hkl1 and hkl2 files do not exist.
    Raises:
        FileNotFoundError: If the hkl_input_file does not exist.
        Exception: If there is an error in creating or submitting the job script.
    Usage:
    >>> hkl_input_file = 'path/to/input.hkl'
    >>> highres = 1.5
    >>> pg = 'P1'
    >>> pdb = 'path/to/unit_cell.pdb'
    >>> nsh = 10    
    >>> suffix = 'offset_value'
    >>> run_partialator(hkl_input_file, highres, pg, pdb, nsh,
    suffix)
    This will create a job script in the same directory as the input hkl file,
    submit it to the SLURM scheduler, and return the path to the CCstar.dat file
    and the error filename to parse.
    If the hkl1 and hkl2 files do not exist, it will print a
    message and return (None, None).
    Note: Ensure that the necessary modules for `compare_hkl` and `check_hkl
    are available in the environment where this script is run.
    """
    
    path = os.path.dirname(os.path.abspath(hkl_input_file))
    os.chdir(path)
    print(f'We are in {os.getcwd()}')
    data = os.path.basename(hkl_input_file).split('.')[0]
    data_output_name = data if len(suffix) == 0 else f"{data}_offset_{suffix.replace('.', '_')}"

    if os.path.exists(f'{data}.hkl1') and os.path.exists(f'{data}.hkl2'):
        
        job_file = os.path.join(path, "%s.sh" % data_output_name)
        
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
            
            # Get the directory where the current script is located
            current_script_dir = os.path.dirname(os.path.abspath(__file__))

            # Script to search for
            script_name = "many_plots-upt-v2.py"

            # Search starting from the current script's directory
            script_path = find_script(script_name, current_script_dir)

            if script_path:
                command = (
                    f"python3 {script_path} -i {data_output_name}_CCstar.dat "
                    f"-x '1/d' -y 'CC*' -o {data_output_name}.png "
                    f"-add_nargs {data_output_name}_Rsplit.dat -yad 'Rsplit/%' "
                    f"-x_lim_dw 1. -x_lim_up {max_dd} -t {data_output_name} "
                    f"-legend {data_output_name} >> output.err\n"
                )
                fh.writelines(command)
            else:
                print(f"Could not find {script_name} under {current_script_dir}. Skipping execution.")

        print(f'The {job_file} is going to be submitted')    
        os.system("sbatch %s" % job_file)
        
        return "%s_CCstar.dat" % os.path.join(path, data_output_name), "%s.err" % os.path.join(path, data_output_name)
    else:
        print(f'You do not have hkl1 and/or hkl2 files for {hkl_input_file}')
        return None, None
