import os
import sys
import glob
import shlex
import subprocess
import time
import csv
from typing import List, Dict
from statistics import mean
import numpy as np
from pathlib import Path

LIMIT_FOR_RESERVED_NODES = 1000
SLEEP_TIME = 15  # seconds

csv_header = ["POSITION", "Unit Cell a","Unit Cell b", "Unit Cell c", "Angle alpha", "Angle beta", "Angle gamma", "Space Group", "CC1/2", "R-meas", "R-pim", "R-merge"]

def are_the_reserved_nodes_overloaded(node_list: str) -> bool:
    """
    Check if the reserved nodes are overloaded by counting running jobs.
    Returns True if the number of jobs exceeds LIMIT_FOR_RESERVED_NODES.
    """
    try:
        jobs_cmd = f'squeue -w {node_list}'
        lines = subprocess.check_output(shlex.split(jobs_cmd)).decode().splitlines()
        # squeue prints a header line; subtract it
        job_count = max(0, len(lines) - 1)
    except subprocess.CalledProcessError:
        job_count = 0
    return job_count > LIMIT_FOR_RESERVED_NODES

def read_gxparm_xds(gxparm_xds_filepath: str) -> Tuple[Optional[int], Optional[Tuple[float, float, float, float, float, float]]]:
    """
    Reads a GXPARM.XDS file and extracts space group number and unit cell parameters.
    Returns (space_group, (a,b,c,alpha,beta,gamma)) or (None, None) if not found.
    """
    with open(gxparm_xds_filepath, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 7:
                try:
                    sg = int(parts[0])
                    a, b, c, alpha, beta, gamma = map(float, parts[1:])
                    # Quick sanity check: positive cell edges
                    if a > 0 and b > 0 and c > 0:
                        return sg, (a, b, c, alpha, beta, gamma)
                except ValueError:
                    continue
    return None, None

def read_correct_lp_file(correct_lp_path: str):
    # Extract CC/2, R-merge, R-meas, R-pim from CORRECT.LP
    CChalf = Rmerge = Rmeas = Rpim = ""
    if os.path.isfile(correct_lp_path):
        with open(correct_lp_path, 'r') as lp:
            lp_lines = lp.readlines()
            for i, line in enumerate(lp_lines):
                if "WILSON STATISTICS" in line:
                    # Using relative line offsets based on Bash script
                    try:
                        CChalf = lp_lines[i-13].split()[10]
                        Rmerge = lp_lines[i-13].split()[4]
                        Rmeas = lp_lines[i-13].split()[12]
                        Rpim = lp_lines[i-21].split()[9]
                    except Exception:
                        print("Error reading CORRECT.LP")
                        return None, None, None, None
                    break
    return CChalf, Rmerge, Rmeas, Rpim

from typing import List, Dict
from statistics import mean

def group_humidity_plateaus(humidity: List[float], temperatures: List[float], positions: List[int], tolerance: float = 0.5) -> Dict[float, Dict]:
    groups: Dict[float, Dict] = {}
    if not humidity:
        return groups
    
    tolerance /= 100
    
    current_positions = [positions[0]]
    current_humidities = [humidity[0]]
    current_temperatures = [temperatures[0]]

    for h, t, p in zip(humidity[1:], temperatures[1:], positions[1:]):
        # Check if this humidity is within tolerance of the current group
        if abs(h - mean(current_humidities)) <= tolerance:
            current_positions.append(p)
            current_humidities.append(h)
            current_temperatures.append(t)
        else:
            # Finalize the current plateau
            avg_hum = round(mean(current_humidities), 2)
            groups[avg_hum] = {
                "avg_temperature": mean(current_temperatures),
                "positions": current_positions
            }
            # Start a new plateau
            current_positions = [p]
            current_humidities = [h]
            current_temperatures = [t]

    # Finalize the last plateau
    avg_hum = round(mean(current_humidities), 2)
    groups[avg_hum] = {
        "avg_temperature": mean(current_temperatures),
        "positions": current_positions
    }

    return groups

def xscale_inp_generating(
    input_path: str,
    current_data_processing_folder: str,
    space_group_number: str = "!SPACE_GROUP_NUMBER",
    unit_cell_constants: str = "!UNIT_CELL_CONSTANTS",
    reference_dataset: str = "!REFERENCE_DATA_SET",
    xds_ascii_files: List[str] = [],
    output_hkl: str = "output.hkl",
) -> Path:
    """
    Generate XSCALE.INP in current_data_processing_folder.
    The parameters must already be full XSCALE lines, e.g.:
    space_group_number="SPACE_GROUP_NUMBER=96"
    unit_cell_constants="UNIT_CELL_CONSTANTS=a b c alpha beta gamma"
    reference_dataset="REFERENCE_DATA_SET=path/to/reference.mtz" or "!REFERENCE_DATA_SET"
    """
    xscale_inp = Path(current_data_processing_folder) / "XSCALE.INP"
    header = (
        f"{space_group_number}\n"
        f"{unit_cell_constants}\n"
        f"{reference_dataset}\n"
        f"OUTPUT_FILE={output_hkl}\n"
        f"FRIEDEL'S_LAW=TRUE"
    )
    with open(xscale_inp, "w") as f:
        f.write(header + "\n")
        for file in xds_ascii_files:
            f.write(f"INPUT_FILE={file}\n")
    return xscale_inp

def read_info_txt(info_txt_file: str) -> Tuple[List[float], List[float], List[int]]:
    pass

def xscale_start(
    current_data_processing_folder: str,
    user: str,
    reserved_nodes: str,
    slurm_partition: str,
    ssh_private_key_path: str,
    ssh_public_key_path: str,  # currently unused; kept for interface compatibility
    login_node: str | None = None,
):
    """Prepare and submit the XSCALE job (locally or via SSH)."""
    os.chdir(current_data_processing_folder)
    job_name = Path(current_data_processing_folder).name
    slurmfile = Path(current_data_processing_folder) / f"{job_name}_XSCALE.sh"
    err_file = Path(current_data_processing_folder) / f"{job_name}_XSCALE.err"
    out_file = Path(current_data_processing_folder) / f"{job_name}_XSCALE.out"

    def get_slurm_header(partition: str, reservation: str | None = None, extras: list[str] | None = None):
        lines = [
            "#!/bin/sh\n",
            f"#SBATCH --job-name={job_name}\n",
            f"#SBATCH --partition={partition}\n",
            "#SBATCH --nodes=1\n",
            f"#SBATCH --output={out_file}\n",
            f"#SBATCH --error={err_file}\n",
        ]
        if reservation:
            lines.append(f"#SBATCH --reservation={reservation}\n")
        if extras:
            lines.extend(extras)
        return lines

    def get_common_xscale_commands():
        return [
            "source /etc/profile.d/modules.sh\n",
            "module load xray ccp4/7.1\n",
            "xscale\n",
        ]

    sbatch_file: list[str] = []
    ssh_command = ""

    is_maxwell = "maxwell" in reserved_nodes
    if is_maxwell:
        slurm_extras = [
            "#SBATCH --time=12:00:00\n",
            "#SBATCH --nice=100\n",
            "#SBATCH --mem=500000\n",
        ]
        sbatch_file += get_slurm_header("allcpu,upex,short", extras=slurm_extras)
        sbatch_file += get_common_xscale_commands()
    else:
        reserved_nodes_overloaded = are_the_reserved_nodes_overloaded(reserved_nodes)
        partition = slurm_partition if not reserved_nodes_overloaded else "allcpu,upex,short"
        reservation = reserved_nodes if not reserved_nodes_overloaded else None
        sbatch_file += get_slurm_header(partition, reservation)
        sbatch_file += get_common_xscale_commands()

        if login_node:
            ssh_command = (
                f"/usr/bin/ssh -o BatchMode=yes -o CheckHostIP=no "
                f"-o StrictHostKeyChecking=no -o GSSAPIAuthentication=no "
                f"-o GSSAPIDelegateCredentials=no -o PasswordAuthentication=no "
                f"-o PubkeyAuthentication=yes -o PreferredAuthentications=publickey "
                f"-o ConnectTimeout=10 -l {user} -i {ssh_private_key_path} {login_node}"
            )

    with open(slurmfile, "w") as fh:
        fh.writelines(sbatch_file)
    os.chmod(slurmfile, 0o755)

    # NOTE: If ssh_command is set, slurmfile must exist at the same path on the remote host (shared FS).
    submit_command = f'{ssh_command} "sbatch {slurmfile}"' if ssh_command else f"sbatch {slurmfile}"
    subprocess.run(submit_command, shell=True, check=True)


def main():
    input_path = sys.argv[1]
    current_data_processing_folder = sys.argv[2]
    output_csv = os.path.join(current_data_processing_folder, "summary.csv")
    
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(csv_header)
    
    reference_dataset = sys.argv[3]
    user = sys.argv[4]
    reserved_nodes = sys.argv[5]
    slurm_partition = sys.argv[6]
    ssh_private_key_path = sys.argv[7]
    ssh_public_key_path = sys.argv[8]
    login_node = sys.argv[9] if len(sys.argv) > 9 else None
    if "maxwell" in reserved_nodes:
        login_node = None  # local submit on maxwell

    while not os.path.exists(input_path):
        time.sleep(SLEEP_TIME)
    os.makedirs(current_data_processing_folder, exist_ok=True)
    os.chmod(current_data_processing_folder, 0o755)
    
    info_txt_file = os.path.join(input_path, 'info.txt')
    
    if not os.path.exists(info_txt_file):
        print("Info txt does not exist, not gonna work for you, dummy")
        exit()
    
    humidity, temperature, position = read_info_txt(info_txt_file)
    
    groups = group_humidity_plateaus(humidity, temperature, position)
    
    csvfile = open(output_csv, 'a', newline='')
    for humidity_level, data in groups.items():
        os.makedirs(os.path.join(current_data_processing_folder, f'humidity_{str(humidity_level).replace(".", "p")}'), exist_ok=True)
        positions = data['positions']
        xds_ascii_files = []
        space_groups = []
        UC_a = []
        UC_b = []
        UC_c = []
        alpha = []
        beta = []
        gamma = []
        for position in positions:
            position_path = os.path.join(input_path, str(position))
            if not os.path.exists(position_path):
                print(f"Position path {position_path} does not exist, skipping")
                continue
            # Derive SG and UC
            gxparm_filepath = os.path.join(position_path, 'GXPARM.XDS')
            if not os.path.exists(gxparm_filepath):
                print("GXPARM file does not exist, skipping this position")
                continue
            space_group, unit_cell = read_gxparm_xds(gxparm_filepath)
            space_groups.append(space_group if space_group else np.nan)
            UC_a.append(unit_cell[0] if unit_cell else np.nan)
            UC_b.append(unit_cell[1] if unit_cell else np.nan)
            UC_c.append(unit_cell[2] if unit_cell else np.nan)
            alpha.append(unit_cell[3] if unit_cell else np.nan)
            beta.append(unit_cell[4] if unit_cell else np.nan)
            gamma.append(unit_cell[5] if unit_cell else np.nan)
            
            correct_lp_path = os.path.join(position_path, 'CORRECT.LP')
            if not os.path.exists(correct_lp_path):
                print("CORRECT.LP file does not exist, skipping this position")
                continue
            CChalf, Rmerge, Rmeas, Rpim = read_correct_lp_file(correct_lp_path)
            
            csvfile_writer = csv.writer(csvfile)
            csvfile_writer.writerow([position, UC_a[-1], UC_b[-1], UC_c[-1], alpha[-1], beta[-1], gamma[-1], space_groups[-1], CChalf, Rmeas, Rpim, Rmerge])
            xds_ascii_file = os.path.join(position_path, 'XDS_ASCII.HKL')
            if not os.path.exists(xds_ascii_file):
                print("XDS_ASCII.HKL file does not exist, skipping this position")
                continue
            xds_ascii_files.append(xds_ascii_file)
        median_space_group = np.median([sg for sg in space_groups if not np.isnan(sg)]) if space_groups else None
        avg_UC_a = mean([a for a in UC_a if not np.isnan(a)]) if UC_a else None
        avg_UC_b = mean([b for b in UC_b if not np.isnan(b)]) if UC_b else None
        avg_UC_c = mean([c for c in UC_c if not np.isnan(c)]) if UC_c else None
        avg_alpha = mean([al for al in alpha if not np.isnan(al)]) if alpha else None
        avg_beta = mean([be for be in beta if not np.isnan(be)]) if beta else None   
        avg_gamma = mean([ga for ga in gamma if not np.isnan(ga)]) if gamma else None
        space_group_line = f"SPACE_GROUP_NUMBER={int(median_space_group)}" if median_space_group else "!SPACE_GROUP_NUMBER"
        unit_cell_line = f"UNIT_CELL_CONSTANTS={avg_UC_a} {avg_UC_b} {avg_UC_c} {avg_alpha} {avg_beta} {avg_gamma}" if all(v is not None for v in [avg_UC_a, avg_UC_b, avg_UC_c, avg_alpha, avg_beta, avg_gamma]) else "!UNIT_CELL_CONSTANTS"
        reference_dataset_line = f"REFERENCE_DATA_SET={reference_dataset}" if reference_dataset.lower() != "none" else "!REFERENCE_DATA_SET"
        xscale_inp = xscale_inp_generating(
            input_path=input_path,
            current_data_processing_folder=os.path.join(current_data_processing_folder, f'humidity_{str(humidity_level).replace(".", "p")}'),
            space_group_number=space_group_line,
            unit_cell_constants=unit_cell_line,
            reference_dataset=reference_dataset_line,
            xds_ascii_files=xds_ascii_files,
            output_hkl=f"humidity_{str(humidity_level).replace('.', 'p')}.ahkl"
        ) 
        print(f"Generated XSCALE.INP at {xscale_inp}")
        # Submit
        xscale_start(
            current_data_processing_folder=os.path.join(current_data_processing_folder, f'humidity_{str(humidity_level).replace(".", "p")}'),
            user=user,
            reserved_nodes=reserved_nodes,
            slurm_partition=slurm_partition,
            ssh_private_key_path=ssh_private_key_path,
            ssh_public_key_path=ssh_public_key_path,
            login_node=login_node,
        )
    csvfile.close()
      
if __name__ == "__main__":
    main()