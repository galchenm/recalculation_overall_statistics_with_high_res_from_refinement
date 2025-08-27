import os
import sys
import glob
import shlex
import subprocess
from collections import Counter
from pathlib import Path
from typing import Optional, Tuple

LIMIT_FOR_RESERVED_NODES = 1000
SLEEP_TIME = 15  # seconds

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


def xscale_inp_generating(
    input_path: str,
    current_data_processing_folder: str,
    space_group_number: str = "!SPACE_GROUP_NUMBER",
    unit_cell_constants: str = "!UNIT_CELL_CONSTANTS",
    reference_dataset: str = "!REFERENCE_DATA_SET",
    output_hkl: str = "output.hkl",
) -> Path:
    """
    Generate XSCALE.INP in current_data_processing_folder.
    The parameters must already be full XSCALE lines, e.g.:
    space_group_number="SPACE_GROUP_NUMBER=96"
    unit_cell_constants="UNIT_CELL_CONSTANTS=a b c alpha beta gamma"
    reference_dataset="REFERENCE_DATA_SET=path/to/reference.mtz" or "!REFERENCE_DATA_SET"
    """
    xds_ascii_files = glob.glob(os.path.join(input_path, "**", "XDS_ASCII.HKL"), recursive=True)
    if not xds_ascii_files:
        raise FileNotFoundError(f"No XDS_ASCII.HKL files found under: {input_path}")

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


def read_gxparm_xds(filepath: str) -> Tuple[Optional[int], Optional[Tuple[float, float, float, float, float, float]]]:
    """
    Reads a GXPARM.XDS file and extracts space group number and unit cell parameters.
    Returns (space_group, (a,b,c,alpha,beta,gamma)) or (None, None) if not found.
    """
    with open(filepath, "r") as f:
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


def read_and_average_gxparm(input_path: str) -> Tuple[Optional[int], Optional[Tuple[float, float, float, float, float, float]]]:
    """
    Reads all GXPARM.XDS files under input_path.
    Returns:
        - most frequent space group number (mode),
        - average unit cell (a,b,c,alpha,beta,gamma) computed ONLY over files with that SG.
    """
    sg_list = []
    uc_list_by_sg = {}
    filepaths = glob.glob(os.path.join(input_path, "**", "GXPARM.XDS"), recursive=True)

    for filepath in filepaths:
        sg, uc = read_gxparm_xds(filepath)
        if sg is None or uc is None:
            continue
        sg_list.append(sg)
        uc_list_by_sg.setdefault(sg, []).append(uc)

    if not sg_list:
        return None, None

    most_common_sg = Counter(sg_list).most_common(1)[0][0]
    cell_list = uc_list_by_sg.get(most_common_sg, [])
    if not cell_list:
        return most_common_sg, None

    n = len(cell_list)
    avg_uc = tuple(sum(vals[i] for vals in cell_list) / n for i in range(6))
    return most_common_sg, avg_uc


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
    if len(sys.argv) < 9:
        raise SystemExit(
            "Usage: script.py <input_path> <proc_folder> <reference_dataset> <user> "
            "<reserved_nodes> <slurm_partition> <ssh_private_key> <ssh_public_key> [login_node]"
        )

    input_path = sys.argv[1]
    current_data_processing_folder = sys.argv[2]
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

    # Derive SG and UC
    space_group, unit_cell = read_and_average_gxparm(input_path)

    # Build XSCALE fields
    sg_line = f"SPACE_GROUP_NUMBER={space_group}" if isinstance(space_group, int) else "!SPACE_GROUP_NUMBER"
    if unit_cell and all(x is not None for x in unit_cell):
        uc_str = " ".join(f"{x:.4f}" for x in unit_cell)
        uc_line = f"UNIT_CELL_CONSTANTS={uc_str}"
    else:
        uc_line = "!UNIT_CELL_CONSTANTS"
    ref_line = f"REFERENCE_DATA_SET={reference_dataset}" if reference_dataset and reference_dataset.lower() != "none" else "!REFERENCE_DATA_SET"

    # Write XSCALE.INP
    xscale_inp_generating(
        input_path=input_path,
        current_data_processing_folder=current_data_processing_folder,
        space_group_number=sg_line,
        unit_cell_constants=uc_line,
        reference_dataset=ref_line,
        output_hkl="output.hkl",
    )

    # Submit
    xscale_start(
        current_data_processing_folder=current_data_processing_folder,
        user=user,
        reserved_nodes=reserved_nodes,
        slurm_partition=slurm_partition,
        ssh_private_key_path=ssh_private_key_path,
        ssh_public_key_path=ssh_public_key_path,
        login_node=login_node,
    )


if __name__ == "__main__":
    main()
