import os
import glob
import sys

def xscale_inp_generating(input_path, output_path):
    XDS_ASCII_files = glob.glob(os.path.join(input_path, "*/**/XDS_ASCII.HKL"), recursive=True)
    if XDS_ASCII_files:
        with open(os.path.join(output_path, "XSCALE.INP"), "w") as f:
            header = f"""{space_group_number}
            {unit_cell_constants}
            {reference_dataset}
            OUTPUT_FILE={output_hkl_name}.ahkl
            FRIEDEL'S_LAW=TRUE"""

            f.write(header + '\n')
            for file in XDS_ASCII_files:
                f.write(f"INPUT_FILE={file}\n")

def main():
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    
    os.makedirs(output_path, exist_ok=True)
    os.chmod(output_path, 0o755)
    os.chdir(output_path)

    xscale_inp_generating(input_path)

if __name__ == "__main__":
    main()

