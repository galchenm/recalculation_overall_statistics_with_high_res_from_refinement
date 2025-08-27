import os
import sys
import time
import glob
import csv
from datetime import datetime
import re

# 1 step
def extract_cbf_files_modification_times(input_path: str, output_csv: str = 'file_mod_times.csv') -> str:
    """
    Reads a list of files from `input_path` and writes a CSV `output_csv` with
    filename and last modification time (HH:MM).
    """
    files = glob.glob(os.path.join(input_path, "**", "*.cbf"), recursive=True)

    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Filename", "ModificationTime"])  # CSV header

        for file_path in files:
            if os.path.isfile(file_path):
                filename = os.path.basename(file_path)
                mod_time = datetime.fromtimestamp(os.path.getmtime(file_path)).strftime("%H:%M")
                writer.writerow([filename, mod_time])
            else:
                print(f"Warning: '{file_path}' does not exist or is not a regular file")

    return output_csv

# 2 step
def extract_shroud_temperatures_times(input_log_file: str, output_file: str):
    """
    Reads a log CSV with datetime and multiple values, extracts HH:MM from the first column,
    keeps only the first 4 value columns, and writes a new CSV with a header.
    """
    if not os.path.isfile(input_log_file):
        raise FileNotFoundError(f"Input file '{input_log_file}' not found.")

    with open(input_log_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Write header
        outfile.write("Time;Value1;Value2;Value3;Value4\n")

        # Skip the first line (header)
        next(infile)

        for line in infile:
            line = line.strip()
            if not line:
                continue  # skip empty lines

            columns = line.split(';')
            if len(columns) < 5:
                continue  # skip lines that don't have enough columns

            datetime_str = columns[0]  # e.g., "2025-05-02 23:38:48"
            hh_mm = ":".join(datetime_str.split(' ')[1].split(':')[:2])

            # Keep only first 4 value columns
            values = columns[1:5]

            # Write to output
            outfile.write(f"{hh_mm};{';'.join(values)}\n")

    return output_file

# 3 step
def filter_unique_times_from_log(input_file: str, output_file: str):
    """
    Reads a CSV with a header and 'HH:MM' times in the first column,
    keeps only the first occurrence of each time, and writes to a new CSV.
    """
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"Input file '{input_file}' not found.")

    seen_times = set()

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Read header
        header = infile.readline()
        outfile.write(header)

        # Process remaining lines
        for line in infile:
            line = line.strip()
            if not line:
                continue  # skip empty lines

            hh_mm = line.split(';')[0]
            if hh_mm not in seen_times:
                seen_times.add(hh_mm)
                outfile.write(line + '\n')

    return output_file

# 4 step
def mapping_files_with_humidity_level_temperature_via_timestamp(file1: str, file2: str, output_file: str):
    """
    Combines two CSV files based on matching times.

    Args:
        file1: CSV with "Filename, ModificationTime"
        file2: CSV with "Time;Value1;Value2;Value3;Value4"
        output_file: CSV to save combined results with "Filename, ModificationTime, Value1"
    """
    if not os.path.isfile(file1):
        raise FileNotFoundError(f"Input file1 '{file1}' not found.")
    if not os.path.isfile(file2):
        raise FileNotFoundError(f"Input file2 '{file2}' not found.")

    # Load file2 into a dictionary: Time -> Value1
    time_to_value1 = {}
    with open(file2, 'r') as f2:
        reader = csv.reader(f2, delimiter=';')
        for row in reader:
            if not row or len(row) < 2:
                continue
            time_to_value1[row[0].strip()] = row[1].strip()

    # Process file1 and write to output
    with open(file1, 'r') as f1, open(output_file, 'w', newline='') as fout:
        writer = csv.writer(fout)
        writer.writerow(["Filename", "ModificationTime", "Value1"])

        reader1 = csv.reader(f1)
        for row in reader1:
            if not row or len(row) < 2:
                continue

            filename = row[0].strip()
            mod_time = row[1].strip()

            value1 = time_to_value1.get(mod_time)
            if value1:
                writer.writerow([filename, mod_time, value1])
            else:
                print(f"Warning: No matching Time for {filename} with ModificationTime {mod_time}")

    return output_file

# 5 step. weird 
def generate_crystal_csv(
    search_directory,
    csv_file,
    search_term1="chito",
    search_term2="lyso"
):
    """
    Searches for GXPARM.XDS files in subdirectories matching search terms,
    extracts unit cell constants and CC/R statistics, and writes to a CSV.
    """
    # CSV header
    header = [
        "Directory", "Unit Cell a", "Unit Cell b", "Unit Cell c",
        "Unit Cell alpha", "Unit Cell beta", "Unit Cell gamma",
        "CC/2", "R-merge", "R-meas", "R-pim", "Subdirectory", "Subdirectory_12"
    ]

    with open(csv_file, 'w', newline='') as f_out:
        writer = csv.writer(f_out)
        writer.writerow(header)

        # Walk through directories to find GXPARM.XDS files
        for root, dirs, files in os.walk(search_directory):
            for file in files:
                if file == "GXPARM.XDS":
                    directory = os.path.abspath(root)
                    if (search_term1 in directory) or (search_term2 and search_term2 in directory):
                        gxparm_path = os.path.join(directory, file)

                        # Extract unit cell constants from GXPARM.XDS (line 4, columns 2-7)
                        try:
                            with open(gxparm_path, 'r') as gx:
                                lines = gx.readlines()
                                line4 = lines[3].split()  # 0-based index
                                unit_cell = line4[1:7]  # columns 2-7
                        except Exception:
                            unit_cell = [""] * 6

                        # Extract CC/2, R-merge, R-meas, R-pim from CORRECT.LP
                        correct_lp_path = os.path.join(directory, "CORRECT.LP")
                        cc2 = r_merge = r_meas = r_pim = ""
                        if os.path.isfile(correct_lp_path):
                            with open(correct_lp_path, 'r') as lp:
                                lp_lines = lp.readlines()
                                for i, line in enumerate(lp_lines):
                                    if "WILSON STATISTICS" in line:
                                        # Using relative line offsets based on Bash script
                                        try:
                                            cc2 = lp_lines[i-13].split()[10]
                                            r_merge = lp_lines[i-13].split()[4]
                                            r_meas = lp_lines[i-13].split()[12]
                                            r_pim = lp_lines[i-21].split()[9]
                                        except Exception:
                                            pass
                                        break

                        # Extract subdirectory names
                        path_parts = directory.split(os.sep)
                        subdirectory = path_parts[-5] if len(path_parts) >= 5 else ""
                        subdirectory_12 = path_parts[-3] if len(path_parts) >= 3 else ""

                        # Write row to CSV
                        writer.writerow([directory, *unit_cell, cc2, r_merge, r_meas, r_pim, subdirectory, subdirectory_12])

    return csv_file

# 6 step
def modify_correct_output(input_file="correct-readout.out", output_file="modify_correct_output.list"):
    """
    Processes a CSV file, extracts POSITION from the first column, selects specific columns,
    and writes a modified CSV with header: POSITION,Unit Cell a,Unit Cell c,CC/2,R-meas.
    """
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"Input file '{input_file}' not found.")

    with open(input_file, 'r') as fin, open(output_file, 'w', newline='') as fout:
        reader = csv.reader(fin)
        writer = csv.writer(fout)

        # Write header
        writer.writerow(["POSITION", "Unit Cell a", "Unit Cell c", "CC/2", "R-meas"])

        # Skip header in input
        next(reader, None)

        for row in reader:
            if not row or len(row) < 10:
                continue  # skip empty or malformed lines

            # Extract POSITION from the last part of the first column
            path_parts = row[0].split('/')
            last_part = path_parts[-1]
            pos_parts = last_part.split("POSITION_")
            position = pos_parts[1] if len(pos_parts) == 2 else ""

            # Select specific columns
            unit_cell_a = row[1]
            unit_cell_c = row[3]
            cc2 = row[7]
            r_meas = row[9]

            writer.writerow([position, unit_cell_a, unit_cell_c, cc2, r_meas])

    return output_file

# 7 step. wtf is that
def merge_xds_data(
    humid_file="05_combine_02and04.list",
    raw_position_file="correct-readout.out",
    output_file="08_merge_xds_data.list"
):
    """
    Merge humidity data and unit cell data by POSITION into a single CSV.
    
    Args:
        humid_file: CSV from combine_02and04.list
        raw_position_file: CSV from correct-readout.out
        output_file: final merged CSV
    """
    if not os.path.isfile(humid_file) or not os.path.isfile(raw_position_file):
        raise FileNotFoundError("One or both input files not found.")

    # Step 1: Preprocess raw position file
    corrected_position_file = "tmp_corrected_position.csv"
    with open(raw_position_file, 'r') as f_in, open(corrected_position_file, 'w', newline='') as f_out:
        reader = csv.reader(f_in)
        writer = csv.writer(f_out)
        writer.writerow(["POSITION","Cell_a","Cell_b","Cell_c","Alpha","Beta","Gamma","CC12","Rmerge","Rmeas","Rpim","Subdir1","Subdir2"])

        next(reader, None)  # skip header
        for row in reader:
            if len(row) < 8:
                continue
            dir_field, unitcell_str, cc2, rmerge, rmeas, rpim, sdir1, sdir2 = row[:8]

            match = re.search(r"POSITION_([0-9]+)", dir_field)
            if not match:
                continue
            pos = match.group(1)

            try:
                a, b, c, alpha, beta, gamma = unitcell_str.split()
            except ValueError:
                a = b = c = alpha = beta = gamma = ""

            writer.writerow([pos, a, b, c, alpha, beta, gamma, cc2, rmerge, rmeas, rpim, sdir1, sdir2])

    # Step 2: Extract run numbers from humidity file
    tmp_humid_file = "tmp_humid.csv"
    with open(humid_file, 'r') as f_in, open(tmp_humid_file, 'w', newline='') as f_out:
        reader = csv.reader(f_in)
        writer = csv.writer(f_out)
        writer.writerow(["POSITION","Filename","ModificationTime","Value1"])
        next(reader, None)  # skip header

        for row in reader:
            if len(row) < 3:
                continue
            filename = row[0]
            mod_time = row[1]
            value1 = row[2]

            match = re.search(r"_([0-9]+)_00001\.cbf", filename)
            if not match:
                continue
            run_num = int(match.group(1))
            writer.writerow([run_num, filename, mod_time, value1])

    # Step 3: Load both CSVs into dictionaries for merging
    position_data = {}
    with open(corrected_position_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            position_data[row["POSITION"]] = row

    humid_data = []
    with open(tmp_humid_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            humid_data.append(row)

    # Step 4: Merge and write final CSV
    header = ["POSITION","Filename","ModificationTime","Value1","Cell_a","Cell_b","Cell_c","Alpha","Beta","Gamma","CC12","Rmerge","Rmeas","Rpim","Subdir1","Subdir2"]
    with open(output_file, 'w', newline='') as f_out:
        writer = csv.writer(f_out)
        writer.writerow(header)

        for row in humid_data:
            pos = str(row["POSITION"])
            if pos in position_data:
                pos_row = position_data[pos]
                merged_row = [
                    pos,
                    row["Filename"],
                    row["ModificationTime"],
                    row["Value1"],
                    pos_row["Cell_a"],
                    pos_row["Cell_b"],
                    pos_row["Cell_c"],
                    pos_row["Alpha"],
                    pos_row["Beta"],
                    pos_row["Gamma"],
                    pos_row["CC12"],
                    pos_row["Rmerge"],
                    pos_row["Rmeas"],
                    pos_row["Rpim"],
                    pos_row["Subdir1"],
                    pos_row["Subdir2"]
                ]
                writer.writerow(merged_row)

    # Clean up temporary files
    os.remove(corrected_position_file)
    os.remove(tmp_humid_file)

    return output_file
