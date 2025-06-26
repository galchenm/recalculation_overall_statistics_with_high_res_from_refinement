import re

def parse_cryst1_from_pdb(file_path):
    """Parse the CRYST1 line from a PDB file to extract unit cell parameters.
    This function reads a PDB file and extracts the unit cell parameters from the CRYST1 line.
    The CRYST1 line contains the following parameters:
    - a: length of the unit cell in the x direction
    - b: length of the unit cell in the y direction
    - c: length of the unit cell in the z direction
    - alpha: angle between b and c in degrees
    - beta: angle between a and c in degrees
    - gamma: angle between a and b in degrees
    The function returns these parameters as a tuple (a, b, c, alpha, beta, gamma).
    Args:
        file_path (str): The path to the PDB file.
    Returns:
        tuple: A tuple containing the unit cell parameters (a, b, c, alpha, beta, gamma).
    Raises:
        ValueError: If the CRYST1 line is not found in the PDB file.
    """
    # Regular expression to match the CRYST1 line and capture the parameters
    pattern = r"CRYST1\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)"
    
    # Open and read the PDB file
    with open(file_path, 'r') as file:
        for line in file:
            # Search for the CRYST1 line
            if line.startswith("CRYST1"):
                # Match and capture the unit cell parameters
                match = re.match(pattern, line.strip())
                if match:
                    # Extract the values from the match
                    a = float(match.group(1))
                    b = float(match.group(2))
                    c = float(match.group(3))
                    alpha = float(match.group(4))
                    beta = float(match.group(5))
                    gamma = float(match.group(6))
                    
                    return a, b, c, alpha, beta, gamma
    # If no CRYST1 line is found
    raise ValueError("No CRYST1 line found in the provided PDB file.")

def parse_UC_file(UC_file):
    """Parse a unit cell file to extract unit cell parameters.
    This function checks the file extension and calls the appropriate parsing function. 
    If the file is a PDB file, it uses the parse_cryst1_from_pdb function.
    If the file is not a PDB file, it uses a regular expression to extract the unit cell parameters.
    Args:
        UC_file (str): The path to the unit cell file.
    Returns:
        tuple: A tuple containing the unit cell parameters (a, b, c, alpha, beta, gamma).
    Raises:
        ValueError: If the unit cell parameters are not found in the file.
    """
    if UC_file.endswith('pdb'):
        return parse_cryst1_from_pdb(UC_file)
    else:
        # Regular expressions to capture the unit cell parameters
        pattern = r"a = ([\d.]+) A.*?b = ([\d.]+) A.*?c = ([\d.]+) A.*?al = ([\d.]+) deg.*?be = ([\d.]+) deg.*?ga = ([\d.]+) deg"
        
        # Read the content of the file
        with open(UC_file, 'r') as file:
            file_content = file.read()

        # Search for the pattern in the file content
        match = re.search(pattern, file_content, re.DOTALL)

        if match:
            # Extract values from the match
            a = float(match.group(1))
            b = float(match.group(2))
            c = float(match.group(3))
            al = float(match.group(4))
            be = float(match.group(5))
            ga = float(match.group(6))
            return a, b, c, al, be, ga
        else:
            raise ValueError("Unit cell parameters not found in the provided file.")
            return None, None, None, None, None, None