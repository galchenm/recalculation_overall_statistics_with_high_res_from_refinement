import numpy as np
import re
import sys
import pylab
from mpl_toolkits.mplot3d import Axes3D
import os

# Load filename from command line arguments
streamFileName = sys.argv[1]
markerSize = 0.5  # Default marker size

# Optionally, get a marker size from command line
if len(sys.argv) >= 3:
    markerSize = float(sys.argv[2])

# Open the stream file safely using 'with' to ensure it closes
with open(streamFileName, 'r') as f:
    stream = f.read()

# Create output directory for plots
output_path = os.path.join(os.path.dirname(os.path.abspath(streamFileName)), 'plots_res')
os.makedirs(output_path, exist_ok=True)

# Define colors and names for plotting
colors = ["r", "g", "b"]
xStarNames = ["astar", "bstar", "cstar"]

# Loop over each name to create a separate plot for each
for i in range(3):
    # Compile regex for each pattern
    pattern = re.compile(rf"{xStarNames[i]} = ([\+\-\d\.]+ [\+\-\d\.]+ [\+\-\d\.]+)")
    xStarStrings = pattern.findall(stream)  # Find all matches

    if not xStarStrings:
        print(f"No matches found for {xStarNames[i]}")
        continue

    # Initialize the array to store coordinates
    xStars = np.zeros((3, len(xStarStrings)), float)

    # Convert matched strings to floats and store in xStars
    for j in range(len(xStarStrings)):
        try:
            xStars[:, j] = np.array([float(s) for s in xStarStrings[j].split()])
        except ValueError as e:
            print(f"Error parsing values for {xStarNames[i]} at index {j}: {e}")
            continue

    # Clear any previous plot and create a new figure
    pylab.clf()
    fig = pylab.figure()
    ax = fig.add_subplot(111, projection='3d')  # Proper 3D initialization
    ax.scatter(xStars[0, :], xStars[1, :], xStars[2, :], marker=".", color=colors[i], s=markerSize)
    pylab.title(f"{xStarNames[i]}s")

    # Define output file path and save the plot
    out = os.path.join(output_path, os.path.basename(streamFileName).split('.')[0] + "_" + xStarNames[i] + '.png')
    pylab.savefig(out)
    print(f"Saved plot to {out}")

# Close all pyplot figures after saving
pylab.close()
