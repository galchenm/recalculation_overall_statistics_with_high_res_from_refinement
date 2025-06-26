# -*- coding: utf-8 -*-

import re
import numpy as np
import pylab
from pathlib import Path


def orientation_plot(stream_file_name, run_name, markerSize=0.5):
    """
    Generate a 3D scatter plot of astar vectors from a CrystFEL stream file.

    Parameters:
        stream_file_name (str or Path): Path to the .stream file.
        run_name (str): A string used to name the output plot image.
        markerSize (float): Size of points in the 3D scatter plot.

    Returns:
        str: Path to the saved PNG plot file.
    """
    output_filename = run_name.replace("/", "_")
    stream_file = Path(stream_file_name)

    output_path = stream_file.parent / 'plots_res'
    output_path.mkdir(exist_ok=True)

    out = output_path / f"{output_filename}_astars_new.png"

    if out.exists():
        return str(out)

    print('Plotting')

    # Extract 'astar' vectors from stream
    with stream_file.open('r') as f:
        stream = f.read()

    pattern = re.compile(r"astar = ([\+\-\d\.eE]+ [\+\-\d\.eE]+ [\+\-\d\.eE]+)")
    aStarStrings = pattern.findall(stream)

    if not aStarStrings:
        print("No astar vectors found.")
        return ""

    aStars = np.array([list(map(float, x.split())) for x in aStarStrings])

    # Plot
    pylab.clf()
    fig = pylab.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(aStars[:, 0], aStars[:, 1], aStars[:, 2],
               marker=".", color="r", s=markerSize)
    ax.set_title("astars")
    ax.set_xlabel("a*_x")
    ax.set_ylabel("a*_y")
    ax.set_zlabel("a*_z")

    pylab.savefig(out)
    pylab.close()

    return str(out)
