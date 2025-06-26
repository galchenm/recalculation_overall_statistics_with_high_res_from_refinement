# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
import logging


def ave_resolution_plot(stream_filename):
    """
    Generate a histogram of diffraction resolution limits from a CrystFEL stream file.

    Parameters:
        stream_filename (str): Path to the .stream file.

    Returns:
        str: Path to the saved plot image file, or empty string on failure.
    """
    logger = logging.getLogger('app')
    path_to_plots = os.path.join(os.path.dirname(stream_filename), 'plots_res')
    os.makedirs(path_to_plots, exist_ok=True)

    output_file = os.path.join(
        path_to_plots,
        os.path.splitext(os.path.basename(stream_filename))[0] + '-ave-resolution.png'
    )

    resolutions = []

    try:
        with open(stream_filename, 'r') as f:
            for line in f:
                if "diffraction_resolution_limit" in line:
                    try:
                        res = float(line.split('= ')[1].split()[0])
                        resolutions.append(res)
                    except (IndexError, ValueError):
                        logger.warning(f"Could not parse resolution in line: {line.strip()}")

        if not resolutions:
            logger.info(f'No resolution data found in {os.path.basename(stream_filename)}')
            return ""

        res_array = np.array(resolutions)
        mean_val = np.mean(res_array)
        max_val = np.max(res_array)
        min_val = np.min(res_array)
        std_val = np.std(res_array)

        print(f"Mean: {mean_val:.2f} nm⁻¹ = {10.0 / mean_val:.2f} Å")
        print(f"Best: {max_val:.2f} nm⁻¹ = {10.0 / max_val:.2f} Å")
        print(f"Worst: {min_val:.2f} nm⁻¹ = {10.0 / min_val:.2f} Å")
        print(f"Std deviation: {std_val:.2f} nm⁻¹")

        # Plot histogram
        plt.figure()
        plt.hist(resolutions, bins=30, color='skyblue', edgecolor='black')
        plt.title('Resolution Based on Indexing Results')
        plt.xlabel('Resolution (nm⁻¹)')
        plt.ylabel('Frequency')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(output_file)
        plt.close()

        return output_file

    except Exception as e:
        logger.error(f'Error generating ave-res plot for {os.path.basename(stream_filename)}: {e}')
        return ""
