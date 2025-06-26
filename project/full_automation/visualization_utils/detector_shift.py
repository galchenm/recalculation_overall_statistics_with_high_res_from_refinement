# -*- coding: utf-8 -*-

import os
import re
import numpy as np
import matplotlib.pyplot as plt


def detector_shift(filename, geom=None, rerun_detector_shift=False):
    shift_xy_pattern = re.compile(r"^predict_refine/det_shift\sx\s=\s([-0-9.]+)\sy\s=\s([-0-9.]+)\smm$")
    shift_z_pattern = re.compile(r"^predict_refine/clen_shift\s=\s([-0-9.]+)\smm$")

    x_shifts, y_shifts, z_shifts = [], [], []

    with open(filename, 'r') as f:
        for line in f:
            match_xy = shift_xy_pattern.match(line)
            if match_xy:
                x_shifts.append(float(match_xy.group(1)))
                y_shifts.append(float(match_xy.group(2)))
                continue

            match_z = shift_z_pattern.match(line)
            if match_z:
                z_shifts.append(float(match_z.group(1)))

    if not x_shifts or not y_shifts:
        print("No detector shift data found.")
        return

    mean_x = sum(x_shifts) / len(x_shifts)
    mean_y = sum(y_shifts) / len(y_shifts)
    print('Mean shifts: dx = {:.2f} mm,  dy = {:.2f} mm'.format(mean_x, mean_y))

    if rerun_detector_shift and geom:
        print("Apply shifts to geometry")
        output_geom = "{}-{}.geom".format(
            os.path.splitext(geom)[0],
            os.path.basename(os.path.dirname(os.path.dirname(filename)))
        )
        print('Applying corrections to {}, output filename {}'.format(geom, output_geom))
        print('Shifts: dx = {:.2f} mm, dy = {:.2f} mm'.format(mean_x, mean_y))

        panel_resolutions = {}
        default_res = 0.0

        with open(geom, 'r') as g, open(output_geom, 'w') as h:
            for line in g:
                match1 = re.match(r"^\s*res\s+=\s+([0-9.]+)", line)
                if match1:
                    default_res = float(match1.group(1))
                    h.write(line)
                    continue

                match2 = re.match(r"^\s*(.*)/res\s+=\s+([0-9.]+)", line)
                if match2:
                    panel = match2.group(1)
                    panel_resolutions[panel] = float(match2.group(2))
                    default_res = float(match2.group(2))
                    h.write(line)
                    continue

                match3 = re.match(r"^\s*(.*)/corner_x\s+=\s+([-0-9.]+)", line)
                if match3:
                    panel = match3.group(1)
                    cx = float(match3.group(2))
                    res = panel_resolutions.get(panel, default_res)
                    shift_x = mean_x * res * 1e-3
                    h.write('%s/corner_x = %f\n' % (panel, cx + shift_x))
                    continue

                match4 = re.match(r"^\s*(.*)/corner_y\s+=\s+([-0-9.]+)", line)
                if match4:
                    panel = match4.group(1)
                    cy = float(match4.group(2))
                    res = panel_resolutions.get(panel, default_res)
                    shift_y = mean_y * res * 1e-3
                    h.write('%s/corner_y = %f\n' % (panel, cy + shift_y))
                    continue

                h.write(line)

        print('Saved new geometry')
    elif rerun_detector_shift and not geom:
        print("Warning: rerun_detector_shift is True but no geometry file provided. Skipping geometry update.")
    else:
        print("Don't apply shifts to geometry")

    def plot_new_centre(x, y):
        circle = plt.Circle((x, y), 0.1, color='r', fill=False)
        plt.gca().add_artist(circle)
        plt.plot(x, y, 'm8')
        plt.grid(True)

    nbins = 200
    H, xedges, yedges = np.histogram2d(x_shifts, y_shifts, bins=nbins)
    H = np.rot90(H)
    H = np.flipud(H)
    Hmasked = np.ma.masked_where(H == 0, H)

    fig2 = plt.figure()
    ax = fig2.add_subplot(111)
    c = ax.pcolormesh(xedges, yedges, Hmasked)
    ax.set_title('Detector shifts according to prediction refinement')
    ax.set_xlabel('x shift / mm')
    ax.set_ylabel('y shift / mm')
    ax.plot(0, 0, 'cH')
    plot_new_centre(mean_x, mean_y)
    fig2.colorbar(c, ax=ax, label='Counts')

    output_dir = os.path.join(os.path.dirname(filename), 'plots_res')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_img = os.path.join(
        output_dir, os.path.basename(filename).split('.')[0] + '-detector-shift.png'
    )
    if os.path.exists(output_img):
        os.remove(output_img)
    fig2.savefig(output_img)
    plt.close(fig2)
    print('Saved detector shift plot to {}'.format(output_img))