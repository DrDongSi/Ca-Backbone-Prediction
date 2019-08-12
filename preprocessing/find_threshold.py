"""Responsible for automatically finding threshold level if it was not provided
by the user

The threshold value is later used to set all values below the threshold to zero
in order to produce a density map which resembles those which were generated
by the mrc2pdb script with which the CNN was trained.

Threshold values are calculated using two different methods. First, we want to
find a threshold value such that the density map has a surface area to volume
ratio of 0.9. Secondly, we want to find a threshold value such that the ratio of
number of values larger than the threshold to number of non-zero values is
0.484. The median of both numbers is then used as the threshold value and stored
in a file at /preprocessing/threshold.
"""

import subprocess
from scipy.optimize import root_scalar
import os
import json
import mrcfile
from copy import deepcopy
import numpy as np

__author__ = 'Jonas Pfab'


def update_paths(paths):
    paths['threshold'] = paths['output'] + 'threshold'


def execute(paths):
    """Finds threshold level and writes it to threshold file if no threshold
    level was provided by the user"""
    if is_threshold_provided(paths):
        return

    map_data = deepcopy(mrcfile.open(paths['input'], mode='r').data).ravel()
    num_non_zero_values = count_values(map_data, 0)

    # Find threshold value such that the surface area to volume ratio is 0.9
    threshold1 = root_scalar(lambda t: 0.9163020188305991 - sav(t, paths),
                             bracket=[0, 10]).root
    # Finding threshold value such that the ratio of number of values larger
    # than the threshold to number of non-zero values is 0.4
    threshold2 = root_scalar(lambda t: 0.4640027954434909 - (count_values(map_data, t) / num_non_zero_values),
                             bracket=[0, 10]).root

    threshold = (threshold1 + threshold2) / 2

    with open(paths['threshold'], 'w') as f:
        f.write(str(threshold))


def count_values(map_data, threshold):
    """Count number of values in 'map_data' larger than given 'threshold'

    Parameters
    ----------
    map_data: array
        Protein density map data
    threshold: float
        Threshold level which values have to exceed in order to be counted

    Returns
    ----------
    num_values: float
        The number of values which are larger than the given 'threshold'
    """
    return len(np.where(map_data > threshold)[0])


def sav(threshold, paths):
    """Calculates surface area to volume ratio for cleaned map with given
    threshold

    Parameters
    ----------
    threshold: float
        Threshold level which is applied to cleaned map before the surface area
        and volume are measured
    paths: dict
        Contains path to cleaned map

    Returns
    ----------
    sav: float
        Surface area to volume ratio
    """
    chimera_script = open(paths['output'] + 'measure.cmd', 'w')
    chimera_script.write('open ' + paths['input'] + '\n'
                         'volume #0 level ' + str(threshold) + '\n'
                         'measure volume #0\n'
                         'measure area #0\n')
    chimera_script.close()

    output = subprocess.check_output(['/usr/local/bin/chimera', '--nogui', chimera_script.name])
    volume, surface_area = parse_sav(output)

    os.remove(chimera_script.name)

    return float('inf') if volume == 0 else surface_area / volume


def parse_sav(output):
    """Parses surface area to volume ratio from given chimera output"""
    volume, area = None, None
    lines = str(output).split('\\n')
    for line in lines:
        if 'area' in line and area is None:
            area = float(line.split(' = ')[-1])
        elif 'volume' in line and volume is None:
            volume = float(line.split(' = ')[-1])

    return volume, area


def is_threshold_provided(paths):
    """Checks if threshold value is already provided by user"""
    if 'thresholds_file' in paths:
        emdb_id = paths['input'].split('/')[-2]
        with open(paths['thresholds_file']) as f:
            thresholds = json.load(f)

        if emdb_id in thresholds:
            return True

    return False
