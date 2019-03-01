import subprocess
from scipy.optimize import root_scalar
import os
import json


def update_paths(paths):
    paths['threshold'] = paths['output'] + 'threshold'


def execute(paths):
    if 'thresholds_file' in paths:
        emdb_id = paths['input'].split('/')[-2]
        with open(paths['thresholds_file']) as f:
            thresholds = json.load(f)

        if emdb_id in thresholds:
            return

    threshold = root_scalar(lambda t: 0.7 - sav(t, paths), bracket=[0, 10]).root

    with open(paths['threshold'], 'w') as f:
        f.write(str(threshold))


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
    chimera_script.write('open ' + paths['cleaned_map'] + '\n'
                         'volume #0 level ' + str(threshold) + '\n'
                         'measure volume #0\n'
                         'measure area #0\n')
    chimera_script.close()

    output = subprocess.check_output(['/usr/local/bin/chimera', '--nogui', chimera_script.name])
    volume, surface_area = parse_chimera_output(output)

    os.remove(chimera_script.name)

    return 0 if volume == 0 else surface_area / volume


def parse_chimera_output(output):
    volume, area = None, None
    lines = str(output).split('\\n')
    for line in lines:
        if 'area' in line and area is None:
            area = float(line.split(' = ')[-1])
        elif 'volume' in line and volume is None:
            volume = float(line.split(' = ')[-1])

    return volume, area
