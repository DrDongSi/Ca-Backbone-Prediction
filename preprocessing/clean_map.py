"""Cleans input map

The module is responsible for cleaning the input map and re-sampling it to a
voxel size of 1. This is accomplished by creating a script of chimera commands
and then pass it to chimera in no GUI mode.

The cleaning can be applied on a mrc file by calling the 'resample_map'
method
"""

import subprocess
import os
import json
import mrcfile
import numpy as np
from copy import deepcopy


def update_paths(paths):
    paths['cleaned_map'] = paths['output'] + 'cleaned_map.mrc'
    paths['blank_ent'] = paths['output'] + 'blank.ent'
    paths['blank_centered_ent'] = paths['output'] + 'blank_centered.ent'
    paths['bounding_box'] = paths['output'] + 'bounding_box.pdb'


def execute(paths):
    """Creates the cleaning and re-sampling script and passes it to chimera

    Parameters
    ----------
    paths: dict
        Contains relevant paths for input and output files for the current
        prediction
    """
    input_map = mrcfile.open(paths['input'])
    origin = input_map.header.origin.item(0)
    input_data = deepcopy(input_map.data)
    input_data[input_data < get_threshold(paths)] = 0

    x_is, y_is, z_is = np.nonzero(input_data)
    x_is, y_is, z_is = sorted(x_is), sorted(y_is), sorted(z_is)

    atoms = [
        [x_is[0], 0, 0],
        [x_is[-1], 0, 0],
        [0, y_is[0], 0],
        [0, y_is[-1], 0],
        [0, 0, z_is[0]],
        [0, 0, z_is[-1]]
    ]

    with open(paths['bounding_box'], 'w') as fp:
        n = 0
        for i in range(min(100, len(atoms))):
            atom = atoms[i]
            n += 1
            fp.write('ATOM      1  CA  GLY A' + str(n).rjust(4) +
                     '    ' + '{0:.3f}'.format(atom[2] + origin[0]).rjust(8) +
                     '{0:.3f}'.format(atom[1] + origin[1]).rjust(8) +
                     '{0:.3f}'.format(atom[0] + origin[2]).rjust(8) +
                     '  1.00  0.00           C  \n')

    chimera_script = open(paths['output'] + 'resample.cmd', 'w')
    chimera_script.write(
        'open %s\n' % paths['input'] +
        'cofr models\n' +
        'cofr fixed \n' +
        'open %s\n' % paths['bounding_box'] +
        'mov cofr mod #1\n' +
        'molmap #1 6 gridSpacing 1\n' +
        'vop resample #0 onGrid #1.1\n' +
        'volume #2 level %f\n' % get_threshold(paths) +
        'volume #2 step 1\n' +
        'sop hideDust #2 size 30\n' +
        'sop invertShown #2\n' +
        'mask #2 #2 invert true\n' +
        'vop gaussian #3 sDev 0.5\n' +
        'volume #4 save %s\n' % paths['cleaned_map']
    )
    chimera_script.close()

    subprocess.run(['/usr/local/bin/chimera', '--nogui', chimera_script.name])
    os.remove(chimera_script.name)


def get_threshold(paths):
    if 'thresholds_file' in paths:
        emdb_id = paths['input'].split('/')[-2]

        with open(paths['thresholds_file']) as f:
            thresholds = json.load(f)

        if emdb_id in thresholds:
            return thresholds[emdb_id]

    with open(paths['threshold']) as f:
        return float(f.readline())
