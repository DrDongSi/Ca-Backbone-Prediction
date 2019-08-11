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
    paths['bounding_box'] = paths['output'] + 'bounding_box.pdb'
    paths['bounding_box_centered'] = paths['output'] + 'bounding_box_centered.pdb'


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

    chimera_run(paths, [
        'open %s' % paths['input'],
        'cofr models',
        'cofr fixed',
        'open %s' % paths['bounding_box'],
        'mov cofr mod #1',
        'write relative #0 #1 %s' % paths['bounding_box_centered'],
        'close #1',
        'open %s' % paths['bounding_box_centered'],
        'molmap #1 6 gridSpacing 1',
        'vop resample #0 onGrid #1.1',
        # 'vop gaussian #2 sDev 0.5',
        'volume #2 save %s' % paths['cleaned_map']
    ])

    # This can fail if the map does not contain any dust. However, the
    # prediction can continue.
    chimera_run(paths, [
        'open %s' % paths['cleaned_map'],
        'volume #0 level %f' % get_threshold(paths),
        'volume #0 step 1',
        'sop hideDust #0 size 30',
        'sop invertShown #0',
        'mask #0 #0 invert true',
        'volume #1 save %s' % paths['cleaned_map']
    ])


def chimera_run(paths, commands):
    with open(paths['output'] + 'clean_map.cmd', 'w') as fp:
        fp.write('\n'.join(commands))

    subprocess.run(['/usr/local/bin/chimera', '--nogui', fp.name])
    os.remove(fp.name)


def get_threshold(paths):
    if 'thresholds_file' in paths:
        emdb_id = paths['input'].split('/')[-2]

        with open(paths['thresholds_file']) as f:
            thresholds = json.load(f)

        if emdb_id in thresholds:
            return thresholds[emdb_id]

    with open(paths['threshold']) as f:
        return float(f.readline())
