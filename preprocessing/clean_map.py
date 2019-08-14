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
from postprocessing.pdb_reader_writer import PDB_Reader_Writer


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
    input_data = deepcopy(input_map.data)

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
        for atom in atoms:
            n += 1
            PDB_Reader_Writer.write_single_pdb(
                file=fp,
                type='ATOM',
                chain='A',
                node=np.array([atom[2], atom[1], atom[0]]),
                seqnum=n
            )


    chimera_script = open(paths['output'] + 'resample.cmd', 'w')

    if 'hidedusts_file' in paths:
        emdb_id = paths['input'].split('/')[-2]
        with open(paths['hidedusts_file']) as f:
            hidedusts = json.load(f)

    if 'hidedusts_file' in paths and emdb_id in hidedusts:
        level, hidedust_size = hidedusts[emdb_id]

        chimera_script.write('open ' + paths['input'] + '\n'
                         'cofr models\n'
                         'cofr fixed\n'
                         'open ' + paths['bounding_box'] + '\n'                         
                         'move cofr mod #1\n'
                         'write relative #0 #1 ' + paths['bounding_box_centered'] + '\n'
                         'open ' + paths['bounding_box_centered'] + '\n'
                         'molmap #2 6 gridSpacing 1\n'                        
                         'volume #0 level ' + str(level) + '\n'
                         'sop hideDust #0 size ' + str(hidedust_size) + '\n'
                         'sel #0\n'
                         'mask sel #0\n'
                         'vop resample #3 onGrid #2.1\n'
                         'volume #4 save ' + paths['cleaned_map'])
    else:
        chimera_script.write('open ' + paths['input'] + '\n'
                         'cofr models\n'
                         'cofr fixed\n'
                         'open ' + paths['bounding_box'] + '\n'                         
                         'move cofr mod #1\n'
                         'write relative #0 #1 ' + paths['bounding_box_centered'] + '\n'
                         'open ' + paths['bounding_box_centered'] + '\n'
                         'molmap #2 6 gridSpacing 1\n'                        
                         'sel #0\n'
                         'vop resample sel onGrid #2.1\n'
                         'volume #3 save ' + paths['cleaned_map'])
                     
    chimera_script.close()

    script_finished = False
    while not script_finished:
        try:
            subprocess.run(['/usr/local/bin/chimera', '--nogui', chimera_script.name])
            script_finished = True
        except FileNotFoundError as error:
            if not create_symbolic_link():
                raise error

    os.remove(chimera_script.name)


def create_symbolic_link():
    """Creates symbolic link to chimera bin in /usr/local/bin if user wants to

    Returns
    -------
    link_created: bool
        Indicates whether or not the symbolic link was created
    """
    print('It looks like there is no link to chimera in /usr/local/bin')
    if input('Do you want to create one? (y/n) ') in ['y', 'yes']:
        chimera_bin = input('Enter path to chimera bin file: ')
        subprocess.run(['ln', '-s', chimera_bin, '/usr/local/bin/chimera'])

        return True
    else:
        return False
