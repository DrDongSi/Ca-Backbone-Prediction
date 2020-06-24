"""Cleans input map

The module is responsible for cleaning the input map and re-sampling it to a
voxel size of 1. This is accomplished by creating a script of chimera commands
and then pass it to chimera in no GUI mode.

The cleaning can be applied on a mrc file by calling the 'resample_map'
method
"""

import subprocess
import os
from shutil import copyfile
import json


def update_paths(paths):
    paths['cleaned_map'] = paths['output'] + 'cleaned_map.mrc'
    paths['bounding_box'] = paths['output'] + 'bounding_box.ent'
    paths['bounding_box_centered'] = paths['output'] + 'bounding_box_centered.ent'


def execute(paths):
    """Creates the cleaning and re-sampling script and passes it to chimera

    Parameters
    ----------
    paths: dict
        Contains relevant paths for input and output files for the current
        prediction
	"""
    # remove the Ca-Backbone-Prediction from the copyfile
    copyfile(os.getcwd() + '/preprocessing/bounding_box.ent', paths['bounding_box'])

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
            # Changed the chimera path to an input
            subprocess.run([paths['chimera_path'], '--nogui', chimera_script.name])
            script_finished = True
        except FileNotFoundError as error:
            raise error
            #if not create_symbolic_link():
                #raise error

    os.remove(chimera_script.name)

# Removing the function, as the chimera link parameter should handle the symbolic link.