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


def update_paths(paths):
    paths['cleaned_map'] = paths['output'] + 'cleaned_map.mrc'
    paths['blank_ent'] = paths['output'] + 'blank.ent'
    paths['blank_centered_ent'] = paths['output'] + 'blank_centered.ent'


def execute(paths):
    """Creates the cleaning and re-sampling script and passes it to chimera

    Parameters
    ----------
    paths: dict
        Contains relevant paths for input and output files for the current
        prediction
    """

    copyfile(os.getcwd() + '/Ca-Backbone-Prediction/preprocessing/blank.ent', paths['blank_ent'])

    chimera_script = open(paths['output'] + 'resample.cmd', 'w')

    chimera_script.write('open ' + paths['input'] + '\n'
                         'cofr models\n'
                         'cofr fixed\n'
                         'open ' + paths['blank_ent'] + '\n'                         
                         'move cofr mod #1\n'
                         'write relative #0 #1 ' + paths['blank_centered_ent'] + '\n'
                         'open ' + paths['blank_centered_ent'] + '\n'
                         'molmap #2 6 gridSpacing 1\n'                        
                         'volume #0 level 0.1\n'
                         'sop hideDust #0 size 1.0\n'
                         'sel #0\n'
                         'open ' + os.getcwd() + '/Ca-Backbone-Prediction/preprocessing/surfinvert.py\n'
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
