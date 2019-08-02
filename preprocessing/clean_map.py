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


def update_paths(paths):
    paths['cleaned_map'] = paths['output'] + 'cleaned_map.mrc'


def execute(paths):
    """Creates the cleaning and re-sampling script and passes it to chimera

    Parameters
    ----------
    paths: dict
        Contains relevant paths for input and output files for the current
        prediction
    """

    chimera_script = open(paths['output'] + 'resample.cmd', 'w')

    chimera_script.write('open %s\n' % paths['input'] +
                         'volume #0 voxelSize 1\n'
                         'volume #0 level %d\n' % get_threshold(paths) +
                         'sop hideDust #0 size 30\n'
                         'vop gaussian #0 sDev 0.5\n'
                         'volume #1 save ' + paths['cleaned_map'])

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


def get_threshold(paths):
    if 'thresholds_file' in paths:
        emdb_id = paths['input'].split('/')[-2]

        with open(paths['thresholds_file']) as f:
            thresholds = json.load(f)

        if emdb_id in thresholds:
            return thresholds[emdb_id]

    with open(paths['threshold']) as f:
        return float(f.readline())


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
