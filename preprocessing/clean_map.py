"""Cleans input map

The module is responsible for cleaning the input map and re-sampling it to a
voxel size of 1. This is accomplished by creating a script of chimera commands
and then pass it to chimera in no GUI mode.

The cleaning can be applied on a mrc file by calling the 'resample_map'
method
"""

import subprocess
import os
from threading import Lock

lock = Lock()

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
    lock.acquire()

    chimera_script = open('resample.cmd', 'w')
    chimera_script.write('open ' + paths['input'] + '\n'
                         'open ' + paths['ground_truth'] + '\n'
                         'molmap #1 6 gridSpacing 1\n'
                         'volume #1.1 level 0.1\n'
                         'mask #0 #1.1\n'
                         'vop resample #2 onGrid #1.1\n'
                         'volume #3 save ' + paths['cleaned_map'])
    chimera_script.close()

    subprocess.run(['/usr/local/bin/chimera', '--nogui', chimera_script.name])
    os.remove('resample.cmd')

    lock.release()
