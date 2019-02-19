"""Coordinates the complete prediction pipeline

All prediction steps have to be added to the prediction pipeline.
"""

import sys
import os
from shutil import copyfile
from multiprocessing import cpu_count
from multiprocessing.dummy import Pool as ThreadPool
from time import time
from evaluation import Evaluator
import preprocessing as pre
import cnn
import postprocessing as post


# List contains every prediction step that is executed in order to produce
# the final prediction
PREDICTION_PIPELINE = [
    pre.clean_map,
    pre.normalize_map,
    cnn.predict_with_module,
    post.build_backbone_trace,
    post.helix_refinement
]


def run_predictions(input_path, output_path, thresholds_file, num_skip):
    """Creates thread pool which will concurrently run the prediction for every
    protein map in the 'input_path'

    Parameters
    ----------
    input_path: str
        Path of the input directory where the different protein directories are
        located

    output_path: str
        Path of the folder where all generated files will be stored

    thresholds_file: str
        Path of the JSON file which contains the threshold values for the input
        files

    num_skip: int
        The number of prediction steps that should be skipped
    """
    evaluator = Evaluator(input_path)

    # Create list of arguments for every prediction
    args_list = [(input_path, output_path, thresholds_file, num_skip, evaluator, emdb_id)
                 for emdb_id in filter(lambda d: os.path.isdir(input_path + d), os.listdir(input_path))]

    start_time = time()
    pool = ThreadPool(cpu_count())
    pool.map(run_prediction, args_list)

    evaluator.create_report(output_path, time() - start_time)

def run_prediction(args):
    """Coordinates the execution of every prediction step in the prediction
    pipeline

    Parameters
    ----------
    args: tuple
        Tuple of arguments required for the prediction. They are unpacked at
        the beginning of the method
    """
    # Unpack parameters
    input_path, output_path, thresholds_file, num_skip, evaluator, emdb_id = args

    mrc_file = get_file(input_path + emdb_id, ['mrc', 'map'])
    gt_file = get_file(input_path + emdb_id, ['pdb', 'ent'])
    # Directory that contains paths to all relevant files. This will be
    # updated with every prediction step
    paths = {
        'input': input_path + emdb_id + '/' + mrc_file,
        'ground_truth': input_path + emdb_id + '/' + gt_file,
        'thresholds_file': thresholds_file
    }

    start_time = time()
    for prediction_step in PREDICTION_PIPELINE:
        paths['output'] = output_path + emdb_id + '/' + prediction_step.__name__.split('.')[0] + '/'
        os.makedirs(paths['output'], exist_ok=True)

        prediction_step.update_paths(paths)
        if num_skip <= 0:
            prediction_step.execute(paths)
        else:
            num_skip -= 1

    evaluator.evaluate(emdb_id, paths['traces_refined'], paths['ground_truth'], time() - start_time)
    copyfile(paths['traces_refined'], output_path + emdb_id + '/' + emdb_id + '.pdb')


def get_file(path, allowed_extensions):
    """Returns file in path with allowed extension"""
    return next(f for f in os.listdir(path) if f.split('.')[-1] in allowed_extensions)


if __name__ == '__main__':
    try:
        run_predictions(sys.argv[1] + ('/' if sys.argv[1][-1] != '/' else ''),
                        sys.argv[2] + ('/' if sys.argv[1][-1] != '/' else ''),
                        sys.argv[3],
                        0 if '-s' not in sys.argv else int(sys.argv[sys.argv.index('-s') + 1]))
    except IndexError:
        print('Missing argument')
    except ValueError:
        print('Invalid argument')
