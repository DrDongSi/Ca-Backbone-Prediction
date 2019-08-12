"""Coordinates the complete prediction pipeline

All prediction steps have to be added to the prediction pipeline.
"""

import os
import sys
from shutil import copyfile
from multiprocessing import cpu_count, Pool, Semaphore
from time import time
import traceback
from .evaluation import Evaluator
import preprocessing as pre
import cnn
import postprocessing as post


# List contains every prediction step that is executed in order to produce
# the final prediction
PREDICTION_PIPELINE = [
    pre.clean_map,
    pre.find_threshold,
    pre.normalize_map,
    cnn.predict_with_module,
    post.build_backbone_trace,
    post.merge_chains,
    post.helix_refinement
]


def run_predictions(input_path, output_path, thresholds_file, num_skip, check_existing, hidedusts_file, debug):
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

    check_existing: bool
        If set prediction steps are only executed if their results are not
        existing in the output path yet
    """
    # Create list of parameters for every prediction
    params_list = [(emdb_id, input_path, output_path, thresholds_file, num_skip, check_existing, hidedusts_file, debug)
                   for emdb_id in filter(lambda d: os.path.isdir(input_path + d), os.listdir(input_path))]

    start_time = time()
    max_processes_allowed_to_access_tensorflow = 4
    semaphore = Semaphore(min(min(cpu_count(), len(params_list)), max_processes_allowed_to_access_tensorflow))
    pool = Pool(min(cpu_count(), len(params_list)), initializer=init_child, initargs=(semaphore,))
    results = pool.map(run_prediction, params_list)

    # Filter 'None' results
    results = filter(lambda r: r is not None, results)

    evaluator = Evaluator(input_path)
    for emdb_id, predicted_file, gt_file, execution_time in results:
        evaluator.evaluate(emdb_id, predicted_file, gt_file, execution_time)

    evaluator.create_report(output_path, time() - start_time)

def init_child(semaphore_):
    global semaphore
    semaphore = semaphore_


def run_prediction(params):
    """Coordinates the execution of every prediction step in the prediction
    pipeline

    Parameters
    ----------
    params: tuple
        Tuple of parameters required for the prediction. They are unpacked at
        the beginning of the method

    Returns
    ----------
    result: tuple
        Result as tuple containing the emdb id, predicted file, ground truth
        file, and execution time respectively
    """
    # Unpack parameters
    emdb_id, input_path, output_path, thresholds_file, num_skip, check_existing, hidedusts_file, debug = params
    paths = make_paths(input_path, emdb_id, thresholds_file, hidedusts_file)

    start_time = time()
    for prediction_step in PREDICTION_PIPELINE:
        paths['output'] = output_path + emdb_id + '/' + prediction_step.__name__.split('.')[0] + '/'
        os.makedirs(paths['output'], exist_ok=True)

        try:
            prediction_step.update_paths(paths)
            if num_skip > 0 or (check_existing and files_exist(paths)):
                num_skip -= 1
            else:
                prediction_step.execute(paths)
        except BaseException:
            exc_info = sys.exc_info()
            traceback.print_exception(*exc_info)

            return None

    if not os.path.isfile(paths['traces_refined']):
        return None

    if 'traces_refined' in paths:
        copyfile(paths['traces_refined'], output_path + emdb_id + '/' + emdb_id + '.pdb')

    if debug is False:
        try:
            os.remove(paths['cleaned_map'])
            os.remove(paths['normalized_map'])
            os.remove(paths['loops_confidence'])
            os.remove(paths['sheet_confidence'])
            os.remove(paths['helix_confidence'])
            os.remove(paths['backbone_confidence'])
            os.remove(paths['ca_confidence'])
        except:
            pass

    return emdb_id, paths['traces_refined'], paths['ground_truth'], time() - start_time


def make_paths(input_path, emdb_id, thresholds_file, hidedusts_file):
    """Creates base paths dictionary with density map, ground truth, and
    optionally the thresholds file"""
    mrc_file = get_file(input_path + emdb_id, ['mrc', 'map'])
    gt_file = get_file(input_path + emdb_id, ['pdb', 'ent'])
    # Directory that contains paths to all relevant files. This will be
    # updated with every prediction step
    paths = {
        'input': input_path + emdb_id + '/' + mrc_file,
        'ground_truth': input_path + emdb_id + '/' + gt_file,
    }

    if thresholds_file is not None:
        paths['thresholds_file'] = thresholds_file

    if hidedusts_file is not None:
        paths['hidedusts_file'] = hidedusts_file

    return paths


def files_exist(paths):
    """Checks if all files specified in the 'paths' dict exist"""
    for path in paths.values():
        if not os.path.isdir(path) and not os.path.isfile(path):
            return False

    return True


def get_file(path, allowed_extensions):
    """Returns file in path with allowed extension"""
    return next(f for f in os.listdir(path) if f.split('.')[-1] in allowed_extensions)
