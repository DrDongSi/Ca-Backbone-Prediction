"""Coordinates the complete prediction pipeline

All prediction steps have to be added to the prediction pipeline.
"""

import sys
import os
from shutil import copyfile
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


def run_prediction(input_path, output_path, thresholds_file, num_skip):
    """Coordinates the execution of the prediction pipeline

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

    # Iterate over every directory in input path
    for emdb_id in filter(lambda d: os.path.isdir(input_path + d), os.listdir(input_path)):
        mrc_file = get_file(input_path + emdb_id, ['mrc', 'map'])
        gt_file = get_file(input_path + emdb_id, ['pdb', 'ent'])
        # Directory that contains paths to all relevant files. This will be
        # updated with every prediction step
        paths = {
            'input': input_path + emdb_id + '/' + mrc_file,
            'ground_truth': input_path + emdb_id + '/' + gt_file,
            'thresholds_file': thresholds_file
        }

        for prediction_step in PREDICTION_PIPELINE:
            paths['output'] = output_path + emdb_id + '/' + prediction_step.__name__.split('.')[0] + '/'
            os.makedirs(paths['output'], exist_ok=True)

            prediction_step.update_paths(paths)
            if num_skip <= 0:
                prediction_step.execute(paths)
            else:
                num_skip -= 1

        evaluator.calculate_accuracy(emdb_id, paths['traces_refined'], paths['ground_truth'])
        copyfile(paths['traces_refined'], output_path + emdb_id + '/' + emdb_id + '.pdb')

    evaluator.create_report(output_path)


def get_file(path, allowed_extensions):
    """Returns file in path with allowed extension"""
    return next(f for f in os.listdir(path) if f.split('.')[-1] in allowed_extensions)


if __name__ == '__main__':
    try:
        run_prediction(sys.argv[1] + ('/' if sys.argv[1][-1] != '/' else ''),
                       sys.argv[2] + ('/' if sys.argv[1][-1] != '/' else ''),
                       sys.argv[3],
                       0 if '-s' not in sys.argv else int(sys.argv[sys.argv.index('-s') + 1]))
    except IndexError:
        print('Missing argument')
    except ValueError:
        print('Invalid argument')
