import argparse
from prediction.prediction import run_predictions
from gui.gui import run_gui

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Ca Backbone Prediction from High Resolution CryoEM Data')
    parser.add_argument('input', type=str, help='Folder containing protein maps')
    parser.add_argument('output', type=str, help='Folder where prediction results will be stored')
    parser.add_argument('-t', '--thresholds', metavar='Thresholds', type=str,
                        help='JSON file which contains the thresholds')
    parser.add_argument('-s', '--skip', metavar='N', type=int, nargs=1, default=[0],
                        help='Number of prediction steps that should be skipped')
    parser.add_argument('-c', '--check_existing', action='store_const', const=True, default=False,
                        help='Check if results already exists and if so skip prediction step')
    parser.add_argument('-d', '--hidedusts', metavar='HideDusts', type=str,
                        help='JSON file which contains the hide dust sizes')
    parser.add_argument('-b', '--debug', action='store_const', const=True, default=False,
                        help='Enter debug mode, where mrc files are kept, otherwise mrc files are deleted at end to save memory')
    parser.add_argument('--nogui', action='store_const', const=True, default=False,
                        help='Run without GUI interface')

    args = parser.parse_args()

    args.input += '/' if args.input[-1] != '/' else ''
    args.output += '/' if args.output[-1] != '/' else ''

    if args.nogui:
        run_predictions(args.input, args.output, args.thresholds, args.skip[0], args.check_existing, args.hidedusts, args.debug)
    else:
        run_gui()
