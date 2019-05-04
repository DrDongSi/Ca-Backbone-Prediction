import argparse
from prediction.prediction import run_predictions

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

    args = parser.parse_args()

    args.input += '/' if args.input[-1] != '/' else ''
    args.output += '/' if args.output[-1] != '/' else ''

    run_predictions(args.input, args.output, args.thresholds, args.skip[0], args.check_existing)
