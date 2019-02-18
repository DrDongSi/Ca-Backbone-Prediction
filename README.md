# Cα Backbone Prediction
Deep Learning for Cα Backbone Prediction from High Resolution CryoEM Data

<img src="https://i.ibb.co/pQbMPrW/sample-prediction.png" alt="sample-prediction" border="0">

## Installing Required Packages
In order to run the backbone prediction we need to install all required Python packages.
This can be done by creating a virtual environment with `python -m venv env` and activating it with `source ./env/bin/activate`. Once the virtual Python environment is activated, the required packages can be installed with pip using `pip install -r requirements.txt `. 

## Usage
The backbone prediction can be run by invoking the `prediction.py` script located in the `prediction` package. It requires three arguments:

* Input path
* Output path
* Threshold file

Additionally, we can pass an optional flag `-s` followed by a number `n` if we want to skip the first `n` prediction steps.
> Skipping prediction steps is only possible if the results of the skipped steps are already available in the output path

An example command to execute the prediction could therefore be the following.

`python prediction/prediction.py INPUT_PATH OUTPUT_PATH THRESHOLD_FILE`