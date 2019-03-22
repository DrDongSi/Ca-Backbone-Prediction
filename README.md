## Citation
Paper citation with more information: Spencer Moritz, Jonas Pfab, Tianqi Wu, Jie Hou, Jianlin Cheng, Renzhi Cao, Liguo Wang, Dong Si. (2019). Cascaded-CNN: Deep Learning to Predict Protein Backbone Structure from High-Resolution Cryo-EM Density Maps. 10.1101/572990.

## Cα Backbone Prediction from Cryo-EM
Deep Learning for Cα Backbone Prediction from High Resolution CryoEM Data

<img src="https://i.ibb.co/pQbMPrW/sample-prediction.png" alt="sample-prediction" border="0">

## Installing Required Packages
In order to run the backbone prediction we need to install all required Python packages (Python version of 3.5 or higher is required). This can be done by creating a virtual environment with `python3 -m venv env` and activating it with `source ./env/bin/activate`. Once the virtual Python environment is activated, the required packages can be installed with pip using `pip install -r requirements.txt`.

Additionally, we need to have Chimera installed on the system and a symbolic link to the chimera binary file in `/usr/local/bin/chimera` must exist.

## Usage
Before we can run the prediction we have to set the source root of the project using the following command.

`export PYTHONPATH="${PYTHONPATH}:/PATH_TO_REPOSITORY/Ca-Backbone-Prediction/`

The backbone prediction can be run by invoking the `prediction.py` script located in the `prediction` package. It requires two positional arguments:

* Input path
* Output path

> Paths must be passed as absolute paths

In addition to the input and output path we can also provide a JSON file containing threshold values for each protein using `-t THRESHOLD_FILE `. If this file is not provided threshold values are determined automatically which can lead to worse prediction results.

An optional flag `-s` followed by a number `n` can be passed as an argument to skip the first `n` prediction steps.
> Skipping prediction steps is only possible if the results of the skipped steps are already available in the output path

Another optional flag `-c` can be set if you don't want to re-predict protein maps for which all/part of the results are already available in the output path. If set only prediction steps for which the results are not there yet are executed.

An example command to execute the prediction could therefore be the following.

`python prediction/prediction.py INPUT_PATH OUTPUT_PATH -t THRESHOLD_FILE`

### Sample Execution
A sample execution of the prediction will be demonstrated for the [5778 protein density map](https://www.emdataresource.org/EMD-5778) shown in the following image.

<img src="https://i.ibb.co/XtbqC9z/density-map.png" alt="density-map" border="0">

The program will require the map as well as the fitted PDB as an input. Both files need to be stored in the same folder which is ideally named after the EMDB id. We also need to create a folder where the results will be stored. The folder structure could then look as following.

<img src="https://i.ibb.co/V9f3tJ4/folder-structure.png" alt="folder-structure" border="0">

Furthermore, we want to provide a file containing the threshold value that should be used for the protein map. Therefore, we create a `thresholds.json` file with the following content:

```
{
	"5778": 4.50
}
```

> Note that the threshold value is read based on the folder name in which the density map is stored in

Now, we can start the prediction using the following command.

`python prediction/prediction.py input output -t thresholds.json`

During the execution all prediction steps are run and the artifacts created by each step are stored in the output folder. In the following diagram we can see the flow of the different steps and the files they create.

<img src="https://i.ibb.co/w0DWq2M/prediction-pipeline.png" alt="prediction-pipeline" border="0">

> Note that the execution can take a significant amount of time

After the execution finished the folder structure should look as following.

<img src="https://i.ibb.co/nfRRgjH/folder-structure2.png" alt="folder-structure2" border="0">

The final prediction is stored in the **5778.pdb** file. The created backbone trace is shown in the following image.

<img src="https://i.ibb.co/hsq9VW5/prediction.png" alt="prediction" border="0">

Additionally, a **results.xls** file is created in the output folder containing metrics about the prediction results. In the case of the 5778 protein map this could e.g. look like this.

<img src="https://i.ibb.co/37HCVtm/results.png" alt="results" border="0">
