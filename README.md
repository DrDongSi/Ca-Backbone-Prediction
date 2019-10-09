## Citation
Paper citation with more information: Spencer Moritz, Jonas Pfab, Tianqi Wu, Jie Hou, Jianlin Cheng, Renzhi Cao, Liguo Wang, Dong Si. (2019). Cascaded-CNN: Deep Learning to Predict Protein Backbone Structure from High-Resolution Cryo-EM Density Maps. 10.1101/572990.

## Cα Backbone Prediction from Cryo-EM
Deep Learning for Cα Backbone Prediction from High Resolution CryoEM Data

## Installing Required Packages
In order to run the backbone prediction we need to install all required Python packages (Python version of 3.5 or higher is required). This can be done by creating a virtual environment with `python3 -m venv env` and activating it with `source ./env/bin/activate`. Once the virtual Python environment is activated, the required packages can be installed with pip using `pip install -r requirements.txt`.

Additionally, we need to have Chimera installed on the system and a symbolic link to the chimera binary file in `/usr/local/bin/chimera` must exist.

## Usage

The backbone prediction can be run by invoking the `main.py` script located in the root of the project. It requires two positional arguments:

* Input path
* Output path

> Paths must be passed as absolute paths

In addition to the input and output path we can also provide a JSON file containing threshold values for each protein using `-t THRESHOLD_FILE `. If this file is not provided threshold values are determined automatically which can lead to worse prediction results.

Another optional JSON file that can be passed is the hide dusts json file using `-d HIDEDUSTS_FILE`. If this file is not provided results can be slightly worse due to dust noise. 

An optional flag `-s` followed by a number `n` can be passed as an argument to skip the first `n` prediction steps.
> Skipping prediction steps is only possible if the results of the skipped steps are already available in the output path

Another optional flag `-c` can be set if you don't want to re-predict protein maps for which all/part of the results are already available in the output path. If set only prediction steps for which the results are not there yet are executed.

Another optional flag `-b` can be set if you want to keep mrc files for debugging purposes. Note: This will take up a lot more memory.

An example command to execute the prediction could therefore be the following.

`python main.py INPUT_PATH OUTPUT_PATH -t THRESHOLD_FILE`

### Sample Execution
A sample execution of the prediction will be demonstrated for the [6272 protein density map](https://www.emdataresource.org/EMD-6272) shown in the following image.

<img src="https://i.ibb.co/N7rS1Pn/6272mrc.png" alt="6272mrc" border="0">

The program will require the map as well as the fitted PDB as an input. Both files need to be stored in the same folder which is ideally named after the EMDB id. We also need to create a folder where the results will be stored. The folder structure could then look as following.

<img src="https://i.ibb.co/mzRg1XT/folder-structure.png" alt="folder-structure" border="0">

Furthermore, we want to provide a file containing the threshold value that should be used for the protein map. Therefore, we create a `thresholds.json` file with the following content:

```
{
	"6272": 0.009
}
```

> Note that the threshold value is read based on the folder name in which the density map is stored in

Now, we can start the prediction using the following command.

`python main.py input output -t thresholds.json`

Also, an option to apply hide dust sizes can be provided in a `hidedusts.json` file. Format is [contour level, hideDust size].

```
{
  "6272": [1.0, 10.0]
}
``` 

Additionally, we can start the prediction using the following command.

`python main.py input output -d hidedusts.json`

During the execution all prediction steps are run and the artifacts created by each step are stored in the output folder. In the following diagram we can see the flow of the different steps and the files they create.

<img src="https://i.ibb.co/w0DWq2M/prediction-pipeline.png" alt="prediction-pipeline" border="0">

> Note that the execution can take a significant amount of time

After the execution finished the folder structure should look as following.

<img src="https://i.ibb.co/XVKX1q7/folder-structure2.png" alt="folder-structure2" border="0">

The final prediction is stored in the **6272.pdb** file. The created backbone trace is shown in the following image.

<img src="https://i.ibb.co/nbnbtkQ/6272pdb.png" alt="6272pdb" border="0">

Additionally, a **results.xls** file is created in the output folder containing metrics about the prediction results.



## Sequence Mapping (Updating)
The code which maps protein sequences into the predicted Cα traces can be found [here](https://github.com/DrJieHou/CaTrace2Seq/). We put the paper's version in the the folder 'sequence_mapping', users can directly setup program in this folder. The latest version can be found at [here](https://github.com/DrJieHou/CaTrace2Seq/).

Test Environment
--------------------------------------------------------------------------------------
Red Hat Enterprise Linux Server release 6.4 (Santiago), perl 5, version 16, subversion 3 (v5.16.3)

Installation Steps
--------------------------------------------------------------------------------------

**(1) Configure CaTrace2Seq (required)**

```
cd sequence_mapping

perl setup_env.pl

(Note: since we use quality assessment tool with non-redundent sequence database, the tool requires around 33G)


```

**(2) Run CaTrace2Seq (required)**

```
sh bin/run_CaTrace2Seq.sh  <path of fasta sequence> <path of Ca trace> <length threshold for fragment> <output-directory> <num of cpus>
```

**(4) Practice the examples** 

```
cd example

sh run_6272.sh


 *****************************************************************************
 *                                 TM-SCORE                                  *
 * A scoring function to assess the quality of protein structure predictions *
 * Based on statistics:                                                      *
 *       0.0 < TM-score < 0.17, Random predictions                           *
 *       0.4 < TM-score < 1.00, Meaningful predictions                       *
 * Reference: Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710     *
 * For comments, please email to: yzhang@ku.edu                              *
 *****************************************************************************

Structure1: /home/jh7x  Length=  397
Structure2: 6272/pdb3j  Length=  397 (by which all scores are normalized)
Number of residues in common=  397
RMSD of  the common residues=    2.150

TM-score    = 0.9314  (d0= 7.20, TM10= 0.9314)
MaxSub-score= 0.7438  (d0= 3.50)
GDT-TS-score= 0.8086 %(d<1)=0.5819 %(d<2)=0.7431 %(d<4)=0.9093 %(d<8)=1.0000
GDT-HA-score= 0.6039 %(d<0.5)=0.1814 %(d<1)=0.5819 %(d<2)=0.7431 %(d<4)=0.9093

 -------- rotation matrix to rotate Chain-1 to Chain-2 ------
 i          t(i)         u(i,1)         u(i,2)         u(i,3)
 1     -0.0638134694   0.9999870448   0.0020705458   0.0046500556
 2      0.0280063468  -0.0020814379   0.9999950989   0.0023387491
 3      0.1425030895  -0.0046451903  -0.0023483976   0.9999864535
 
```

