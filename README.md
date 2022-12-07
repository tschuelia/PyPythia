# PyPythia: Phylogenetic Difficulty Prediction Library

Lightweight python library to predict the difficulty of Multiple Sequence Alignments (MSA). 

## Installation and Requirements
In order to use this difficulty prediction, you need RAxML-NG installed somewhere on your system. 
You can find the install instructions [here](https://github.com/amkozlov/raxml-ng).

To install Pythia, run the following steps:

1. Clone this repository: `git clone https://github.com/tschuelia/PyPythia.git`  
Note: You can clone the repo to any location on your system, it does not need to be in the same directory as RAxML-NG 

2. `cd` into the directory:   
`cd PyPythia`

3. Checkout the latest release:
`git checkout tags/1.0.1`

4. Install the python package by running.
`pip install .`  

5. Verify the correct installation by running `pythia -h`. If you are having trouble running `pythia`, you can also replace `pythia` with `python pypythia/prediction.py`. 
For this to work you need to be in the PyPythia directory (which you should still be in at this point :-)).


Note: as of version 1.0.1, PyPythia includes a Python script that allows predictions from code without installing Pythia. See Section `Usage Without Installation` for details.


### Predictor
Per default this library uses the trained LightGBM boosted tree predictor `predictors/predictor_lgb_v1.0.0.pckl`. 
We will regularly retrain and update this predictor. You can also use an old version of the predictor by explicitly setting the path to the predictor to use (see usage below).
In future versions the set of features might change, so we do not guarantee backwards compatibility for old predictor versions unless explicitly stated. 

The predictor predicts the difficulty on a scale of 0.0 to 1.0. An MSA with a difficulty of 0.0 is predicted to be easy to analyze, 1.0 means the MSA is difficult. 

## Usage
This library can be used in two ways: either directly as command line tool, or the prediction can be called from other python code.

### Command Line Tool
If you only want to predict the difficulty for a single MSA, you can query the predictor using the command line interface, for example like this:
```commandline
pythia --msa examples/example.phy --raxmlng /path/to/raxml-ng
```
The output will be something like `The predicted difficulty for MSA examples/example.phy is: 0.08.`, telling us that example.phy is an easy dataset. In fact, this dataset exhibits a single likelihood peak.  
*Note that Pythia can also handle FASTA input files, see section Input Data below.*

The following options are available:
```commandline
PyPythia version 1.0.1 released by The Exelixis Lab
Developed by: Julia Haag
Latest version: https://github.com/tschuelia/PyPythia
Questions/problems/suggestions? Please open an issue on GitHub.

usage: pythia [-h] -m MSA -r RAXMLNG [-p PREDICTOR] [-o OUTPUT] [-prec PRECISION] [-sT] [--removeDuplicates] [-v] [-b]
              [-q]

Parser for optional config file setting.

options:
  -h, --help            show this help message and exit
  -m MSA, --msa MSA     Multiple Sequence Alignment to predict the difficulty for. Must be in either phylip or fasta
                        format.
  -r RAXMLNG, --raxmlng RAXMLNG
                        Path to the binary of RAxML-NG. For install instructions see
                        https://github.com/amkozlov/raxml-ng.
  -p PREDICTOR, --predictor PREDICTOR
                        Filepath of the predictor to use. If not set, assume it is
                        'predictors/predictor_lgb_v1.0.0.pckl' in the project directory.
  -o OUTPUT, --output OUTPUT
                        Option to specify a filepath where the result will be written to. The file will contain a
                        single line with only the difficulty.
  -prec PRECISION, --precision PRECISION
                        Set the number of decimals the difficulty should be rounded to. Recommended and default is 2.
  -sT, --storeTrees     If set, stores the parsimony trees as '{msa_name}.parsimony.trees' file.
  --removeDuplicates    Pythia refuses to predict the difficulty for MSAs containing duplicate sequences. If this
                        option is set, PyPythia removes the duplicate sequences, stores the reduced MSA as
                        '{msa_name}.{phy/fasta}.pythia.reduced' and predicts the difficulty for the reduced alignment.
  -v, --verbose         If set, additionally prints the MSA features.
  -b, --benchmark       If set, time the runtime of the prediction.
  -q, --quiet           If set, Pythia does not print progress updates and only prints the predicted difficulty.
```


### From Code
You can also use the library as a regular python library by installing it in your current environment with 
`pip install -e .` 
Then you can query the prediction like this:

```python
from pypythia.predictor import DifficultyPredictor
from pypythia.prediction import get_all_features
from pypythia.raxmlng import RAxMLNG
from pypythia.msa import MSA

predictor = DifficultyPredictor(open("pypythia/predictors/predictor_lgb_v1.0.0.pckl", "rb"))
raxmlng = RAxMLNG("/path/to/raxml-ng")
msa = MSA("examples/example.phy")

msa_features = get_all_features(raxmlng, msa)
difficulty = predictor.predict(msa_features)
print(difficulty)
```
*Note that Pythia can also handle FASTA input files, see section Input Data below.*

##### Using Python multiprocessing
There are reported issues with multiprocessing in Python and LightGBM based predictors (see for example the [LightGBM FAQ](https://lightgbm.readthedocs.io/en/latest/FAQ.html#lightgbm-hangs-when-multithreading-openmp-and-using-forking-in-linux-at-the-same-time)). 
We added a type check in the `predictor.py` prediction code that sets the number of threads to 1 for the prediction (`num_threads=1`) if the predictor is a LightGBM predictor. 
This should not affect the previous Pythia versions using the scikit-learn predictors. Since the multithreading issues do not occur consistently, this issue is hard to debug. 
If you encounter any issues with Python multiprocessing and Pythia please open a GitHub issue.

### Usage Without Installation
As of version 1.0.1, PyPythia includes a script `prediction_no_install.py` in the root directory. This script contains the single function `predict_difficulty`. 
Provided a path to an MSA, a path to a trained difficulty predictor (e.g. `pypythia/predictors/predictor_lgb_v1.0.0.pckl`), and a path to an executable of RAxML-NG, this fucntion
returns the predicted difficulty without requiring an installation of PyPythia. Note that this script can only be called from PyPythia's root directory.

### Input data
#### Supported file types
The input for Pythia is an MSA file in either Phylip or FASTA format. 

#### Supported  data types
Pythia supports DNA, AA, and morphological data.

Please note that by morphological data we refer to biological data. According to our analyses, the attributes of biological morphological data are similar to the attributes of DNA and AA data.
However, when analyzing language data (cognate, sound-class, and morphosyntactic data) we observed substantially distinct attributes and concluded that morphological language data is not 
comparable to DNA, AA, or biological morphological data. Thus, at the moment Pythia is not able to reliably predict the difficulty for language alignments.

#### Taxon names
Make sure that the MSA only contains RAxML-NG compatible taxon names. 
In particular, taxon labels with spaces, tabs, newlines, commas, colons, semicolons and parenthesis are invalid.

#### MSAs with duplicate sequences
As of version 1.0.0 Pythia refuses to predict the difficulty for MSAs containing multiple exactly identical sequences (duplicate sequences).
The reason for this is that duplicate sequences can have a substantial impact on the resulting topologies during the maximum parsimony tree inference 
and therefore on the topological distance measures. 

If you set the command line option `--removeDuplicates`, Pythia will create a reduced alignment with all duplicates removed and predict the difficulty for this reduced alignment.
For duplicate sequences, the first occurrence of the sequence is kept.
WARNING: The resulting predicted difficulty is only applicable to the reduced MSA! We recommend to only use the created reduced alignment for your subsequent analyses. 

## C Library
The same functionality is also available as C library [here](https://github.com/tschuelia/difficulty_prediction). 
Since the C library depends on [Coraxlib](https://codeberg.org/Exelixis-Lab/coraxlib) it is not as easy and fast to use as this python library.
If you are only interested in the difficulty of your MSA, we recommend using this Python library. 
If you want to incorporate the difficulty prediction in a phylogenetic tool, we recommend using the faster C library.


## Retraining
To continuously and automatically improve the prediction accuracy of Pythia, we regularly extend the training data set and subsequently retrain the predictor. 
We extend the training data using the anonymized MSAs that we continuously obtain during our RAxML Grove database updates. 
Note that these MSAs are only available internally in RAxML Grove and are not publicly available.

Older versions of the trained predictors are available in the `predictors` directory.


## Support
If you encounter any trouble using Pythia, have a question, or you find a bug, please feel free to open an issue here.


## Publication
The paper explaining the details of Pythia is published in MBE:    
Haag, J., Höhler, D., Bettisworth, B., & Stamatakis, A. (2022). **From Easy to Hopeless - Predicting the Difficulty of Phylogenetic Analyses.** *Molecular Biology and Evolution*, 39(12). [https://doi.org/10.1093/molbev/msac254](https://doi.org/10.1093/molbev/msac254)

## References
* A. M. Kozlov, D. Darriba, T. Flouri, B. Morel, and A. Stamatakis (2019) 
**RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference** 
*Bioinformatics*, 35(21): 4453–4455. 
[https://doi.org/10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305)

* D. Höhler, W. Pfeiffer, V. Ioannidis, H. Stockinger, A. Stamatakis (2022)
**RAxML Grove: an empirical phylogenetic tree database**
*Bioinformatics*, 38(6):1741–1742.
[https://doi.org/10.1093/bioinformatics/btab863](https://doi.org/10.1093/bioinformatics/btab863)
