# PyPythia: Phylogenetic Difficulty Prediction Library

Lightweight python library to predict the difficulty of Multiple Sequence Alignments (MSA). 

## Installation and Requirements
In order to use this difficulty prediction, you need RAxML-NG installed somewhere on your system. 
You can find the install instructions [here](https://github.com/amkozlov/raxml-ng).

Install the python package by running
    ```
    pip install -e .
    ```

### Predictor
Per default this library uses the trained scikit-learn random forest predictor `predictor.pckl`. 
We will regularly retrain and update this predictor. You can also use an old version of the predictor by explicitly setting the path to the predictor to use (see usage below).
In future versions the set of features might change, so we do not guarantee backwards compatibility for old predictor versions. 

The predictor predicts the difficulty on a scale of 0.0 to 1.0. An MSA with a difficulty of 0.0 is predicted to be easy to analyze, 1.0 means the MSA is difficult. 

## Usage
This library can be used in two ways: either directly as command line tool, or the prediction can be called from other python code.

### Command Line Tool
If you only want to predict the difficulty for a single MSA, you can query the predictor using the command line interface, for example like this:
```commandline
pythia --msa examples/exmple.phy --raxmlng /path/to/raxml-ng
```
The output will be something like `The predicted difficulty for MSA examples/example.phy is: 0.12.`, telling us that example.phy is an easy dataset. In fact, this dataset exhibits a single likelihood peak.

The following options are available:
```commandline
usage: pythia [-h] --msa MSA [--model MODEL] --raxmlng RAXMLNG [--predictor PREDICTOR] [--storeTrees] [--verbose]

Parser for optional config file setting.

options:
  -h, --help            show this help message and exit
  --msa MSA             Multiple Sequence Alignment to predict the difficulty for. Must be in either phylip or fasta format.
  --model MODEL         Model to use for the prediction. This can be either a model string (e.g. GTR+G) or a path to a partition file.If not set the data type is automatically inferred, and the model is set to GTR+G for DNA MSAs and
                        to LG+G for Protein MSAs.
  --raxmlng RAXMLNG     Path to the binary of RAxML-NG. For install instructions see https://github.com/amkozlov/raxml-ng.
  --predictor PREDICTOR
                         Filepath of the predictor to use. If not set, assume it is 'predictor.pckl' in the project directory.
  --storeTrees          If set, stores the parsimony trees as '{msa_name}.parsimony.trees' file
  --verbose             If set, prints the MSA features
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

predictor = DifficultyPredictor(open("pypythia/predictor.pckl", "rb"))
raxmlng = RAxMLNG("/path/to/raxml-ng")
msa = MSA("examples/example.phy")
model = "GTR+G"

msa_features = get_all_features(raxmlng, msa, model)
difficulty = predictor.predict(msa_features)
print(difficulty)
```


## C Library
The same functionality is also available as C library [here](https://github.com/tschuelia/difficulty_prediction). 
Since the C library depends on [Coraxlib](https://codeberg.org/Exelixis-Lab/coraxlib) it is not as easy and fast to use as this python library.
If you are only interested in the difficulty of your MSA, we recommend using this Python library. 
If you want to incorporate the difficulty prediction in a phylogenetic tool, we recommend using the faster C library.

## References
* A. M. Kozlov, D. Darriba, T. Flouri, B. Morel, and A. Stamatakis (2019) 
**RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference** 
*Bioinformatics*, 35(21): 4453–4455. 
doi: [10.1093/bioinformatics/btz305](http://dx.doi.org/10.1093/bioinformatics/btz305)
