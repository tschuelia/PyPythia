This library can be used in two ways: either directly as command line tool, or the prediction can be called from other python code.

## Command Line Tool
If you only want to predict the difficulty for a single MSA, you can query the predictor using the command line interface, for example like this:
```commandline
pythia --msa examples/example.phy --raxmlng /path/to/raxml-ng
```
Note that when you installed PyPythia using conda, you will have to download the `example.phy` and adjust the path accordingly.

The output will be something like `The predicted difficulty for MSA examples/example.phy is: 0.16.`, telling us that example.phy is an easy dataset. In fact, this dataset exhibits a single likelihood peak. Depending on the predictor version you are using, the actual value might slightly differ. This is expected and nothing to worry about 🙂

*Note that Pythia can also handle FASTA input files, see section Input Data below.*

The following options are available:
```commandline
PyPythia version 1.2.0 released by The Exelixis Lab
Developed by: Julia Haag
Latest version: https://github.com/tschuelia/PyPythia
Questions/problems/suggestions? Please open an issue on GitHub.

usage: pythia [-h] -m MSA -r RAXMLNG [-t THREADS] [-p PREDICTOR] [-o OUTPUT] [-prec PRECISION] [-sT] [--removeDuplicates] [--forceDuplicates]
              [--shap] [-v] [-b] [-q]

Parser for Pythia command line options.

options:
  -h, --help            show this help message and exit
  -m MSA, --msa MSA     Multiple Sequence Alignment to predict the difficulty for. Must be in either phylip or fasta format.
  -r RAXMLNG, --raxmlng RAXMLNG
                        Path to the binary of RAxML-NG. For install instructions see https://github.com/amkozlov/raxml-ng.
  -t THREADS, --threads THREADS
                        Number of threads to use for parallel parsimony tree inference. If none is set, Pythia uses the parallelization scheme
                        of RAxML-NG that automatically detects the optimal number of threads for your machine.
  -p PREDICTOR, --predictor PREDICTOR
                        Filepath of the predictor to use. If not set, assume it is 'predictors/latest.pckl' in the project directory.
  -o OUTPUT, --output OUTPUT
                        Option to specify a filepath where the result will be written to. The file will contain a single line with only the
                        difficulty.
  -prec PRECISION, --precision PRECISION
                        Set the number of decimals the difficulty should be rounded to. Recommended and default is 2.
  -sT, --storeTrees     If set, stores the parsimony trees as '{msa_name}.parsimony.trees' file.
  --removeDuplicates    Pythia refuses to predict the difficulty for MSAs containing duplicate sequences. If this option is set, PyPythia
                        removes the duplicate sequences, stores the reduced MSA as '{msa_name}.{phy/fasta}.pythia.reduced' and predicts the
                        difficulty for the reduced alignment.
  --forceDuplicates     Per default, Pythia refuses to predict the difficulty for MSAs containing duplicate sequences. Set this option if you
                        are absolutely sure that you want to predict the difficulty for this MSA.
  --shap                If set, computes the shapley values of the prediction as waterfall plot in '{msa_name}.shap.pdf'. When using this
                        option, make sure you understand what shapley values are and how to interpret this plot.For details on shapley values
                        refer to the wiki: https://github.com/tschuelia/PyPythia/wiki/Usage#shapley-values.
  -v, --verbose         If set, additionally prints the MSA features.
  -b, --benchmark       If set, time the runtime of the prediction.
  -q, --quiet           If set, Pythia does not print progress updates and only prints the predicted difficulty.
```


## From Code
You can also use the library as a regular python library by installing it in your current environment.
Then you can query the prediction like this:

```python
from pypythia.predictor import DifficultyPredictor
from pypythia.prediction import get_all_features
from pypythia.raxmlng import RAxMLNG
from pypythia.msa import MSA

predictor = DifficultyPredictor(open("pypythia/predictors/latest.pckl", "rb"))
raxmlng = RAxMLNG("/path/to/raxml-ng")
msa = MSA("examples/example.phy")

msa_features = get_all_features(raxmlng, msa)
difficulty = predictor.predict(msa_features)
print(difficulty)
```
*Note that Pythia can also handle FASTA input files, see section Input Data below.*

#### Using Python multiprocessing
There are reported issues with multiprocessing in Python and LightGBM based predictors (see for example the [LightGBM FAQ](https://lightgbm.readthedocs.io/en/latest/FAQ.html#lightgbm-hangs-when-multithreading-openmp-and-using-forking-in-linux-at-the-same-time)).
We added a type check in the `predictor.py` prediction code that sets the number of threads to 1 for the prediction (`num_threads=1`) if the predictor is a LightGBM predictor.
This should not affect the previous Pythia versions using the scikit-learn predictors. Since the multithreading issues do not occur consistently, this issue is hard to debug.
If you encounter any issues with Python multiprocessing and Pythia please open a GitHub issue.

## Usage Without Installation
As of version 1.0.1, PyPythia includes a script `prediction_no_install.py` in the root directory. This script contains the single function `predict_difficulty`.
Provided a path to an MSA, a path to a trained difficulty predictor (e.g. `pypythia/predictors/latest.pckl`), and a path to an executable of RAxML-NG, this fucntion
returns the predicted difficulty without requiring an installation of PyPythia. Note that this script can only be called from PyPythia's root directory.

To use this script, open it using your favorite text editor / python IDE and add the following at the end:
```python
msa_file = "path/to/your/msa"  # the file path of your MSA, can be either relative or absolute
raxmlng_exe_path = "path/to/raxml-ng/bin/raxml-ng"  # path pointing to the RAxML-NG executable on your system
predictor_path = "pypythia/predictors/latest.pckl"
predict_difficulty(msa_file, predictor_path, raxmlng_exe_path)
```

# Input data
### Supported file types
The input for Pythia is an MSA file in either Phylip or FASTA format.

### Supported  data types
Pythia supports DNA, AA, and morphological data.

Please note that by morphological data we refer to biological data. According to our analyses, the attributes of biological morphological data are similar to the attributes of DNA and AA data.
However, when analyzing language data (cognate, sound-class, and morphosyntactic data) we observed substantially distinct attributes and concluded that morphological language data is not
comparable to DNA, AA, or biological morphological data. Thus, at the moment Pythia is not able to reliably predict the difficulty for language alignments.

### Taxon names
Make sure that the MSA only contains RAxML-NG compatible taxon names.
In particular, taxon labels with spaces, tabs, newlines, commas, colons, semicolons and parenthesis are invalid.

### MSAs with duplicate sequences
As of version 1.0.0 Pythia refuses to predict the difficulty for MSAs containing multiple exactly identical sequences (duplicate sequences).
The reason for this is that duplicate sequences can have a substantial impact on the resulting topologies during the maximum parsimony tree inference
and therefore on the topological distance measures.

If you set the command line option `--removeDuplicates`, Pythia will create a reduced alignment with all duplicates removed and predict the difficulty for this reduced alignment.
For duplicate sequences, the first occurrence of the sequence is kept.
WARNING: The resulting predicted difficulty is only applicable to the reduced MSA! We recommend to only use the created reduced alignment for your subsequent analyses.


# Predictors
To continuously and automatically improve the prediction accuracy of Pythia, we regularly extend the training data set and subsequently retrain the predictor.
We extend the training data using the anonymized MSAs that we continuously obtain during our RAxML Grove database updates.
Note that these MSAs are only available internally in RAxML Grove and are not publicly available.
As per default, PyPythia uses the lastest predictor `predictors/latest.pckl`. Older versions of the trained predictors are available in the `predictors` directory and can be passed to Pythia (see Usage instructions above). All predictors of versions >= 1.0.0 are trained using DNA, AA, and morphological MSAs.

Note that the predictions for the same MSA can be different when using different versions of Pythia.


# Shapley Values
As of version 1.1.0, Pythia includes an option to plot Shapley values for a prediction. The interpretation of Shapley values is not straight-forward and we emphasize the importance of learning about these values before drawing conclusions based on the resulting plot!
We provide the Shapley values as waterfall plot.
In the following, we briefly describe what Shapley values are, what a waterfall plot is, and how you can interpret this plot.
It is important to note that Shapley values are not the same as feature importances. Predicting the difficulty of two distinct MSAs will lead to two distinct waterfall plots.

## Shapley Values
Based on the training data, our difficulty predictor Pythia has learned a base line difficulty. This base line difficulty is the expected value for every new prediction. Starting off this base line, Pythia adjusts its prediction using the features of the MSA. To determine how much each feature contributes to this change, ultimately leading to the final prediction is estimated by Shapley values. Since Pythia is a tree-based regressor, computing the Shapley values requires some advanced mathematics that I won't go into detail about here. If you are interested in this check out the links in the More Details section below. Due to the calculation of Shapley values, the value for one feature is NOT the difference in prediction when removing this feature. The Shapley value for one feature can only be interpreted considering all feature values together for a specific set of feature values.

## Waterfall plot
The following figure shows an exemplary waterfall plot output for the MSA `example/example.py` and Pythia version 1.1.0.



The x-axis depicts the difficulty and the y-axis the features alongside the respective feature value. The features are sorted by their Shapley value with the highest contribution on top. You can read the plot as follows. The base line difficulty that Pythia v1.1.0 learned is 0.35, as indicated by the `E[f(x)] = 0.35` on the x-axis. The `proportion_invariant` feature contributed to the overall prediction with a shift towards `1.0` (more difficult) of `0.01`, so *in combination with the other features*, a `proportion_invariant` of `0.341` indicates that the MSA is slightly more difficult than the average difficulty in the training set. We emphasize that the *combination with the other features* part, since the same value for `proportion_invariant` with a different MSA and different feature values for the remaining features might lead to a substantially different contribution to the overall prediction.
The feature with the highest impact for this example is the patterns-over-taxa ratio (`num_patterns/num_taxa`). The overall contribution is 0.23 towards `0.0`, meaning it shifts the overall prediction towards `easy`.

<img src="https://github.com/tschuelia/PyPythia/blob/master/examples/example.phy.shap.png" width="700">

## More Details
For further information please refer to [this great book on interpretable ML](https://christophm.github.io/interpretable-ml-book/shapley.html), the [documentation of the `shap` package](https://shap.readthedocs.io/en/latest/index.html), especially [their notes on the interpretability of Shapley values](https://shap.readthedocs.io/en/latest/example_notebooks/overviews/Be%20careful%20when%20interpreting%20predictive%20models%20in%20search%20of%20causal%C2%A0insights.html#Be-careful-when-interpreting-predictive-models-in-search-of-causal%C2%A0insights).
