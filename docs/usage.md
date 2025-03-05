# Using PyPythia

This library can be used in two ways: either directly as command line tool, or the prediction can be called from other
python code.

## Command Line Interface

If you only want to predict the difficulty for a single MSA, you can query the predictor using the command line
interface, for example like this:

```commandline
pythia --msa examples/example.phy --raxmlng /path/to/raxml-ng
```

Note that when you installed PyPythia using conda, you will have to download the `example.phy` and adjust the path
accordingly.

The output will be something like `The predicted difficulty for MSA examples/example.phy is: 0.02.`, telling us that
example.phy is an easy dataset. In fact, this dataset exhibits a single likelihood peak. Depending on the predictor
version or operating system you are using, the actual value might slightly differ.
This is expected and nothing to worry about 🙂

*Note that Pythia can also handle FASTA input files, see section Input Data below.*

The following options are available:

```commandline
PyPythia version 2.0.0 released by The Exelixis Lab
Developed by: Julia Haag
Latest version: https://github.com/tschuelia/PyPythia
Questions/problems/suggestions? Please open an issue on GitHub.

usage: pythia [-h] -m MSA -r RAXMLNG [-t THREADS] [-s SEED] [-p PREFIX]
              [--predictor PREDICTOR] [--shap] [--forceDuplicates] [--forceFullGaps]
              [--nofiles] [-V]

Parser for Pythia command line options.

options:
  -h, --help            show this help message and exit
  -m MSA, --msa MSA     Multiple Sequence Alignment to predict the difficulty for. Must be
                        in either phylip or fasta format.
  -r RAXMLNG, --raxmlng RAXMLNG
                        Path to the binary of RAxML-NG. For install instructions see
                        https://github.com/amkozlov/raxml-ng.(default: 'raxml-ng' if in
                        $PATH, otherwise this option is mandatory).
  -t THREADS, --threads THREADS
                        Number of threads to use for parallel parsimony tree inference
                        (default: RAxML-NG autoconfig).
  -s SEED, --seed SEED  Seed for the RAxML-NG parsimony tree inference (default: 0).
  -p PREFIX, --prefix PREFIX
                        Prefix of the PyPythia log and result file (default: MSA file name).
  --predictor PREDICTOR
                        Filepath of the alternative predictor to use (default: latest
                        Pythia).
  --shap                If set, computes the shapley values of the prediction as waterfall
                        plot in '{prefix}.shap.pdf'. When using this option, make sure you
                        understand what shapley values are and how to interpret this
                        plot.For details on shapley values refer to the documentation:
                        https://tschuelia.github.io/PyPythia/latest/usage/#shap-waterfall-
                        plot (default: False).
  --forceDuplicates     Per default, Pythia refuses to predict the difficulty for MSAs
                        containing duplicate sequences,and removes duplicate sequences prior
                        to predicting the difficulty. Only set this option if you are
                        absolutely sure that you want to predict the difficulty for this MSA
                        (default: False).
  --forceFullGaps       Per default, Pythia refuses to predict the difficulty for MSAs
                        containing sequences with only gaps,and removes full-gap sequences
                        prior to predicting the difficulty. Only set this option if you are
                        absolutely sure that you want to predict the difficulty for this MSA
                        (default: False).
  --nofiles             Prevent Pythia from writing any files and only print logs/results to
                        the terminal (default: False). WARNING: in this case and if your MSA
                        contains duplicate/full-gap sequences the reduced MSA will not be
                        stored.
  -V, --version         Print the version number and exit.
```

### Result files

Pythia will write the following files:
- A logfile containing the same information as printed to the terminal: `{result_prefix}.pythia.log`
- The reduced MSA file in case the input MSA contained duplicate/full-gap sequences (and the reduction was not disabled): `{result_prefix}.reduced.phy`
- The inferred parsimony trees in Newick format: `{result_prefix}.pythia.trees`
- The shapley values as waterfall plot (if --shap is set): `{result_prefix}.shap.pdf`
- The features and predicted difficulty as CSV file: `{result_prefix}.pythia.csv`

The result_prefix can be set using the `--prefix` command line option. If not set, Pythia uses the MSA file as prefix. You can prevent Pythia from writing any files via the flag `--nofiles`.

## From Code

You can also use the library as a regular python library by installing it in your current environment.
The following code snippet shows how to predict the difficulty for an MSA using PyPythia:

```python
from pypythia.prediction import predict_difficulty
import pathlib

msa = pathlib.Path("examples/example.phy")
difficulty = predict_difficulty(msa)
print(f"The predicted difficulty for MSA {msa} is: {round(difficulty, 2)}.")
```

And the output will be the same as for the CLI: `The predicted difficulty for MSA examples/example.phy is: 0.02.`.

If you want to get all features, or do more specific analyses of your MSA, see the API Reference for further details on
all available classes and methods.

## Input data

### Supported file types

The input for Pythia is an MSA file in either Phylip or FASTA format.

### Supported  data types

Pythia supports DNA, AA, and morphological data.

Please note that by morphological data we refer to biological data. According to our analyses, the attributes of
biological morphological data are similar to the attributes of DNA and AA data.
However, when analyzing language data (cognate, sound-class, and morphosyntactic data) we observed substantially
distinct attributes and concluded that morphological language data is not
comparable to DNA, AA, or biological morphological data. Thus, at the moment Pythia is not able to reliably predict the
difficulty for language alignments.

### Taxon names

Make sure that the MSA only contains RAxML-NG compatible taxon names.
In particular, taxon labels with spaces, tabs, newlines, commas, colons, semicolons and parenthesis are invalid.

### MSAs with duplicate sequences

Pythia refuses to predict the difficulty for MSAs containing duplicate sequences or MSAs containing sequences containing
only gaps.
As of version 2.0.0, Pythia removes duplicates and full-gap sequences per default and predicts the difficulty for this
reduced MSA.
If you absolutely want to predict the difficulty for the original MSA, set the command line flags `--forceDuplicates`
and `--forceFullGaps`.

As of version 1.0.0 Pythia refuses to predict the difficulty for MSAs containing multiple exactly identical sequences (
duplicate sequences).
The reason for this is that duplicate sequences can have a substantial impact on the resulting topologies during the
maximum parsimony tree inference
and therefore on the topological distance measures.

If you set the command line option `--removeDuplicates`, Pythia will create a reduced alignment with all duplicates
removed and predict the difficulty for this reduced alignment.
For duplicate sequences, the first occurrence of the sequence is kept.
WARNING: The resulting predicted difficulty is only applicable to the reduced MSA! We recommend to only use the created
reduced alignment for your subsequent analyses.

## Predictors

To continuously and automatically improve the prediction accuracy of Pythia, we regularly extend the training data set
and subsequently retrain the predictor.
We extend the training data using the anonymized MSAs that we continuously obtain during our RAxML Grove database
updates.
Note that these MSAs are only available internally in RAxML Grove and are not publicly available.
As per default, PyPythia uses the latest predictor `predictors/latest.txt`. If you want to use an older version of Pythia,
please install the respective PyPythia version.
You can also pass a custom predictor file using the `--predictor` option. However, this will only work if the passed
file contains a LightGBM Booster model.

Note that the predictions for the same MSA can be different when using different versions of Pythia.

## SHAP Waterfall Plot

As of version 1.1.0, Pythia includes an option to plot Shapley values for a prediction. The interpretation of Shapley
values is not straight-forward, and we emphasize the importance of learning about these values before drawing conclusions
based on the resulting plot!
We provide the Shapley values as waterfall plot.
In the following, we briefly describe what Shapley values are, what a waterfall plot is, and how you can interpret this
plot.
It is important to note that Shapley values are not the same as feature importances. Predicting the difficulty of two
distinct MSAs will lead to two distinct waterfall plots.

### Shapley Values

Based on the training data, our difficulty predictor Pythia has learned a base line difficulty. This base line
difficulty is the expected value for every new prediction. Starting off this base line, Pythia adjusts its prediction
using the features of the MSA. To determine how much each feature contributes to this change, ultimately leading to the
final prediction is estimated by Shapley values. Since Pythia is a tree-based regressor, computing the Shapley values
requires some advanced mathematics that I won't go into detail about here. If you are interested in this check out the
links in the More Details section below. Due to the calculation of Shapley values, the value for one feature is NOT the
difference in prediction when removing this feature. The Shapley value for one feature can only be interpreted
considering all feature values together for a specific set of feature values.

### Waterfall plot

The following figure shows an exemplary waterfall plot output for the MSA `example/example.py` and Pythia version 1.1.0.

The x-axis depicts the difficulty and the y-axis the features alongside the respective feature value. The features are
sorted by their Shapley value with the most important feature on top. You can read the plot as follows. The baseline
difficulty that Pythia v1.1.0 learned is 0.35, as indicated by the `E[f(x)] = 0.35` on the x-axis. The
`proportion_invariant` feature contributed to the overall prediction with a shift towards `1.0` (more difficult) of
`0.01`, so *in combination with the other features*, a `proportion_invariant` of `0.341` indicates that the MSA is
slightly more difficult than the average difficulty in the training set. We emphasize that the *combination with the
other features* part, since the same value for `proportion_invariant` with a different MSA and different feature values
for the remaining features might lead to a substantially different contribution to the overall prediction.
The feature with the highest impact for this example is the patterns-over-taxa ratio (`num_patterns/num_taxa`). The
overall contribution is 0.23 towards `0.0`, meaning it shifts the overall prediction towards `easy`.

<img src="../img/example.phy.shap.pdf" width="700">

### More Details

For further information please refer
to [this great book on interpretable ML](https://christophm.github.io/interpretable-ml-book/shapley.html),
the [documentation of the `shap` package](https://shap.readthedocs.io/en/latest/index.html),
especially [their notes on the interpretability of Shapley values](https://shap.readthedocs.io/en/latest/example_notebooks/overviews/Be%20careful%20when%20interpreting%20predictive%20models%20in%20search%20of%20causal%C2%A0insights.html#Be-careful-when-interpreting-predictive-models-in-search-of-causal%C2%A0insights).
