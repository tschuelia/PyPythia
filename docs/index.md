# Home

Welcome to the PyPythia documentation.

Pythia is a lightweight python library to predict the difficulty of Multiple Sequence Alignments (MSA). Phylogenetic
analyzes under the Maximum-Likelihood (ML) model are time and resource intensive. To adequately capture the vastness of
tree space, one needs to infer multiple independent trees. On some datasets, multiple tree inferences converge to
similar tree topologies, on others to multiple, topologically highly distinct yet statistically indistinguishable
topologies. Pythia predicts the degree of difficulty of analyzing a dataset prior to initiating ML-based tree
inferences. Predicting the difficulty using Pythia is substantially faster than inferring multiple ML trees using
RAxML-NG. Pythia can be used to increase user awareness with respect to the amount of signal and uncertainty to be
expected in phylogenetic analyzes, and hence inform an appropriate (post-)analysis setup. Further, it can be used to
select appropriate search algorithms for easy-, intermediate-, and hard-to-analyze datasets. Pythia supports DNA, AA,
and morphological data in Phylip and FASTA format.

### Support

If you encounter any trouble using Pythia, have a question, or you find a bug, please feel free to open an
issue [here](https://github.com/tschuelia/PyPythia/issues).

### Publication

The paper explaining the details of Pythia is published in MBE:
Haag, J., Höhler, D., Bettisworth, B., & Stamatakis, A. (2022). **From Easy to Hopeless - Predicting the Difficulty of
Phylogenetic Analyses.** *Molecular Biology and Evolution*, 39(12). [https://doi.org/10.1093/molbev/msac254](https://doi.org/10.1093/molbev/msac254)

> [!WARNING]
> Since this publication, we made some considerable changes to Pythia.
> The most important change is that we switched from using a Random Forest Regressor to using a LightGBM Gradient
> Boosted Tree Regressor.
> This affects all Pythia versions >= 1. If you use Pythia in your work, please state the correct learning algorithm. If
> you are unsure, feel free to reach out to me 🙂
>
> We will soon publish a pre-print that explains the changes in detail, stay tuned!


### CPythia

The same functionality is also available as C library [here](https://github.com/tschuelia/difficulty_prediction).
Since the C library depends on [Coraxlib](https://codeberg.org/Exelixis-Lab/coraxlib) it is not as easy and fast to use
as this python library.
If you are only interested in the difficulty of your MSA, we recommend using this Python library.
If you want to incorporate the difficulty prediction in a phylogenetic tool, we recommend using the faster C library.
