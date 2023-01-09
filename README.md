# PyPythia: Phylogenetic Difficulty Prediction Library
Pythia is a lightweight python library to predict the difficulty of Multiple Sequence Alignments (MSA). 
Pythia supports DNA, AA, and morphological data in phylip and FASTA format.

### Documentation
For detailed instructions on how to install and use Pythia see [the wiki](https://github.com/tschuelia/PyPythia/wiki).

### Support
If you encounter any trouble using Pythia, have a question, or you find a bug, please feel free to open an issue here.


### Publication
The paper explaining the details of Pythia is published in MBE:    
Haag, J., Höhler, D., Bettisworth, B., & Stamatakis, A. (2022). **From Easy to Hopeless - Predicting the Difficulty of Phylogenetic Analyses.** *Molecular Biology and Evolution*, 39(12). [https://doi.org/10.1093/molbev/msac254](https://doi.org/10.1093/molbev/msac254)

### C Library
The same functionality is also available as C library [here](https://github.com/tschuelia/difficulty_prediction). 
Since the C library depends on [Coraxlib](https://codeberg.org/Exelixis-Lab/coraxlib) it is not as easy and fast to use as this python library.
If you are only interested in the difficulty of your MSA, we recommend using this Python library. 
If you want to incorporate the difficulty prediction in a phylogenetic tool, we recommend using the faster C library.

### References
* A. M. Kozlov, D. Darriba, T. Flouri, B. Morel, and A. Stamatakis (2019) 
**RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference** 
*Bioinformatics*, 35(21): 4453–4455. 
[https://doi.org/10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305)

* D. Höhler, W. Pfeiffer, V. Ioannidis, H. Stockinger, A. Stamatakis (2022)
**RAxML Grove: an empirical phylogenetic tree database**
*Bioinformatics*, 38(6):1741–1742.
[https://doi.org/10.1093/bioinformatics/btab863](https://doi.org/10.1093/bioinformatics/btab863)
