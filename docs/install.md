# Installing PyPythia

## Requirements

In order to use this difficulty prediction, you need RAxML-NG installed somewhere on your system. You can find the
installation instructions [here](https://github.com/amkozlov/raxml-ng).

## Install using conda (recommended)

The easiest (and recommended) way to install PyPythia is by using conda:

```
conda install pypythia -c conda-forge -c nodefaults
```

## Install using pip

You can also install Pythia using the python package manager pip:

```
pip install pythiaphylopredictor
```

Please note that this can lead to issues with package versions and dependencies when installing in an existing (conda)
environment.

Verify the correct installation by running `pythia -h`.

### Installing a specific tag/version

You can again use pip for this and simply specify the tag you wish to install, e.g. for version `0.0.1` run:

```
pip install git+https://github.com/tschuelia/PyPythia.git@0.0.1
```

## Installation from source

You can install Pythia from source if you want to explore the code or get the lastest development version.
To do so run the following steps:

```
git clone https://github.com/tschuelia/PyPythia.git
cd PyPythia
pip install .
```

Verify the correct installation by running `pythia -h`.

## Troubleshooting

Most issues when installing Pythia seem to arise from broken or non-working LightGBM installations. If you encounter any
such problem, and none of the following options help, please refer to the
LightGBM [installation instructions](https://github.com/microsoft/LightGBM/tree/master/python-package) for your
operating system and install LightGBM manually _before_ repeating the Pythia installation as described above.

### Python version

Since Pythia version 1.1.0 we provide the option to output the Shapley values for your prediction. Currently, the `shap`
package does not support Python Version 3.11. The requirements should take care of the correct Python version, but if
you encounter any issues, please first check that the Python version is <3.11. You can do so by typing
`python --version` in your terminal and checking the output.

### Installing on M1 chips

Installing on MacBooks with M1 chips caused some trouble for some users that seem to be caused by LightGBM's
multiprocessing support. If you encounter any errors with the log pointing to LightGBM, the first thing you could try is
to install LightGBM using [homebrew](https://brew.sh/index):

```
brew install lightgbm
```

This might take a few minutes to finish.
Once this ran successfully you can try to rerun the install instructions above.

If this does not solve your problem, you can try to install LightGBM manually using pip and disabling the
multiprocessing:

```
pip install lightgbm --install-option=--nomp
```

and then rerun the installation of PyPythia. Thanks [@willbour](https://github.com/willbour) for finding the fix for
this!

### Using a clean conda environment

When using conda and installing PyPythia using pip in an existing environment, you might encounter dependency or version
related issues. To check whether this is the case or you have a general issue with Pythia please try to create a new,
clean conda environment:

1. Use the provided environment file `etc/environment.yml` and create a new conda environment:

```
conda env create --file etc/environment.yml
```

If you want to install a different version of Pythia, you can add the git tag by appending `@[version]` (e.g. for
version 1.1.0 append `@1.1.0`) after `.git` in the `etc/environment.yml` file.

2. Activate the conda environment: `conda activate pythia`
3. Try to (re)run Pythia.

### Running Pythia

If you are having trouble running pythia, you can also replace `pythia` with `python pypythia/main.py`. For this
to work you need to install Pythia from source and you need to be in the PyPythia directory (which you should be after
the installation).
