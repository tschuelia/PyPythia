# Installing PyPythia

## Requirements

PyPythia uses RAxML-NG to calculate alignment features. Conda and pip installations therefore require a working
RAxML-NG executable on the system; see the upstream
[installation instructions](https://github.com/amkozlov/raxml-ng).

## Install using conda (recommended)

The easiest way to install a released version of PyPythia is from conda-forge:

```shell
conda install pythiaphylopredictor -c conda-forge -c nodefaults
```

## Install using pip

PyPythia is also available from PyPI:

```shell
pip install pythiaphylopredictor
```

Installing with pip into an existing environment can lead to dependency conflicts. Verify either installation with
`pythia -h`.

### Installing a specific tag/version

Install a tagged source version by including the tag in the Git URL, for example:

```shell
pip install git+https://github.com/tschuelia/PyPythia.git@0.0.1
```

For older conda releases, search for `pypythia` instead of `pythiaphylopredictor`.

## Installation from source

Pixi is the supported way to create a development environment. It installs the Python dependencies, development tools,
and RAxML-NG into a project-local, reproducible environment:

```shell
git clone https://github.com/tschuelia/PyPythia.git
cd PyPythia
pixi install
pixi run cli-smoke
```

Use `pixi shell` to work in the environment interactively. The default environment uses Python 3.14; the `py311` and
`py314` environments reproduce the minimum and maximum Python versions tested in CI.

## Troubleshooting

### LightGBM

Many pip installation problems are caused by a broken LightGBM installation. Refer to the
[LightGBM installation instructions](https://github.com/microsoft/LightGBM/tree/master/python-package) and install its
system prerequisites before retrying PyPythia. Pixi and conda installations obtain LightGBM and its native dependencies
from conda-forge.

### Python version

PyPythia supports Python 3.11 through Python 3.14. Select a locked Pixi environment explicitly with
`pixi run --environment py311 ...` or `pixi run --environment py314 ...`.

### Reproducing the development environment

The committed `pixi.lock` contains resolved dependencies for Linux, Intel macOS, and Apple Silicon macOS. Run
`pixi install --frozen` to reproduce it without changing dependency versions. Use `pixi update` deliberately when the
lock file should be refreshed, and commit the resulting lock-file change together with the manifest change.

### Running PyPythia

When working from source, run commands through Pixi, for example `pixi run pythia --help`. Pip and conda-forge
installations expose the `pythia` command directly.
