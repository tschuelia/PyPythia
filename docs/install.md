# Installing PyPythia

## Requirements

PyPythia uses RAxML-NG to calculate alignment features. Installations of a released PyPythia package with global
Pixi, Conda, or pip therefore require a working RAxML-NG executable on the system; see the upstream
[RAxML-NG installation instructions](https://github.com/amkozlov/raxml-ng).

## Pixi (recommended)

First [install Pixi](https://pixi.prefix.dev/latest/installation/), then install PyPythia globally from conda-forge
and verify that its command-line interface is available:

```shell
pixi global install pythiaphylopredictor --channel conda-forge
pythia -h
```

Remember to install RAxML-NG separately, as described in [Requirements](#requirements).

## Conda (alternative)

PyPythia is also available from conda-forge:

```shell
conda install pythiaphylopredictor -c conda-forge -c nodefaults
pythia -h
```

Remember to install RAxML-NG separately, as described in [Requirements](#requirements).

## pip (alternative)

PyPythia is available from PyPI:

```shell
pip install pythiaphylopredictor
pythia -h
```

Installing with pip into an existing environment can lead to dependency conflicts. Remember to install RAxML-NG
separately, as described in [Requirements](#requirements).

### Installing a specific tag/version

Install a tagged source version by including the tag in the Git URL, for example:

```shell
pip install git+https://github.com/tschuelia/PyPythia.git@0.0.1
```

For older conda releases, search for `pypythia` instead of `pythiaphylopredictor`.

## Local development from source

Pixi is the supported way to create a development environment. It installs the Python dependencies, development
tools, and the pinned RAxML-NG version into project-local, reproducible environments.

Clone the repository, install the default environment, and verify the CLI:

```shell
git clone https://github.com/tschuelia/PyPythia.git
cd PyPythia
pixi install
pixi run cli-smoke
```

To work in the environment interactively, start a Pixi shell. Commands such as `pythia` and `pytest` are then
available directly:

```shell
pixi shell
pythia -h
```

Alternatively, you can prepend every command with `pixi run`. For example `pixi run pythia -h`.

Run the test suite through its Pixi task:

```shell
pixi run test
```

The default environment uses Python 3.14. The `py311` and `py314` environments reproduce the minimum and maximum
Python versions tested in CI:

```shell
pixi run --environment py311 test
pixi run --environment py314 test
```

Install the pre-commit hooks once, then run all checks on demand using the dedicated environment:

```shell
pixi run --environment pre-commit pre-commit-install
pixi run --environment pre-commit pre-commit-run
```

Build or serve the documentation using the documentation environment:

```shell
pixi run --environment docs docs-build
pixi run --environment docs docs-serve
```

### Reproducing and updating the development environment

The committed `pixi.lock` contains resolved dependencies for Linux, Intel macOS, and Apple Silicon macOS. Run
`pixi install --frozen` to reproduce it without changing dependency versions. Use `pixi update` deliberately when the
lock file should be refreshed, and commit the resulting lock-file change together with the manifest change.

## Troubleshooting

### LightGBM

Many pip installation problems are caused by a broken LightGBM installation. Refer to the
[LightGBM installation instructions](https://github.com/microsoft/LightGBM/tree/master/python-package) and install its
system prerequisites before retrying PyPythia. Pixi and conda installations obtain LightGBM and its native dependencies
from conda-forge.

### Python version

PyPythia supports Python 3.11 through Python 3.14. For source development, select a locked Pixi environment explicitly
with `pixi run --environment py311 ...` or `pixi run --environment py314 ...`.

### Running PyPythia

When working from source, run commands through Pixi, for example `pixi run cli-smoke`, or start an interactive
environment with `pixi shell`. Global Pixi, Conda, and pip installations expose the `pythia` command directly.
