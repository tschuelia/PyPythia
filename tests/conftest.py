import pathlib

import pandas as pd
import pytest

from pypythia.predictor import DifficultyPredictor
from pypythia.raxmlng import RAxMLNG

from .test_config import RAXMLNG_COMMAND

df = pd.read_parquet("tests/data/msa_test_data.parquet")
df["msa_file"] = df["msa_file"].apply(
    lambda x: pathlib.Path.cwd() / "tests" / "data" / x
)


@pytest.fixture
def msa_test_data():
    return df


@pytest.fixture(
    params=df.iterrows(),
    ids=lambda x: f"{x[1].msa_file.parent.name}/{x[1].msa_file.name}",
)
def msa_test_data_row(request):
    idx, row = request.param
    return row


@pytest.fixture
def phylip_msa_file():
    return pathlib.Path.cwd() / "tests" / "data" / "DNA" / "0.phy"


@pytest.fixture
def small_msa_file():
    return pathlib.Path.cwd() / "tests" / "data" / "DNA" / "small.fasta"


@pytest.fixture
def predictor():
    return DifficultyPredictor(pathlib.Path("pypythia/predictors/latest.txt"))


@pytest.fixture
def raxmlng_command():
    return pathlib.Path(RAXMLNG_COMMAND)


@pytest.fixture
def raxmlng(raxmlng_command):
    return RAxMLNG(raxmlng_command)


@pytest.fixture
def data_dir():
    return pathlib.Path.cwd() / "tests" / "data"


@pytest.fixture
def multiple_trees_path():
    return pathlib.Path.cwd() / "tests" / "data" / "trees" / "many.trees"


@pytest.fixture
def raxmlng_rfdistance_log():
    return pathlib.Path.cwd() / "tests" / "data" / "logs" / "raxml.rfdistance.log"


@pytest.fixture
def raxmlng_inference_log():
    return pathlib.Path.cwd() / "tests" / "data" / "logs" / "raxml.inference.log"
