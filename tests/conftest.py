import pathlib

import pandas as pd
import pytest

from pypythia.msa import MSA, parse
from pypythia.predictor import DifficultyPredictor
from pypythia.raxmlng import RAxMLNG

from .test_config import RAXMLNG_COMMAND


@pytest.fixture
def msa_test_data():
    df = pd.read_parquet("tests/data/msa_test_data.parquet")
    df["msa_file"] = df["msa_file"].apply(
        lambda x: pathlib.Path.cwd() / "tests" / "data" / x
    )
    return df


@pytest.fixture
def phylip_msa_file():
    return pathlib.Path.cwd() / "tests" / "data" / "DNA" / "0.phy"


@pytest.fixture
def fasta_msa_file():
    return pathlib.Path.cwd() / "tests" / "data" / "DNA" / "0.fasta"


@pytest.fixture
def small_msa_file():
    return pathlib.Path.cwd() / "tests" / "data" / "DNA" / "small.fasta"


@pytest.fixture
def msa_with_duplicates_and_full_gap_sequences():
    return pathlib.Path.cwd() / "tests" / "data" / "DNA" / "5.phy"


@pytest.fixture
def dna_phylip_msa(phylip_msa_file):
    return parse(phylip_msa_file)


@pytest.fixture
def dna_fasta_msa(fasta_msa_file):
    return parse(fasta_msa_file)


@pytest.fixture
def small_msa(small_msa_file):
    return parse(small_msa_file)


@pytest.fixture
def small_msa_with_signal():
    pth = pathlib.Path.cwd() / "tests" / "data" / "DNA" / "3.phy"
    return parse(pth)


@pytest.fixture
def msa_with_duplicate_sequences(dna_phylip_msa):
    return dna_phylip_msa


@pytest.fixture
def msa_without_duplicate_sequences(small_msa_with_signal):
    return small_msa_with_signal


@pytest.fixture
def all_msa_files_with_model():
    morph_models = {"0.phy": "MULTI3_GTR", "1.phy": "MULTI2_GTR"}

    files_and_models = []

    for data_type in ["DNA", "AA", "MORPH"]:
        base_dir = pathlib.Path.cwd() / "tests" / "data" / data_type
        for f in base_dir.iterdir():
            model = ""
            if data_type == "DNA":
                model = "GTR+G"
            elif data_type == "AA":
                model = "LG+G"
            elif data_type == "MORPH":
                model = morph_models[f.name]

            files_and_models.append((f, model))

    return files_and_models


@pytest.fixture
def predictor():
    return DifficultyPredictor(open("pypythia/predictors/latest.pckl", "rb"))


@pytest.fixture
def sklearn_predictor():
    return DifficultyPredictor(
        open("pypythia/predictors/predictor_sklearn_rf_v0.0.1.pckl", "rb")
    )


@pytest.fixture
def raxmlng_command():
    return pathlib.Path(RAXMLNG_COMMAND)


@pytest.fixture
def raxmlng(raxmlng_command):
    return RAxMLNG(raxmlng_command)


@pytest.fixture
def multiple_trees_path():
    return pathlib.Path.cwd() / "tests" / "data" / "trees" / "many.trees"


@pytest.fixture
def raxmlng_rfdistance_log():
    return pathlib.Path.cwd() / "tests" / "data" / "logs" / "raxml.rfdistance.log"


@pytest.fixture
def raxmlng_inference_log():
    return pathlib.Path.cwd() / "tests" / "data" / "logs" / "raxml.inference.log"
