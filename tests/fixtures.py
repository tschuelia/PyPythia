import pathlib

import pytest

from .test_config import RAXMLNG_COMMAND

from pypythia.msa import MSA
from pypythia.raxmlng import RAxMLNG
from pypythia.predictor import DifficultyPredictor


@pytest.fixture
def example_msa_path():
    return pathlib.Path.cwd() / "tests" / "data" / "DNA" / "0.phy"


@pytest.fixture
def dna_phylip_msa(example_msa_path):
    return MSA(example_msa_path)


@pytest.fixture
def dna_fasta_msa():
    pth = pathlib.Path.cwd() / "tests" / "data" / "DNA" / "0.fasta"
    return MSA(pth)


@pytest.fixture
def small_msa():
    pth = pathlib.Path.cwd() / "tests" / "data" / "DNA" / "small.fasta"
    return MSA(pth)


@pytest.fixture
def small_msa_with_signal():
    pth = pathlib.Path.cwd() / "tests" / "data" / "DNA" / "3.phy"
    return MSA(pth)


@pytest.fixture
def msa_with_duplicate_sequences(dna_phylip_msa):
    return dna_phylip_msa


@pytest.fixture
def msa_without_duplicate_sequences(small_msa_with_signal):
    return small_msa_with_signal


@pytest.fixture
def all_msa_files_with_model():
    morph_models = {
        "0.phy": "MULTI3_GTR",
        "1.phy": "MULTI2_GTR"
    }

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
    return DifficultyPredictor(open("pypythia/predictors/predictor_sklearn_rf_v0.0.1.pckl", "rb"))


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
