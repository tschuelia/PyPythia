import pytest
import os

from pypythia.msa import MSA
from pypythia.raxmlng import RAxMLNG
from pypythia.predictor import DifficultyPredictor


@pytest.fixture
def example_msa_path():
    cwd = os.getcwd()
    return os.path.join(cwd, "tests", "data", "DNA", "0.phy")


@pytest.fixture
def dna_phylip_msa(example_msa_path):
    return MSA(example_msa_path)


@pytest.fixture
def dna_fasta_msa():
    cwd = os.getcwd()
    return MSA(os.path.join(cwd, "tests", "data", "DNA", "0.fasta"))


@pytest.fixture
def small_msa():
    cwd = os.getcwd()
    return MSA(os.path.join(cwd, "tests", "data", "DNA", "small.fasta"))


@pytest.fixture
def small_msa_with_signal():
    cwd = os.getcwd()
    return MSA(os.path.join(cwd, "tests", "data", "DNA", "3.phy"))


@pytest.fixture
def all_msa_files_with_model():
    cwd = os.getcwd()

    morph_models = {
        "0.phy": "MULTI3_GTR",
        "1.phy": "MULTI2_GTR"
    }

    files_and_models = []

    for data_type in ["DNA", "AA", "MORPH"]:
        base_dir = os.path.join(cwd, "tests", "data", data_type)
        for f in os.listdir(base_dir):
            msa_file = os.path.join(base_dir, f)
            model = ""
            if data_type == "DNA":
                model = "GTR+G"
            elif data_type == "AA":
                model = "LG+G"
            elif data_type == "MORPH":
                model = morph_models[f]

            files_and_models.append((msa_file, model))

    return files_and_models


@pytest.fixture
def predictor():
    return DifficultyPredictor(open("pypythia/predictor.pckl", "rb"))


@pytest.fixture
def raxmlng_command():
    return "/Users/julia/Desktop/Promotion/software/raxml-ng/bin/raxml-ng"


@pytest.fixture
def raxmlng(raxmlng_command):
    return RAxMLNG(raxmlng_command)


@pytest.fixture
def multiple_trees_path():
    cwd = os.getcwd()
    return os.path.join(cwd, "tests", "data", "trees", "many.trees")


@pytest.fixture
def raxmlng_rfdistance_log():
    cwd = os.getcwd()
    return os.path.join(cwd, "tests", "data", "logs", "raxml.rfdistance.log")


@pytest.fixture
def raxmlng_inference_log():
    cwd = os.getcwd()
    return os.path.join(cwd, "tests", "data", "logs", "raxml.inference.log")
