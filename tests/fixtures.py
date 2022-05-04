import pytest
import os

from msa import MSA
from raxmlng import RAxMLNG


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
