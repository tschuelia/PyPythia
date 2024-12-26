import pathlib
import tempfile
from tempfile import TemporaryDirectory

import pytest

from pypythia.custom_errors import RAxMLNGError
from pypythia.raxmlng import _get_raxmlng_rfdist_results, run_raxmlng_command


def test_get_raxmlng_rfdist_results(raxmlng_rfdistance_log):
    num_topos, rel_rfdist, abs_rfdist = _get_raxmlng_rfdist_results(
        raxmlng_rfdistance_log
    )

    assert num_topos == 2
    assert rel_rfdist == pytest.approx(1 / 3, abs=0.01)
    assert abs_rfdist == pytest.approx(2)


def test_get_raxmlng_rfdist_results_raises_value_error(raxmlng_inference_log):
    with pytest.raises(ValueError):
        test_get_raxmlng_rfdist_results(raxmlng_inference_log)


def test_infer_parsimony_trees(raxmlng, phylip_msa_file):
    with TemporaryDirectory() as tmpdir:
        file_path = raxmlng.infer_parsimony_trees(
            msa_file=phylip_msa_file,
            model="GTR+G",
            prefix=pathlib.Path(tmpdir),
            n_trees=10,
        )

        trees = file_path.open().readlines()
        trees = [t.strip() for t in trees if t]

        assert len(trees) == 10


def test_get_rfdistance_results(raxmlng, multiple_trees_path):
    num_topos, rel_rfdist, abs_rfdist = raxmlng.get_rfdistance_results(
        multiple_trees_path
    )
    assert num_topos == 6
    assert rel_rfdist == pytest.approx(0.114, abs=0.01)
    assert abs_rfdist == pytest.approx(22.269, abs=0.1)


def test_run_raxmlng_command(raxmlng_command, phylip_msa_file):
    with tempfile.TemporaryDirectory() as tmpdir:
        prefix = pathlib.Path(tmpdir) / "test"
        raxmlng_example = [
            raxmlng_command,
            "--parse",
            "--msa",
            phylip_msa_file,
            "--model",
            "GTR+G",
            "--prefix",
            prefix,
        ]
        run_raxmlng_command(raxmlng_example)


def test_run_raxmlng_command_raises_raxmlng_error_with_error_output(raxmlng_command):
    with pytest.raises(
        RAxMLNGError, match="ERROR: Alignment file not found: this_does_no_exist.phy"
    ):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = pathlib.Path(tmpdir) / "test"
            raxmlng_failure = [
                str(raxmlng_command.absolute()),
                "--parse",
                "--msa",
                "this_does_no_exist.phy",
                "--model",
                "GTR+G",
                "--prefix",
                str(prefix.absolute()),
            ]
            run_raxmlng_command(raxmlng_failure)


def test_run_raxmlng_raises_runtime_error():
    with pytest.raises(RuntimeError, match="Running RAxML-NG command failed."):
        run_raxmlng_command(["this", "is", "invalid"])
