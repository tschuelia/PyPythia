import os
import pathlib
import tempfile
from tempfile import TemporaryDirectory

import pytest

import pypythia.raxmlng as raxmlng_module
from pypythia.custom_errors import RAxMLNGError
from pypythia.raxmlng import RAxMLNG, get_raxmlng_rfdist_results, run_raxmlng_command


def _raxmlng_with_version(tmp_path, monkeypatch, version):
    executable = tmp_path / "raxml-ng"
    executable.touch()
    monkeypatch.setattr(
        raxmlng_module.subprocess,
        "check_output",
        lambda *args, **kwargs: f"RAxML-NG v. {version} released for testing.",
    )
    return RAxMLNG(executable)


def test_raxmlng_init(raxmlng_command):
    assert raxmlng_command.exists()
    # Should work without any problems
    RAxMLNG(raxmlng_command)


@pytest.mark.parametrize(
    ("version", "expected_major_version"),
    [("1.2.2", 1), ("2.0.2", 2), ("2.0-beta3", 2)],
)
def test_raxmlng_init_parses_major_version(
    tmp_path, monkeypatch, version, expected_major_version
):
    raxmlng = _raxmlng_with_version(tmp_path, monkeypatch, version)

    assert raxmlng._major_version == expected_major_version


def test_raxmlng_init_fails_for_unparseable_version(tmp_path, monkeypatch):
    executable = tmp_path / "raxml-ng"
    executable.touch()
    monkeypatch.setattr(
        raxmlng_module.subprocess,
        "check_output",
        lambda *args, **kwargs: "RAxML-NG development build",
    )

    with pytest.raises(RuntimeError, match="Could not determine the RAxML-NG version"):
        RAxMLNG(executable)


def test_raxmlng_init_fails_non_existing_exe():
    with pytest.raises(FileNotFoundError, match="RAxML-NG executable not found."):
        RAxMLNG(pathlib.Path("this_does_not_exist"))


def test_raxmlng_init_fails_wrong_exe(raxmlng_command):
    with pytest.raises(
        RuntimeError, match="Your RAxML-NG executable does not seem to work."
    ):
        tmpfile = tempfile.NamedTemporaryFile("wb", delete=False)
        # Manually break the RAxML-NG file to trigger the executable-broken error
        tmpfile.write(b"NonSense" + raxmlng_command.read_bytes())
        tmpfile.close()
        os.chmod(tmpfile.name, 0o777)
        RAxMLNG(pathlib.Path(tmpfile.name))
        os.unlink(tmpfile.name)


def test_raxmlng_init_fails_non_raxmlng_exe():
    with pytest.raises(
        RuntimeError,
        match="The given executable `.*` does not seem to be a RAxML-NG executable.",
    ):
        tmpfile = tempfile.NamedTemporaryFile("w", delete=False)
        tmpfile.write("#!/bin/bash\n")
        tmpfile.write("echo test\n")
        tmpfile.close()
        os.chmod(tmpfile.name, 0o777)
        RAxMLNG(pathlib.Path(tmpfile.name))
        os.unlink(tmpfile.name)


def test_get_raxmlng_rfdist_results(raxmlng_rfdistance_log):
    num_topos, rel_rfdist = get_raxmlng_rfdist_results(raxmlng_rfdistance_log)

    assert num_topos == 2
    assert rel_rfdist == pytest.approx(1 / 3, abs=0.01)


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


@pytest.mark.parametrize(
    ("version", "expected_adaptive_setting"),
    [("1.2.2", None), ("2.0.2", "off"), ("3.0.0", "off")],
)
def test_infer_parsimony_trees_sets_adaptive_by_raxmlng_version(
    tmp_path, monkeypatch, version, expected_adaptive_setting
):
    raxmlng = _raxmlng_with_version(tmp_path, monkeypatch, version)
    commands = []
    monkeypatch.setattr(raxmlng_module, "run_raxmlng_command", commands.append)

    raxmlng.infer_parsimony_trees(
        msa_file=tmp_path / "msa.phy",
        model="GTR+G",
        prefix=tmp_path / "inference",
    )

    command = commands[0]
    if expected_adaptive_setting is None:
        assert "--adaptive" not in command
    else:
        adaptive_index = command.index("--adaptive")
        assert command[adaptive_index + 1] == expected_adaptive_setting
        assert command.count("--adaptive") == 1


@pytest.mark.parametrize(
    ("version", "expected_adaptive_setting"), [("1.2.2", None), ("2.0.2", "off")]
)
def test_infer_parsimony_trees_overrides_adaptive_setting(
    tmp_path, monkeypatch, version, expected_adaptive_setting
):
    raxmlng = _raxmlng_with_version(tmp_path, monkeypatch, version)
    commands = []
    monkeypatch.setattr(raxmlng_module, "run_raxmlng_command", commands.append)

    raxmlng.infer_parsimony_trees(
        msa_file=tmp_path / "msa.phy",
        model="GTR+G",
        prefix=tmp_path / "inference",
        adaptive="on",
    )

    command = commands[0]
    if expected_adaptive_setting is None:
        assert "--adaptive" not in command
    else:
        adaptive_index = command.index("--adaptive")
        assert command[adaptive_index + 1] == expected_adaptive_setting
        assert command.count("--adaptive") == 1


def test_get_rfdistance_results(raxmlng, multiple_trees_path):
    num_topos, rel_rfdist = raxmlng.get_rfdistance_results(multiple_trees_path)
    assert num_topos == 6
    assert rel_rfdist == pytest.approx(0.114, abs=0.01)


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
