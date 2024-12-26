import pathlib
import shutil
import subprocess
from tempfile import TemporaryDirectory
from typing import Optional, Union

from pypythia.custom_errors import RAxMLNGError


def run_raxmlng_command(cmd: list[str]) -> None:
    try:
        subprocess.check_output(cmd, encoding="utf-8")
    except subprocess.CalledProcessError as e:
        raise RAxMLNGError(subprocess_exception=e)
    except Exception as e:
        raise RuntimeError("Running RAxML-NG command failed.") from e


_raxmlng = shutil.which("raxml-ng")
DEFAULT_RAXMLNG_EXE = pathlib.Path(_raxmlng) if _raxmlng else None


def _get_value_from_line(line: str, search_string: str) -> float:
    line = line.strip()
    if search_string in line:
        _, value = line.rsplit(" ", 1)
        return float(value)

    raise ValueError(
        f"The given line '{line}' does not contain the search string '{search_string}'."
    )


def _get_raxmlng_rfdist_results(log_file: pathlib.Path) -> tuple[float, float, float]:
    """Method that parses the number of unique topologies, relative RF-Distance, and absolute RF-Distance in the given log file.

    Args:
        log_file (str): Filepath of a RAxML-NG log file.

    Returns:
        num_topos (int): Number of unique topologies of the given set of trees.
        rel_rfdist (float): Relative RF-Distance of the given set of trees. Computed as average over all pairwise RF-Distances. Value between 0.0 and 1.0.
        abs_rfdist (float): Absolute RF-Distance of the given set of trees.

    Raises:
        ValueError: If the given log file does not contain the unique topologies, relative RF-Distance, or absolute RF-Distance.
    """
    abs_rfdist = None
    rel_rfdist = None
    num_topos = None

    for line in log_file.open().readlines():
        line = line.strip()

        if "Average absolute RF distance in this tree set:" in line:
            abs_rfdist = _get_value_from_line(
                line, "Average absolute RF distance in this tree set:"
            )
        elif "Average relative RF distance in this tree set:" in line:
            rel_rfdist = _get_value_from_line(
                line, "Average relative RF distance in this tree set:"
            )
        elif "Number of unique topologies in this tree set:" in line:
            num_topos = _get_value_from_line(
                line, "Number of unique topologies in this tree set:"
            )

    if abs_rfdist is None or rel_rfdist is None or num_topos is None:
        raise ValueError("Error parsing raxml-ng log.")

    return num_topos, rel_rfdist, abs_rfdist


class RAxMLNG:
    """Class structure for features computed using RAxML-NG.

    This class provides methods for computing MSA attributes using RAxML-NG.

    Args:
        exe_path (Executable): Path to an executable of RAxML-NG. See https://github.com/amkozlov/raxml-ng for install instructions.

    Attributes:
        exe_path (Executable): Path to an executable of RAxML-NG.
    """

    def __init__(self, exe_path: Optional[pathlib.Path] = DEFAULT_RAXMLNG_EXE):
        self.exe_path = exe_path

    def _base_cmd(
        self, msa_file: pathlib.Path, model: str, prefix: pathlib.Path, **kwargs
    ) -> list[str]:
        additional_settings = []
        for key, value in kwargs.items():
            if value is None:
                additional_settings += [f"--{key}"]
            else:
                additional_settings += [f"--{key}", str(value)]

        return [
            str(self.exe_path.absolute()),
            "--msa",
            str(msa_file.absolute()),
            "--model",
            model,
            "--prefix",
            str(prefix.absolute()),
            *additional_settings,
        ]

    def _run_alignment_parse(
        self, msa_file: pathlib.Path, model: str, prefix: pathlib.Path, **kwargs
    ) -> None:
        cmd = self._base_cmd(msa_file, model, prefix, parse=None, **kwargs)
        run_raxmlng_command(cmd)

    def _run_rfdist(
        self, trees_file: pathlib.Path, prefix: pathlib.Path, **kwargs
    ) -> None:
        additional_settings = []
        for key, value in kwargs.items():
            if value is None:
                additional_settings += [f"--{key}"]
            else:
                additional_settings += [f"--{key}", str(value)]
        cmd = [
            str(self.exe_path.absolute()),
            "--rfdist",
            str(trees_file.absolute()),
            "--prefix",
            str(prefix.absolute()),
            *additional_settings,
        ]
        run_raxmlng_command(cmd)

    def infer_parsimony_trees(
        self,
        msa_file: pathlib.Path,
        model: str,
        prefix: pathlib.Path,
        n_trees: int = 24,
        **kwargs,
    ) -> pathlib.Path:
        """Method that infers n_trees using the RAxML-NG implementation of maximum parsimony.

        Args:
            msa_file (str): Filepath of the MSA to compute the parsimony trees for.
            model (str): String representation of the substitution model to use. Needs to be a valid RAxML-NG model. For example "GTR+G" for DNA data or "LG+G" for protein data.
            prefix (str): Prefix of where to store the RAxML-NG results.
            n_trees (int): Number of parsimony trees to compute.
            **kwargs: Optional additional RAxML-NG settings.
                The name of the kwarg needs to be a valid RAxML-NG flag.
                For flags with a value pass it like this: "flag=value", for flags without a value pass it like this: "flag=None".
                See https://github.com/amkozlov/raxml-ng for all options.

        Returns:
            output_trees_file (str): Filepath pointing to the computed trees.

        """
        cmd = self._base_cmd(
            msa_file, model, prefix, start=None, tree=f"pars{{{n_trees}}}", **kwargs
        )
        run_raxmlng_command(cmd)
        return pathlib.Path(f"{prefix}.raxml.startTree")

    def get_rfdistance_results(
        self, trees_file: pathlib.Path, prefix: pathlib.Path = None, **kwargs
    ) -> tuple[float, float, float]:
        """Method that computes the number of unique topologies, relative RF-Distance, and absolute RF-Distance for the given set of trees.

        Args:
            trees_file: Filepath of a file containing > 1 Newick strings.
            prefix (str): Optional prefix to use when running RAxML-NG

        Returns:
            num_topos (float): Number of unique topologies of the given set of trees.
            rel_rfdist (float): Relative RF-Distance of the given set of trees. Computed as average over all pairwise RF-Distances. Value between 0.0 and 1.0.
            abs_rfdist (float): Absolute RF-Distance of the given set of trees.
        """
        with TemporaryDirectory() as tmpdir:
            tmpdir = pathlib.Path(tmpdir)
            if not prefix:
                prefix = tmpdir / "rfdist"
            self._run_rfdist(trees_file, prefix, **kwargs)
            log_file = pathlib.Path(f"{prefix}.raxml.log")
            return _get_raxmlng_rfdist_results(log_file)
