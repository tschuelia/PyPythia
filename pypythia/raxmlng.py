from pypythia.utils import run_cmd
from pypythia.custom_types import *
from pypythia.raxmlng_parser import *

from tempfile import TemporaryDirectory


class RAxMLNG:
    """Class structure for features computed using RAxML-NG.

    This class provides methods for computing MSA attributes using RAxML-NG.

    Args:
        exe_path (str): Path to an executable of RAxML-NG. See https://github.com/amkozlov/raxml-ng for install instructions.

    Attributes:
        exe_path (str): Path to an executable of RAxML-NG.
    """

    def __init__(self, exe_path: Executable):
        self.exe_path = exe_path

    def _base_cmd(
            self, msa_file: FilePath, model: Model, prefix: str, **kwargs
    ) -> Command:
        additional_settings = []
        for key, value in kwargs.items():
            if value is None:
                additional_settings += [f"--{key}"]
            else:
                additional_settings += [f"--{key}", str(value)]

        return [
            self.exe_path,
            "--msa",
            msa_file,
            "--model",
            model,
            "--prefix",
            prefix,
            *additional_settings,
        ]

    def _run_alignment_parse(
            self, msa_file: FilePath, model: Model, prefix: str, **kwargs
    ) -> None:
        cmd = self._base_cmd(msa_file, model, prefix, parse=None, **kwargs)
        run_cmd(cmd)

    def _run_rfdist(self, trees_file: FilePath, prefix: str, **kwargs) -> None:
        additional_settings = []
        for key, value in kwargs.items():
            if value is None:
                additional_settings += [f"--{key}"]
            else:
                additional_settings += [f"--{key}", str(value)]
        cmd = [
            self.exe_path,
            "--rfdist",
            trees_file,
            "--prefix",
            prefix,
            *additional_settings,
        ]
        run_cmd(cmd)

    def infer_parsimony_trees(
            self,
            msa_file: FilePath,
            model: Model,
            prefix: str,
            n_trees: int = 100,
            **kwargs,
    ) -> FilePath:
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
        run_cmd(cmd)
        return prefix + ".raxml.startTree"

    def get_rfdistance_results(
            self, trees_file: FilePath, prefix: str = None, **kwargs
    ) -> Tuple[float, float, float]:
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
            if not prefix:
                prefix = tmpdir
            self._run_rfdist(trees_file, prefix, **kwargs)
            log_file = prefix + ".raxml.log"
            return get_raxmlng_rfdist_results(log_file)

    def get_patterns_gaps_invariant(
            self, msa_file: FilePath, model: Model, prefix: str = None
    ) -> Tuple[int, float, float]:
        """Method that obtains the number of patterns, proportion of gaps, and proportion of invariant sites in the given MSA.

        Args:
            msa_file (str): Filepath of the MSA to compute the parsimony trees for.
            model (str): String representation of the substitution model to use. Needs to be a valid RAxML-NG model. For example "GTR+G" for DNA data or "LG+G" for protein data.
            prefix (str): Optional prefix to use when running RAxML-NG

        Returns:
            n_patterns (int): Number of unique patterns in the given MSA.
            prop_gaps (float): Proportion of gaps in the given MSA.
            prop_inv (float): Proportion of invariant sites in the given MSA.
        """
        with TemporaryDirectory() as tmpdir:
            if not prefix:
                prefix = tmpdir + "/parse"
            self._run_alignment_parse(msa_file, model, prefix)
            return get_patterns_gaps_invariant(f"{prefix}.raxml.log")
