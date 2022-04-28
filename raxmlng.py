from utils import run_cmd
from raxmlng_parser import *

from tempfile import TemporaryDirectory


class RAxMLNG:
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

    def run_alignment_parse(
        self, msa_file: FilePath, model: Model, prefix: str, **kwargs
    ) -> None:
        cmd = self._base_cmd(msa_file, model, prefix, parse=None, **kwargs)
        run_cmd(cmd)

    def run_rfdist(self, trees_file: FilePath, **kwargs) -> None:
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
    ) -> None:
        cmd = self._base_cmd(
            msa_file, model, prefix, start=None, tree=f"pars{{{n_trees}}}", **kwargs
        )
        run_cmd(cmd)

    def get_rfdistance_results(
        self, trees_file: FilePath, **kwargs
    ) -> Tuple[float, float, float]:
        self.run_rfdist(trees_file, **kwargs)
        log_file = trees_file + ".raxml.log"
        num_topos = get_raxmlng_num_unique_topos(log_file)
        rel_rfdist = get_raxmlng_rel_rf_distance(log_file)
        abs_rfdist = get_raxmlng_abs_rf_distance(log_file)

        return num_topos, rel_rfdist, abs_rfdist

    def get_patterns_gaps_invariant(
        self, msa_file: FilePath, model: Model
    ) -> Tuple[int, float, float]:
        with TemporaryDirectory() as tmpdir:
            prefix = tmpdir + "/parse"
            self.run_alignment_parse(msa_file, model, prefix)
            return get_patterns_gaps_invariant(f"{prefix}.raxml.log")
