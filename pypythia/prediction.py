import pathlib
import shutil
from tempfile import TemporaryDirectory
from typing import Optional

import numpy as np
import pandas as pd

from build.lib.pypythia.custom_errors import PyPythiaException
from pypythia.logger import log_runtime_information, logger
from pypythia.msa import MSA, parse
from pypythia.predictor import DEFAULT_MODEL_FILE, DifficultyPredictor
from pypythia.raxmlng import DEFAULT_RAXMLNG_EXE, RAxMLNG


def predict_difficulty(
    msa_file: pathlib.Path,
    model_file: Optional[pathlib.Path] = DEFAULT_MODEL_FILE,
    raxmlng: Optional[pathlib.Path] = DEFAULT_RAXMLNG_EXE,
    threads: int = None,
    seed: int = 0,
) -> np.float64:
    """
    Predicts the difficulty of an MSA using the given difficulty predictor.

    Args:
        msa_file (FilePath): Path to the MSA file the difficulty should be predicted for. The file must be either in "fasta" or "phylip" format.
        model_file (FilePath): Path to a trained difficulty predictor.
        raxmlng (Executable): Path to an executable of RAxML-NG. See https://github.com/amkozlov/raxml-ng for install instructions.
        threads (int, optional): The number of threads to use for parallel parsimony tree inference. Uses the RAxML-NG auto parallelization scheme if none is set.
        seed (int, optional): Seed for the RAxML-NG parsimony tree inference. Default is 0.

    Returns:
        difficulty (float): The predicted difficulty for the given MSA.

    Raises:
        ValueError: If the file format of the given MSA is not FASTA or PHYLIP.
        ValueError: If the data type of the given MSA cannot be inferred.
        PyPythiaException: If the provided difficulty predictor was trained with a subset incompatible to Pythia.
    """

    predictor = DifficultyPredictor(model_file=model_file)

    if raxmlng is None:
        raise PyPythiaException(
            "Path to the RAxML-NG executable is required if 'raxml-ng' is not in $PATH."
        )

    raxmlng = RAxMLNG(**{"exe_path": raxmlng} if raxmlng else {})
    msa = parse(msa_file)
    msa_features = collect_features(
        msa, msa_file, raxmlng, log_info=False, threads=threads, seed=seed
    )
    difficulty = predictor.predict(msa_features)

    return difficulty[0]


def collect_features(
    msa: MSA,
    msa_file: pathlib.Path,
    raxmlng: RAxMLNG,
    pars_trees_file: Optional[pathlib.Path] = None,
    log_info: bool = True,
    threads: int = None,
    seed: int = 0,
) -> pd.DataFrame:
    """Helper function to collect all features required for predicting the difficulty of the MSA.

    Args:
        msa (MSA): MSA object corresponding to the MSA file to compute the features for.
        raxmlng (RAxMLNG): Initialized RAxMLNG object.
        store_trees (bool, optional): If True, store the inferred parsimony trees as "{msa_name}.parsimony.trees" file in the current workdir.
        log_info (bool, optional): If True, log intermediate progress information using the default logger.
        threads (int, optional): The number of threads to use for parallel parsimony tree inference. Uses the RAxML-NG auto parallelization scheme if none is set.
    Returns:
        all_features (Dict): Dictionary containing all features required for predicting the difficulty of the MSA. The keys correspond to the feature names the predictor was trained with.
    """
    if not log_info:
        logger.remove()

    with TemporaryDirectory() as tmpdir:
        msa_file = msa_file
        model = msa.get_raxmlng_model()

        log_runtime_information("Retrieving num_taxa, num_sites.", log_runtime=True)

        n_pars_trees = 24
        log_runtime_information(
            f"Inferring {n_pars_trees} parsimony trees with random seed {seed}.",
            log_runtime=True,
        )
        trees = raxmlng.infer_parsimony_trees(
            msa_file,
            model,
            pathlib.Path(tmpdir) / "pars",
            redo=None,
            seed=seed,
            n_trees=n_pars_trees,
            **dict(threads=threads) if threads else {},
        )
        if pars_trees_file is not None:
            log_runtime_information(
                f"Storing the inferred parsimony trees in the file {pars_trees_file}."
            )
            shutil.copy(trees, pars_trees_file)

        log_runtime_information(
            "Computing the RF-Distance for the parsimony trees.", log_runtime=True
        )
        num_topos, rel_rfdist, _ = raxmlng.get_rfdistance_results(trees, redo=None)

        features = {
            "num_taxa": msa.n_taxa,
            "num_sites": msa.n_sites,
            "num_patterns": msa.n_patterns,
            "num_patterns/num_taxa": msa.n_patterns / msa.n_taxa,
            "num_sites/num_taxa": msa.n_sites / msa.n_taxa,
            "num_patterns/num_sites": msa.n_patterns / msa.n_sites,
            "proportion_gaps": msa.percentage_gaps,
            "proportion_invariant": msa.percentage_invariant,
            "entropy": msa.entropy(),
            "bollback": msa.bollback_multinomial(),
            "pattern_entropy": msa.pattern_entropy(),
            "avg_rfdist_parsimony": rel_rfdist,
            "proportion_unique_topos_parsimony": num_topos / n_pars_trees,
        }
        return pd.DataFrame(features, index=[0])
