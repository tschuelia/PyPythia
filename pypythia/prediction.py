import pathlib
import shutil
from tempfile import TemporaryDirectory
from typing import Optional

import numpy as np
import pandas as pd

from pypythia.config import DEFAULT_MODEL_FILE, DEFAULT_RAXMLNG_EXE
from pypythia.custom_errors import PyPythiaException
from pypythia.custom_types import DataType, FileFormat
from pypythia.logger import log_runtime_information, logger
from pypythia.msa import MSA, deduplicate_sequences, parse, remove_full_gap_sequences
from pypythia.predictor import DifficultyPredictor
from pypythia.raxmlng import RAxMLNG


def predict_difficulty(
    msa_file: pathlib.Path,
    model_file: Optional[pathlib.Path] = DEFAULT_MODEL_FILE,
    raxmlng: Optional[pathlib.Path] = DEFAULT_RAXMLNG_EXE,
    threads: int = None,
    seed: int = 0,
    file_format: Optional[FileFormat] = None,
    data_type: Optional[DataType] = None,
    deduplicate: bool = True,
    remove_full_gaps: bool = True,
    reduced_msa_file: Optional[pathlib.Path] = None,
) -> np.float64:
    """Predict the difficulty of an MSA using the PyPythia difficulty predictor.

    Per default, the MSA is deduplicated and full gap sequences are removed before the difficulty is predicted.

    Args:
        msa_file (pathlib.Path): Path to the MSA file. Note that the MSA file must be in either FASTA or PHYLIP format.
        model_file (pathlib.Path, optional): Path to the trained difficulty predictor model.
            Defaults to the latest model shipped with PyPythia.
        raxmlng (pathlib.Path, optional): Path to the RAxML-NG executable.
            If not set, uses the RAxML-NG binary found in the PATH environment variable.
        threads (int, optional): Number of threads to use for parallel parsimony tree inference. If not set, uses the
            RAxML-NG auto parallelization scheme.
        seed (int, optional): Random seed to use for the parsimony tree inference. Defaults to 0.
        file_format (FileFormat, optional): File format of the MSA file. Defaults to None. In this case, the file format
            is inferred based on the file content. See `pypythia.msa.parse` for information on when this is required.
        data_type (DataType, optional): Data type of the MSA sequences. Defaults to None. In this case, the data type
            is inferred based on the file content. See `pypythia.msa.parse` for information on when this is required.
        deduplicate (bool, optional): If True, remove duplicate sequences from the MSA. Defaults to True.
        remove_full_gaps (bool, optional): If True, remove full gap sequences from the MSA. Defaults to True.
        reduced_msa_file (pathlib.Path, optional): Path to store the reduced MSA after deduplication and removal of full gap sequences.

    Returns:
        np.float64: Predicted difficulty of the MSA.
    """

    predictor = DifficultyPredictor(model_file=model_file)

    if raxmlng is None:
        raise PyPythiaException(
            "Path to the RAxML-NG executable is required if 'raxml-ng' is not in $PATH."
        )

    raxmlng = RAxMLNG(**{"exe_path": raxmlng} if raxmlng else {})
    msa = parse(msa_file, file_format=file_format, data_type=data_type)

    if deduplicate and msa.contains_duplicate_sequences():
        msa = deduplicate_sequences(msa)
    if remove_full_gaps and msa.contains_full_gap_sequences():
        msa = remove_full_gap_sequences(msa)

    if reduced_msa_file:
        msa.write(reduced_msa_file)

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
        pars_trees_file (pathlib.Path, optional): Path to store the inferred parsimony trees. Defaults to None.
            In this case, the trees are not stored.
        log_info (bool, optional): If True, log intermediate progress information using the default logger.
        threads (int, optional): The number of threads to use for parallel parsimony tree inference. Defaults to None.
            Uses the RAxML-NG auto parallelization scheme if none is set.
        seed (int, optional): Random seed to use for the parsimony tree inference. Defaults to 0.
    Returns:
        Dataframe containing a single row with all features required for predicting the difficulty of the MSA.
        The columns correspond to the feature names the predictor was trained with.
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
            "proportion_gaps": msa.proportion_gaps,
            "proportion_invariant": msa.proportion_invariant,
            "entropy": msa.entropy(),
            "bollback": msa.bollback_multinomial(),
            "pattern_entropy": msa.pattern_entropy(),
            "avg_rfdist_parsimony": rel_rfdist,
            "proportion_unique_topos_parsimony": num_topos / n_pars_trees,
        }
        return pd.DataFrame(features, index=[0])
