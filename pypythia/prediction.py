import pathlib
import shutil
import tempfile
from tempfile import TemporaryDirectory
from typing import Optional

import numpy as np
import pandas as pd

from pypythia.config import DEFAULT_MODEL_FILE, DEFAULT_RAXMLNG_EXE
from pypythia.custom_errors import PyPythiaException
from pypythia.custom_types import DataType, FileFormat
from pypythia.logger import log_runtime_information, logger
from pypythia.msa import (
    MSA,
    deduplicate_sequences,
    parse_msa,
    remove_full_gap_sequences,
)
from pypythia.predictor import DifficultyPredictor
from pypythia.raxmlng import RAxMLNG


def _handle_duplicates(msa: MSA, deduplicate: bool, log_info: bool = False) -> MSA:
    contains_duplicates = msa.contains_duplicate_sequences()
    if contains_duplicates and deduplicate:
        log_info and log_runtime_information(
            "The input MSA contains duplicate sequences. Removing duplicates before predicting the difficulty."
        )
        return deduplicate_sequences(msa)
    elif contains_duplicates:
        log_info and logger.warning(
            "WARNING: The provided MSA contains duplicate sequences, but deduplication is disabled. "
            "Pythia will predict the difficulty for the MSA with duplicates."
        )
        return msa
    else:
        return msa


def _handle_full_gap_sequences(
    msa: MSA, remove_full_gaps: bool, log_info: bool = False
) -> MSA:
    contains_full_gaps = msa.contains_full_gap_sequences()
    if contains_full_gaps and remove_full_gaps:
        log_info and log_runtime_information(
            "The input MSA contains sequences with only gaps. Removing full gap sequences before predicting the difficulty."
        )
        return remove_full_gap_sequences(msa)
    elif contains_full_gaps:
        log_info and log_runtime_information(
            "WARNING: The provided MSA contains sequences with only gaps, but gap removal is disabled. "
            "Pythia will predict the difficulty for the MSA with full gap sequences."
        )
        return msa
    else:
        return msa


def collect_features(
    msa: MSA,
    msa_file: pathlib.Path,
    raxmlng: RAxMLNG,
    pars_trees_file: Optional[pathlib.Path] = None,
    log_info: bool = False,
    threads: int = None,
    seed: int = 0,
) -> pd.DataFrame:
    """Helper function to collect all features required for predicting the difficulty of the MSA.

    Args:
        msa (MSA): MSA object corresponding to the MSA file to compute the features for.
        raxmlng (RAxMLNG): Initialized RAxMLNG object.
        pars_trees_file (pathlib.Path, optional): Path to store the inferred parsimony trees. Defaults to None.
            In this case, the trees are not stored.
        log_info (bool, optional): If True, log intermediate progress information using the default logger. Defaults to False.
        threads (int, optional): The number of threads to use for parallel parsimony tree inference. Defaults to None.
            Uses the RAxML-NG auto parallelization scheme if none is set.
        seed (int, optional): Random seed to use for the parsimony tree inference. Defaults to 0.
    Returns:
        Dataframe containing a single row with all features required for predicting the difficulty of the MSA.
        The columns correspond to the feature names the predictor was trained with.
    """
    # If the MSA contains less than 4 sequences, RAxML-NG will fail as there is only a single possible
    # tree topology for this MSA. In this case, any phylogenetic inference is meaningless and we raise a
    # PyPythia exception to inform the user.
    if msa.n_taxa < 4:
        raise PyPythiaException(
            "The MSA contains less than 4 sequences. "
            "Phylogenetic inference is not meaningful for such small MSAs as there exists only a single possible tree topology. "
        )

    with TemporaryDirectory() as tmpdir:
        n_pars_trees = 24
        log_info and log_runtime_information(
            f"Inferring {n_pars_trees} parsimony trees with random seed {seed}.",
        )
        trees = raxmlng.infer_parsimony_trees(
            msa_file,
            msa.get_raxmlng_model(),
            pathlib.Path(tmpdir) / "pars",
            redo=None,
            seed=seed,
            n_trees=n_pars_trees,
            **dict(threads=threads) if threads else {},
        )
        if pars_trees_file is not None:
            log_info and log_runtime_information(
                f"Storing the inferred parsimony trees in the file {pars_trees_file}."
            )
            shutil.copy(trees, pars_trees_file)

        log_info and log_runtime_information(
            "Computing the RF-Distance for the parsimony trees."
        )
        num_topos, rel_rfdist = raxmlng.get_rfdistance_results(trees, redo=None)

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


def predict_difficulty(
    msa_file: pathlib.Path,
    raxmlng: Optional[pathlib.Path] = DEFAULT_RAXMLNG_EXE,
    threads: int = None,
    seed: int = 0,
    file_format: Optional[FileFormat] = None,
    data_type: Optional[DataType] = None,
    deduplicate: bool = True,
    remove_full_gaps: bool = True,
    result_prefix: Optional[pathlib.Path] = None,
    store_results: bool = True,
    plot_shap: bool = False,
    model_file: pathlib.Path = DEFAULT_MODEL_FILE,
    log_info: bool = False,
) -> np.float64:
    """Predict the difficulty of an MSA using the PyPythia difficulty predictor.

    Per default, the MSA is deduplicated and full gap sequences are removed before the difficulty is predicted.

    Args:
        msa_file (pathlib.Path): Path to the MSA file. Note that the MSA file must be in either FASTA or PHYLIP format.
        raxmlng (pathlib.Path, optional): Path to the RAxML-NG executable.
            If not set, uses the RAxML-NG binary found in the PATH environment variable.
        threads (int, optional): Number of threads to use for parallel parsimony tree inference. If not set, uses the
            RAxML-NG auto parallelization scheme.
        seed (int, optional): Random seed to use for the parsimony tree inference. Defaults to 0.
        file_format (FileFormat, optional): File format of the MSA file. Defaults to None. In this case, the file format
            is inferred based on the file content. See `pypythia.msa.parse_msa` for information on when this is required.
        data_type (DataType, optional): Data type of the MSA sequences. Defaults to None. In this case, the data type
            is inferred based on the file content. See `pypythia.msa.parse_msa` for information on when this is required.
        deduplicate (bool, optional): If True, remove duplicate sequences from the MSA. Defaults to True.
        remove_full_gaps (bool, optional): If True, remove full gap sequences from the MSA. Defaults to True.
        result_prefix (pathlib.Path, optional): Prefix for the result files. Defaults to None. In this case, the prefix
            is set to the MSA file name.
        store_results (bool, optional): If True, store intermediate results as file. Defaults to True.
            In this case, the following files are stored:
            - The reduced MSA in PHYLIP format (if duplicates or full gap sequences were removed) in `{result_prefix}.reduced.phy`
            - The inferred parsimony trees in Newick format in `{result_prefix}.pythia.trees`
            - The shapley values as waterfall plot in `{result_prefix}.shap.pdf` (if plot_shap=True)
            - The features and predicted difficulty as CSV file in `{result_prefix}.pythia.csv`
        plot_shap (bool, optional): If True, plot the shapley values as waterfall plot. Defaults to False.
        model_file (pathlib.Path): Path to the trained difficulty predictor model.
            Defaults to the latest model shipped with PyPythia.
        log_info (bool, optional): If True, log intermediate progress information using the default logger. Defaults to False.

    Returns:
        np.float64: Predicted difficulty of the MSA.
    """
    if not msa_file.exists():
        raise PyPythiaException(f"The given MSA {msa_file} file does not exist.")

    if raxmlng is None:
        raise PyPythiaException(
            "Path to the RAxML-NG executable is required if 'raxml-ng' is not in $PATH."
        )

    result_prefix = pathlib.Path(result_prefix) if result_prefix else msa_file

    pars_trees_file = pathlib.Path(f"{result_prefix}.pythia.trees")
    shap_file = pathlib.Path(f"{result_prefix}.shap.pdf")
    results_file = pathlib.Path(f"{result_prefix}.pythia.csv")

    # We definitely need to store the reduced MSA somewhere for RAxML-NG
    if store_results:
        # If the user wants to keep the results, use the result_prefix
        reduced_msa_file = pathlib.Path(f"{result_prefix}.reduced.phy")
        _tmpfile = None
    else:
        # Else, use a temporary file
        _tmpfile = tempfile.NamedTemporaryFile(mode="w", suffix=".phy")
        reduced_msa_file = pathlib.Path(_tmpfile.name)

    log_info and log_runtime_information(
        message=f"Starting prediction for MSA {msa_file}."
    )

    # Init RAxML-NG
    try:
        raxmlng = RAxMLNG(**{"exe_path": raxmlng} if raxmlng else {})
    except Exception as e:
        raise PyPythiaException("Initializing RAxML-NG failed.") from e

    # Init the prediction model
    log_info and log_runtime_information(message=f"Loading predictor {model_file.name}")
    try:
        predictor = DifficultyPredictor(model_file=model_file)
    except Exception as e:
        raise PyPythiaException("Initializing the difficulty predictor failed.") from e

    # Load the MSA
    log_info and log_runtime_information(message="Loading MSA")
    msa = parse_msa(msa_file, file_format=file_format, data_type=data_type)

    # Deduplicate the MSA if necessary
    reduced_msa = _handle_duplicates(msa, deduplicate)

    # Remove full gap sequences if necessary
    reduced_msa = _handle_full_gap_sequences(reduced_msa, remove_full_gaps)

    # Check if the reduced MSA is different from the original MSA
    is_reduced = msa != reduced_msa
    if is_reduced:
        if reduced_msa.n_taxa < 4:
            raise PyPythiaException(
                "During preprocessing, Pythia reduced the input MSA by removing duplicate sequences and/or "
                "sequences containing only gaps leading to an MSA with less than 4 sequences. "
                "RAxML-NG refuses to infer trees for such small MSAs as there exists only a single possible tree topology. "
                "You can rerun the prediction and disable deduplication and gap removal to use the original MSA. "
            )

        # If the reduced MSA is different from the original MSA, proceed with the reduced MSA
        msa = reduced_msa

        log_info and log_runtime_information(
            "The input MSA contained duplicate sequences and/or sequences containing only gaps. "
            "WARNING: This predicted difficulty is only applicable to the reduced MSA (duplicate sequences removed). ",
        )

        # Save the reduced MSA
        msa_file = reduced_msa_file
        msa.write(msa_file)

        log_info and log_runtime_information(
            f"Saving a reduced alignment as {reduced_msa_file}.\n"
            f"We recommend to only use the reduced alignment {reduced_msa_file} for your subsequent analyses.\n",
        )

    # Compute the MSA Features
    log_info and log_runtime_information(
        f"Starting to compute MSA features for MSA {msa_file}"
    )

    log_info and log_runtime_information(
        "Number of threads not specified, using RAxML-NG autoconfig."
        if threads is None
        else f"Using {threads} threads for parallel parsimony tree inference."
    )

    msa_features = collect_features(
        msa=msa,
        msa_file=msa_file,
        raxmlng=raxmlng,
        pars_trees_file=pars_trees_file if store_results else None,
        log_info=log_info,
        threads=threads,
        seed=seed,
    )

    # Predict the difficulty
    log_info and log_runtime_information("Predicting the difficulty")
    difficulty = predictor.predict(msa_features)

    if plot_shap and store_results:
        # Plot shapley values
        # this only makes sense if store_results=True, otherwise the figure would be lost
        fig = predictor.plot_shapley_values(msa_features)
        fig.tight_layout()
        fig.savefig(fname=shap_file)

    log_info and log_runtime_information("Done")

    # Log the feature values
    if log_info:
        logger.info("─" * 20)
        logger.info("FEATURES: ")
        for feat, val in msa_features.items():
            logger.info(f"{feat}: {round(val[0], 2)}")

    if store_results:
        # Write the features + difficulty
        msa_features["difficulty"] = difficulty
        msa_features["msa_file"] = str(msa_file)

        msa_features.to_csv(results_file, index=False)

    if log_info:
        logger.info("")
        logger.info(f"Results: {results_file}.")
        is_reduced and logger.info(f"Reduced MSA: {reduced_msa_file}.")
        logger.info(f"Inferred parsimony trees: {pars_trees_file}.")
        logger.info(f"SHAP waterfall plot: {shap_file}.")
        logger.warning(
            "WARNING: When using shap plots, make sure you understand what shapley values are and how you can interpret"
            " this plot. For details refer to the wiki: https://github.com/tschuelia/PyPythia/wiki/Usage#shapley-values"
        )

    if _tmpfile is not None:
        # store_results was false, so we stored the reduced MSA in a temporary file, which we need to clean up
        _tmpfile.close()

    return difficulty[0]
