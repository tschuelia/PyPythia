import argparse
import pathlib
import sys
import time

from pypythia.config import DEFAULT_MODEL_FILE, DEFAULT_RAXMLNG_EXE
from pypythia.logger import get_header, log_runtime_information, logger
from pypythia.msa import MSA, deduplicate_sequences, parse, remove_full_gap_sequences
from pypythia.prediction import collect_features
from pypythia.predictor import DifficultyPredictor
from pypythia.raxmlng import RAxMLNG


def _parse_cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Parser for Pythia command line options."
    )

    parser.add_argument(
        "-m",
        "--msa",
        type=str,
        required=True,
        help="Multiple Sequence Alignment to predict the difficulty for. Must be in either phylip or fasta format.",
    )

    parser.add_argument(
        "-r",
        "--raxmlng",
        type=str,
        default=DEFAULT_RAXMLNG_EXE,
        required=DEFAULT_RAXMLNG_EXE is None,
        help="Path to the binary of RAxML-NG. For install instructions see https://github.com/amkozlov/raxml-ng."
        "(default: 'raxml-ng' if in $PATH, otherwise this option is mandatory).",
    )

    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=False,
        help="Number of threads to use for parallel parsimony tree inference (default: RAxML-NG autoconfig).",
    )

    parser.add_argument(
        "-s",
        "--seed",
        type=int,
        default=0,
        required=False,
        help="Seed for the RAxML-NG parsimony tree inference (default: 0).",
    )

    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        required=False,
        help="Prefix of the PyPythia log and result file (default: MSA file name).",
    )

    parser.add_argument(
        "--predictor",
        type=str,
        default=DEFAULT_MODEL_FILE,
        required=False,
        help="Filepath of the alternative predictor to use (default: latest Pythia).",
    )

    parser.add_argument(
        "-prec",
        "--precision",
        type=int,
        default=2,
        required=False,
        help="Set the number of decimals the difficulty should be rounded to (default: 2).",
    )

    parser.add_argument(
        "-sT",
        "--storeTrees",
        help="If set, stores the parsimony trees as '{prefix}.pythia.trees' file (default: False).",
        action="store_true",
    )

    parser.add_argument(
        "--forceDuplicates",
        help="Per default, Pythia refuses to predict the difficulty for MSAs containing duplicate sequences. "
        "Only set this option if you are absolutely sure that you want to predict the difficulty "
        "for this MSA (default: False). ",
        action="store_true",
    )

    parser.add_argument(
        "--forceFullGaps",
        help="Per default, Pythia refuses to predict the difficulty for MSAs containing sequences with only gaps. "
        "Only set this option if you are absolutely sure that you want to predict the difficulty "
        "for this MSA (default: False). ",
        action="store_true",
    )

    parser.add_argument(
        "--shap",
        help="If set, computes the shapley values of the prediction as waterfall plot in '{prefix}.shap.pdf'. "
        "When using this option, make sure you understand what shapley values are and how to interpret this plot."
        "For details on shapley values refer to the wiki: "
        "https://github.com/tschuelia/PyPythia/wiki/Usage#shapley-values (default: False).",
        action="store_true",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        help="If set, additionally prints the MSA features (default: False).",
        action="store_true",
    )

    return parser.parse_args()


def _handle_duplicates(msa: MSA, force_duplicates: bool) -> MSA:
    if msa.contains_duplicate_sequences() and force_duplicates:
        logger.warning(
            "WARNING: The provided MSA contains duplicate sequences. "
            "The setting 'forceDuplicates' is set, Pythia will predict the difficulty for the MSA with duplicates."
        )
        return msa

    if msa.contains_duplicate_sequences():
        log_runtime_information(
            "The input MSA contains duplicate sequences. Removing duplicates before predicting the difficulty."
        )
        return deduplicate_sequences(msa)
    else:
        return msa


def _handle_full_gap_sequences(msa: MSA, force_full_gaps: bool) -> MSA:
    if msa.contains_full_gap_sequences() and force_full_gaps:
        log_runtime_information(
            "WARNING: The provided MSA contains sequences with only gaps. "
            "The setting 'forceFullGaps' is set, Pythia will predict the difficulty for the MSA with full gap sequences."
        )
        return msa

    if msa.contains_full_gap_sequences():
        log_runtime_information(
            "The input MSA contains sequences with only gaps. Removing full gap sequences before predicting the difficulty."
        )
        return remove_full_gap_sequences(msa)
    else:
        return msa


def main():
    logger.info(get_header())
    args = _parse_cli()

    # Format all paths to pathlib.Path objects and set a default value if not provided
    msa_file = pathlib.Path(args.msa)
    raxmlng_executable = pathlib.Path(args.raxmlng)
    prefix = pathlib.Path(args.prefix) if args.prefix else msa_file
    predictor_file = pathlib.Path(args.predictor)

    # Setup all result files, log the Pythia header
    log_file = pathlib.Path(f"{prefix}.pythia.log")
    logger.add(log_file, format="{message}")
    log_file.write_text(get_header() + "\n")

    results_file = pathlib.Path(f"{prefix}.pythia.csv")
    pars_trees_file = pathlib.Path(f"{prefix}.pythia.trees")
    shap_file = pathlib.Path(f"{prefix}.shap.pdf")
    reduced_msa_file = pathlib.Path(f"{prefix}.reduced.phy")

    # Start the actual prediction
    SCRIPT_START = time.perf_counter()

    logger.info(
        f"Pythia was called at {time.strftime('%d-%b-%Y %H:%M:%S')} as follows:\n"
    )
    logger.info(" ".join(sys.argv))
    logger.info("")

    log_runtime_information(
        message=f"Starting prediction for MSA {msa_file}.", log_runtime=True
    )

    raxmlng = RAxMLNG(raxmlng_executable)

    log_runtime_information(
        message=f"Loading predictor {predictor_file.name}", log_runtime=True
    )
    predictor = DifficultyPredictor(predictor_file)

    log_runtime_information(message="Loading MSA", log_runtime=True)

    msa = parse(msa_file)

    # First, deduplicate the MSA if necessary
    reduced_msa = _handle_duplicates(msa, args.forceDuplicates)

    # Second, remove full gap sequences if necessary
    reduced_msa = _handle_full_gap_sequences(reduced_msa, args.forceFullGaps)

    # check if the reduced MSA is different from the original MSA
    is_reduced = msa != reduced_msa
    if is_reduced:
        reduced_msa.write(reduced_msa_file)
        msa = reduced_msa
        msa_file = reduced_msa_file

        log_runtime_information(
            "The input MSA contained duplicate sequences and/or sequences containing only gaps. "
            f"Saving a reduced alignment as {reduced_msa_file}.\n"
            "WARNING: This predicted difficulty is only applicable to the reduced MSA (duplicate sequences removed). "
            f"We recommend to only use the reduced alignment {reduced_msa_file} for your subsequent analyses.\n",
            log_runtime=True,
        )

    log_runtime_information(
        f"Starting to compute MSA features for MSA {msa_file}", log_runtime=True
    )

    if args.threads is None:
        log_runtime_information(
            "Number of threads not specified, using RAxML-NG autoconfig.",
            log_runtime=True,
        )
    else:
        log_runtime_information(
            f"Using {args.threads} threads for parallel parsimony tree computation.",
            log_runtime=True,
        )

    features_start = time.perf_counter()
    msa_features = collect_features(
        msa=msa,
        msa_file=msa_file,
        raxmlng=raxmlng,
        pars_trees_file=pars_trees_file if args.storeTrees else None,
        log_info=True,
        threads=args.threads,
        seed=args.seed,
    )
    features_end = time.perf_counter()

    log_runtime_information("Predicting the difficulty", log_runtime=True)

    prediction_start = time.perf_counter()
    difficulty = predictor.predict(msa_features)[0]
    prediction_end = time.perf_counter()

    script_end = time.perf_counter()

    log_runtime_information("Done")

    if args.shap:
        fig = predictor.plot_shapley_values(msa_features)
        fig.tight_layout()
        fig.savefig(fname=shap_file)

    if args.verbose:
        logger.info("─" * 20)
        logger.info("FEATURES: ")
        for feat, val in msa_features.items():
            logger.info(f"{feat}: {round(val[0], args.precision)}")

    msa_features["difficulty"] = difficulty
    msa_features["msa_file"] = str(msa_file)

    msa_features.to_csv(results_file, index=False)
    logger.info("")
    logger.info(f"Results: {results_file}.")

    if is_reduced:
        logger.info(f"Reduced MSA: {reduced_msa_file}.")

    if args.storeTrees:
        logger.info(f"Inferred parsimony trees: {pars_trees_file}.")

    if args.shap:
        logger.info(f"SHAP waterfall plot: {shap_file}.")
        logger.warning(
            "WARNING: When using shap plots, make sure you understand what shapley values are and how you can interpret"
            " this plot. For details refer to the wiki: https://github.com/tschuelia/PyPythia/wiki/Usage#shapley-values"
        )

    logger.info("")
    total_runtime = script_end - SCRIPT_START
    hours, remainder = divmod(total_runtime, 3600)
    minutes, seconds = divmod(remainder, 60)

    if hours > 0:
        logger.info(
            f"Total runtime: {int(hours):02d}:{int(minutes):02d}:{seconds:02d} hours ({round(total_runtime)} seconds)."
        )
    elif minutes > 0:
        logger.info(
            f"Total runtime: {int(minutes):02d}:{int(seconds):02d} minutes ({round(total_runtime)} seconds)."
        )
    else:
        logger.info(f"Total runtime: {seconds:.2f} seconds.")

    logger.info(
        f"\nThe predicted difficulty for MSA {msa_file} is: {round(difficulty, args.precision)}\n"
    )


if __name__ == "__main__":
    main()
