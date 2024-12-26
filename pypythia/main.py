import argparse
import pathlib
import time

from pypythia.custom_errors import PyPythiaException
from pypythia.logger import SCRIPT_START, get_header, log_runtime_information, logger
from pypythia.msa import parse
from pypythia.prediction import collect_features
from pypythia.predictor import DEFAULT_MODEL_FILE, DifficultyPredictor
from pypythia.raxmlng import DEFAULT_RAXMLNG_EXE, RAxMLNG


def _setup_argparse() -> argparse.ArgumentParser:
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
        help="Path to the binary of RAxML-NG. For install instructions see https://github.com/amkozlov/raxml-ng.",
    )

    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=False,
        help="Number of threads to use for parallel parsimony tree inference. "
        "If none is set, Pythia uses the parallelization scheme of RAxML-NG "
        "that automatically detects the optimal number of threads for your machine.",
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
        help="Filepath of the alternative predictor to use. Uses the latest Pythia predictor per default.",
    )

    parser.add_argument(
        "-prec",
        "--precision",
        type=int,
        default=2,
        required=False,
        help="Set the number of decimals the difficulty should be rounded to. Recommended and default is 2.",
    )

    parser.add_argument(
        "-sT",
        "--storeTrees",
        help="If set, stores the parsimony trees as '{msa_name}.parsimony.trees' file.",
        action="store_true",
    )

    parser.add_argument(
        "--removeDuplicates",
        help="Pythia refuses to predict the difficulty for MSAs containing duplicate sequences. "
        "If this option is set, PyPythia removes the duplicate sequences, "
        "stores the reduced MSA as '{msa_name}.{phy/fasta}.pythia.reduced' "
        "and predicts the difficulty for the reduced alignment.",
        action="store_true",
    )

    parser.add_argument(
        "--forceDuplicates",
        help="Per default, Pythia refuses to predict the difficulty for MSAs containing duplicate sequences. "
        "Set this option if you are absolutely sure that you want to predict the difficulty for this MSA. ",
        action="store_true",
    )

    parser.add_argument(
        "--shap",
        help="If set, computes the shapley values of the prediction as waterfall plot in '{msa_name}.shap.pdf'. "
        "When using this option, make sure you understand what shapley values are and how to interpret this plot."
        "For details on shapley values refer to the wiki: https://github.com/tschuelia/PyPythia/wiki/Usage#shapley-values.",
        action="store_true",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        help="If set, additionally prints the MSA features.",
        action="store_true",
    )

    parser.add_argument(
        "-b",
        "--benchmark",
        help="If set, time the runtime of the prediction.",
        action="store_true",
    )

    return parser


def main():
    logger.info(get_header())
    parser = _setup_argparse()
    args = parser.parse_args()

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

    log_runtime_information(message="Starting prediction.", log_runtime=True)

    raxmlng = RAxMLNG(raxmlng_executable)

    log_runtime_information(
        message=f"Loading predictor {predictor_file.name}", log_runtime=True
    )
    predictor = DifficultyPredictor(predictor_file)

    log_runtime_information(message="Checking MSA", log_runtime=True)

    msa = parse(msa_file)
    final_warning_string = None

    if msa.contains_duplicate_sequences() and not (
        args.removeDuplicates or args.forceDuplicates
    ):
        raise PyPythiaException(
            "The provided MSA contains sequences that are exactly identical (duplicate sequences). "
            "Duplicate sequences influence the topological distances and distort the difficulty. "
            "If you are absolutely sure that you want to predict the difficulty for this MSA, "
            "set the option --forceDuplicates."
        )

    if not msa.contains_duplicate_sequences() and args.removeDuplicates:
        logger.warning(
            "WARNING: The provided MSA does not contain duplicate sequences. "
            "The setting 'removeDuplicates' has no effect."
        )

    if msa.contains_duplicate_sequences() and args.removeDuplicates:
        reduced_msa = msa_file + ".pythia.reduced"
        log_runtime_information(
            f"The input alignment {msa_file} contains duplicate sequences: "
            f"saving a reduced alignment as {reduced_msa}\n",
            log_runtime=True,
        )
        # TODO: save reduced MSA
        # msa.save_reduced_alignment(reduced_msa_file=reduced_msa, replace_original=True)
        # msa_file = reduced_msa

        final_warning_string = (
            f"WARNING: This predicted difficulty is only applicable to the reduced MSA (duplicate sequences removed). "
            f"We recommend to only use the reduced alignment {msa_file} for your subsequent analyses.\n"
        )

    if msa.contains_duplicate_sequences() and args.forceDuplicates:
        logger.warning(
            "WARNING: The provided MSA contains duplicate sequences. "
            "The setting 'forceDuplicates' is set, Pythia will predict the difficulty for the MSA with duplicates."
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

    if final_warning_string:
        logger.warning(final_warning_string)

    if args.shap:
        fig = predictor.plot_shapley_values(msa_features)
        fig.tight_layout()
        fig.savefig(fname=shap_file)

    if args.verbose:
        logger.info("─" * 20)
        logger.info("FEATURES: ")
        for feat, val in msa_features.items():
            logger.info(f"{feat}: {round(val[0], args.precision)}")

    if args.benchmark:
        feature_time = round(features_end - features_start, 3)
        prediction_time = round(prediction_end - prediction_start, 3)
        runtime_script = round(script_end - SCRIPT_START, 3)

        logger.info(
            f"{'─' * 20}\n"
            f"RUNTIME SUMMARY:\n"
            f"Feature computation runtime:\t{feature_time} seconds\n"
            f"Prediction:\t\t\t{prediction_time} seconds\n"
            f"----\n"
            f"Total script runtime:\t\t{runtime_script} seconds\n"
            f"Note that all runtimes include the overhead for python and python subprocess calls.\n"
            f"For a more accurate and fine grained benchmark, call the respective feature computations from code."
        )

        msa_features["feature_computation_time"] = feature_time
        msa_features["prediction_time"] = prediction_time
        msa_features["script_runtime"] = runtime_script

    msa_features["difficulty"] = difficulty
    msa_features["msa_file"] = str(msa_file)

    msa_features.to_csv(results_file, index=False)
    logger.info("")
    logger.info(f"Results saved to {results_file}.")

    if args.storeTrees:
        logger.info(f"Inferred parsimony trees saved to {pars_trees_file}")

    if args.shap:
        logger.info(f"Waterfall plot of shapley values saved to {shap_file}")
        logger.warning(
            "WARNING: When using shap plots, make sure you understand what shapley values are and how you can interpret"
            " this plot. For details refer to the wiki: https://github.com/tschuelia/PyPythia/wiki/Usage#shapley-values"
        )

    logger.info(
        f"\nThe predicted difficulty for MSA {msa_file} is: {round(difficulty, args.precision)}\n"
    )


if __name__ == "__main__":
    main()
