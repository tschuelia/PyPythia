import argparse
import os.path
import shutil
import time
import logging
from tempfile import TemporaryDirectory

from pypythia.custom_types import *
from pypythia.custom_errors import PyPythiaException
from pypythia.msa import MSA
from pypythia.predictor import DifficultyPredictor
from pypythia.raxmlng import RAxMLNG
from pypythia import __version__

logging.basicConfig(level=logging.WARNING, format="%(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

SCRIPT_START = time.perf_counter()


def get_all_features(
    raxmlng: RAxMLNG,
    msa: MSA,
    store_trees: bool = False,
    log_info: bool = True,
    threads: int = None,
) -> Dict:
    """Helper function to collect all features required for predicting the difficulty of the MSA.

    Args:
        raxmlng (RAxMLNG): Initialized RAxMLNG object.
        msa (MSA): MSA object corresponding to the MSA file to compute the features for.
        store_trees (bool, optional): If True, store the inferred parsimony trees as "{msa_name}.parsimony.trees" file in the current workdir.
        log_info (bool, optional): If True, log intermediate progress information using the default logger.
        threads (int, optional): The number of threads to use for parallel parsimony tree inference. Uses the RAxML-NG auto parallelization scheme if none is set.
    Returns:
        all_features (Dict): Dictionary containing all features required for predicting the difficulty of the MSA. The keys correspond to the feature names the predictor was trained with.
    """
    with TemporaryDirectory() as tmpdir:
        msa_file = msa.msa_file
        model = msa.get_raxmlng_model()

        if log_info:
            log_runtime_information(
                "Retrieving num_patterns, percentage_gaps, percentage_invariant",
                log_runtime=True,
            )
        patterns, gaps, invariant = raxmlng.get_patterns_gaps_invariant(msa_file, model)

        if log_info:
            log_runtime_information("Retrieving num_taxa, num_sites", log_runtime=True)

        ntaxa = msa.number_of_taxa()
        nsites = msa.number_of_sites()

        n_pars_trees = 100
        if log_info:
            log_runtime_information(
                f"Inferring {n_pars_trees} parsimony trees", log_runtime=True
            )
        trees = raxmlng.infer_parsimony_trees(
            msa_file,
            model,
            os.path.join(tmpdir, "pars"),
            redo=None,
            seed=0,
            n_trees=n_pars_trees,
            **dict(threads=threads) if threads else {}
        )
        if store_trees:
            fn = f"{msa.msa_name}.parsimony.trees"
            log_runtime_information(
                f"Storing the inferred parsimony trees in the file {fn}"
            )
            shutil.copy(trees, fn)

        if log_info:
            log_runtime_information(
                "Computing the RF-Distance for the parsimony trees", log_runtime=True
            )
        num_topos, rel_rfdist, _ = raxmlng.get_rfdistance_results(trees, redo=None)

        return {
            "num_taxa": ntaxa,
            "num_sites": nsites,
            "num_patterns": patterns,
            "num_patterns/num_taxa": patterns / ntaxa,
            "num_sites/num_taxa": nsites / ntaxa,
            "num_patterns/num_sites": patterns / nsites,
            "proportion_gaps": gaps,
            "proportion_invariant": invariant,
            "entropy": msa.entropy(),
            "bollback": msa.bollback_multinomial(),
            "pattern_entropy": msa.pattern_entropy(),
            "avg_rfdist_parsimony": rel_rfdist,
            "proportion_unique_topos_parsimony": num_topos / n_pars_trees,
        }


def print_header():
    logger.info(
        f"PyPythia version {__version__} released by The Exelixis Lab\n"
        f"Developed by: Julia Haag\n"
        f"Latest version: https://github.com/tschuelia/PyPythia\n"
        f"Questions/problems/suggestions? Please open an issue on GitHub.\n",
    )


def log_runtime_information(message, log_runtime=True):
    if log_runtime:
        seconds = time.perf_counter() - SCRIPT_START
        fmt_time = time.strftime("%H:%M:%S", time.gmtime(seconds))
        time_string = f"[{fmt_time}] "
    else:
        time_string = ""
    logger.info(f"{time_string}{message}")


def main():
    print_header()

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
        required=True,
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
        "-p",
        "--predictor",
        type=argparse.FileType("rb"),
        default=os.path.join(
            os.path.dirname(__file__), "predictors/latest.pckl"
        ),
        required=False,
        help="Filepath of the predictor to use. If not set, "
        "assume it is 'predictors/latest.pckl' in the project directory.",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=argparse.FileType("w"),
        required=False,
        help="Option to specify a filepath where the result will be written to. "
        "The file will contain a single line with only the difficulty.",
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
        "--shap",
        help="If set, computes the shapley values of the prediction as waterfall plot in '{msa_name}.shap.pdf'. "
             "When using this option, make sure you understand what shapley values are and how to interpret this plot."
             "For details on shapley values refer to the wiki: https://github.com/tschuelia/PyPythia/wiki/Usage#shapley-values.",
        action="store_true"
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

    parser.add_argument(
        "-q",
        "--quiet",
        help="If set, Pythia does not print progress updates and only prints the predicted difficulty.",
        action="store_true",
    )

    args = parser.parse_args()

    if args.quiet:
        logger.setLevel(logging.WARNING)

    log_runtime_information(message="Starting prediction.", log_runtime=True)

    raxmlng_executable = args.raxmlng
    raxmlng = RAxMLNG(raxmlng_executable)

    log_runtime_information(
        message=f"Loading predictor {args.predictor.name}", log_runtime=True
    )
    predictor = DifficultyPredictor(args.predictor)

    log_runtime_information(message="Checking MSA", log_runtime=True)

    msa_file = args.msa
    msa = MSA(msa_file)
    final_warning_string = None

    if msa.contains_duplicate_sequences() and not args.removeDuplicates:
        raise PyPythiaException(
            "The provided MSA contains sequences that are exactly identical (duplicate sequences). "
            "Duplicate sequences influence the topological distances and distort the difficulty."
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
        msa.save_reduced_alignment(reduced_msa_file=reduced_msa, replace_original=True)
        msa_file = reduced_msa

        final_warning_string = (
            f"WARNING: This predicted difficulty is only applicable to the reduced MSA (duplicate sequences removed). "
            f"We recommend to only use the reduced alignment {msa_file} for your subsequent analyses.\n"
        )

    log_runtime_information(
        f"Starting to compute MSA features for MSA {msa_file}", log_runtime=True
    )

    if args.threads is not None:
        log_runtime_information(f"Using {args.threads} threads for parallel parsimony tree computation.", log_runtime=True)
    else:
        log_runtime_information(f"Number of threads not specified, using RAxML-NG autoconfig.", log_runtime=True)

    features_start = time.perf_counter()
    msa_features = get_all_features(
        raxmlng=raxmlng,
        msa=msa,
        store_trees=args.storeTrees,
        log_info=True,
        threads=args.threads,
    )
    features_end = time.perf_counter()

    log_runtime_information("Predicting the difficulty", log_runtime=True)

    prediction_start = time.perf_counter()
    difficulty = predictor.predict(msa_features)
    prediction_end = time.perf_counter()

    script_end = time.perf_counter()

    log_runtime_information("Done")

    logger.info(
        f"\nThe predicted difficulty for MSA {msa_file} is: {round(difficulty, args.precision)}\n"
    )

    if final_warning_string:
        logger.warning(final_warning_string)

    if args.shap:
        fig = predictor.plot_shapley_values(msa_features)
        fig.tight_layout()
        fig.savefig(fname=f"{msa.msa_name}.shap.pdf")
        logger.info(
            f"Waterfall plot of shapley values saved to {msa.msa_name}.shap.pdf"
        )
        logger.warning("WARNING: When using this plot make sure you understand what shapley values are and how you can interpret"
                       " this plot. For details refer to the wiki: https://github.com/tschuelia/PyPythia/wiki/Usage#shapley-values")

    if args.storeTrees:
        logger.info("─" * 20)
        logger.info(
            f"Inferred parsimony trees saved to {msa.msa_name}.parsimony.trees"
        )

    if args.verbose:
        logger.info("─" * 20)
        logger.info("FEATURES: ")
        for feat, val in msa_features.items():
            logger.info(f"{feat}: {round(val, args.precision)}")

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

    if args.output:
        args.output.write(str(round(difficulty, args.precision)))


if __name__ == "__main__":
    main()
