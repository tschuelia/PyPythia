import argparse
import pathlib
import sys
import time

from pypythia import __version__
from pypythia.config import DEFAULT_MODEL_FILE, DEFAULT_RAXMLNG_EXE
from pypythia.logger import get_header, logger
from pypythia.prediction import predict_difficulty


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
        "--shap",
        help="If set, computes the shapley values of the prediction as waterfall plot in '{prefix}.shap.pdf'. "
        "When using this option, make sure you understand what shapley values are and how to interpret this plot."
        "For details on shapley values refer to the documentation: "
        "https://tschuelia.github.io/PyPythia/latest/usage/#shap-waterfall-plot (default: False).",
        action="store_true",
    )

    parser.add_argument(
        "--forceDuplicates",
        help="Per default, Pythia refuses to predict the difficulty for MSAs containing duplicate sequences,"
        "and removes duplicate sequences prior to predicting the difficulty. "
        "Only set this option if you are absolutely sure that you want to predict the difficulty "
        "for this MSA (default: False). ",
        action="store_true",
    )

    parser.add_argument(
        "--forceFullGaps",
        help="Per default, Pythia refuses to predict the difficulty for MSAs containing sequences with only gaps,"
        "and removes full-gap sequences prior to predicting the difficulty. "
        "Only set this option if you are absolutely sure that you want to predict the difficulty "
        "for this MSA (default: False). ",
        action="store_true",
    )

    parser.add_argument(
        "--nofiles",
        help="Prevent Pythia from writing any files and only print logs/results to the terminal (default: False). "
        "WARNING: in this case and if your MSA contains duplicate/full-gap sequences the reduced MSA will not be stored.",
        action="store_true",
    )

    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=__version__,
        help="Print the version number and exit.",
    )

    return parser.parse_args()


def main():
    logger.info(get_header())
    args = _parse_cli()

    # Format all paths to pathlib.Path objects and set a default value if not provided
    msa_file = pathlib.Path(args.msa)
    prefix = pathlib.Path(args.prefix) if args.prefix else msa_file

    store_results = not args.nofiles

    if store_results:
        # Setup the logfile and log the Pythia header
        log_file = pathlib.Path(f"{prefix}.pythia.log")
        logger.add(log_file, format="{message}")
        log_file.write_text(get_header() + "\n")

    logger.info(
        f"Pythia was called at {time.strftime('%d-%b-%Y %H:%M:%S')} as follows:\n"
    )
    logger.info(" ".join(sys.argv))
    logger.info("")

    # Start the actual prediction
    script_start = time.perf_counter()

    difficulty = predict_difficulty(
        msa_file=msa_file,
        model_file=pathlib.Path(args.predictor),
        raxmlng=pathlib.Path(args.raxmlng),
        threads=args.threads,
        seed=args.seed,
        deduplicate=not args.forceDuplicates,
        remove_full_gaps=not args.forceFullGaps,
        result_prefix=prefix,
        store_results=store_results,
        plot_shap=args.shap,
        log_info=True,
    )

    script_end = time.perf_counter()

    logger.info("")
    total_runtime = script_end - script_start
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
        f"\nThe predicted difficulty for MSA {msa_file} is: {round(difficulty, 2)}\n"
    )


if __name__ == "__main__":
    main()
