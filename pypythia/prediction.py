import argparse
import os.path
import shutil
import textwrap
from tempfile import TemporaryDirectory
import time

from pypythia.custom_types import *
from pypythia.msa import MSA
from pypythia.predictor import DifficultyPredictor
from pypythia.raxmlng import RAxMLNG


def get_all_features(
    raxmlng: RAxMLNG,
    msa: MSA,
    store_trees: bool = False,
) -> Dict:
    """Helper function to collect all features required for predicting the difficulty of the MSA.

    Args:
        raxmlng (RAxMLNG): Initialized RAxMLNG object.
        msa (MSA): MSA object corresponding to the MSA file to compute the features for.
        store_trees (bool): If True, store the inferred parsimony trees as "{msa_name}.parsimony.trees" file in the current workdir.
    Returns:
        all_features (Dict): Dictionary containing all features required for predicting the difficulty of the MSA. The keys correspond to the feature names the predictor was trained with.
    """
    with TemporaryDirectory() as tmpdir:
        msa_file = msa.msa_file
        model = msa.get_raxmlng_model()
        patterns, gaps, invariant = raxmlng.get_patterns_gaps_invariant(msa_file, model)

        ntaxa = msa.number_of_taxa()
        nsites = msa.number_of_sites()

        n_pars_trees = 100
        trees = raxmlng.infer_parsimony_trees(
            msa_file, model, tmpdir, redo=None, seed=0, n_trees=n_pars_trees
        )
        num_topos, rel_rfdist, _ = raxmlng.get_rfdistance_results(trees, redo=None)

        if store_trees:
            shutil.copy(trees, f"{msa.msa_name}.parsimony.trees")

        return {
            "num_taxa": ntaxa,
            "num_sites": nsites,
            "num_patterns": patterns,
            "num_patterns/num_taxa": patterns / ntaxa,
            "num_sites/num_taxa": nsites / ntaxa,
            "proportion_gaps": gaps,
            "proportion_invariant": invariant,
            "entropy": msa.entropy(),
            "bollback": msa.bollback_multinomial(),
            "avg_rfdist_parsimony": rel_rfdist,
            "proportion_unique_topos_parsimony": num_topos / n_pars_trees,
        }


def main():
    parser = argparse.ArgumentParser(
        description="Parser for optional config file setting."
    )
    parser.add_argument(
        "--msa",
        type=str,
        required=True,
        help="Multiple Sequence Alignment to predict the difficulty for. Must be in either phylip or fasta format.",
    )

    parser.add_argument(
        "--raxmlng",
        type=str,
        required=True,
        help="Path to the binary of RAxML-NG. For install instructions see https://github.com/amkozlov/raxml-ng.",
    )

    parser.add_argument(
        "--predictor",
        type=argparse.FileType("rb"),
        default=os.path.join(os.path.dirname(__file__), "predictor.pckl"),
        required=False,
        help="Filepath of the predictor to use. If not set, assume it is 'predictor.pckl' in the project directory.",
    )

    parser.add_argument(
        "--storeTrees",
        help="If set, stores the parsimony trees as '{msa_name}.parsimony.trees' file",
        action="store_true",
    )

    parser.add_argument(
        "--verbose",
        help="If set, prints the MSA features",
        action="store_true",
    )

    parser.add_argument(
        "--benchmark",
        help="If set, time the runtime of the prediction",
        action="store_true",
    )

    args = parser.parse_args()

    script_start = time.perf_counter()

    raxmlng_executable = args.raxmlng
    raxmlng = RAxMLNG(raxmlng_executable)
    msa_file = args.msa
    predictor = DifficultyPredictor(args.predictor)

    try:
        msa = MSA(msa_file)
    except Exception as e:
        raise RuntimeError("Error reading the provided MSA: ", msa_file) from e

    features_start = time.perf_counter()
    msa_features = get_all_features(raxmlng, msa, args.storeTrees)
    features_end = time.perf_counter()

    prediction_start = time.perf_counter()
    difficulty = predictor.predict(msa_features)
    prediction_end = time.perf_counter()

    script_end = time.perf_counter()

    print(f"The predicted difficulty for MSA {msa_file} is: {round(difficulty, 2)}.")

    if args.benchmark:
        feature_time = round(features_end - features_start, 3)
        prediction_time = round(prediction_end - prediction_start, 3)
        runtime_script = round(script_end - script_start, 3)

        print(textwrap.dedent(
                f"""
                === RUNTIME SUMMARY:
                Feature computation runtime:\t{feature_time} seconds
                Prediction:\t\t\t{prediction_time} seconds
                ----
                Total script runtime:\t\t{runtime_script} seconds
                
                Note that all runtimes include the overhead for python and python subprocess calls. 
                For a more accurate and fine grained benchmark, call the respective feature computations from code.
                """
            ))
