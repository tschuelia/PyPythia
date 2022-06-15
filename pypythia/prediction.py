import argparse
import os.path
import shutil
from tempfile import TemporaryDirectory

from pypythia.custom_types import *
from pypythia.msa import MSA
from pypythia.predictor import DifficultyPredictor
from pypythia.raxmlng import RAxMLNG


def get_all_features(raxmlng: RAxMLNG, msa: MSA, model: str, store_trees: bool = False) -> Dict:
    """Helper function to collect all features required for predicting the difficulty of the MSA.

    Args:
        raxmlng (RAxMLNG): Initialized RAxMLNG object.
        msa (MSA): MSA object corresponding to the MSA file to compute the features for.
        model (str): String representation of the substitution model to use. Needs to be a valid RAxML-NG model. For example "GTR+G" for DNA data or "LG+G" for protein data.
        store_trees (bool): If True, store the inferred parsimony trees as "{msa_name}.parsimony.trees" file in the current workdir.
    Returns:
        all_features (Dict): Dictionary containing all features required for predicting the difficulty of the MSA. The keys correspond to the feature names the predictor was trained with.
    """
    with TemporaryDirectory() as tmpdir:
        msa_file = msa.msa_file
        patterns, gaps, invariant = raxmlng.get_patterns_gaps_invariant(msa_file, model)

        ntaxa = msa.number_of_taxa()
        nsites = msa.number_of_sites()

        n_pars_trees = 100
        trees = raxmlng.infer_parsimony_trees(msa_file, model, tmpdir, redo=None, seed=0, n_trees=n_pars_trees)
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
            "proportion_unique_topos_parsimony": num_topos / n_pars_trees
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
        "--model",
        type=str,
        required=False,
        help="Model to use for the prediction. This can be either a model string (e.g. GTR+G) or a path to a partition file."
             "If not set the data type is automatically inferred, and the model is set to GTR+G for DNA MSAs and to LG+G for Protein MSAs.",
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

    args = parser.parse_args()

    raxmlng_executable = args.raxmlng
    raxmlng = RAxMLNG(raxmlng_executable)
    msa_file = args.msa
    msa = MSA(msa_file)

    if args.model:
        model = args.model
    else:
        data_type = msa.guess_data_type()
        model = "GTR+G" if data_type == "DNA" else "LG+G"

    predictor = DifficultyPredictor(args.predictor)

    msa_features = get_all_features(raxmlng, msa, model, args.storeTrees)

    if args.verbose:
        print("FEATURES: ")
        for feat, val in msa_features.items():
            print(f"{feat}: {round(val, 2)}")
        print("---------")

    if args.storeTrees:
        print(f"Inferred parsimony trees saved to {msa.msa_name}.parsimony.trees")
        print("---------")

    difficulty = predictor.predict(msa_features)
    print(
        f"The predicted difficulty for MSA {msa_file} is: {round(difficulty, 2)}."
    )
