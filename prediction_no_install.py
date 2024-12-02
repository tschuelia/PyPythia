import pathlib

from pypythia.msa import MSA
from pypythia.prediction import get_all_features
from pypythia.predictor import DifficultyPredictor
from pypythia.raxmlng import RAxMLNG


def predict_difficulty(
    msa_file: pathlib.Path, predictor_path: pathlib.Path, raxmlng_exe_path: pathlib.Path
) -> float:
    """
    Predicts the difficulty of an MSA using the given difficulty predictor.

    Args:
        msa_file (FilePath): Path to the MSA file the difficulty should be predicted for. The file must be either in "fasta" or "phylip" format.
        predictor_path (FilePath): Path to a trained difficulty predictor.
        raxmlng_exe_path (Executable): Path to an executable of RAxML-NG. See https://github.com/amkozlov/raxml-ng for install instructions.

    Returns:
        difficulty (float): The predicted difficulty for the given MSA.

    Raises:
        ValueError: If the file format of the given MSA is not FASTA or PHYLIP.
        ValueError: If the data type of the given MSA cannot be inferred.
        PyPythiaException: If the provided difficulty predictor was trained with a subset incompatible to Pythia.
    """
    predictor = DifficultyPredictor(open(predictor_path, "rb"))
    raxmlng = RAxMLNG(raxmlng_exe_path)
    msa = MSA(msa_file)
    msa_features = get_all_features(raxmlng, msa, log_info=False)
    difficulty = predictor.predict(msa_features)
    print(f"The difficulty of your MSA {msa_file} is {round(difficulty, 2)}")
    return difficulty
