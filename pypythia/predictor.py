import pickle

import pandas as pd

from pypythia.custom_types import *


class DifficultyPredictor:
    """Class structure for the trained difficulty predictor.

    This class provides methods for predicting the difficulty of an MSA.

    Args:
        predictor_path (file object): Open file object for the trained predictor.

    Attributes:
        predictor: Loaded prediction algorithm.
        features: Names of the features the predictor was trained with.
    """

    def __init__(self, predictor_path) -> None:
        self.predictor = pickle.load(predictor_path)
        self.features = [
            "num_patterns/num_taxa",
            "num_sites/num_taxa",
            "proportion_gaps",
            "proportion_invariant",
            "entropy",
            "bollback",
            "avg_rfdist_parsimony",
            "proportion_unique_topos_parsimony"
        ]

    def predict(self, query: Dict) -> float:
        """ Predicts the difficulty for the given set of MSA features.

        Args:
            query (Dict): Dict containing the features of the MSA to predict the difficulty for.
                query needs to contain at least the features the predictor was trained with.
                rYou can check this using the DifficultyPredictor.features attribute

        Returns:
            difficulty (float): The predicted difficulty for the given set of MSA features.

        Raises:
            ValueError: If not all features the predictor was trained with are present in the given query.
        """
        df = pd.DataFrame()
        for feature in self.features:
            value = query.get(feature)
            if value is None:
                raise ValueError(
                    f"The value for feature {feature} is not present in the query. Make sure to pass the correct set of features.")
            df[feature] = [value]

        return self.predictor.predict(df)[0]
