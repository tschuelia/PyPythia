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
            "proportion_unique_topos_parsimony",
        ]

    def predict(self, query: Dict) -> float:
        """Predicts the difficulty for the given set of MSA features.

        Args:
            query (Dict): Dict containing the features of the MSA to predict the difficulty for.
                query needs to contain at least the features the predictor was trained with.
                You can check this using the DifficultyPredictor.features attribute

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
                    f"The value for feature {feature} is not present in the query. "
                    f"Make sure to pass the correct set of features."
                    f"The required set of features is: {self.features}."
                )
            elif isinstance(value, list):
                if len(value) != 1:
                    raise ValueError(f"The value for feature {feature} is a list of length {len(value)}. "
                                     f"Either provide a single value or a list with a single value only.")
                else:
                    value = value[0]
            df[feature] = [value]

        try:
            prediction = self.predictor.predict(df)[0]
        except Exception as e:
            raise RuntimeError(
                "An error occurred predicting the difficulty for the provided set of MSA features."
            ) from e

        return prediction
