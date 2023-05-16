import pickle

from lightgbm import LGBMRegressor
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor

from pypythia.custom_types import *
from pypythia.custom_errors import PyPythiaException


class DifficultyPredictor:
    """Class structure for the trained difficulty predictor.

    This class provides methods for predicting the difficulty of an MSA.

    Args:
        predictor_handle (file object): Open file handle for the trained predictor.
        features (optional list[string]):
            Names of the features the passed predictor was trained with.
            If you are using a LightGBM based predictor, the order of the features needs to be the same
            order as the order the predictor was trained with!

    Attributes:
        predictor: Loaded trained predictor.
        features: Names of the features the predictor was trained with.
    """

    def __init__(self, predictor_handle, features=None) -> None:
        self.predictor = pickle.load(predictor_handle)

        # Pythia version < 1.0.0 was a scikit-learn based predictor
        # starting from version 1.0.0 we have a lightGBM based predictor
        # so we need to distinguish the predictor type to obtain the feature names
        if isinstance(self.predictor, RandomForestRegressor):
            self.features = list(self.predictor.feature_names_in_)
        elif isinstance(self.predictor, LGBMRegressor):
            self.features = list(self.predictor.feature_name_)
        else:
            # hard-code the set of features we can use based on prediction.py
            self.features = [
                "num_patterns/num_taxa",
                "num_sites/num_taxa",
                "num_patterns/num_sites",
                "proportion_gaps",
                "proportion_invariant",
                "entropy",
                "pattern_entropy",
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
            PyPythiaException: If not all features the predictor was trained with are present in the given query.
        """
        df = pd.DataFrame()
        for feature in self.features:
            value = query.get(feature)
            if value is None:
                raise PyPythiaException(
                    f"The value for feature {feature} is not present in the query. "
                    f"Make sure to pass the correct set of features."
                    f"The required set of features is: {self.features}."
                )
            elif isinstance(value, list):
                if len(value) != 1:
                    raise PyPythiaException(
                        f"The value for feature {feature} is a list of length {len(value)}. "
                        f"Either provide a single value or a list with a single value only."
                    )
                else:
                    value = value[0]
            df[feature] = [value]

        df = df.reindex(columns=self.features)
        if np.all(np.isinf(df)):
            raise PyPythiaException("All features in this set are infinite. "
                                    "Something went wrong during feature computation.")

        try:
            if isinstance(self.predictor, LGBMRegressor):
                prediction = self.predictor.predict(df, num_threads=1)
            else:
                prediction = self.predictor.predict(df)

            prediction = prediction.clip(min=0.0, max=1.0)[0]
        except Exception as e:
            raise PyPythiaException(
                "An error occurred predicting the difficulty for the provided set of MSA features."
            ) from e

        return prediction