import pickle
import shap

from lightgbm import LGBMRegressor
from matplotlib.figure import Figure
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor

from pypythia.custom_types import *
from pypythia.custom_errors import PyPythiaException


class DifficultyPredictor:
    """Class structure for the trained difficulty predictor.

    This class provides methods for predicting the difficulty of an MSA.

    Args:
        predictor_handle (file object):
            Open file handle for the trained predictor. We do not guarantee the functionality of this class
            for predictors other than lightGBM Regressors and scikit-learn RandomForestRegressors
        features (optional list[string]):
            Names of the features the passed predictor was trained with.
            If you are using a LightGBM based predictor, the order of the features needs to be the same
            order as the order the predictor was trained with!
            If no list is passed and the predictor is either a lightGBM Regressor or scikit-learn RandomForestRegressor,
            the features will be automatically determined.
            For any other predictor type features cannot be None.

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
            if features is None:
                raise PyPythiaException(
                    "If passing a predictor other than a lightGMB Regressor or scikit-learn RandomForestRegressor, "
                    "you also need to pass a list of features the predictor was trained with."
                )
            self.features = features

    def _check_and_pack_query(self, query: Dict) -> pd.DataFrame:
        """
        Checks whether the given query is in correct format and packs the features as pandas Dataframe
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

        return df

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
        df = self._check_and_pack_query(query)

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

    def plot_shapley_values(self, query: Dict) -> Figure:
        explainer = shap.TreeExplainer(self.predictor)
        df = self._check_and_pack_query(query)
        shap_values = explainer.shap_values(df)
        base_values = explainer.expected_value

        if isinstance(self.predictor, RandomForestRegressor):
            base_values = base_values[0]

        return shap.plots.waterfall(
            shap.Explanation(
                values=shap_values[0],
                base_values=base_values,
                data=df.iloc[0]
            ),
            show=False
        )


