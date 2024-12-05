import pathlib
import warnings

import lightgbm as lgb
import pandas as pd
from matplotlib.figure import Figure

from pypythia.custom_errors import PyPythiaException


class DifficultyPredictor:
    """Class structure for the trained difficulty predictor.

    This class provides methods for predicting the difficulty of an MSA.

    Args:
        predictor_file (file object):
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

    def __init__(self, model_file: pathlib.Path, features: list[str] = None) -> None:
        self.predictor = lgb.Booster(model_file=model_file)
        self.features = self.predictor.feature_name() if features is None else features

    def predict(self, query: pd.DataFrame) -> float:
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
        if not set(self.features).issubset(query.columns):
            raise PyPythiaException(
                "The provided query does not contain all features the predictor was trained with."
            )

        try:
            prediction = self.predictor.predict(query[self.features])
            prediction = prediction.clip(min=0.0, max=1.0)
        except Exception as e:
            raise PyPythiaException(
                "An error occurred predicting the difficulty for the provided set of MSA features."
            ) from e

        return prediction

    def plot_shapley_values(self, query: pd.DataFrame) -> Figure:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            import shap

        if not set(self.features).issubset(query.columns):
            raise PyPythiaException(
                "The provided query does not contain all features the predictor was trained with."
            )

        df = query[self.features]

        explainer = shap.TreeExplainer(self.predictor)
        shap_values = explainer.shap_values(df)
        base_values = explainer.expected_value

        return shap.plots.waterfall(
            shap.Explanation(
                values=shap_values[0], base_values=base_values, data=df.iloc[0]
            ),
            show=False,
        ).figure
