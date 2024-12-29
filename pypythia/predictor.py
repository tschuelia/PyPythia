import pathlib
import warnings
from typing import Optional

import lightgbm as lgb
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from numpy import typing as npt

from pypythia.config import DEFAULT_MODEL_FILE
from pypythia.custom_errors import PyPythiaException


class DifficultyPredictor:
    """Class structure for the trained difficulty predictor.

    This class provides methods for predicting the difficulty and plot the shapley values for an MSA.

    Args:
        model_file (pathlib.Path, optional): Path to the trained difficulty predictor model.
            Defaults to the latest model shipped with PyPythia.
            Note that this model file must be in the LightGBM .txt format.
        features (list[str], optional): Names of the features the predictor was trained with.
            Defaults to None. In this case, the features are inferred from the model file.

    Attributes:
        predictor (lgb.Booster): The trained LightGBM model used for predicting the difficulty.
        features (list[str]): Names of the features the predictor was trained with.
    """

    def __init__(
        self,
        model_file: Optional[pathlib.Path] = DEFAULT_MODEL_FILE,
        features: list[str] = None,
    ) -> None:
        self.model_file = model_file
        self.predictor = lgb.Booster(model_file=model_file)
        self.features = self.predictor.feature_name() if features is None else features

    def __str__(self):
        return f"DifficultyPredictor(model_file={self.model_file}, features={self.features})"

    def __repr__(self):
        return self.__str__()

    def _check_query(self, query: pd.DataFrame):
        if not set(self.features).issubset(query.columns):
            missing_features = set(self.features) - set(query.columns)
            raise PyPythiaException(
                "The provided query does not contain all features the predictor was trained with. "
                "Missing features: " + ", ".join(missing_features)
            )

    def predict(self, query: pd.DataFrame) -> npt.NDArray[np.float64]:
        """Predict the difficulty for a set of MSAs defined by rows in the given query dataframe.

        Args:
            query (pd.DataFrame): DataFrame containing the features for which to predict the difficulty.
                Each row in the DataFrame corresponds to a single MSA and the columns correspond to the features.

        Returns:
            A numpy array of predicted difficulties for the provided set of MSAs in float64 format.
            The difficulties are values in the range [0, 1] where higher values indicate higher difficulty.

        """
        self._check_query(query)

        try:
            prediction = self.predictor.predict(query[self.features])
            prediction = prediction.clip(min=0.0, max=1.0)
            return prediction
        except Exception as e:
            raise PyPythiaException(
                "An error occurred predicting the difficulty for the provided set of MSA features."
            ) from e

    def plot_shapley_values(self, query: pd.DataFrame) -> Figure:
        """Plot the shapley values for the **first** MSA in the given query dataframe.

        Please read our notes on SHAP values in the documentation to understand the plot.

        Args:
            query (pd.DataFrame): DataFrame containing the features for which to plot the shapley values.

        Returns:
            A matplotlib Figure object containing the waterfall plot of the shapley values for the first MSA in the query.
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            import shap

        self._check_query(query)

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
