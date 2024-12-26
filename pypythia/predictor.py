import pathlib
import warnings
from typing import Optional

import lightgbm as lgb
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from numpy import typing as npt

from pypythia.custom_errors import PyPythiaException

DEFAULT_MODEL_FILE = pathlib.Path(__file__).parent / "predictors/latest.txt"


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
        """Predicts the difficulty for the given set of MSA features.
        TODO: adjust documentation -> also allows batch prediction!

        Args:
            query (Dict): Dict containing the features of the MSA to predict the difficulty for.
                query needs to contain at least the features the predictor was trained with.
                You can check this using the DifficultyPredictor.features attribute

        Returns:
            difficulty (float): The predicted difficulty for the given set of MSA features.

        Raises:
            PyPythiaException: If not all features the predictor was trained with are present in the given query.
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
