import matplotlib.figure
import numpy as np
import pandas as pd
import pytest

from pypythia.custom_errors import PyPythiaException


class TestPredictor:
    def test_str_and_repr(self, predictor):
        expected_str = (
            "DifficultyPredictor("
            "model_file=pypythia/predictors/latest.txt, features=['num_patterns/num_taxa', "
            "'num_sites/num_taxa', 'proportion_gaps', 'proportion_invariant', 'entropy', "
            "'bollback', 'num_patterns/num_sites', 'pattern_entropy', 'avg_rfdist_parsimony', "
            "'proportion_unique_topos_parsimony'])"
        )
        assert str(predictor) == expected_str
        assert repr(predictor) == expected_str

    def test_predict(self, predictor):
        query = pd.DataFrame(
            {
                "num_patterns/num_taxa": 0.0,
                "num_sites/num_taxa": 0.0,
                "num_patterns/num_sites": 0.0,
                "proportion_gaps": 0.0,
                "proportion_invariant": 0.0,
                "entropy": 0.0,
                "pattern_entropy": 0.0,
                "bollback": 0.0,
                "avg_rfdist_parsimony": 0.0,
                "proportion_unique_topos_parsimony": 0.0,
            },
            index=[0],
        )

        prediction = predictor.predict(query)[0]

        # don't test against an actual value, as this might change once the predictor changes
        # but the prediction should be a float between 0.0 and 1.0
        assert isinstance(prediction, float)
        assert 0.0 <= prediction <= 1.0

    def test_predict_with_multiple_queries(self, predictor):
        query = pd.DataFrame(
            {
                "num_patterns/num_taxa": [0.0, 1.0],
                "num_sites/num_taxa": [0.0, 1.0],
                "num_patterns/num_sites": [0.0, 1.0],
                "proportion_gaps": [0.0, 1.0],
                "proportion_invariant": [0.0, 1.0],
                "entropy": [0.0, 1.0],
                "pattern_entropy": [0.0, 1.0],
                "bollback": [0.0, 1.0],
                "avg_rfdist_parsimony": [0.0, 1.0],
                "proportion_unique_topos_parsimony": [0.0, 1.0],
            }
        )

        prediction = predictor.predict(query)

        assert prediction.shape == (2,)
        # don't test against an actual value, as this might change once the predictor changes
        # but the prediction should be a float between 0.0 and 1.0
        assert np.all(0.0 <= prediction) and np.all(prediction <= 1.0)

    def test_predict_with_missing_features_raises_pypythia_exception(self, predictor):
        query = pd.DataFrame(
            {
                "num_patterns/num_taxa": 0.0,
                "num_sites/num_taxa": 0.0,
            },
            index=[0],
        )

        with pytest.raises(
            PyPythiaException,
            match="The provided query does not contain all features the predictor was trained with.",
        ):
            predictor.predict(query)

    def test_predict_force_error(self, predictor):
        query = pd.DataFrame(
            {
                "num_patterns/num_taxa": "foo",  # The predictor can't handle string features
                "num_sites/num_taxa": 0.0,
                "num_patterns/num_sites": 0.0,
                "proportion_gaps": 0.0,
                "proportion_invariant": 0.0,
                "entropy": 0.0,
                "pattern_entropy": 0.0,
                "bollback": 0.0,
                "avg_rfdist_parsimony": 0.0,
                "proportion_unique_topos_parsimony": 0.0,
            },
            index=[0],
        )

        with pytest.raises(
            PyPythiaException,
            match="An error occurred predicting the difficulty for the provided set of MSA features.",
        ):
            predictor.predict(query)

    def test_plot_shapley_values(self, predictor):
        query = pd.DataFrame(
            {
                "num_patterns/num_taxa": 0.0,
                "num_sites/num_taxa": 0.0,
                "num_patterns/num_sites": 0.0,
                "proportion_gaps": 0.0,
                "proportion_invariant": 0.0,
                "entropy": 0.0,
                "pattern_entropy": 0.0,
                "bollback": 0.0,
                "avg_rfdist_parsimony": 0.0,
                "proportion_unique_topos_parsimony": 0.0,
            },
            index=[0],
        )

        fig = predictor.plot_shapley_values(query)
        assert isinstance(fig, matplotlib.figure.Figure)

    def test_plot_shapley_values_with_missing_features_raises_pypythia_exception(
        self, predictor
    ):
        query = pd.DataFrame(
            {
                "num_patterns/num_taxa": 0.0,
            },
            index=[0],
        )

        with pytest.raises(
            PyPythiaException,
            match="The provided query does not contain all features the predictor was trained with.",
        ):
            predictor.plot_shapley_values(query)
