import matplotlib.figure

from tests.fixtures import *
import numpy as np

from pypythia.custom_errors import PyPythiaException


class TestPredictor:
    def test_predict(self, predictor):
        query = {
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
        }

        prediction = predictor.predict(query)

        # don't test against an actual value, as this might change once the predictor changes
        # but the prediction should be a float between 0.0 and 1.0
        assert isinstance(prediction, float)
        assert 0.0 <= prediction <= 1.0

    def test_predict_with_list_values(self, predictor):
        query = {
            "num_patterns/num_taxa": [0.0],
            "num_sites/num_taxa": [0.0],
            "num_patterns/num_sites": [0.0],
            "proportion_gaps": [0.0],
            "proportion_invariant": 0.0,
            "entropy": 0.0,
            "pattern_entropy": 0.0,
            "bollback": 0.0,
            "avg_rfdist_parsimony": 0.0,
            "proportion_unique_topos_parsimony": 0.0,
        }

        prediction = predictor.predict(query)

        # don't test against an actual value, as this might change once the predictor changes
        # but the prediction should be a float between 0.0 and 1.0
        assert isinstance(prediction, float)
        assert 0.0 <= prediction <= 1.0

    def test_predict_with_missing_features_raises_pypythia_exception(self, predictor):
        query = {
            "num_patterns/num_taxa": 0.0,
            "num_sites/num_taxa": 0.0,
        }

        with pytest.raises(PyPythiaException):
            predictor.predict(query)

    def test_predict_with_list_of_features_raises_pypythia_exception(self, predictor):
        query = {
            "num_patterns/num_taxa": [0.0, 1.0],
            "num_sites/num_taxa": 0.0,
        }

        with pytest.raises(PyPythiaException):
            predictor.predict(query)

    def test_prediction_fails_with_pypythia_exception(self, predictor):
        query = {
            "num_patterns/num_taxa": -np.inf,
            "num_sites/num_taxa": -np.inf,
            "num_patterns/num_sites": -np.inf,
            "proportion_gaps": -np.inf,
            "proportion_invariant": -np.inf,
            "entropy": -np.inf,
            "bollback": -np.inf,
            "pattern_entropy": -np.inf,
            "avg_rfdist_parsimony": -np.inf,
            "proportion_unique_topos_parsimony": -np.inf,
        }

        with pytest.raises(PyPythiaException):
            predictor.predict(query)

    def test_plot_shapley_values(self, predictor):
        query = {
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
        }

        fig = predictor.plot_shapley_values(query)
        assert isinstance(fig, matplotlib.figure.Figure)


class TestBackwardsSklearnCompatibility:
    def test_predict(self, sklearn_predictor):
        query = {
            "num_patterns/num_taxa": 0.0,
            "num_sites/num_taxa": 0.0,
            "proportion_gaps": 0.0,
            "proportion_invariant": 0.0,
            "entropy": 0.0,
            "bollback": 0.0,
            "avg_rfdist_parsimony": 0.0,
            "proportion_unique_topos_parsimony": 0.0,
        }

        prediction = sklearn_predictor.predict(query)

        # don't test against an actual value, as this might change once the predictor changes
        # but the prediction should be a float between 0.0 and 1.0
        assert isinstance(prediction, float)
        assert 0.0 <= prediction <= 1.0

    def test_plot_shapley_values(self, sklearn_predictor):
        query = {
            "num_patterns/num_taxa": 0.0,
            "num_sites/num_taxa": 0.0,
            "proportion_gaps": 0.0,
            "proportion_invariant": 0.0,
            "entropy": 0.0,
            "bollback": 0.0,
            "avg_rfdist_parsimony": 0.0,
            "proportion_unique_topos_parsimony": 0.0,
        }

        fig = sklearn_predictor.plot_shapley_values(query)
        assert isinstance(fig, matplotlib.figure.Figure)