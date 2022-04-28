import warnings

from custom_types import *
import pickle
import pandas as pd


class DifficultyPredictor:
    def __init__(self, pickle_path) -> None:
        self.pickle_path = pickle_path
        self.predictor = pickle.load(self.pickle_path)
        self.features = [
            "num_patterns/num_taxa",
            "num_sites/num_taxa",
            "proportion_gaps",
            "proportion_invariant",
            "entropy",
            "bollback",
            "avg_rfdist_parsimony",
        ]

    def predict(self, query: Dict) -> float:
        df = pd.DataFrame()
        for feature in self.features:
            value = query.get(feature)
            if value is None:
                warnings.warn(f"The value for feature {feature} is not present in the query. Make sure to pass the correct set of features.")
            df[feature] = [value]

        return self.predictor.predict(df)[0]
