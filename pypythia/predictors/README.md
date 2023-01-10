# Available trained predictors
You can use any of the provided trained predictors to query Pythia with. 
You can set the predictor to use by passing the respective path via the `-p` switch.
We recommend you to only use this option for reproducing older results as the default predictor will yield the most accurate results.

- latest.pckl
  - equals the predictor `predictor_lgb_v1.0.0.pckl`

- predictor_sklearn_rf_v0.0.1.pckl
  - trained random forest predictor used as default predictor in Pythia versions <= 0.0.1
  - trained using `sklearn.ensemble.RandomForestRegressor`
  - predictor used to produce the results as presented in our Pythia publication (Preprint)
  - trained on 3250 empirical datasets obtained from TreeBase
  - MAE = 0.09
  - MAPE = 2.9%

- predictor_lgb_v1.0.0.pckl
  - trained boosted tree predictor used as default predictor in Pythia version 1.0.0
  - trained using `lightgbm.LGBMRegressor`
  - trained on 4262 empirical datasets
    - 3250 DNA + AA datasets obtained from TreeBase
    - 538 DNA + AA datasets obtained from our RAxML-Grove
    - 474 morphological datasets obtained from TreeBase
  - MAE = 0.09
  - MAPE = 2.7%
  - The model was trained using the following set of parameters (optimal set of parameters determined using the optuna optimizer):
    - learning_rate: 0.14337340122167785
    - max_depth: 8
    - lambda_l1: 5.5476436252362316e-08
    - lambda_l2: 0.001223580762903873
    - num_leaves: 26
    - bagging_fraction: 0.7529716269905308
    - bagging_freq: 6,
    - min_child_samples: 16

#### Abbreviations
- MAE = mean absolute error (prediction error on test data not seen during training)
- MAPE = mean absolute percentage error (prediction error on test data not seen during training)
