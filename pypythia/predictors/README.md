# Available trained predictors
You can use any of the provided trained predictors to query Pythia with. 
You can set the predictor to use by passing the respective path via the `-p` switch.
We recommend you to only use this option for reproducing older results as the default predictor will yield the most accurate results.

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
  - trained on 3250 empirical datasets obtained from TreeBase + 538 empirical datasets obtained via our RAxML-Grove
  - MAE = 0.09
  - MAPE = 2.5%



#### Abbreviations
- MAE = mean absolute error (prediction error on test data not seen during training)
- MAPE = mean absolute percentage error (prediction error on test data not seen during training)
