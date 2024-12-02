# Available trained predictors
You can use any of the provided trained predictors to query Pythia with.
You can set the predictor to use by passing the respective path via the `-p` switch.
We recommend you to only use this option for reproducing older results as the default predictor will yield the most accurate results.

- latest.pckl
  - equals the predictor `predictor_lgb_v1.2.0.pckl`

- predictor_sklearn_rf_v0.0.1.pckl
  - trained random forest predictor used as default predictor in Pythia versions <= 0.0.1
  - trained using `sklearn.ensemble.RandomForestRegressor`
  - predictor used to produce the results as presented in our Pythia publication (Preprint)
  - trained on 3250 empirical datasets obtained from TreeBase
  - MAE = 0.09
  - MAPE = 2.9%

- predictor_lgb_v1.0.0.pckl
  - trained boosted tree predictor used as default predictor in Pythia version 1.0.*
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
    - bagging_freq: 6
    - min_child_samples: 16

- predictor_lgb_v1.1.0.pckl
  - trained boosted tree predictor used as default predictor in Pythia version 1.1.*
  - trained using `lightgbm.LGBMRegressor`
  - trained on 12 547 empirical datasets
    - 11 108 DNA MSAs
    - 979 AA MSAs
    - 460 morphological MSAs
  - MAE = 0.07
  - MAPE = 1.7%
  - The model was trained using the following set of parameters (optimal set of parameters determined using the optuna optimizer):
    - learning_rate: 0.08145759700381745
    - max_depth: 10
    - lambda_l1: 0.0002391390668080122
    - lambda_l2: 1.593402265932751e-08
    - num_leaves: 70
    - bagging_fraction: 0.9660277131424135
    - bagging_freq: 7
    - min_child_samples: 13

- predictor_lgb_v1.2.0.pckl
  - trained boosted tree predictor used as default predictor in Pythia version 1.2.*
  - trained using `lightgbm.LGBMRegressor`
  - trained on 12 547 empirical datasets
    - 11 108 DNA MSAs
    - 979 AA MSAs
    - 460 morphological MSAs
  - MAE = 0.06
  - MAPE =
  - The model was trained using the following set of parameters (optimal set of parameters determined using the optuna optimizer):
    - learning_rate: 0.1149633834134289
    - max_depth: 10
    - lambda_l1: 0.23830191835616107
    - lambda_l2: 1.3926298656531715
    - num_leaves: 241
    - bagging_fraction: 0.7763797800146144
    - bagging_freq: 1,
    - min_child_samples: 19

#### Abbreviations
- MAE = mean absolute error (prediction error on test data not seen during training)
- MAPE = mean absolute percentage error (prediction error on test data not seen during training)
