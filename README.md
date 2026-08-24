# Joint Longitudinal Survival Models: Predictive Value of Repeated Measures for Graft and Patient Failure

> **Note:** This analysis was done in R (RStudio), not Python/Colab, so the scripts should run as-is given the data files, aside from updating the hardcoded file paths to your own machine.

## What this is

A coursework assignment on 838 kidney transplant recipients, evaluating whether *repeated* post-transplant biomarker measurements (kidney function, creatinine, blood pressure) improve prediction of graft failure beyond a standard baseline-only survival model, comparing a time-dependent Cox model against full joint longitudinal-survival models.

## What I did

- **Missing data imputation**: used MICE (`mice` package), restricted to the 21 repeated biomarker measurements (GFR, serum creatinine, MAP at 7 follow-up visits), generating 10 imputed datasets.
- **Time-dependent Cox regression**: built a counting-process (start/stop) dataset (`tmerge`) to model graft failure as a function of both fixed baseline covariates and the time-varying biomarker values, pooling results across the 10 imputed datasets with Rubin's rule.
- **Longitudinal (LME) sub-models**: for each biomarker (GFR, creatinine, MAP), compared three mixed-effects model structures (random intercept, random intercept + slope, random intercept + spline) via AIC to find the best-fitting trajectory shape, then checked residual normality and homoscedasticity.
- **Joint Cox-LME models**: fit three separate joint models (via the `JM` package, Weibull proportional-hazards submodel), each linking a biomarker's longitudinal trajectory to the hazard of graft failure, and checked model assumptions and convergence.
- **Model comparison**: compared effect estimates from the time-dependent Cox model against the three joint models, finding consistent direction and significance across all three biomarkers.
- **Dynamic, patient-specific prediction**: used the fitted joint models to generate individualized, updating graft-survival predictions for a specific patient as new biomarker measurements accumulate over time.

## Key finding

Higher GFR was associated with lower graft-failure risk, while higher creatinine and higher blood pressure (MAP) were associated with higher risk, consistent across both modeling approaches. The joint models additionally captured each biomarker's full trajectory shape and produced dynamic, updating survival predictions for individual patients, which a standard time-dependent Cox model cannot do on its own.

## Repo contents

- **`mice_imputation_cox_models.R`**: MICE imputation and pooled Cox models (baseline and time-dependent).
- **`joint_lme_survival_models.R`**: LME model selection per biomarker and the three joint Cox-LME models.
- **`repeated_measures_graft_patient_failure_presentation.pptx`**: the final presentation reporting all results, assumption checks, and dynamic prediction plots.
- **`data/d.csv`**: baseline (wide-format) patient data.
- **`data/dlong.csv`**: long-format repeated biomarker measurements.
- **`data/kidney_transplant.RData`**: R data file with the prepared objects used across both scripts.

## Tools

R · mice · survival · mitools · nlme · splines · JM · dplyr · tidyr · ggplot2

## Notes

This was a coursework exercise for a biostatistics module (MAM03), not a production project. It's shared here to demonstrate handling repeated-measures clinical data end-to-end: multiple imputation restricted to the right variables, comparing a standard time-varying-covariate approach against joint modeling, and producing individualized dynamic predictions rather than only population-level hazard ratios.
