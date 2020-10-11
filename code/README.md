# Code

## `glm-ftest.R`

Contains the function `glm.ftest.v2` to conduct an analysis of deviance and *F*-tests on the appropriate error term for generalized linear models. This function was used to produce Tables S1-S3 and S5 in the Supplementary Material.

## `plot-feasibility-domain.R`

Contains code for calculating normalized angles from critical boundaries, as well as plotting feasibility domains (as in Fig. 4 and S1 of the paper).

## `prep-time-series.R`

This code processes `data/arabidopsis_clean_df.csv` and `data/ExperimentPlantBiomass.csv` into `output/timeseries_df.csv`, which is then used for the Bayesian multivariate autoregressive models and structural stability analysis.

## `simulate-community-dynamics.R`

This code is used for the non-equilibrium simulations presented in the Supplementary Material.

## `temperature-structural-stability.R`

This code reproduces the underlying figure for Fig. S1, which was then annotated using (Gimp)[https://www.gimp.org/].

## `AOP2-LYER-Ptoid-persistence.R`

This code allows me to check whether the Baysian multivariate autoregressive models were able to reproduce the positive effect of *AOP2*$-$ on food-chain persistence. This was important in guiding my model selection (final model given by equation 2 in the Supplementary Material).
