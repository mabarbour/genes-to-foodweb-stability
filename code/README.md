# Code

## `glm-ftest.R`

Contains the function `glm.ftest.v2` to conduct an analysis of deviance and *F*-tests on the appropriate error term for generalized linear models. This function was used to produce Tables S1-S4 in the Supplementary Material, as well as the mulinomial ANOVA GLMs presented in the main text of the paper.

## `plot-feasibility-domain.R`

Contains code for calculating normalized angles from critical boundaries, as well as plotting feasibility domains (as in Fig. 4 and S1 of the paper).

## `prep-time-series.R`

This code processes `data/arabidopsis_clean_df.csv` and `data/ExperimentPlantBiomass.csv` into `output/timeseries_df.csv`, which is then used for the Bayesian multivariate autoregressive models and structural stability analysis.

## `simulate-community-dynamics.R`

This code is used for the non-equilibrium simulations presented in the Supplementary Material.

## `temperature-structural-stability.R`

This code reproduces the underlying figure for Fig. S1, which was then annotated using (Gimp)[https://www.gimp.org/].
