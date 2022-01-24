# Output

The folder contains processed data files. Descriptions are given below:

## `.rds` files

Contain saved versions of each model analyzed in `structural-stability.Rmd`. I saved each version because it considerably sped up generating `structural-stability.html`. Most of these files were too large to use version control on GitHub; therefore, I used the [piggyback package](https://docs.ropensci.org/piggyback/articles/intro.html) to attach it to v3.0 (see Assets in Release v3.0)

## `time-series-data.RData`

Organized time-series data from `InsectAbundanceSurvival.csv` for Bayesian multivariate autoregressive models and structural stability analysis.
