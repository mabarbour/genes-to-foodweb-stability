# general function for evaluating what percentage of posterior is aop2 > AOP2 for the LYER-Ptoid boundary
# this function was important for guiding my model selection
pp_aop2_LP_persist <- function(brm_model,
                               temp.cond = 1.5, # condition on imaginary intermediate temp between treatments
                               aop2.cond = 2){
  # get posterior predictions
  pp.aop2_model <- posterior_samples(brm_model, pars = "^b")

  # calculate structural stability with aop2
  aop2_stability <- apply(
    pp.aop2_model,
    MARGIN = 1,
    FUN = function(x) {
      tmp.mat = matrix(c(ifelse(is.na(x["b_log1pLYERt1_log1pLYER_t"]) == FALSE,
                                x["b_log1pLYERt1_log1pLYER_t"],
                                0) +
                           ifelse(is.na(x["b_log1pLYERt1_log1pLYER_t:aop2_genotypes"]) == FALSE,
                                  x["b_log1pLYERt1_log1pLYER_t:aop2_genotypes"]*aop2.cond,
                                  0) +
                           ifelse(is.na(x["b_log1pLYERt1_log1pLYER_t:temp"]) == FALSE,
                                  x["b_log1pLYERt1_log1pLYER_t:temp"]*temp.cond,
                                  0),
                         ifelse(is.na(x["b_log1pLYERt1_log1pPtoid_t"]) == FALSE,
                                x["b_log1pLYERt1_log1pPtoid_t"],
                                0) +
                           ifelse(is.na(x["b_log1pLYERt1_log1pPtoid_t:aop2_genotypes"]) == FALSE,
                                  x["b_log1pLYERt1_log1pPtoid_t:aop2_genotypes"]*aop2.cond,
                                  0) +
                           ifelse(is.na(x["b_log1pLYERt1_log1pPtoid_t:temp"]) == FALSE,
                                  x["b_log1pLYERt1_log1pPtoid_t:temp"]*temp.cond,
                                  0),
                         ifelse(is.na(x["b_log1pPtoidt1_log1pLYER_t"]) == FALSE,
                                x["b_log1pPtoidt1_log1pLYER_t"],
                                0) +
                           ifelse(is.na(x["b_log1pPtoidt1_log1pLYER_t:aop2_genotypes"]) == FALSE,
                                  x["b_log1pPtoidt1_log1pLYER_t:aop2_genotypes"]*aop2.cond,
                                  0) +
                           ifelse(is.na(x["b_log1pPtoidt1_log1pLYER_t:temp"]) == FALSE,
                                  x["b_log1pPtoidt1_log1pLYER_t:temp"]*temp.cond,
                                  0),
                         ifelse(is.na(x["b_log1pPtoidt1_log1pPtoid_t"]) == FALSE,
                                x["b_log1pPtoidt1_log1pPtoid_t"],
                                0) +
                           ifelse(is.na(x["b_log1pPtoidt1_log1pPtoid_t:aop2_genotypes"]) == FALSE,
                                  x["b_log1pPtoidt1_log1pPtoid_t:aop2_genotypes"]*aop2.cond,
                                  0) +
                           ifelse(is.na(x["b_log1pPtoidt1_log1pPtoid_t:temp"]) == FALSE,
                                  x["b_log1pPtoidt1_log1pPtoid_t:temp"]*temp.cond,
                                  0)),
                       ncol = 2, byrow = TRUE)
      tmp.r = matrix(c(x["b_log1pLYERt1_intercept"] +
                         ifelse(is.na(x["b_log1pLYERt1_aop2_genotypes"]) == FALSE,
                                x["b_log1pLYERt1_aop2_genotypes"]*aop2.cond,
                                0) +
                         ifelse(is.na(x["b_log1pLYERt1_temp"]) == FALSE,
                                x["b_log1pLYERt1_temp"]*temp.cond,
                                0),
                       x["b_log1pPtoidt1_intercept"] +
                         ifelse(is.na(x["b_log1pPtoidt1_aop2_genotypes"]) == FALSE,
                                x["b_log1pPtoidt1_aop2_genotypes"]*aop2.cond,
                                0) +
                         ifelse(is.na(x["b_log1pPtoidt1_temp"]) == FALSE,
                                x["b_log1pPtoidt1_temp"]*temp.cond,
                                0)),
                     ncol = 1)
      FeasibilityBoundaryLYER.Ptoid = as.numeric(BoundaryLYER.Ptoid(tmp.mat, tmp.r)["boundary"])
      c(FeasibilityBoundaryLYER.Ptoid = FeasibilityBoundaryLYER.Ptoid)
    })
  aop2_stability.df <- data.frame(
    allele = "aop2",
    posterior_sample = 1:nrow(pp.aop2_model),
    FeasibilityBoundaryLYER.Ptoid = aop2_stability
  )

  # calculate structural stability with AOP2
  AOP2_stability <- apply(
    pp.aop2_model,
    MARGIN = 1,
    FUN = function(x) {
      tmp.mat = matrix(c(ifelse(is.na(x["b_log1pLYERt1_log1pLYER_t"]) == FALSE,
                                x["b_log1pLYERt1_log1pLYER_t"],
                                0) +
                           ifelse(is.na(x["b_log1pLYERt1_log1pLYER_t:AOP2_genotypes"]) == FALSE,
                                  x["b_log1pLYERt1_log1pLYER_t:AOP2_genotypes"]*aop2.cond,
                                  0) +
                           ifelse(is.na(x["b_log1pLYERt1_log1pLYER_t:temp"]) == FALSE,
                                  x["b_log1pLYERt1_log1pLYER_t:temp"]*temp.cond,
                                  0),
                         ifelse(is.na(x["b_log1pLYERt1_log1pPtoid_t"]) == FALSE,
                                x["b_log1pLYERt1_log1pPtoid_t"],
                                0) +
                           ifelse(is.na(x["b_log1pLYERt1_log1pPtoid_t:AOP2_genotypes"]) == FALSE,
                                  x["b_log1pLYERt1_log1pPtoid_t:AOP2_genotypes"]*aop2.cond,
                                  0) +
                           ifelse(is.na(x["b_log1pLYERt1_log1pPtoid_t:temp"]) == FALSE,
                                  x["b_log1pLYERt1_log1pPtoid_t:temp"]*temp.cond,
                                  0),
                         ifelse(is.na(x["b_log1pPtoidt1_log1pLYER_t"]) == FALSE,
                                x["b_log1pPtoidt1_log1pLYER_t"],
                                0) +
                           ifelse(is.na(x["b_log1pPtoidt1_log1pLYER_t:AOP2_genotypes"]) == FALSE,
                                  x["b_log1pPtoidt1_log1pLYER_t:AOP2_genotypes"]*aop2.cond,
                                  0) +
                           ifelse(is.na(x["b_log1pPtoidt1_log1pLYER_t:temp"]) == FALSE,
                                  x["b_log1pPtoidt1_log1pLYER_t:temp"]*temp.cond,
                                  0),
                         ifelse(is.na(x["b_log1pPtoidt1_log1pPtoid_t"]) == FALSE,
                                x["b_log1pPtoidt1_log1pPtoid_t"],
                                0) +
                           ifelse(is.na(x["b_log1pPtoidt1_log1pPtoid_t:AOP2_genotypes"]) == FALSE,
                                  x["b_log1pPtoidt1_log1pPtoid_t:AOP2_genotypes"]*aop2.cond,
                                  0) +
                           ifelse(is.na(x["b_log1pPtoidt1_log1pPtoid_t:temp"]) == FALSE,
                                  x["b_log1pPtoidt1_log1pPtoid_t:temp"]*temp.cond,
                                  0)),
                       ncol = 2, byrow = TRUE)
      tmp.r = matrix(c(x["b_log1pLYERt1_intercept"] +
                         ifelse(is.na(x["b_log1pLYERt1_AOP2_genotypes"]) == FALSE,
                                x["b_log1pLYERt1_AOP2_genotypes"]*aop2.cond,
                                0) +
                         ifelse(is.na(x["b_log1pLYERt1_temp"]) == FALSE,
                                x["b_log1pLYERt1_temp"]*temp.cond,
                                0),
                       x["b_log1pPtoidt1_intercept"] +
                         ifelse(is.na(x["b_log1pPtoidt1_AOP2_genotypes"]) == FALSE,
                                x["b_log1pPtoidt1_AOP2_genotypes"]*aop2.cond,
                                0) +
                         ifelse(is.na(x["b_log1pPtoidt1_temp"]) == FALSE,
                                x["b_log1pPtoidt1_temp"]*temp.cond,
                                0)),
                     ncol = 1)
      FeasibilityBoundaryLYER.Ptoid = as.numeric(BoundaryLYER.Ptoid(tmp.mat, tmp.r)["boundary"])
      c(FeasibilityBoundaryLYER.Ptoid = FeasibilityBoundaryLYER.Ptoid)
    })
  AOP2_stability.df <- data.frame(
    allele = "AOP2",
    posterior_sample = 1:nrow(pp.aop2_model),
    FeasibilityBoundaryLYER.Ptoid = AOP2_stability
  )
  # combine data
  all.aop2_stability.df <- bind_rows(aop2_stability.df, AOP2_stability.df)

  # organize data
  aop2_LPbound <- all.aop2_stability.df %>%
    select(allele, posterior_sample, FeasibilityBoundaryLYER.Ptoid) %>%
    spread(allele, FeasibilityBoundaryLYER.Ptoid) %>%
    mutate(allele_effect = `aop2` - `AOP2`)

  # calculate percentage where aop2 > AOP2 effect
  aop2_LPbound_BayesP <- mean(aop2_LPbound$allele_effect > 0)

  # calculate aop2 > AOP2 effect size
  aop2_LPbound_effect <- mean(aop2_LPbound$allele_effect)

  return(c(aop2_LPbound_BayesP = aop2_LPbound_BayesP, aop2_LPbound_effect = aop2_LPbound_effect))
}

# function for evaluating mean aop2 matrices and mean LYER-Ptoid boundary
aop2_LP_persist <- function(brm_model,
                            temp.cond = 1.5, # condition on imaginary intermediate temp between treatments
                            aop2.cond = 2){
  # get mean coefficients
  coef.aop2 <- fixef(brm_model)[,"Estimate"]

  # interaction matrix for all levels of aop2
  aop2.mat <- matrix(c(ifelse(is.na(coef.aop2["log1pLYERt1_log1pLYER_t"]) == FALSE,
                              coef.aop2["log1pLYERt1_log1pLYER_t"],
                              0) +
                         ifelse(is.na(coef.aop2["log1pLYERt1_log1pLYER_t:aop2_genotypes"]) == FALSE,
                                coef.aop2["log1pLYERt1_log1pLYER_t:aop2_genotypes"]*aop2.cond,
                                0) +
                         ifelse(is.na(coef.aop2["log1pLYERt1_log1pLYER_t:temp"]) == FALSE,
                                coef.aop2["log1pLYERt1_log1pLYER_t:temp"]*temp.cond,
                                0),
                       ifelse(is.na(coef.aop2["log1pLYERt1_log1pPtoid_t"]) == FALSE,
                              coef.aop2["log1pLYERt1_log1pPtoid_t"],
                              0) +
                         ifelse(is.na(coef.aop2["log1pLYERt1_log1pPtoid_t:aop2_genotypes"]) == FALSE,
                                coef.aop2["log1pLYERt1_log1pPtoid_t:aop2_genotypes"]*aop2.cond,
                                0) +
                         ifelse(is.na(coef.aop2["log1pLYERt1_log1pPtoid_t:temp"]) == FALSE,
                                coef.aop2["log1pLYERt1_log1pPtoid_t:temp"]*temp.cond,
                                0),
                       ifelse(is.na(coef.aop2["log1pPtoidt1_log1pLYER_t"]) == FALSE,
                              coef.aop2["log1pPtoidt1_log1pLYER_t"],
                              0) +
                         ifelse(is.na(coef.aop2["log1pPtoidt1_log1pLYER_t:aop2_genotypes"]) == FALSE,
                                coef.aop2["log1pPtoidt1_log1pLYER_t:aop2_genotypes"]*aop2.cond,
                                0) +
                         ifelse(is.na(coef.aop2["log1pPtoidt1_log1pLYER_t:temp"]) == FALSE,
                                coef.aop2["log1pPtoidt1_log1pLYER_t:temp"]*temp.cond,
                                0),
                       ifelse(is.na(coef.aop2["log1pPtoidt1_log1pPtoid_t"]) == FALSE,
                              coef.aop2["log1pPtoidt1_log1pPtoid_t"],
                              0) +
                         ifelse(is.na(coef.aop2["log1pPtoidt1_log1pPtoid_t:aop2_genotypes"]) == FALSE,
                                coef.aop2["log1pPtoidt1_log1pPtoid_t:aop2_genotypes"]*aop2.cond,
                                0) +
                         ifelse(is.na(coef.aop2["log1pPtoidt1_log1pPtoid_t:temp"]) == FALSE,
                                coef.aop2["log1pPtoidt1_log1pPtoid_t:temp"]*temp.cond,
                                0)),
                     ncol = 2, byrow = TRUE)

  AOP2.mat <- matrix(c(ifelse(is.na(coef.aop2["log1pLYERt1_log1pLYER_t"]) == FALSE,
                              coef.aop2["log1pLYERt1_log1pLYER_t"],
                              0) +
                         ifelse(is.na(coef.aop2["log1pLYERt1_log1pLYER_t:AOP2_genotypes"]) == FALSE,
                                coef.aop2["log1pLYERt1_log1pLYER_t:AOP2_genotypes"]*aop2.cond,
                                0) +
                         ifelse(is.na(coef.aop2["log1pLYERt1_log1pLYER_t:temp"]) == FALSE,
                                coef.aop2["log1pLYERt1_log1pLYER_t:temp"]*temp.cond,
                                0),
                       ifelse(is.na(coef.aop2["log1pLYERt1_log1pPtoid_t"]) == FALSE,
                              coef.aop2["log1pLYERt1_log1pPtoid_t"],
                              0) +
                         ifelse(is.na(coef.aop2["log1pLYERt1_log1pPtoid_t:AOP2_genotypes"]) == FALSE,
                                coef.aop2["log1pLYERt1_log1pPtoid_t:AOP2_genotypes"]*aop2.cond,
                                0) +
                         ifelse(is.na(coef.aop2["log1pLYERt1_log1pPtoid_t:temp"]) == FALSE,
                                coef.aop2["log1pLYERt1_log1pPtoid_t:temp"]*temp.cond,
                                0),
                       ifelse(is.na(coef.aop2["log1pPtoidt1_log1pLYER_t"]) == FALSE,
                              coef.aop2["log1pPtoidt1_log1pLYER_t"],
                              0) +
                         ifelse(is.na(coef.aop2["log1pPtoidt1_log1pLYER_t:AOP2_genotypes"]) == FALSE,
                                coef.aop2["log1pPtoidt1_log1pLYER_t:AOP2_genotypes"]*aop2.cond,
                                0) +
                         ifelse(is.na(coef.aop2["log1pPtoidt1_log1pLYER_t:temp"]) == FALSE,
                                coef.aop2["log1pPtoidt1_log1pLYER_t:temp"]*temp.cond,
                                0),
                       ifelse(is.na(coef.aop2["log1pPtoidt1_log1pPtoid_t"]) == FALSE,
                              coef.aop2["log1pPtoidt1_log1pPtoid_t"],
                              0) +
                         ifelse(is.na(coef.aop2["log1pPtoidt1_log1pPtoid_t:AOP2_genotypes"]) == FALSE,
                                coef.aop2["log1pPtoidt1_log1pPtoid_t:AOP2_genotypes"]*aop2.cond,
                                0) +
                         ifelse(is.na(coef.aop2["log1pPtoidt1_log1pPtoid_t:temp"]) == FALSE,
                                coef.aop2["log1pPtoidt1_log1pPtoid_t:temp"]*temp.cond,
                                0)),
                     ncol = 2, byrow = TRUE)

  # growth rates with aop2 allele (both genotypes)
  aop2.IGR <- matrix(c(coef.aop2["log1pLYERt1_intercept"] +
                         ifelse(is.na(coef.aop2["log1pLYERt1_aop2_genotypes"]) == FALSE,
                                coef.aop2["log1pLYERt1_aop2_genotypes"]*aop2.cond,
                                0) +
                         ifelse(is.na(coef.aop2["log1pLYERt1_temp"]) == FALSE,
                                coef.aop2["log1pLYERt1_temp"]*temp.cond,
                                0),
                       coef.aop2["log1pPtoidt1_intercept"] +
                         ifelse(is.na(coef.aop2["log1pPtoidt1_aop2_genotypes"]) == FALSE,
                                coef.aop2["log1pPtoidt1_aop2_genotypes"]*aop2.cond,
                                0) +
                         ifelse(is.na(coef.aop2["log1pPtoidt1_temp"]) == FALSE,
                                coef.aop2["log1pPtoidt1_temp"]*temp.cond,
                                0)),
                     ncol = 1)

  # growth rates with AOP2 allele (both genotypes)
  AOP2.IGR <- matrix(c(coef.aop2["log1pLYERt1_intercept"] +
                         ifelse(is.na(coef.aop2["log1pLYERt1_AOP2_genotypes"]) == FALSE,
                                coef.aop2["log1pLYERt1_AOP2_genotypes"]*aop2.cond,
                                0) +
                         ifelse(is.na(coef.aop2["log1pLYERt1_temp"]) == FALSE,
                                coef.aop2["log1pLYERt1_temp"]*temp.cond,
                                0),
                       coef.aop2["log1pPtoidt1_intercept"] +
                         ifelse(is.na(coef.aop2["log1pPtoidt1_AOP2_genotypes"]) == FALSE,
                                coef.aop2["log1pPtoidt1_AOP2_genotypes"]*aop2.cond,
                                0) +
                         ifelse(is.na(coef.aop2["log1pPtoidt1_temp"]) == FALSE,
                                coef.aop2["log1pPtoidt1_temp"]*temp.cond,
                                0)),
                     ncol = 1)

  return(list(aop2.mat = aop2.mat, AOP2.mat = AOP2.mat, aop2.IGR = aop2.IGR, AOP2.IGR = AOP2.IGR,
              aop2_feasibility = -inv(aop2.mat) %*% aop2.IGR,
              AOP2_feasibility = -inv(AOP2.mat) %*% AOP2.IGR,
              aop2_LP_boundary = as.numeric(BoundaryLYER.Ptoid(aop2.mat, aop2.IGR)["boundary"]),
              AOP2_LP_boundary = as.numeric(BoundaryLYER.Ptoid(AOP2.mat, AOP2.IGR)["boundary"])))
}
