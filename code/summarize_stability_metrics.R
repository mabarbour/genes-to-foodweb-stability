
## Median effects ----

aop2_vs_AOP2_median_effects_unadj <- function(brms_model, n.geno = 1, temp.value = 0, logbiomass.value = 0){

  # focus on median estimates
  coef.aop2_vs_AOP2 <- fixef(brms_model)[,"Estimate"]

  # baseline interaction matrix
  base.mat <- matrix(c(coef.aop2_vs_AOP2["logBRBRt1_melogBRBR_tse_logBRBRt"],
                       coef.aop2_vs_AOP2["logBRBRt1_melogLYER_tse_logLYERt"],
                       coef.aop2_vs_AOP2["logBRBRt1_melogPtoid_tse_logPtoidt"],
                       coef.aop2_vs_AOP2["logLYERt1_melogBRBR_tse_logBRBRt"],
                       coef.aop2_vs_AOP2["logLYERt1_melogLYER_tse_logLYERt"],
                       coef.aop2_vs_AOP2["logLYERt1_melogPtoid_tse_logPtoidt"],
                       coef.aop2_vs_AOP2["logPtoidt1_melogBRBR_tse_logBRBRt"],
                       coef.aop2_vs_AOP2["logPtoidt1_melogLYER_tse_logLYERt"],
                       coef.aop2_vs_AOP2["logPtoidt1_melogPtoid_tse_logPtoidt"]),
                     ncol = 3, byrow = TRUE)
  base.mat[is.na(base.mat)] <- 0 # NA interactions become zeros

  # baseline intrinsic growth rates
  base.IGR <- matrix(c(coef.aop2_vs_AOP2["logBRBRt1_Intercept"],
                       coef.aop2_vs_AOP2["logLYERt1_Intercept"],
                       coef.aop2_vs_AOP2["logPtoidt1_Intercept"]),
                     ncol = 1)
  base.IGR[is.na(base.IGR)] <- 0 # NA r become zeros

  ## aop2 effects: I will add this to the base.mat to get effect of aop2 on matrix
  # interaction matrix
  delta.aop2.mat <- matrix(c(coef.aop2_vs_AOP2["logBRBRt1_melogBRBR_tse_logBRBRt:aop2_genotypes"],
                             coef.aop2_vs_AOP2["logBRBRt1_melogLYER_tse_logLYERt:aop2_genotypes"],
                             coef.aop2_vs_AOP2["logBRBRt1_melogPtoid_tse_logPtoidt:aop2_genotypes"],
                             coef.aop2_vs_AOP2["logLYERt1_melogBRBR_tse_logBRBRt:aop2_genotypes"],
                             coef.aop2_vs_AOP2["logLYERt1_melogLYER_tse_logLYERt:aop2_genotypes"],
                             coef.aop2_vs_AOP2["logLYERt1_melogPtoid_tse_logPtoidt:aop2_genotypes"],
                             coef.aop2_vs_AOP2["logPtoidt1_melogBRBR_tse_logBRBRt:aop2_genotypes"],
                             coef.aop2_vs_AOP2["logPtoidt1_melogLYER_tse_logLYERt:aop2_genotypes"],
                             coef.aop2_vs_AOP2["logPtoidt1_melogPtoid_tse_logPtoidt:aop2_genotypes"]),
                           ncol = 3, byrow = TRUE)
  delta.aop2.mat[is.na(delta.aop2.mat)] <- 0 # no aop2 effect if NA
  # intrinsic growth rates
  delta.aop2.IGR <- matrix(c(coef.aop2_vs_AOP2["logBRBRt1_aop2_genotypes"],
                             coef.aop2_vs_AOP2["logLYERt1_aop2_genotypes"],
                             coef.aop2_vs_AOP2["logPtoidt1_aop2_genotypes"]),
                           ncol = 1)
  delta.aop2.IGR[is.na(delta.aop2.IGR)] <- 0 # no aop2 effect if NA

  ## AOP2 effects
  # interaction matrix
  delta.AOP2.mat <- matrix(c(coef.aop2_vs_AOP2["logBRBRt1_melogBRBR_tse_logBRBRt:AOP2_genotypes"],
                             coef.aop2_vs_AOP2["logBRBRt1_melogLYER_tse_logLYERt:AOP2_genotypes"],
                             coef.aop2_vs_AOP2["logBRBRt1_melogPtoid_tse_logPtoidt:AOP2_genotypes"],
                             coef.aop2_vs_AOP2["logLYERt1_melogBRBR_tse_logBRBRt:AOP2_genotypes"],
                             coef.aop2_vs_AOP2["logLYERt1_melogLYER_tse_logLYERt:AOP2_genotypes"],
                             coef.aop2_vs_AOP2["logLYERt1_melogPtoid_tse_logPtoidt:AOP2_genotypes"],
                             coef.aop2_vs_AOP2["logPtoidt1_melogBRBR_tse_logBRBRt:AOP2_genotypes"],
                             coef.aop2_vs_AOP2["logPtoidt1_melogLYER_tse_logLYERt:AOP2_genotypes"],
                             coef.aop2_vs_AOP2["logPtoidt1_melogPtoid_tse_logPtoidt:AOP2_genotypes"]),
                           ncol = 3, byrow = TRUE)
  delta.AOP2.mat[is.na(delta.AOP2.mat)] <- 0 # no aop2 effect if NA
  # intrinsic growth rates
  delta.AOP2.IGR <- matrix(c(coef.aop2_vs_AOP2["logBRBRt1_AOP2_genotypes"],
                             coef.aop2_vs_AOP2["logLYERt1_AOP2_genotypes"],
                             coef.aop2_vs_AOP2["logPtoidt1_AOP2_genotypes"]),
                           ncol = 1)
  delta.AOP2.IGR[is.na(delta.AOP2.IGR)] <- 0 # no aop2 effect if NA

  # condition on temperature for intrinsic growth rates
  delta.temp.r = matrix(c(coef.aop2_vs_AOP2["logBRBRt1_temp"],
                          coef.aop2_vs_AOP2["logLYERt1_temp"],
                          coef.aop2_vs_AOP2["logPtoidt1_temp"]),
                        ncol = 1)
  delta.temp.r[is.na(delta.temp.r)] <- 0 # NA r become zeros

  # condition on plant biomass value for intrinsic growth rates
  delta.logbiomass.r = matrix(c(coef.aop2_vs_AOP2["logBRBRt1_logBiomass_g_t1"],
                                coef.aop2_vs_AOP2["logLYERt1_logBiomass_g_t1"],
                                coef.aop2_vs_AOP2["logPtoidt1_logBiomass_g_t1"]),
                              ncol = 1)
  delta.logbiomass.r[is.na(delta.logbiomass.r)] <- 0 # NA r become zeros

  # make aop2 and AOP2 matrices and IGR
  aop2.mat <- base.mat + delta.aop2.mat*n.geno - diag(3) # - diag(3) puts on 'continuous time' scale
  AOP2.mat <- base.mat + delta.AOP2.mat*n.geno - diag(3) # - diag(3) puts on 'continuous time' scale
  aop2.IGR <- base.IGR + delta.aop2.IGR*n.geno + delta.temp.r*temp.value + delta.logbiomass.r*logbiomass.value
  AOP2.IGR <- base.IGR + delta.AOP2.IGR*n.geno + delta.temp.r*temp.value + delta.logbiomass.r*logbiomass.value

  # 2 genotypes of each
  # aop2aop2.mat <- base.mat + delta.aop2.mat  + delta.aop2.mat - diag(3) # - diag(3) puts on 'continuous time' scale
  # AOP2AOP2.mat <- base.mat + delta.AOP2.mat + delta.AOP2.mat- diag(3) # - diag(3) puts on 'continuous time' scale
  # aop2aop2.IGR <- base.IGR + delta.aop2.IGR + delta.aop2.IGR
  # AOP2AOP2.IGR <- base.IGR + delta.AOP2.IGR + delta.AOP2.IGR

  # output
  return(list(base.mat = base.mat - diag(3),
              base.IGR = base.IGR,
              aop2.mat = aop2.mat,
              AOP2.mat = AOP2.mat,
              aop2.IGR = aop2.IGR,
              AOP2.IGR = AOP2.IGR))
  # aop2aop2.mat = aop2aop2.mat,
  # AOP2AOP2.mat = AOP2AOP2.mat,
  # aop2aop2.IGR = aop2aop2.IGR,
  # AOP2AOP2.IGR = AOP2AOP2.IGR))
}


aop2_vs_AOP2_median_effects_adj <- function(brms_model, n.geno = 1, temp.value = 0, logbiomass.value = 0){

  # focus on median estimates
  coef.aop2_vs_AOP2 <- fixef(brms_model)[,"Estimate"]

  # baseline interaction matrix
  base.mat <- matrix(c(coef.aop2_vs_AOP2["logBRBRt1_melogBRBR_tse_logBRBRt"],
                       coef.aop2_vs_AOP2["logBRBRt1_melogLYER_t_adjse_logLYERt"],
                       coef.aop2_vs_AOP2["logBRBRt1_melogPtoid_t_adjse_logPtoidt"],
                       coef.aop2_vs_AOP2["logLYERt1adj_melogBRBR_tse_logBRBRt"],
                       coef.aop2_vs_AOP2["logLYERt1adj_melogLYER_t_adjse_logLYERt"],
                       coef.aop2_vs_AOP2["logLYERt1adj_melogPtoid_t_adjse_logPtoidt"],
                       coef.aop2_vs_AOP2["logPtoidt1adj_melogBRBR_tse_logBRBRt"],
                       coef.aop2_vs_AOP2["logPtoidt1adj_melogLYER_t_adjse_logLYERt"],
                       coef.aop2_vs_AOP2["logPtoidt1adj_melogPtoid_t_adjse_logPtoidt"]),
                     ncol = 3, byrow = TRUE)
  base.mat[is.na(base.mat)] <- 0 # NA interactions become zeros

  # baseline intrinsic growth rates
  base.IGR <- matrix(c(coef.aop2_vs_AOP2["logBRBRt1_Intercept"],
                       coef.aop2_vs_AOP2["logLYERt1adj_Intercept"],
                       coef.aop2_vs_AOP2["logPtoidt1adj_Intercept"]),
                     ncol = 1)
  base.IGR[is.na(base.IGR)] <- 0 # NA r become zeros

  ## aop2 effects: I will add this to the base.mat to get effect of aop2 on matrix
  # interaction matrix
  delta.aop2.mat <- matrix(c(coef.aop2_vs_AOP2["logBRBRt1_melogBRBR_tse_logBRBRt:aop2_genotypes"],
                             coef.aop2_vs_AOP2["logBRBRt1_melogLYER_t_adjse_logLYERt:aop2_genotypes"],
                             coef.aop2_vs_AOP2["logBRBRt1_melogPtoid_t_adjse_logPtoidt:aop2_genotypes"],
                             coef.aop2_vs_AOP2["logLYERt1adj_melogBRBR_tse_logBRBRt:aop2_genotypes"],
                             coef.aop2_vs_AOP2["logLYERt1adj_melogLYER_t_adjse_logLYERt:aop2_genotypes"],
                             coef.aop2_vs_AOP2["logLYERt1adj_melogPtoid_t_adjse_logPtoidt:aop2_genotypes"],
                             coef.aop2_vs_AOP2["logPtoidt1adj_melogBRBR_tse_logBRBRt:aop2_genotypes"],
                             coef.aop2_vs_AOP2["logPtoidt1adj_melogLYER_t_adjse_logLYERt:aop2_genotypes"],
                             coef.aop2_vs_AOP2["logPtoidt1adj_melogPtoid_t_adjse_logPtoidt:aop2_genotypes"]),
                           ncol = 3, byrow = TRUE)
  delta.aop2.mat[is.na(delta.aop2.mat)] <- 0 # no aop2 effect if NA
  # intrinsic growth rates
  delta.aop2.IGR <- matrix(c(coef.aop2_vs_AOP2["logBRBRt1_aop2_genotypes"],
                             coef.aop2_vs_AOP2["logLYERt1adj_aop2_genotypes"],
                             coef.aop2_vs_AOP2["logPtoidt1adj_aop2_genotypes"]),
                           ncol = 1)
  delta.aop2.IGR[is.na(delta.aop2.IGR)] <- 0 # no aop2 effect if NA

  ## AOP2 effects
  # interaction matrix
  delta.AOP2.mat <- matrix(c(coef.aop2_vs_AOP2["logBRBRt1_melogBRBR_tse_logBRBRt:AOP2_genotypes"],
                             coef.aop2_vs_AOP2["logBRBRt1_melogLYER_t_adjse_logLYERt:AOP2_genotypes"],
                             coef.aop2_vs_AOP2["logBRBRt1_melogPtoid_t_adjse_logPtoidt:AOP2_genotypes"],
                             coef.aop2_vs_AOP2["logLYERt1adj_melogBRBR_tse_logBRBRt:AOP2_genotypes"],
                             coef.aop2_vs_AOP2["logLYERt1adj_melogLYER_t_adjse_logLYERt:AOP2_genotypes"],
                             coef.aop2_vs_AOP2["logLYERt1adj_melogPtoid_t_adjse_logPtoidt:AOP2_genotypes"],
                             coef.aop2_vs_AOP2["logPtoidt1adj_melogBRBR_tse_logBRBRt:AOP2_genotypes"],
                             coef.aop2_vs_AOP2["logPtoidt1adj_melogLYER_t_adjse_logLYERt:AOP2_genotypes"],
                             coef.aop2_vs_AOP2["logPtoidt1adj_melogPtoid_t_adjse_logPtoidt:AOP2_genotypes"]),
                           ncol = 3, byrow = TRUE)
  delta.AOP2.mat[is.na(delta.AOP2.mat)] <- 0 # no aop2 effect if NA
  # intrinsic growth rates
  delta.AOP2.IGR <- matrix(c(coef.aop2_vs_AOP2["logBRBRt1_AOP2_genotypes"],
                             coef.aop2_vs_AOP2["logLYERt1adj_AOP2_genotypes"],
                             coef.aop2_vs_AOP2["logPtoidt1adj_AOP2_genotypes"]),
                           ncol = 1)
  delta.AOP2.IGR[is.na(delta.AOP2.IGR)] <- 0 # no aop2 effect if NA

  # condition on temperature for intrinsic growth rates
  delta.temp.r = matrix(c(coef.aop2_vs_AOP2["logBRBRt1_temp"],
                          coef.aop2_vs_AOP2["logLYERt1_temp"],
                          coef.aop2_vs_AOP2["logPtoidt1_temp"]),
                        ncol = 1)
  delta.temp.r[is.na(delta.temp.r)] <- 0 # NA r become zeros

  # condition on plant biomass value for intrinsic growth rates
  delta.logbiomass.r = matrix(c(coef.aop2_vs_AOP2["logBRBRt1_logBiomass_g_t1"],
                                coef.aop2_vs_AOP2["logLYERt1_logBiomass_g_t1"],
                                coef.aop2_vs_AOP2["logPtoidt1_logBiomass_g_t1"]),
                              ncol = 1)
  delta.logbiomass.r[is.na(delta.logbiomass.r)] <- 0 # NA r become zeros

  # make aop2 and AOP2 matrices and IGR
  aop2.mat <- base.mat + delta.aop2.mat*n.geno - diag(3) # - diag(3) puts on 'continuous time' scale
  AOP2.mat <- base.mat + delta.AOP2.mat*n.geno - diag(3) # - diag(3) puts on 'continuous time' scale
  aop2.IGR <- base.IGR + delta.aop2.IGR*n.geno + delta.temp.r*temp.value + delta.logbiomass.r*logbiomass.value
  AOP2.IGR <- base.IGR + delta.AOP2.IGR*n.geno + delta.temp.r*temp.value + delta.logbiomass.r*logbiomass.value

  # 2 genotypes of each
  # aop2aop2.mat <- base.mat + delta.aop2.mat  + delta.aop2.mat - diag(3) # - diag(3) puts on 'continuous time' scale
  # AOP2AOP2.mat <- base.mat + delta.AOP2.mat + delta.AOP2.mat- diag(3) # - diag(3) puts on 'continuous time' scale
  # aop2aop2.IGR <- base.IGR + delta.aop2.IGR + delta.aop2.IGR
  # AOP2AOP2.IGR <- base.IGR + delta.AOP2.IGR + delta.AOP2.IGR

  # output
  return(list(base.mat = base.mat - diag(3),
              base.IGR = base.IGR,
              aop2.mat = aop2.mat,
              AOP2.mat = AOP2.mat,
              aop2.IGR = aop2.IGR,
              AOP2.IGR = AOP2.IGR))
              # aop2aop2.mat = aop2aop2.mat,
              # AOP2AOP2.mat = AOP2AOP2.mat,
              # aop2aop2.IGR = aop2aop2.IGR,
              # AOP2AOP2.IGR = AOP2AOP2.IGR))
}

## Full posterior ----

aop2_vs_AOP2_posterior_samples_unadj <- function(brms_model, n.geno = 1, temp.value = 0, logbiomass.value = 0){

  # get posterior predictions
  pp.aop2_vs_AOP2_model <- posterior_samples(brms_model, pars = "^b")

  # get posterior samples of structural stability
  base_stability <- apply(
    pp.aop2_vs_AOP2_model,
    MARGIN = 1,
    FUN = function(x) {
      base.mat = matrix(c(x["bsp_logBRBRt1_melogBRBR_tse_logBRBRt"],
                          x["bsp_logBRBRt1_melogLYER_tse_logLYERt"],
                          x["bsp_logBRBRt1_melogPtoid_tse_logPtoidt"],
                          x["bsp_logLYERt1_melogBRBR_tse_logBRBRt"],
                          x["bsp_logLYERt1_melogLYER_tse_logLYERt"],
                          x["bsp_logLYERt1_melogPtoid_tse_logPtoidt"],
                          x["bsp_logPtoidt1_melogBRBR_tse_logBRBRt"],
                          x["bsp_logPtoidt1_melogLYER_tse_logLYERt"],
                          x["bsp_logPtoidt1_melogPtoid_tse_logPtoidt"]),
                        ncol = 3, byrow = TRUE)
      base.mat[is.na(base.mat)] <- 0 # NA interactions become zeros

      base.r = matrix(c(x["b_logBRBRt1_Intercept"],
                        x["b_logLYERt1_Intercept"],
                        x["b_logPtoidt1_Intercept"]),
                      ncol = 1)
      base.r[is.na(base.r)] <- 0 # NA r become zeros

      # condition on temperature for intrinsic growth rates
      delta.temp.r = matrix(c(x["b_logBRBRt1_temp"],
                              x["b_logLYERt1_temp"],
                              x["b_logPtoidt1_temp"]),
                            ncol = 1)
      delta.temp.r[is.na(delta.temp.r)] <- 0 # NA r become zeros

      # condition on plant biomass value for intrinsic growth rates
      delta.logbiomass.r = matrix(c(x["b_logBRBRt1_logBiomass_g_t1"],
                                    x["b_logLYERt1_logBiomass_g_t1"],
                                    x["b_logPtoidt1_logBiomass_g_t1"]),
                                  ncol = 1)
      delta.logbiomass.r[is.na(delta.logbiomass.r)] <- 0 # NA r become zeros

      # add to get aop2 allele effect
      tmp.mat <- base.mat - diag(3) # put on 'continuous time scale'
      tmp.r <- base.r + delta.temp.r*temp.value + delta.logbiomass.r*logbiomass.value

      Feasibility = as.numeric(FeasibilityBoundary(tmp.mat, tmp.r)["feasibility"])
      FeasibilityBoundary23 = as.numeric(FeasibilityBoundary(tmp.mat, tmp.r)["alpha.A23"])
      Extinctions = paste0(which(-1*inv(tmp.mat) %*% tmp.r < 0), collapse = ",")
      Resilience = max(Re(eigen((tmp.mat + diag(3)))$values)) # on discrete time scale
      ResilienceLP = max(Re(eigen((tmp.mat[2:3,2:3] + diag(2)))$values)) # on discrete time scale
      FeasibilityBoundaryLYER.Ptoid = as.numeric(BoundaryLYER.Ptoid(tmp.mat[2:3,2:3], tmp.r[2:3])["boundary"])
      FeasibilityLYER.Ptoid = as.numeric(BoundaryLYER.Ptoid(tmp.mat[2:3,2:3], tmp.r[2:3])["feasibility"])
      r_LYER = tmp.r[2]
      r_Ptoid = tmp.r[3]
      list(tmp.mat = tmp.mat,
           tmp.r = tmp.r,
           Feasibility = Feasibility,
           FeasibilityBoundary23 = FeasibilityBoundary23,
           Extinctions = Extinctions,
           Resilience = Resilience,
           ResilienceLP = ResilienceLP,
           FeasibilityBoundaryLYER.Ptoid = FeasibilityBoundaryLYER.Ptoid,
           FeasibilityLYER.Ptoid = FeasibilityLYER.Ptoid,
           r_LYER = r_LYER,
           r_Ptoid = r_Ptoid)
    })
  base_stability.df <- data.frame(
    aop2_vs_AOP2 = 0,
    posterior_sample = 1:nrow(pp.aop2_vs_AOP2_model),
    Feasibility = unlist(lapply(base_stability, FUN = function(x) x$Feasibility)),
    FeasibilityBoundary23 = unlist(lapply(base_stability, FUN = function(x) x$FeasibilityBoundary23)),
    Extinctions = unlist(lapply(base_stability, FUN = function(x) x$Extinctions)),
    Resilience = unlist(lapply(base_stability, FUN = function(x) x$Resilience)),
    ResilienceLP = unlist(lapply(base_stability, FUN = function(x) x$ResilienceLP)),
    FeasibilityBoundaryLYER.Ptoid = unlist(lapply(base_stability, FUN = function(x) x$FeasibilityBoundaryLYER.Ptoid)),
    FeasibilityLYER.Ptoid = unlist(lapply(base_stability, FUN = function(x) x$FeasibilityLYER.Ptoid)),
    r_LYER = unlist(lapply(base_stability, FUN = function(x) x$r_LYER)),
    r_Ptoid = unlist(lapply(base_stability, FUN = function(x) x$r_Ptoid))
  )

  aop2_stability <- apply(
    pp.aop2_vs_AOP2_model,
    MARGIN = 1,
    FUN = function(x) {
      base.mat = matrix(c(x["bsp_logBRBRt1_melogBRBR_tse_logBRBRt"],
                          x["bsp_logBRBRt1_melogLYER_tse_logLYERt"],
                          x["bsp_logBRBRt1_melogPtoid_tse_logPtoidt"],
                          x["bsp_logLYERt1_melogBRBR_tse_logBRBRt"],
                          x["bsp_logLYERt1_melogLYER_tse_logLYERt"],
                          x["bsp_logLYERt1_melogPtoid_tse_logPtoidt"],
                          x["bsp_logPtoidt1_melogBRBR_tse_logBRBRt"],
                          x["bsp_logPtoidt1_melogLYER_tse_logLYERt"],
                          x["bsp_logPtoidt1_melogPtoid_tse_logPtoidt"]),
                        ncol = 3, byrow = TRUE)
      base.mat[is.na(base.mat)] <- 0 # NA interactions become zeros

      base.r = matrix(c(x["b_logBRBRt1_Intercept"],
                        x["b_logLYERt1_Intercept"],
                        x["b_logPtoidt1_Intercept"]),
                      ncol = 1)
      base.r[is.na(base.r)] <- 0 # NA r become zeros

      # aop2 effect
      delta.aop2.mat <- matrix(c(x["bsp_logBRBRt1_melogBRBR_tse_logBRBRt:aop2_genotypes"],
                                 x["bsp_logBRBRt1_melogLYER_tse_logLYERt:aop2_genotypes"],
                                 x["bsp_logBRBRt1_melogPtoid_tse_logPtoidt:aop2_genotypes"],
                                 x["bsp_logLYERt1_melogBRBR_tse_logBRBRt:aop2_genotypes"],
                                 x["bsp_logLYERt1_melogLYER_tse_logLYERt:aop2_genotypes"],
                                 x["bsp_logLYERt1_melogPtoid_tse_logPtoidt:aop2_genotypes"],
                                 x["bsp_logPtoidt1_melogBRBR_tse_logBRBRt:aop2_genotypes"],
                                 x["bsp_logPtoidt1_melogLYER_tse_logLYERt:aop2_genotypes"],
                                 x["bsp_logPtoidt1_melogPtoid_tse_logPtoidt:aop2_genotypes"]),
                               ncol = 3, byrow = TRUE)
      delta.aop2.mat[is.na(delta.aop2.mat)] <- 0 # no aop2 effect if NA

      delta.aop2.r <- matrix(c(x["b_logBRBRt1_aop2_genotypes"],
                               x["b_logLYERt1_aop2_genotypes"],
                               x["b_logPtoidt1_aop2_genotypes"]),
                             ncol = 1)
      delta.aop2.r[is.na(delta.aop2.r)] <- 0 # no aop2 effect if NA

      # condition on temperature for intrinsic growth rates
      delta.temp.r = matrix(c(x["b_logBRBRt1_temp"],
                              x["b_logLYERt1_temp"],
                              x["b_logPtoidt1_temp"]),
                      ncol = 1)
      delta.temp.r[is.na(delta.temp.r)] <- 0 # NA r become zeros

      # condition on plant biomass value for intrinsic growth rates
      delta.logbiomass.r = matrix(c(x["b_logBRBRt1_logBiomass_g_t1"],
                                    x["b_logLYERt1_logBiomass_g_t1"],
                                    x["b_logPtoidt1_logBiomass_g_t1"]),
                            ncol = 1)
      delta.logbiomass.r[is.na(delta.logbiomass.r)] <- 0 # NA r become zeros

      # add to get aop2 allele effect
      tmp.mat <- base.mat + delta.aop2.mat*n.geno - diag(3) # put on 'continuous time scale'
      tmp.r <- base.r + delta.aop2.r*n.geno + delta.temp.r*temp.value + delta.logbiomass.r*logbiomass.value

      Feasibility = as.numeric(FeasibilityBoundary(tmp.mat, tmp.r)["feasibility"])
      FeasibilityBoundary23 = as.numeric(FeasibilityBoundary(tmp.mat, tmp.r)["alpha.A23"])
      Extinctions = paste0(which(-1*inv(tmp.mat) %*% tmp.r < 0), collapse = ",")
      Resilience = max(Re(eigen((tmp.mat + diag(3)))$values)) # on discrete time scale
      ResilienceLP = max(Re(eigen((tmp.mat[2:3,2:3] + diag(2)))$values)) # on discrete time scale
      FeasibilityBoundaryLYER.Ptoid = as.numeric(BoundaryLYER.Ptoid(tmp.mat[2:3,2:3], tmp.r[2:3])["boundary"])
      FeasibilityLYER.Ptoid = as.numeric(BoundaryLYER.Ptoid(tmp.mat[2:3,2:3], tmp.r[2:3])["feasibility"])
      r_LYER = tmp.r[2]
      r_Ptoid = tmp.r[3]
      list(tmp.mat = tmp.mat,
           tmp.r = tmp.r,
           Feasibility = Feasibility,
           FeasibilityBoundary23 = FeasibilityBoundary23,
           Extinctions = Extinctions,
           Resilience = Resilience,
           ResilienceLP = ResilienceLP,
           FeasibilityBoundaryLYER.Ptoid = FeasibilityBoundaryLYER.Ptoid,
           FeasibilityLYER.Ptoid = FeasibilityLYER.Ptoid,
           r_LYER = r_LYER,
           r_Ptoid = r_Ptoid)
    })
  aop2_stability.df <- data.frame(
    aop2_vs_AOP2 = 1,
    posterior_sample = 1:nrow(pp.aop2_vs_AOP2_model),
    Feasibility = unlist(lapply(aop2_stability, FUN = function(x) x$Feasibility)),
    FeasibilityBoundary23 = unlist(lapply(aop2_stability, FUN = function(x) x$FeasibilityBoundary23)),
    Extinctions = unlist(lapply(aop2_stability, FUN = function(x) x$Extinctions)),
    Resilience = unlist(lapply(aop2_stability, FUN = function(x) x$Resilience)),
    ResilienceLP = unlist(lapply(aop2_stability, FUN = function(x) x$ResilienceLP)),
    FeasibilityBoundaryLYER.Ptoid = unlist(lapply(aop2_stability, FUN = function(x) x$FeasibilityBoundaryLYER.Ptoid)),
    FeasibilityLYER.Ptoid = unlist(lapply(aop2_stability, FUN = function(x) x$FeasibilityLYER.Ptoid)),
    r_LYER = unlist(lapply(aop2_stability, FUN = function(x) x$r_LYER)),
    r_Ptoid = unlist(lapply(aop2_stability, FUN = function(x) x$r_Ptoid))
  )

  AOP2_stability <- apply(
    pp.aop2_vs_AOP2_model,
    MARGIN = 1,
    FUN = function(x) {
      base.mat = matrix(c(x["bsp_logBRBRt1_melogBRBR_tse_logBRBRt"],
                          x["bsp_logBRBRt1_melogLYER_tse_logLYERt"],
                          x["bsp_logBRBRt1_melogPtoid_tse_logPtoidt"],
                          x["bsp_logLYERt1_melogBRBR_tse_logBRBRt"],
                          x["bsp_logLYERt1_melogLYER_tse_logLYERt"],
                          x["bsp_logLYERt1_melogPtoid_tse_logPtoidt"],
                          x["bsp_logPtoidt1_melogBRBR_tse_logBRBRt"],
                          x["bsp_logPtoidt1_melogLYER_tse_logLYERt"],
                          x["bsp_logPtoidt1_melogPtoid_tse_logPtoidt"]),
                        ncol = 3, byrow = TRUE)
      base.mat[is.na(base.mat)] <- 0 # NA interactions become zeros

      base.r = matrix(c(x["b_logBRBRt1_Intercept"],
                        x["b_logLYERt1_Intercept"],
                        x["b_logPtoidt1_Intercept"]),
                      ncol = 1)
      base.r[is.na(base.r)] <- 0 # NA r become zeros

      # AOP2 effect
      delta.AOP2.mat <- matrix(c(x["bsp_logBRBRt1_melogBRBR_tse_logBRBRt:AOP2_genotypes"],
                                 x["bsp_logBRBRt1_melogLYER_tse_logLYERt:AOP2_genotypes"],
                                 x["bsp_logBRBRt1_melogPtoid_tse_logPtoidt:AOP2_genotypes"],
                                 x["bsp_logLYERt1_melogBRBR_tse_logBRBRt:AOP2_genotypes"],
                                 x["bsp_logLYERt1_melogLYER_tse_logLYERt:AOP2_genotypes"],
                                 x["bsp_logLYERt1_melogPtoid_tse_logPtoidt:AOP2_genotypes"],
                                 x["bsp_logPtoidt1_melogBRBR_tse_logBRBRt:AOP2_genotypes"],
                                 x["bsp_logPtoidt1_melogLYER_tse_logLYERt:AOP2_genotypes"],
                                 x["bsp_logPtoidt1_melogPtoid_tse_logPtoidt:AOP2_genotypes"]),
                               ncol = 3, byrow = TRUE)
      delta.AOP2.mat[is.na(delta.AOP2.mat)] <- 0 # no aop2 effect if NA

      delta.AOP2.r <- matrix(c(x["b_logBRBRt1_AOP2_genotypes"],
                               x["b_logLYERt1_AOP2_genotypes"],
                               x["b_logPtoidt1_AOP2_genotypes"]),
                             ncol = 1)
      delta.AOP2.r[is.na(delta.AOP2.r)] <- 0 # no aop2 effect if NA

      # condition on temperature for intrinsic growth rates
      delta.temp.r = matrix(c(x["b_logBRBRt1_temp"],
                              x["b_logLYERt1_temp"],
                              x["b_logPtoidt1_temp"]),
                            ncol = 1)
      delta.temp.r[is.na(delta.temp.r)] <- 0 # NA r become zeros

      # condition on plant biomass value for intrinsic growth rates
      delta.logbiomass.r = matrix(c(x["b_logBRBRt1_logBiomass_g_t1"],
                                    x["b_logLYERt1_logBiomass_g_t1"],
                                    x["b_logPtoidt1_logBiomass_g_t1"]),
                                  ncol = 1)
      delta.logbiomass.r[is.na(delta.logbiomass.r)] <- 0 # NA r become zeros

      # add to get AOP2 allele
      tmp.mat <- base.mat + delta.AOP2.mat*n.geno - diag(3) # put on 'continuous time scale'
      tmp.r <- base.r + delta.AOP2.r*n.geno + delta.temp.r*temp.value + delta.logbiomass.r*logbiomass.value

      Feasibility = as.numeric(FeasibilityBoundary(tmp.mat, tmp.r)["feasibility"])
      FeasibilityBoundary23 = as.numeric(FeasibilityBoundary(tmp.mat, tmp.r)["alpha.A23"])
      Extinctions = paste0(which(-1*inv(tmp.mat) %*% tmp.r < 0), collapse = ",")
      Resilience = max(Re(eigen((tmp.mat + diag(3)))$values)) # on discrete time scale
      ResilienceLP = max(Re(eigen((tmp.mat[2:3,2:3] + diag(2)))$values)) # on discrete time scale
      FeasibilityBoundaryLYER.Ptoid = as.numeric(BoundaryLYER.Ptoid(tmp.mat[2:3,2:3], tmp.r[2:3])["boundary"])
      FeasibilityLYER.Ptoid = as.numeric(BoundaryLYER.Ptoid(tmp.mat[2:3,2:3], tmp.r[2:3])["feasibility"])
      r_LYER = tmp.r[2]
      r_Ptoid = tmp.r[3]
      list(tmp.mat = tmp.mat,
           tmp.r = tmp.r,
           Feasibility = Feasibility,
           FeasibilityBoundary23 = FeasibilityBoundary23,
           Extinctions = Extinctions,
           Resilience = Resilience,
           ResilienceLP = ResilienceLP,
           FeasibilityBoundaryLYER.Ptoid = FeasibilityBoundaryLYER.Ptoid,
           FeasibilityLYER.Ptoid = FeasibilityLYER.Ptoid,
           r_LYER = r_LYER,
           r_Ptoid = r_Ptoid)
    })
  AOP2_stability.df <- data.frame(
    aop2_vs_AOP2 = -1,
    posterior_sample = 1:nrow(pp.aop2_vs_AOP2_model),
    Feasibility = unlist(lapply(AOP2_stability, FUN = function(x) x$Feasibility)),
    FeasibilityBoundary23 = unlist(lapply(AOP2_stability, FUN = function(x) x$FeasibilityBoundary23)),
    Extinctions = unlist(lapply(AOP2_stability, FUN = function(x) x$Extinctions)),
    Resilience = unlist(lapply(AOP2_stability, FUN = function(x) x$Resilience)),
    ResilienceLP = unlist(lapply(AOP2_stability, FUN = function(x) x$ResilienceLP)),
    FeasibilityBoundaryLYER.Ptoid = unlist(lapply(AOP2_stability, FUN = function(x) x$FeasibilityBoundaryLYER.Ptoid)),
    FeasibilityLYER.Ptoid = unlist(lapply(AOP2_stability, FUN = function(x) x$FeasibilityLYER.Ptoid)),
    r_LYER = unlist(lapply(AOP2_stability, FUN = function(x) x$r_LYER)),
    r_Ptoid = unlist(lapply(AOP2_stability, FUN = function(x) x$r_Ptoid))
  )

  # combine data
  all.aop2_vs_AOP2_stability.df <- bind_rows(aop2_stability.df, AOP2_stability.df)

  ## Aphid and parasitoid intrinsic growth rates ----
  r_LYER_df <- all.aop2_vs_AOP2_stability.df %>%
    select(aop2_vs_AOP2, posterior_sample, r_LYER) %>%
    spread(key = aop2_vs_AOP2, value = r_LYER) %>%
    mutate(allele_effect = `1` - `-1`)
  aop2_r_LYER_effect <- mean(r_LYER_df$allele_effect)
  aop2_r_LYER_95CI <- quantile(r_LYER_df$allele_effect, probs = c(0.025, 0.975))

  r_Ptoid_df <- all.aop2_vs_AOP2_stability.df %>%
    select(aop2_vs_AOP2, posterior_sample, r_Ptoid) %>%
    spread(key = aop2_vs_AOP2, value = r_Ptoid) %>%
    mutate(allele_effect = `1` - `-1`)
  aop2_r_Ptoid_effect <- mean(r_Ptoid_df$allele_effect)
  aop2_r_Ptoid_95CI <- quantile(r_Ptoid_df$allele_effect, probs = c(0.025, 0.975))

  ## Full community structural stability ----
  # organize data
  FeasibilityBoundary23_df <- all.aop2_vs_AOP2_stability.df %>%
    select(aop2_vs_AOP2, posterior_sample, FeasibilityBoundary23) %>%
    spread(key = aop2_vs_AOP2, value = FeasibilityBoundary23) %>%
    mutate(allele_effect = `1` - `-1`)

  # calculate percentage where aop2 > AOP2 effect
  aop2_SS_full_BayesP <- mean(FeasibilityBoundary23_df$allele_effect > 0)

  # calculate aop2 > AOP2 effect size
  aop2_SS_full_effect <- mean(FeasibilityBoundary23_df$allele_effect)

  ## LYER-Ptoid Structural stability ----
  # organize data
  aop2_LPbound_df <- all.aop2_vs_AOP2_stability.df %>%
    select(aop2_vs_AOP2, posterior_sample, FeasibilityBoundaryLYER.Ptoid) %>%
    spread(aop2_vs_AOP2, FeasibilityBoundaryLYER.Ptoid) %>%
    mutate(allele_effect = `1` - `-1`)

  # calculate percentage where aop2 > AOP2 effect
  aop2_SS_LP_BayesP <- mean(aop2_LPbound_df$allele_effect > 0)

  # calculate aop2 > AOP2 effect size
  aop2_SS_LP_effect <- mean(aop2_LPbound_df$allele_effect)
  aop2_SS_LP_95CI <- quantile(aop2_LPbound_df$allele_effect, probs = c(0.025, 0.975))

  ## Full community resilience ----
  # organize data
  aop2_Resilience_df <- all.aop2_vs_AOP2_stability.df %>%
    select(aop2_vs_AOP2, posterior_sample, Resilience) %>%
    spread(aop2_vs_AOP2, Resilience) %>%
    mutate(allele_effect = abs(`-1`) - abs(`1`))
  # on discrete time, resilience varies from -1 to 1, with values closer to zero being more stable.\
  # therefore, I take the absolute value first. Then I substract AOP2 (-1) from aop2 (1). Whenever this difference
  # is greater than zero, it means that aop2 has a stabilizing effect.

  # calculate percentage where aop2 > AOP2 effect
  aop2_Resil_full_BayesP <- mean(aop2_Resilience_df$allele_effect > 0)

  # calculate aop2 > AOP2 effect size
  aop2_Resil_full_effect <- mean(aop2_Resilience_df$allele_effect)

  ## LYER-Ptoid resilience ----
  # organize data
  aop2_ResilienceLP_df <- all.aop2_vs_AOP2_stability.df %>%
    select(aop2_vs_AOP2, posterior_sample, ResilienceLP) %>%
    spread(aop2_vs_AOP2, ResilienceLP) %>%
    mutate(allele_effect = abs(`-1`) - abs(`1`))
  # see note for Full community resilience calculation

  # calculate percentage where aop2 > AOP2 effect
  aop2_Resil_LP_BayesP <- mean(aop2_ResilienceLP_df$allele_effect > 0)

  # calculate aop2 > AOP2 effect size
  aop2_Resil_LP_effect <- mean(aop2_ResilienceLP_df$allele_effect)

  ## Baseline comparisons
  # aop2
  all.base_vs_aop2_stability.df <- bind_rows(aop2_stability.df, base_stability.df)

  # Full community structural stability
  # organize data
  baseaop2_FeasibilityBoundary23_df <- all.base_vs_aop2_stability.df %>%
    select(aop2_vs_AOP2, posterior_sample, FeasibilityBoundary23) %>%
    spread(key = aop2_vs_AOP2, value = FeasibilityBoundary23) %>%
    mutate(allele_effect = `1` - `0`)

  # calculate percentage where aop2 > AOP2 effect
  baseaop2_SS_full_BayesP <- mean(baseaop2_FeasibilityBoundary23_df$allele_effect > 0)

  # calculate aop2 > AOP2 effect size
  baseaop2_SS_full_effect <- mean(baseaop2_FeasibilityBoundary23_df$allele_effect)

  # LYER-Ptoid Structural stability
  # organize data
  baseaop2_LPbound_df <- all.base_vs_aop2_stability.df %>%
    select(aop2_vs_AOP2, posterior_sample, FeasibilityBoundaryLYER.Ptoid) %>%
    spread(aop2_vs_AOP2, FeasibilityBoundaryLYER.Ptoid) %>%
    mutate(allele_effect = `1` - `0`)

  # calculate percentage where aop2 > AOP2 effect
  baseaop2_SS_LP_BayesP <- mean(baseaop2_LPbound_df$allele_effect > 0)

  # calculate aop2 > AOP2 effect size
  baseaop2_SS_LP_effect <- mean(baseaop2_LPbound_df$allele_effect)

  # Full community resilience
  # organize data
  baseaop2_Resilience_df <- all.base_vs_aop2_stability.df %>%
    select(aop2_vs_AOP2, posterior_sample, Resilience) %>%
    spread(aop2_vs_AOP2, Resilience) %>%
    mutate(allele_effect = abs(`0`) - abs(`1`))
  # on discrete time, resilience varies from -1 to 1, with values closer to zero being more stable.\
  # therefore, I take the absolute value first. Then I substract AOP2 (-1) from aop2 (1). Whenever this difference
  # is greater than zero, it means that aop2 has a stabilizing effect.

  # calculate percentage where aop2 > AOP2 effect
  baseaop2_Resil_full_BayesP <- mean(baseaop2_Resilience_df$allele_effect > 0)

  # calculate aop2 > AOP2 effect size
  baseaop2_Resil_full_effect <- mean(baseaop2_Resilience_df$allele_effect)

  # LYER-Ptoid resilience
  # organize data
  baseaop2_ResilienceLP_df <- all.base_vs_aop2_stability.df %>%
    select(aop2_vs_AOP2, posterior_sample, ResilienceLP) %>%
    spread(aop2_vs_AOP2, ResilienceLP) %>%
    mutate(allele_effect = abs(`0`) - abs(`1`))
  # see note for Full community resilience calculation

  # calculate percentage where aop2 > AOP2 effect
  baseaop2_Resil_LP_BayesP <- mean(baseaop2_ResilienceLP_df$allele_effect > 0)

  # calculate aop2 > AOP2 effect size
  baseaop2_Resil_LP_effect <- mean(baseaop2_ResilienceLP_df$allele_effect)

  # output
  return(list(aop2_r_LYER_effect = aop2_r_LYER_effect,
              aop2_r_LYER_95CI = aop2_r_LYER_95CI,
              aop2_r_Ptoid_effect = aop2_r_Ptoid_effect,
              aop2_r_Ptoid_95CI = aop2_r_Ptoid_95CI,
              aop2_SS_full_BayesP = aop2_SS_full_BayesP,
              aop2_SS_full_effect = aop2_SS_full_effect,
              aop2_SS_LP_BayesP = aop2_SS_LP_BayesP,
              aop2_SS_LP_effect = aop2_SS_LP_effect,
              aop2_SS_LP_95CI = aop2_SS_LP_95CI,
              aop2_Resil_full_BayesP = aop2_Resil_full_BayesP,
              aop2_Resil_full_effect = aop2_Resil_full_effect,
              aop2_Resil_LP_BayesP = aop2_Resil_LP_BayesP,
              aop2_Resil_LP_effect = aop2_Resil_LP_effect,
              baseaop2_SS_full_BayesP = baseaop2_SS_full_BayesP,
              baseaop2_SS_full_effect = baseaop2_SS_full_effect,
              baseaop2_SS_LP_BayesP = baseaop2_SS_LP_BayesP,
              baseaop2_SS_LP_effect = baseaop2_SS_LP_effect,
              baseaop2_Resil_full_BayesP = baseaop2_Resil_full_BayesP,
              baseaop2_Resil_full_effect = baseaop2_Resil_full_effect,
              baseaop2_Resil_LP_BayesP = baseaop2_Resil_LP_BayesP,
              baseaop2_Resil_LP_effect = baseaop2_Resil_LP_effect,
              all.aop2_vs_AOP2_stability.df = all.aop2_vs_AOP2_stability.df,
              base_stability.df = base_stability.df))
}

aop2_vs_AOP2_posterior_samples_adj <- function(brms_model, n.geno = 1, temp.value = 0, logbiomass.value = 0){

  # get posterior predictions
  pp.aop2_vs_AOP2_model <- posterior_samples(brms_model, pars = "^b")

  # get posterior samples of structural stability
  base_stability <- apply(
    pp.aop2_vs_AOP2_model,
    MARGIN = 1,
    FUN = function(x) {
      base.mat = matrix(c(x["bsp_logBRBRt1_melogBRBR_tse_logBRBRt"],
                          x["bsp_logBRBRt1_melogLYER_t_adjse_logLYERt"],
                          x["bsp_logBRBRt1_melogPtoid_t_adjse_logPtoidt"],
                          x["bsp_logLYERt1adj_melogBRBR_tse_logBRBRt"],
                          x["bsp_logLYERt1adj_melogLYER_t_adjse_logLYERt"],
                          x["bsp_logLYERt1adj_melogPtoid_t_adjse_logPtoidt"],
                          x["bsp_logPtoidt1adj_melogBRBR_tse_logBRBRt"],
                          x["bsp_logPtoidt1adj_melogLYER_t_adjse_logLYERt"],
                          x["bsp_logPtoidt1adj_melogPtoid_t_adjse_logPtoidt"]),
                        ncol = 3, byrow = TRUE)
      base.mat[is.na(base.mat)] <- 0 # NA interactions become zeros

      base.r = matrix(c(x["b_logBRBRt1_Intercept"],
                        x["b_logLYERt1adj_Intercept"],
                        x["b_logPtoidt1adj_Intercept"]),
                      ncol = 1)
      base.r[is.na(base.r)] <- 0 # NA r become zeros

      # condition on temperature for intrinsic growth rates
      delta.temp.r = matrix(c(x["b_logBRBRt1_temp"],
                              x["b_logLYERt1_temp"],
                              x["b_logPtoidt1_temp"]),
                            ncol = 1)
      delta.temp.r[is.na(delta.temp.r)] <- 0 # NA r become zeros

      # condition on plant biomass value for intrinsic growth rates
      delta.logbiomass.r = matrix(c(x["b_logBRBRt1_logBiomass_g_t1"],
                                    x["b_logLYERt1_logBiomass_g_t1"],
                                    x["b_logPtoidt1_logBiomass_g_t1"]),
                                  ncol = 1)
      delta.logbiomass.r[is.na(delta.logbiomass.r)] <- 0 # NA r become zeros

      # add to get aop2 allele effect
      tmp.mat <- base.mat - diag(3) # put on 'continuous time scale'
      tmp.r <- base.r + delta.temp.r*temp.value + delta.logbiomass.r*logbiomass.value

      Feasibility = as.numeric(FeasibilityBoundary(tmp.mat, tmp.r)["feasibility"])
      FeasibilityBoundary23 = as.numeric(FeasibilityBoundary(tmp.mat, tmp.r)["alpha.A23"])
      Extinctions = paste0(which(-1*inv(tmp.mat) %*% tmp.r < 0), collapse = ",")
      Resilience = max(Re(eigen((tmp.mat + diag(3)))$values)) # on discrete time scale
      ResilienceLP = max(Re(eigen((tmp.mat[2:3,2:3] + diag(2)))$values)) # on discrete time scale
      FeasibilityBoundaryLYER.Ptoid = as.numeric(BoundaryLYER.Ptoid(tmp.mat[2:3,2:3], tmp.r[2:3])["boundary"])
      r_LYER = tmp.r[2]
      r_Ptoid = tmp.r[3]
      list(tmp.mat = tmp.mat,
           tmp.r = tmp.r,
           Feasibility = Feasibility,
           FeasibilityBoundary23 = FeasibilityBoundary23,
           Extinctions = Extinctions,
           Resilience = Resilience,
           ResilienceLP = ResilienceLP,
           FeasibilityBoundaryLYER.Ptoid = FeasibilityBoundaryLYER.Ptoid,
           r_LYER = r_LYER,
           r_Ptoid = r_Ptoid)
    })
  base_stability.df <- data.frame(
    aop2_vs_AOP2 = 0,
    posterior_sample = 1:nrow(pp.aop2_vs_AOP2_model),
    Feasibility = unlist(lapply(base_stability, FUN = function(x) x$Feasibility)),
    FeasibilityBoundary23 = unlist(lapply(base_stability, FUN = function(x) x$FeasibilityBoundary23)),
    Extinctions = unlist(lapply(base_stability, FUN = function(x) x$Extinctions)),
    Resilience = unlist(lapply(base_stability, FUN = function(x) x$Resilience)),
    ResilienceLP = unlist(lapply(base_stability, FUN = function(x) x$ResilienceLP)),
    FeasibilityBoundaryLYER.Ptoid = unlist(lapply(base_stability, FUN = function(x) x$FeasibilityBoundaryLYER.Ptoid)),
    r_LYER = unlist(lapply(base_stability, FUN = function(x) x$r_LYER)),
    r_Ptoid = unlist(lapply(base_stability, FUN = function(x) x$r_Ptoid))
  )

  aop2_stability <- apply(
    pp.aop2_vs_AOP2_model,
    MARGIN = 1,
    FUN = function(x) {
      base.mat = matrix(c(x["bsp_logBRBRt1_melogBRBR_tse_logBRBRt"],
                          x["bsp_logBRBRt1_melogLYER_t_adjse_logLYERt"],
                          x["bsp_logBRBRt1_melogPtoid_t_adjse_logPtoidt"],
                          x["bsp_logLYERt1adj_melogBRBR_tse_logBRBRt"],
                          x["bsp_logLYERt1adj_melogLYER_t_adjse_logLYERt"],
                          x["bsp_logLYERt1adj_melogPtoid_t_adjse_logPtoidt"],
                          x["bsp_logPtoidt1adj_melogBRBR_tse_logBRBRt"],
                          x["bsp_logPtoidt1adj_melogLYER_t_adjse_logLYERt"],
                          x["bsp_logPtoidt1adj_melogPtoid_t_adjse_logPtoidt"]),
                        ncol = 3, byrow = TRUE)
      base.mat[is.na(base.mat)] <- 0 # NA interactions become zeros

      base.r = matrix(c(x["b_logBRBRt1_Intercept"],
                        x["b_logLYERt1adj_Intercept"],
                        x["b_logPtoidt1adj_Intercept"]),
                      ncol = 1)
      base.r[is.na(base.r)] <- 0 # NA r become zeros

      # aop2 effect
      delta.aop2.mat <- matrix(c(x["bsp_logBRBRt1_melogBRBR_tse_logBRBRt:aop2_genotypes"],
                                 x["bsp_logBRBRt1_melogLYER_t_adjse_logLYERt:aop2_genotypes"],
                                 x["bsp_logBRBRt1_melogPtoid_t_adjse_logPtoidt:aop2_genotypes"],
                                 x["bsp_logLYERt1adj_melogBRBR_tse_logBRBRt:aop2_genotypes"],
                                 x["bsp_logLYERt1adj_melogLYER_t_adjse_logLYERt:aop2_genotypes"],
                                 x["bsp_logLYERt1adj_melogPtoid_t_adjse_logPtoidt:aop2_genotypes"],
                                 x["bsp_logPtoidt1adj_melogBRBR_tse_logBRBRt:aop2_genotypes"],
                                 x["bsp_logPtoidt1adj_melogLYER_t_adjse_logLYERt:aop2_genotypes"],
                                 x["bsp_logPtoidt1adj_melogPtoid_t_adjse_logPtoidt:aop2_genotypes"]),
                               ncol = 3, byrow = TRUE)
      delta.aop2.mat[is.na(delta.aop2.mat)] <- 0 # no aop2 effect if NA

      delta.aop2.r <- matrix(c(x["b_logBRBRt1_aop2_genotypes"],
                               x["b_logLYERt1adj_aop2_genotypes"],
                               x["b_logPtoidt1adj_aop2_genotypes"]),
                             ncol = 1)
      delta.aop2.r[is.na(delta.aop2.r)] <- 0 # no aop2 effect if NA

      # condition on temperature for intrinsic growth rates
      delta.temp.r = matrix(c(x["b_logBRBRt1_temp"],
                              x["b_logLYERt1_temp"],
                              x["b_logPtoidt1_temp"]),
                            ncol = 1)
      delta.temp.r[is.na(delta.temp.r)] <- 0 # NA r become zeros

      # condition on plant biomass value for intrinsic growth rates
      delta.logbiomass.r = matrix(c(x["b_logBRBRt1_logBiomass_g_t1"],
                                    x["b_logLYERt1_logBiomass_g_t1"],
                                    x["b_logPtoidt1_logBiomass_g_t1"]),
                                  ncol = 1)
      delta.logbiomass.r[is.na(delta.logbiomass.r)] <- 0 # NA r become zeros

      # add to get aop2 allele effect
      tmp.mat <- base.mat + delta.aop2.mat*n.geno - diag(3) # put on 'continuous time scale'
      tmp.r <- base.r + delta.aop2.r*n.geno + delta.temp.r*temp.value + delta.logbiomass.r*logbiomass.value

      Feasibility = as.numeric(FeasibilityBoundary(tmp.mat, tmp.r)["feasibility"])
      FeasibilityBoundary23 = as.numeric(FeasibilityBoundary(tmp.mat, tmp.r)["alpha.A23"])
      Extinctions = paste0(which(-1*inv(tmp.mat) %*% tmp.r < 0), collapse = ",")
      Resilience = max(Re(eigen((tmp.mat + diag(3)))$values)) # on discrete time scale
      ResilienceLP = max(Re(eigen((tmp.mat[2:3,2:3] + diag(2)))$values)) # on discrete time scale
      FeasibilityBoundaryLYER.Ptoid = as.numeric(BoundaryLYER.Ptoid(tmp.mat[2:3,2:3], tmp.r[2:3])["boundary"])
      r_LYER = tmp.r[2]
      r_Ptoid = tmp.r[3]
      list(tmp.mat = tmp.mat,
           tmp.r = tmp.r,
           Feasibility = Feasibility,
           FeasibilityBoundary23 = FeasibilityBoundary23,
           Extinctions = Extinctions,
           Resilience = Resilience,
           ResilienceLP = ResilienceLP,
           FeasibilityBoundaryLYER.Ptoid = FeasibilityBoundaryLYER.Ptoid,
           r_LYER = r_LYER,
           r_Ptoid = r_Ptoid)
    })
  aop2_stability.df <- data.frame(
    aop2_vs_AOP2 = 1,
    posterior_sample = 1:nrow(pp.aop2_vs_AOP2_model),
    Feasibility = unlist(lapply(aop2_stability, FUN = function(x) x$Feasibility)),
    FeasibilityBoundary23 = unlist(lapply(aop2_stability, FUN = function(x) x$FeasibilityBoundary23)),
    Extinctions = unlist(lapply(aop2_stability, FUN = function(x) x$Extinctions)),
    Resilience = unlist(lapply(aop2_stability, FUN = function(x) x$Resilience)),
    ResilienceLP = unlist(lapply(aop2_stability, FUN = function(x) x$ResilienceLP)),
    FeasibilityBoundaryLYER.Ptoid = unlist(lapply(aop2_stability, FUN = function(x) x$FeasibilityBoundaryLYER.Ptoid)),
    r_LYER = unlist(lapply(aop2_stability, FUN = function(x) x$r_LYER)),
    r_Ptoid = unlist(lapply(aop2_stability, FUN = function(x) x$r_Ptoid))
  )

  AOP2_stability <- apply(
    pp.aop2_vs_AOP2_model,
    MARGIN = 1,
    FUN = function(x) {
      base.mat = matrix(c(x["bsp_logBRBRt1_melogBRBR_tse_logBRBRt"],
                          x["bsp_logBRBRt1_melogLYER_t_adjse_logLYERt"],
                          x["bsp_logBRBRt1_melogPtoid_t_adjse_logPtoidt"],
                          x["bsp_logLYERt1adj_melogBRBR_tse_logBRBRt"],
                          x["bsp_logLYERt1adj_melogLYER_t_adjse_logLYERt"],
                          x["bsp_logLYERt1adj_melogPtoid_t_adjse_logPtoidt"],
                          x["bsp_logPtoidt1adj_melogBRBR_tse_logBRBRt"],
                          x["bsp_logPtoidt1adj_melogLYER_t_adjse_logLYERt"],
                          x["bsp_logPtoidt1adj_melogPtoid_t_adjse_logPtoidt"]),
                        ncol = 3, byrow = TRUE)
      base.mat[is.na(base.mat)] <- 0 # NA interactions become zeros

      base.r = matrix(c(x["b_logBRBRt1_Intercept"],
                        x["b_logLYERt1adj_Intercept"],
                        x["b_logPtoidt1adj_Intercept"]),
                      ncol = 1)
      base.r[is.na(base.r)] <- 0 # NA r become zeros

      # AOP2 effect
      delta.AOP2.mat <- matrix(c(x["bsp_logBRBRt1_melogBRBR_tse_logBRBRt:AOP2_genotypes"],
                                 x["bsp_logBRBRt1_melogLYER_t_adjse_logLYERt:AOP2_genotypes"],
                                 x["bsp_logBRBRt1_melogPtoid_t_adjse_logPtoidt:AOP2_genotypes"],
                                 x["bsp_logLYERt1adj_melogBRBR_tse_logBRBRt:AOP2_genotypes"],
                                 x["bsp_logLYERt1adj_melogLYER_t_adjse_logLYERt:AOP2_genotypes"],
                                 x["bsp_logLYERt1adj_melogPtoid_t_adjse_logPtoidt:AOP2_genotypes"],
                                 x["bsp_logPtoidt1adj_melogBRBR_tse_logBRBRt:AOP2_genotypes"],
                                 x["bsp_logPtoidt1adj_melogLYER_t_adjse_logLYERt:AOP2_genotypes"],
                                 x["bsp_logPtoidt1adj_melogPtoid_t_adjse_logPtoidt:AOP2_genotypes"]),
                               ncol = 3, byrow = TRUE)
      delta.AOP2.mat[is.na(delta.AOP2.mat)] <- 0 # no aop2 effect if NA

      delta.AOP2.r <- matrix(c(x["b_logBRBRt1_AOP2_genotypes"],
                               x["b_logLYERt1adj_AOP2_genotypes"],
                               x["b_logPtoidt1adj_AOP2_genotypes"]),
                             ncol = 1)
      delta.AOP2.r[is.na(delta.AOP2.r)] <- 0 # no aop2 effect if NA

      # condition on temperature for intrinsic growth rates
      delta.temp.r = matrix(c(x["b_logBRBRt1_temp"],
                              x["b_logLYERt1_temp"],
                              x["b_logPtoidt1_temp"]),
                            ncol = 1)
      delta.temp.r[is.na(delta.temp.r)] <- 0 # NA r become zeros

      # condition on plant biomass value for intrinsic growth rates
      delta.logbiomass.r = matrix(c(x["b_logBRBRt1_logBiomass_g_t1"],
                                    x["b_logLYERt1_logBiomass_g_t1"],
                                    x["b_logPtoidt1_logBiomass_g_t1"]),
                                  ncol = 1)
      delta.logbiomass.r[is.na(delta.logbiomass.r)] <- 0 # NA r become zeros

      # add to get AOP2 allele
      tmp.mat <- base.mat + delta.AOP2.mat*n.geno - diag(3) # put on 'continuous time scale'
      tmp.r <- base.r + delta.AOP2.r*n.geno + delta.temp.r*temp.value + delta.logbiomass.r*logbiomass.value

      Feasibility = as.numeric(FeasibilityBoundary(tmp.mat, tmp.r)["feasibility"])
      FeasibilityBoundary23 = as.numeric(FeasibilityBoundary(tmp.mat, tmp.r)["alpha.A23"])
      Extinctions = paste0(which(-1*inv(tmp.mat) %*% tmp.r < 0), collapse = ",")
      Resilience = max(Re(eigen((tmp.mat + diag(3)))$values)) # on discrete time scale
      ResilienceLP = max(Re(eigen((tmp.mat[2:3,2:3] + diag(2)))$values)) # on discrete time scale
      FeasibilityBoundaryLYER.Ptoid = as.numeric(BoundaryLYER.Ptoid(tmp.mat[2:3,2:3], tmp.r[2:3])["boundary"])
      r_LYER = tmp.r[2]
      r_Ptoid = tmp.r[3]
      list(tmp.mat = tmp.mat,
           tmp.r = tmp.r,
           Feasibility = Feasibility,
           FeasibilityBoundary23 = FeasibilityBoundary23,
           Extinctions = Extinctions,
           Resilience = Resilience,
           ResilienceLP = ResilienceLP,
           FeasibilityBoundaryLYER.Ptoid = FeasibilityBoundaryLYER.Ptoid,
           r_LYER = r_LYER,
           r_Ptoid = r_Ptoid)
    })
  AOP2_stability.df <- data.frame(
    aop2_vs_AOP2 = -1,
    posterior_sample = 1:nrow(pp.aop2_vs_AOP2_model),
    Feasibility = unlist(lapply(AOP2_stability, FUN = function(x) x$Feasibility)),
    FeasibilityBoundary23 = unlist(lapply(AOP2_stability, FUN = function(x) x$FeasibilityBoundary23)),
    Extinctions = unlist(lapply(AOP2_stability, FUN = function(x) x$Extinctions)),
    Resilience = unlist(lapply(AOP2_stability, FUN = function(x) x$Resilience)),
    ResilienceLP = unlist(lapply(AOP2_stability, FUN = function(x) x$ResilienceLP)),
    FeasibilityBoundaryLYER.Ptoid = unlist(lapply(AOP2_stability, FUN = function(x) x$FeasibilityBoundaryLYER.Ptoid)),
    r_LYER = unlist(lapply(AOP2_stability, FUN = function(x) x$r_LYER)),
    r_Ptoid = unlist(lapply(AOP2_stability, FUN = function(x) x$r_Ptoid))
  )

  # combine data
  all.aop2_vs_AOP2_stability.df <- bind_rows(aop2_stability.df, AOP2_stability.df)

  ## Aphid and parasitoid intrinsic growth rates ----
  r_LYER_df <- all.aop2_vs_AOP2_stability.df %>%
    select(aop2_vs_AOP2, posterior_sample, r_LYER) %>%
    spread(key = aop2_vs_AOP2, value = r_LYER) %>%
    mutate(allele_effect = `1` - `-1`)
  aop2_r_LYER_effect <- mean(r_LYER_df$allele_effect)
  aop2_r_LYER_95CI <- quantile(r_LYER_df$allele_effect, probs = c(0.025, 0.975))

  r_Ptoid_df <- all.aop2_vs_AOP2_stability.df %>%
    select(aop2_vs_AOP2, posterior_sample, r_Ptoid) %>%
    spread(key = aop2_vs_AOP2, value = r_Ptoid) %>%
    mutate(allele_effect = `1` - `-1`)
  aop2_r_Ptoid_effect <- mean(r_Ptoid_df$allele_effect)
  aop2_r_Ptoid_95CI <- quantile(r_Ptoid_df$allele_effect, probs = c(0.025, 0.975))

  ## Full community structural stability ----
  # organize data
  FeasibilityBoundary23_df <- all.aop2_vs_AOP2_stability.df %>%
    select(aop2_vs_AOP2, posterior_sample, FeasibilityBoundary23) %>%
    spread(key = aop2_vs_AOP2, value = FeasibilityBoundary23) %>%
    mutate(allele_effect = `1` - `-1`)

  # calculate percentage where aop2 > AOP2 effect
  aop2_SS_full_BayesP <- mean(FeasibilityBoundary23_df$allele_effect > 0)

  # calculate aop2 > AOP2 effect size
  aop2_SS_full_effect <- mean(FeasibilityBoundary23_df$allele_effect)

  ## LYER-Ptoid Structural stability ----
  # organize data
  aop2_LPbound_df <- all.aop2_vs_AOP2_stability.df %>%
    select(aop2_vs_AOP2, posterior_sample, FeasibilityBoundaryLYER.Ptoid) %>%
    spread(aop2_vs_AOP2, FeasibilityBoundaryLYER.Ptoid) %>%
    mutate(allele_effect = `1` - `-1`)

  # calculate percentage where aop2 > AOP2 effect
  aop2_SS_LP_BayesP <- mean(aop2_LPbound_df$allele_effect > 0)

  # calculate aop2 > AOP2 effect size
  aop2_SS_LP_effect <- mean(aop2_LPbound_df$allele_effect)
  aop2_SS_LP_95CI <- quantile(aop2_LPbound_df$allele_effect, probs = c(0.025, 0.975))

  ## Full community resilience ----
  # organize data
  aop2_Resilience_df <- all.aop2_vs_AOP2_stability.df %>%
    select(aop2_vs_AOP2, posterior_sample, Resilience) %>%
    spread(aop2_vs_AOP2, Resilience) %>%
    mutate(allele_effect = abs(`-1`) - abs(`1`))
  # on discrete time, resilience varies from -1 to 1, with values closer to zero being more stable.\
  # therefore, I take the absolute value first. Then I substract AOP2 (-1) from aop2 (1). Whenever this difference
  # is greater than zero, it means that aop2 has a stabilizing effect.

  # calculate percentage where aop2 > AOP2 effect
  aop2_Resil_full_BayesP <- mean(aop2_Resilience_df$allele_effect > 0)

  # calculate aop2 > AOP2 effect size
  aop2_Resil_full_effect <- mean(aop2_Resilience_df$allele_effect)

  ## LYER-Ptoid resilience ----
  # organize data
  aop2_ResilienceLP_df <- all.aop2_vs_AOP2_stability.df %>%
    select(aop2_vs_AOP2, posterior_sample, ResilienceLP) %>%
    spread(aop2_vs_AOP2, ResilienceLP) %>%
    mutate(allele_effect = abs(`-1`) - abs(`1`))
  # see note for Full community resilience calculation

  # calculate percentage where aop2 > AOP2 effect
  aop2_Resil_LP_BayesP <- mean(aop2_ResilienceLP_df$allele_effect > 0)

  # calculate aop2 > AOP2 effect size
  aop2_Resil_LP_effect <- mean(aop2_ResilienceLP_df$allele_effect)

  ## Baseline comparisons
  # aop2
  all.base_vs_aop2_stability.df <- bind_rows(aop2_stability.df, base_stability.df)

  # Full community structural stability
  # organize data
  baseaop2_FeasibilityBoundary23_df <- all.base_vs_aop2_stability.df %>%
    select(aop2_vs_AOP2, posterior_sample, FeasibilityBoundary23) %>%
    spread(key = aop2_vs_AOP2, value = FeasibilityBoundary23) %>%
    mutate(allele_effect = `1` - `0`)

  # calculate percentage where aop2 > AOP2 effect
  baseaop2_SS_full_BayesP <- mean(baseaop2_FeasibilityBoundary23_df$allele_effect > 0)

  # calculate aop2 > AOP2 effect size
  baseaop2_SS_full_effect <- mean(baseaop2_FeasibilityBoundary23_df$allele_effect)

  # LYER-Ptoid Structural stability
  # organize data
  baseaop2_LPbound_df <- all.base_vs_aop2_stability.df %>%
    select(aop2_vs_AOP2, posterior_sample, FeasibilityBoundaryLYER.Ptoid) %>%
    spread(aop2_vs_AOP2, FeasibilityBoundaryLYER.Ptoid) %>%
    mutate(allele_effect = `1` - `0`)

  # calculate percentage where aop2 > AOP2 effect
  baseaop2_SS_LP_BayesP <- mean(baseaop2_LPbound_df$allele_effect > 0)

  # calculate aop2 > AOP2 effect size
  baseaop2_SS_LP_effect <- mean(baseaop2_LPbound_df$allele_effect)

  # Full community resilience
  # organize data
  baseaop2_Resilience_df <- all.base_vs_aop2_stability.df %>%
    select(aop2_vs_AOP2, posterior_sample, Resilience) %>%
    spread(aop2_vs_AOP2, Resilience) %>%
    mutate(allele_effect = abs(`0`) - abs(`1`))
  # on discrete time, resilience varies from -1 to 1, with values closer to zero being more stable.\
  # therefore, I take the absolute value first. Then I substract AOP2 (-1) from aop2 (1). Whenever this difference
  # is greater than zero, it means that aop2 has a stabilizing effect.

  # calculate percentage where aop2 > AOP2 effect
  baseaop2_Resil_full_BayesP <- mean(baseaop2_Resilience_df$allele_effect > 0)

  # calculate aop2 > AOP2 effect size
  baseaop2_Resil_full_effect <- mean(baseaop2_Resilience_df$allele_effect)

  # LYER-Ptoid resilience
  # organize data
  baseaop2_ResilienceLP_df <- all.base_vs_aop2_stability.df %>%
    select(aop2_vs_AOP2, posterior_sample, ResilienceLP) %>%
    spread(aop2_vs_AOP2, ResilienceLP) %>%
    mutate(allele_effect = abs(`0`) - abs(`1`))
  # see note for Full community resilience calculation

  # calculate percentage where aop2 > AOP2 effect
  baseaop2_Resil_LP_BayesP <- mean(baseaop2_ResilienceLP_df$allele_effect > 0)

  # calculate aop2 > AOP2 effect size
  baseaop2_Resil_LP_effect <- mean(baseaop2_ResilienceLP_df$allele_effect)

  # output
  return(list(aop2_r_LYER_effect = aop2_r_LYER_effect,
              aop2_r_LYER_95CI = aop2_r_LYER_95CI,
              aop2_r_Ptoid_effect = aop2_r_Ptoid_effect,
              aop2_r_Ptoid_95CI = aop2_r_Ptoid_95CI,
              aop2_SS_full_BayesP = aop2_SS_full_BayesP,
              aop2_SS_full_effect = aop2_SS_full_effect,
              aop2_SS_LP_BayesP = aop2_SS_LP_BayesP,
              aop2_SS_LP_effect = aop2_SS_LP_effect,
              aop2_SS_LP_95CI = aop2_SS_LP_95CI,
              aop2_Resil_full_BayesP = aop2_Resil_full_BayesP,
              aop2_Resil_full_effect = aop2_Resil_full_effect,
              aop2_Resil_LP_BayesP = aop2_Resil_LP_BayesP,
              aop2_Resil_LP_effect = aop2_Resil_LP_effect,
              baseaop2_SS_full_BayesP = baseaop2_SS_full_BayesP,
              baseaop2_SS_full_effect = baseaop2_SS_full_effect,
              baseaop2_SS_LP_BayesP = baseaop2_SS_LP_BayesP,
              baseaop2_SS_LP_effect = baseaop2_SS_LP_effect,
              baseaop2_Resil_full_BayesP = baseaop2_Resil_full_BayesP,
              baseaop2_Resil_full_effect = baseaop2_Resil_full_effect,
              baseaop2_Resil_LP_BayesP = baseaop2_Resil_LP_BayesP,
              baseaop2_Resil_LP_effect = baseaop2_Resil_LP_effect,
              all.aop2_vs_AOP2_stability.df = all.aop2_vs_AOP2_stability.df,
              base_stability.df = base_stability.df))
}
