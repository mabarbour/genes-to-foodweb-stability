require(tidyverse)

# GLM function for manual F-tests
glm.ftest.v2 <- function(model, test.formula, VC.formula = NULL, harmonic.means = NULL){
  # tidy and add mean_deviance to anova table
  original_F <- anova(model, test="F") %>%
    .[-1, ] %>% # remove NULL
    mutate(term = rownames(.)) %>%
    add_row(term = "Residuals", Df = last(.$`Resid. Df`), Deviance = last(.$`Resid. Dev`)) %>%
    select(term, df = Df, dev = Deviance, F, P = `Pr(>F)`) %>%
    mutate(mean_dev = round(dev / .$df, 3),
           dev = round(dev, 3),
           F = round(F, 2),
           P = round(P, 4)) %>%
    select(term, df, dev, mean_dev, F, P)

  manual_F <- anova(model, test="F") %>%
    .[-1, ] %>% # remove NULL
    mutate(term = rownames(.)) %>%
    add_row(term = "Residuals", Df = last(.$`Resid. Df`), Deviance = last(.$`Resid. Dev`)) %>%
    select(term, df = Df, dev = Deviance) %>%
    mutate(mean_dev = dev / df,
           F = NA,
           P = NA)
  rownames(manual_F) <- manual_F$term

  for(i in 1:length(test.formula)){
    # organize formulas
    lhs.side <- test.formula[[i]][1]
    rhs.side <- test.formula[[i]][2]

    manual_F[lhs.side,"F"] <- manual_F[lhs.side,"mean_dev"] / manual_F[rhs.side,"mean_dev"]
    manual_F[lhs.side,"P"] <- pf(manual_F[lhs.side,"F"],
                                 manual_F[lhs.side,"df"],
                                 manual_F[rhs.side,"df"],
                                 lower.tail = FALSE)

  }

  summary.manual_F <- data.frame(treatment = rep(NA, length(test.formula)),
                                 num_df = rep(NA, length(test.formula)),
                                 den_df = rep(NA, length(test.formula)),
                                 deviance = rep(NA, length(test.formula)),
                                 mean_deviance = rep(NA, length(test.formula)),
                                 F = rep(NA, length(test.formula)),
                                 P = rep(NA, length(test.formula)),
                                 error = rep(NA, length(test.formula)))

  for(j in 1:length(test.formula)){
    # organize formulas
    lhs.side <- test.formula[[j]][1]
    rhs.side <- test.formula[[j]][2]

    summary.manual_F[j,"treatment"] <- lhs.side
    summary.manual_F[j,"num_df"] <- manual_F[lhs.side,"df"]
    summary.manual_F[j,"den_df"] <- manual_F[rhs.side,"df"]
    summary.manual_F[j,"deviance"] <- round(manual_F[lhs.side,"dev"],2)
    summary.manual_F[j,"mean_deviance"] <- round(manual_F[lhs.side,"mean_dev"],2)
    summary.manual_F[j,"F"] <- round(manual_F[lhs.side,"F"],3)
    summary.manual_F[j,"P"] <- round(manual_F[lhs.side,"P"],3)
    summary.manual_F[j,"error"] <- rhs.side
    summary.manual_F[j,"sig"] <- ifelse(manual_F[lhs.side,"P"] < 0.001, "***",
                                        ifelse(manual_F[lhs.side,"P"] < 0.01, "**",
                                               ifelse(manual_F[lhs.side,"P"] < 0.05, "*",
                                                      ifelse(manual_F[lhs.side,"P"] > 0.05 & manual_F[lhs.side,"P"] < 0.1, ".",""))))
  }

  if(is.null(VC.formula) == F){
    summary.manual_VC <- data.frame(error = rep(NA, length(VC.formula)),
                                    deviance_component = rep(NA, length(VC.formula)))

    for(j in 1:length(VC.formula)){

      lhs.side <- VC.formula[[j]][1]
      rhs.side <- VC.formula[[j]][2]

      summary.manual_VC[j,"error"] <- lhs.side
      summary.manual_VC[j,"deviance_component"] <- round((manual_F[lhs.side,"mean_dev"] - manual_F[rhs.side,"mean_dev"])/harmonic.means[j],2)
    }

    summary.manual_VC <- summary.manual_VC %>%
      add_row(error = "Residuals", deviance_component = round(manual_F["Residuals","mean_dev"], 2))
  } else{
    summary.manual_VC <- NULL
  }


  manual_F <- manual_F %>%
    # to improve readibility
    mutate(dev = round(dev, 3),
           mean_dev = round(mean_dev, 3),
           F = round(F, 3), P = round(P, 4))

  summary.manual_F <- summary.manual_F %>%
    mutate(P = ifelse(P < 0.001, "<0.001", P))

  return(list(original_F, manual_F, summary.manual_F, summary.manual_VC))
}
