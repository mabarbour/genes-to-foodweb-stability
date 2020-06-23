## Functions used for non-equilibrium simulations ----

# Author: Matthew A. Barbour

## Function to simulate community dynamics from different initial states ----
GenericCommunityDynamics <- function(Initial.States, IGR.Vector, Interaction.Matrix, Duration){

  # initial addition of herbivores into the experiment, 4 individuals of each species
  initial_state <- Initial.States

  # for loop of temporal dynamics
  Nt1 <- list() # abundance at next time step
  for(i in 1:Duration){
    if(i == 1){
      Nt1[[i]] <- initial_state
    }
    if(i > 1){
      Nt1[[i]] <- IGR.Vector + Interaction.Matrix %*% Nt1[[i-1]]
    }
  }
  Nt1.df <- t(as.data.frame(Nt1)) %>%
    as_tibble() %>%
    mutate(Week = 0:(length(Nt1)-1))
}

## Simulate community dynamics from initial aphid additions and later addition of parasitoids ----
SimulateCommunityDynamics <- function(IGR.Vector, Interaction.Matrix, Duration){

  # initial addition of herbivores into the experiment, 4 individuals of each species
  initial_state <- matrix(c(log(4), log(4)), nrow = 2, dimnames = list(c("BRBR","LYER"),"r"))

  # for loop of temporal dynamics
  Nt1 <- list() # abundance at next time step
  for(i in 1:Duration){
    if(i == 1){
      Nt1[[i]] <- initial_state
    }
    if(i == 2){
      Nt1[[i]] <- IGR.Vector[1:2,] + Interaction.Matrix[1:2,1:2] %*% Nt1[[i-1]]
    }
    if(i == 3){
      tmp.Nt1 <- IGR.Vector[1:2,] + Interaction.Matrix[1:2,1:2] %*% Nt1[[i-1]]
      Nt1[[i]] <- rbind(tmp.Nt1, Ptoid = log(2))
    }
    if(i == 4){
      tmp.Nt1 <- IGR.Vector + Interaction.Matrix %*% Nt1[[i-1]]
      tmp.Nt1[3,] <- tmp.Nt1[3,] + log(2)
      Nt1[[i]] <- tmp.Nt1
    }
    if(i > 3){
      Nt1[[i]] <- IGR.Vector + Interaction.Matrix %*% Nt1[[i-1]]
    }
  }
  Nt1[[1]] <- rbind(Nt1[[1]], Ptoid = 0)
  Nt1[[2]] <- rbind(Nt1[[2]], Ptoid = 0)
  Nt1.df <- t(as.data.frame(Nt1)) %>%
    as_tibble() %>%
    mutate(Week = 0:(length(Nt1)-1))
}

## Generic plot for simulated dynamics ----
PlotCommunityDynamics <- function(x){ # SimulateCommunityDynamics object
  x %>%
    gather(Species, value = "Abundance", -Week) %>%
    ggplot(., aes(x = Week, y = Abundance, color = Species)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    scale_y_continuous(breaks = 0:100) +
    scale_x_continuous(breaks = 0:100)
}

## Determine which and at what time the first species went extinct ----
first_extinction <- function(IGR.Vector, Interaction.Matrix, Duration = 17){
  sim_df <- SimulateCommunityDynamics(IGR.Vector, Interaction.Matrix, Duration = Duration)
  # get first week when each species goes extinct
  extinct_df <- data.frame(
    species = c("BRBR","LYER","Ptoid"),
    week = c(filter(sim_df, BRBR < 0)$Week[1], filter(sim_df, LYER < 0)$Week[1], filter(sim_df, Ptoid < 0)$Week[1])
  )
  first_df <- slice(extinct_df, which.min(week))
  if(nrow(first_df) == 0){
    first_df <- data.frame(species = NA, week = NA)
  }
  return(first_df)
}

## Determine which and at what time LYER or Ptoid go extinct in the remaining food chain ----
first_extinction_2sp <- function(Initial.States, IGR.Vector, Interaction.Matrix, Duration = 10, threshold = 0){
  sim_df <- GenericCommunityDynamics(Initial.States, IGR.Vector, Interaction.Matrix, Duration = Duration)
  # get first week when each species goes extinct
  extinct_df <- data.frame(
    species = c("LYER","Ptoid"),
    week = c(filter(sim_df, LYER < threshold)$Week[1], filter(sim_df, Ptoid < threshold)$Week[1])
  )
  first_df <- slice(extinct_df, which.min(week))
  if(nrow(first_df) == 0){
    first_df <- data.frame(species = NA, week = NA)
  }
  return(first_df)
}
