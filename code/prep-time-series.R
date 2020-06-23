## Author: Matthew A. Barbour

## Prep ----

# load required libraries
library(tidyverse)

# load and manage insect time-series data
df <- read_csv("data/arabidopsis_clean_df.csv") %>%
  rename(com = Composition,
         temp = Temperature,
         rich = Richness) %>%
  mutate(frich = as.factor(rich),
         Cage = as.character(Cage))

# load and manage plant biomass time-series data
tmp.biomass <- read_csv("data/ExperimentPlantBiomass.csv") %>%
  select(Cage, Week_match = Week, No_Plants, Biomass_g) %>% # No_Plants not necessary because we're using Biomass_g as a covariate rather than a response
  mutate(Cage = as.character(Cage))
  # note that there is no biomass data on first week, because plants were left in for 2 weeks
  # we are going to interpolate this biomass based on the plants growth at
  # week 2 of the experiment

get_biomass_t1 <- tmp.biomass %>%
  mutate(Week_match = Week_match - 1) %>%
  select(Cage, Week_match, Biomass_g_t1 = Biomass_g)

merge.biomass.df <- left_join(get_biomass_t1, tmp.biomass) %>%
  # assuming a linear relationship in plant growth between weeks, since
  # they had not hit the saturation point on their growth curve
  mutate(Biomass_g = ifelse(Week_match == 1, Biomass_g_t1 - Biomass_g_t1/6, Biomass_g),
         No_Plants = ifelse(Week_match == 1, 8, No_Plants)) %>%
  arrange(Cage, Week_match)

biomass.0.df <- merge.biomass.df %>%
  mutate(Week_match = Week_match - 1) %>%
  select(Cage, Week_match, Biomass_g_t1 = Biomass_g) %>%
  filter(Week_match == 0) %>%
  mutate(Biomass_g = Biomass_g_t1 - Biomass_g_t1/5,
         No_Plants = 8)

biomass.df <- bind_rows(filter(merge.biomass.df, Week_match > 0), biomass.0.df) %>%
  arrange(Cage, Week_match)


# Useful functions
mean_integer <- function(x) as.integer(mean(x, na.rm=T)) # remove NA, otherwise, the count will be set to zero if only 1 counter counted.

replace_NA_zero <- function(x) ifelse(is.na(x)==TRUE, 0, x) # for this analysis, replacing NA with zero is meaningful, since there weren't any individuals of the species in the cage.

# For loop to get lagged insect time-series data for each cage
cage_i_fits <- list()
for(i in 1:60){
  cage_id <- i

  # Get lagged data for each cage
  time_t <- df %>%
    arrange(Cage, Week) %>%
    filter(Cage == i) %>% # grab first cage for testing
    # add a small constant, different for aphids and parasitoids, whenever the insects were labelled as surviving.
    # Note that this also accurately changes things for the parasitoids in the first two weeks of the experiment.
    mutate(BRBR_t = ifelse(BRBR_Survival == 1, BRBR + 5, 0),
           LYER_t = ifelse(LYER_Survival == 1, LYER + 5, 0),
           Ptoid_t = ifelse(Mummy_Ptoids_Survival == 1, Mummy_Ptoids + 1, 0)) %>%
    select(temp, rich, Col, gsm1, AOP2.gsoh, AOP2, com, Cage, Week_match = Week, BRBR_t, LYER_t, Ptoid_t) %>%
    group_by(temp, rich, Col, gsm1, AOP2.gsoh, AOP2, com, Cage, Week_match) %>%
    summarise_all(funs(mean_integer)) %>% # average across counters
    ungroup() %>%
    add_row(Cage = i,
            temp=first(.$temp), rich=first(.$rich), Col=first(.$Col), gsm1=first(.$gsm1), AOP2.gsoh=first(.$AOP2.gsoh), AOP2=first(.$AOP2), com=first(.$com),
            Week_match = 0, BRBR_t = 4, LYER_t = 4, Ptoid_t = NA, .before = 1) # adding initial values of each insect to data

  time_t1 <- filter(time_t, Week_match > 0) %>%
    select(temp, rich, Col, gsm1, AOP2.gsoh, AOP2, com, Cage, Week_match, BRBR_t1 = BRBR_t, LYER_t1 = LYER_t, Ptoid_t1 = Ptoid_t) %>%
    mutate(Week_match = Week_match - 1) # subtract 1 week so everything matches up

  lagged_df <- left_join(time_t, time_t1) %>%
    arrange(Cage, Week_match) %>%
    mutate_at(vars(BRBR_t:Ptoid_t1), funs(replace_NA_zero)) # replace NA values with zeros
    # note that this creates zeros for all insects at t1 for the last week of the experiment
    # therefore, we will only keep data for prediction through week 16

  lagged_df$Ptoid_t[which(lagged_df$Week_match == 2)] <- 2 # added two initial parasitoid females to experiment
  lagged_df$Ptoid_t[which(lagged_df$Week_match == 3)] <- lagged_df$Ptoid_t[which(lagged_df$Week_match == 3)] + 2 # added two additional parasitoid females
  # before we were ignoring the inital increase of parasitoids from week 2 to week 3, since they are at zeros in the dataset, although I added 2 individuals after week 2 and 3 counts.
  # I added the two lines above to incorporate this info

  cage_i_fits[[i]] <- lagged_df %>%
    mutate(fweek = factor(Week_match))
}
cage_fits_df <- plyr::ldply(cage_i_fits) %>%
  left_join(., biomass.df) %>% # add biomass data
  filter(Week_match < 17) # ensures that we only keep data for prediction through week 16

# write out time series data for further analysis
write_csv(cage_fits_df, "output/timeseries_df.csv")
