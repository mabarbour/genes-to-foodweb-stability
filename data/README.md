# Data

## `InsectAbundanceSurvival.csv`

Contains data on experimental treatments, insect abundances, and observed extinctions. Metadata for each column is given below:

- `Cage`: Identity of cage (1--60), coded as a character for analysis.
- `Counter`: Initials of person that counted insects; DTV = Daniel Trujillo-Villegas, MAB = Matthew A. Barbour. For cages that were independently counted by DTV and MAB on the same date, we averaged across them to get insect abundances for further analysis.
- `Sample_Order`: Order in which cages were counted (1--60). Not used in any analysis, but we checked in order to ensure that counting fatigue did not bias our estimates.
- `Date`: Date insects were counted.
- `Week`: Week of experiment (1--17).
- `Temperature`: Treatment temperature of 20 C or 23 C. For analyses, this variables was renamed to `temp`.
- `Richness`: Number of plant genotypes (1, 2, or 4). For analyses, this variables was renamed to `rich`.
- `Composition`: Genetic composition of each plant population (n = 11). Genotypes are separated by an underscore, but the 4-genotype mixture is called 'Poly'. For analyses, this variable was renamed to `com`.
- `AOP2`: Presence (1) or absence (0) of genotype AOP2 in plant population.
- `AOP2.gsoh`: Presence (1) or absence (0) of AOP2/gsoh genotype in plant population.
- `Col`: Presence (1) or absence (0) of Col genotype in plant population.
- `gsm1`: Presence (1) or absence (0) of gsm1 genotype in plant population.
- `No_Plants`: Total number of plants, across both pots, placed into the cage. This was 8 the vast majority of the time, but occassionally 7 due to mortality while growing in the greenhouse.
- `BRBR`: Number of Brevicoryne brassicae individuals. Aphid individuals were counted to a resolution of 1 individual in first 2 weeks, and to 5 individuals in subsequent weeks to allow for efficient counting of often high aphid abundances in all cages within one sampling date.
- `LYER`: Number of Lipaphis erysimi individuals. Aphid individuals were counted to a resolution of 1 individual in first 2 weeks, and to 5 individuals in subsequent weeks to allow for efficient counting of often high aphid abundances in all cages within one sampling date.
- `Mummy_raw`: Number of mummified aphids (individual level).
- `Ptoids_raw`: Number of adult parasitoids (individual level).
- `Mummy_Ptoids`: Combined number of mummified aphids and adult parasitoid individuals. These counts were always to the individual level.
- `*INSECT*_Survival`: Insect population had positive abundance (1) or went extinct (0). Samples after extinction were coded 'NA', except for Mummy_Ptoids, where NA also includes the first two weeks before they were added to the experiment.

## `insect_abundance_data_2018-09-26_underneathleaf.csv`

Contains data on *Lipaphis erysimi* (`LYER`) and mummified aphids (`Mummy`) counted underneath the basal rosettes of plants in the last week of the experiment. This data was combined with the usual abundance data to adjust *L. erysimi* and *D. rapae* counts and assess the robustness of our inferences. It also includes a count of adult parasitoids (`Ptoids_2018_10_02_biomass_collection`) noticed in some cages after the experiment (plants and aphids removed) on the final biomass collection date (2019-10-02, Y-M-D). Columns `Cage`, `BRBR`, and `Counter` are same as `InsectAbundanceSurvival.csv`. `BRBR` is always zero because they never survived until the end of the experiment.

## `ExperimentPlantBiomass.csv`

Contains data on plant biomass over the course of the experiment. Most of the data columns are the same as `InsectAbundanceSurvival.csv`. The one new data column is:

- `Biomass_g`: Oven-dried plant biomass in grams. 'NA' values in `Week` = 1, reflects the fact that plants grew in the experimental cages for 2 weeks at the initiation of the experiment. Therefore, we were not able to quantify their biomass in the first week.

## `PreExperimentNoInsectsPlantBiomass.csv`

Contains data on plant biomass in each experimental treatment, but in the absence of insects. These growing conditions reflected the start of the experiment, where plants were kept in cages for two weeks. Most of the data columns are the same as `InsectAbundanceSurvival.csv` and `ExperimentPlantBiomass.csv`. The new data columns are:

- `Pot`: Identity of each pot (coded 17 or 18) within each cage. Note that we averaged across pots for this analysis.
- `NumberPlants`: Number of plants within each pot.


