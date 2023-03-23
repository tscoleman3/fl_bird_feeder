## --------------- HEADER ------------------------------------------------------
## Script name: 2g_Experiment-2-clean-seeds.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2023-01-23
## Date Last Modified: 2023-02-10
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script cleans the seed data

## We import the raw seed rain matrix. Then, we trade some of the
## placeholder names with an unknown identifier, remove seeds that were
## potential contaminants from the feeder, and arrange the data in 
## alphabetical order. We export this clean data.

## We also filter seeds not likely dispersed by birds (i.e., those from the 
## canopy tree species, wind-dispersed seeds, and epizoochorous seeds 
## dispersed by mammals. We export this filtered data for supplementary models.

## --------------- SET—UP WORKSPACE --------------------------------------------

rm(list=ls())

# Bring in the data
seeds <- read.csv(file = "data/Experiment-2-raw-seeds.csv",
                  header = TRUE, stringsAsFactors = FALSE)

# Fix site label
seeds[147,3] <- 4

# Replace placeholder names with unknown titles

colnames(seeds)[14] <- "UNK16"
colnames(seeds)[16] <- "UNK17"
colnames(seeds)[17] <- "UNK18"
colnames(seeds)[22] <- "UNK19" # Plantago1 
colnames(seeds)[24] <- "UNK20" # Plantago2

# Drop seeds that may have come from the feeder
seeds <- seeds[, -c(9,12,19,20,23)]

# Rearrange columns
seeds <- seeds |>
  dplyr::select(DATE, BLOCK, TRAP, TREATMENT, ACORN, ASTER1, ASTER2, 
                BRASS, DESMO, CONVULV, PARTH, PINE,RHUS, 
                RUBIAC, SMILAX, UNK6, UNK7, UNK8, UNK9, UNK10, UNK11, UNK12, 
                UNK13, UNK14, UNK15,UNK16, UNK17, UNK18,
                UNK19, UNK20, NOTES)

write.csv(seeds, "data/Experiment-2-clean-seeds.csv",
          row.names = FALSE)

# Drop seeds that are likely wind-dispersed or falling from canopy
seeds <- seeds[, -c(5,6,7,9,12)]
# Acorn, Aster, Aster, Desmodium, Pine

# Rearrange columns
seeds <- seeds |>
  dplyr::select(DATE, BLOCK, TRAP, TREATMENT, BRASS, CONVULV, PARTH, RHUS, 
                RUBIAC, SMILAX, UNK6, UNK7, UNK8, UNK9, UNK10, UNK11, UNK12, 
                UNK13, UNK14, UNK15,UNK16, UNK17, UNK18, UNK19, UNK20, NOTES)

write.csv(seeds, "data/Experiment-2-clean-seeds-filt.csv",
          row.names = FALSE)

