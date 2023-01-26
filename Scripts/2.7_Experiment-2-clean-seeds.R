## --------------- HEADER ------------------------------------------------------
## Script name: 4_clean-seeds.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliaton: University of Florida
## Date Created: 2023-01-23
## Date Last Modified: 2023-01-23
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script cleans the seed data

## --------------- SET—UP WORKSPACE --------------------------------------------

rm(list=ls())

# Bring in the data
seeds <- read.csv(file = "data/raw_seed_traps.csv",
                  header = TRUE, stringsAsFactors = FALSE)

# Replace placeholder names with unknown titles
colnames(seeds)[14] <- "UNK16"
colnames(seeds)[16] <- "UNK17"
colnames(seeds)[17] <- "UNK18"

# Drop seeds that came from the feeder
seeds <- seeds[, -c(12,19,20,23)]

# Rearrange columns
seeds <- seeds |>
  dplyr::select(DATE, BLOCK, TRAP, TREATMENT, ACORN, ASTER1, ASTER2, BRASS, 
                DESMO, DICO, CONVULV, PARTH, PINE, PLANTAG, PLANTAG2, RHUS, 
                RUBIAC, SMILAX, UNK6, UNK7, UNK8, UNK9, UNK10, UNK11, UNK12, 
                UNK13, UNK14, UNK15,UNK16, UNK17, UNK18, NOTES)

write.csv(seeds, "data/clean_seed_traps.csv",
          row.names = FALSE)
