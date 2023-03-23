## --------------- HEADER ------------------------------------------------------
## Script name: 2j_Experiment-2-summarize-seed-detect.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2023-01-23
## Date Last Modified: 2023-02-10
## Copyright (c) David S. Mason, 2023
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script summarizes bird data from experiment 2

library(tidyverse)
library(lubridate)

# Clear the decks
rm(list=ls())

# Bring in the full data
seeds <- read.csv(file = "data/Experiment-2-clean-seeds.csv",
                        header = TRUE)
head(seeds)

# Calculate the total number of seeds and morphotypes
row.sum <- rowSums(seeds[5:30])
sum(row.sum) # 123 seeds
# 31 species columns - 5 non-species columns = 26 morphotypes

# Bring in the filtered data
seeds.filt <- read.csv(file = "data/Experiment-2-clean-seeds-filt.csv",
                  header = TRUE)
head(seeds.filt)

# Calculate the total number of seeds and morphotypes
row.sum <- rowSums(seeds.filt[5:25])
sum(row.sum) # 68 seeds
# 68 seeds are wind-dispersed, canopy, or epizoochorous

# Calculate a different way (check the method)
sum(seeds[,5:30]) # 123
sum(seeds.filt[5:25]) # 68

round((68/123) * 100)
# 55% of seeds in filtered data set 
# 45% of seeds are filtered out (epizoo, canopy, wind)

# Calculate the totals for each species
colSums(seeds[,5:30])

# Rhus + Smilax
(12 + 5) / 123 * 100
