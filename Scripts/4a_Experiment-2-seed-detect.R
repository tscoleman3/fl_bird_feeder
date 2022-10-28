## --------------- HEADER ------------------------------------------------------
## Script name: 4a_Experiment-2-seed-detect.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliaton: University of Florida
## Date Created: 2022-10-21
## Date Last Modified: 2022-10-21
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script analyses the seed detection data by
## sampling date in experiment 2.

## --------------- SET—UP WORKSPACE --------------------------------------------

# Load packages
library(tidyverse)
library(lme4)
library(lubridate)
library(lmerTest)
library(emmeans)
library(car)
library(MuMIn)
library(agridat)
library(DHARMa)
library(glmmTMB)
library(fitdistrplus)
library(mvabund)

# Clear the decks
rm(list=ls())

# Bring in the data
seeds <- read.csv(file = "data/seed_traps2.0.csv",
                  header = TRUE, stringsAsFactors = FALSE)

## --------------- MANYGLM -----------------------------------------------------

# Separate the community data from the predictor variables
comm.mat <- seeds[,5:35]
comm.mat.mv <- mvabund(seeds[,5:35])
predictors <- seeds[,1:4]

# Run simple models
m1 <- manyglm(comm.mat.mv ~ predictors$TREATMENT,
              family = "negative_binomial")

m2 <- manyglm(comm.mat.mv ~ predictors$TREATMENT + (1|predictors$BLOCK),
              family = "negative_binomial")

# Convert date to a dummy variable
predictors$DATE.ORD <- NA
for(i in 1:nrow(predictors)){
  if(predictors$DATE[i] == "2020-11-27"){
    predictors$DATE.ORD[i] <- 1
  }
  if(predictors$DATE[i] == "2020-12-04"){
    predictors$DATE.ORD[i] <- 2
  }
  if(predictors$DATE[i] == "2020-12-11"){
    predictors$DATE.ORD[i] <- 3
  }
  if(predictors$DATE[i] == "2021-01-02"){
    predictors$DATE.ORD[i] <- 4
  }
  if(predictors$DATE[i] == "2021-02-03"){
    predictors$DATE.ORD[i] <- 5
  }
} 

# Run nested model
m4 <- manyglm(comm.mat.mv ~ predictors$TREATMENT + predictors$DATE.ORD + (1|predictors$BLOCK), 
              family = "negative_binomial")

plot(m4)
anova.manyglm(m4, p.uni = "adjusted")



seeds %>% 
  dplyr::select(-NOTES) %>% 
  pivot_longer(5:35, names_to = "SPECIES", values_to = "DETECTIONS") %>% 
  filter(SPECIES == "RHUS") %>% 
  group_by(DATE) %>% 
  summarize(DETECT = mean(DETECTIONS),
            SD = sd(DETECTIONS),
            N = n(),
            SE = SD/sqrt(N))
