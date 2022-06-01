## --------------- HEADER ------------------------------------------------------
## Script name: script4.0-seed-time-series.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliaton: University of Florida
## Date Created: 2021-05-03
## Date Last modified: 2022-05-03
## Copyright (c) 2022, David S. Mason, Tyler S. Coleman, James P. Holdgrafer
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script analyzes the seed trap data from the bird
## feeder project as a time series

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

# fall burning ,frequent burning, herbivore
# go in andnkill with herbicide, should do a good job controlling it.
# absolutely correct, you can burn very frequerntly in the traditional 
# dormant nd early growing season and not reeduce abundance in perpetuity.

# historically part of community. not aware of it being more or less dominant. check bartram.

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
m4 <- manyglm(comm.mat.mv ~ predictors$TREATMENT + (1|predictors$BLOCK/predictors$DATE.ORD), 
              family = "negative_binomial")

plot(m4)
anova(m4, p.uni = "adjusted")

## --------------- RICHNESS ----------------------------------------------------

# Convert species matrix to presence absence
library(vegan)
comm.mat.pa <- decostand(comm.mat, "pa")

# Bring back the predictor variables
seeds <- cbind(predictors[,1:4], comm.mat.pa)

# Calculate richness 
seeds$RICH <- NA
for(i in 1:nrow(seeds)){
  seeds$RICH[i] <- sum(seeds[i,5:35])
}  

# Convert treatment to an ordinal variable
seeds$TREATMENT.NUMB <- NA
for(i in 1:nrow(seeds)){
  if(seeds$TREATMENT[i] == "Control"){
    seeds$TREATMENT.NUMB[i] <- 0
  }
  if(seeds$TREATMENT[i] == "Low"){
    seeds$TREATMENT.NUMB[i] <- 4
  }
  if(seeds$TREATMENT[i] == "Medium"){
    seeds$TREATMENT.NUMB[i] <- 8
  }
  if(seeds$TREATMENT[i] == "High"){
    seeds$TREATMENT.NUMB[i] <- 12
  }
}

is.numeric(seeds$TREATMENT.NUMB)

# Convert date to a dummy variable
seeds$DATE.ORD <- NA
for(i in 1:nrow(seeds)){
  if(seeds$DATE[i] == "2020-11-27"){
    seeds$DATE.ORD[i] <- 1
  }
  if(seeds$DATE[i] == "2020-12-04"){
    seeds$DATE.ORD[i] <- 2
  }
  if(seeds$DATE[i] == "2020-12-11"){
    seeds$DATE.ORD[i] <- 3
  }
  if(seeds$DATE[i] == "2021-01-02"){
    seeds$DATE.ORD[i] <- 4
  }
  if(seeds$DATE[i] == "2021-02-03"){
    seeds$DATE.ORD[i] <- 5
  }
}

seeds$DATE.ORD <- factor(seeds$DATE.ORD, order = TRUE,
                     levels = c(1,2,3,4,5))

# Plot distribution
plotdist(seeds$RICH, histo = TRUE, demp = TRUE)
descdist(seeds$RICH, discrete=TRUE, boot=500) # NB

# Calculate 
m6 <- glmmTMB(RICH ~ TREATMENT.NUMB + (1|BLOCK/DATE.ORD),
              data = seeds, family = nbinom2(link = "log"))
sim_m6 <- simulateResiduals(fittedModel = m6, n = 250)
plot(sim_m6)

car::Anova(m6)

# Calculate summary stats

seeds %>% 
  group_by(TREATMENT) %>% 
  summarize(RICHNESS = mean(RICH),
            SD = sd(RICH))

emtrends(m6, ~ TREATMENT.NUMB, var = "TREATMENT.NUMB", transform = "response")

emmeans(m6, ~TREATMENT,
        transform = "response"
)

## --------------- TOTAL DETECTIONS --------------------------------------------

# Reconstitute the seed dataframe
seeds <- cbind(predictors[,1:4], comm.mat)
seeds$DETECTIONS <- NA

# Calculate detections 
for(i in 1:nrow(seeds)){
  seeds$DETECTIONS[i] <- sum(seeds[i,5:35])
}  

# Convert treatment to an ordinal variable
seeds$TREATMENT.NUMB <- NA
for(i in 1:nrow(seeds)){
  if(seeds$TREATMENT[i] == "Control"){
    seeds$TREATMENT.NUMB[i] <- 0
  }
  if(seeds$TREATMENT[i] == "Low"){
    seeds$TREATMENT.NUMB[i] <- 4
  }
  if(seeds$TREATMENT[i] == "Medium"){
    seeds$TREATMENT.NUMB[i] <- 8
  }
  if(seeds$TREATMENT[i] == "High"){
    seeds$TREATMENT.NUMB[i] <- 12
  }
}

# Convert date to a dummy variable
seeds$DATE.ORD <- NA
for(i in 1:nrow(seeds)){
  if(seeds$DATE[i] == "2020-11-27"){
    seeds$DATE.ORD[i] <- 1
  }
  if(seeds$DATE[i] == "2020-12-04"){
    seeds$DATE.ORD[i] <- 2
  }
  if(seeds$DATE[i] == "2020-12-11"){
    seeds$DATE.ORD[i] <- 3
  }
  if(seeds$DATE[i] == "2021-01-02"){
    seeds$DATE.ORD[i] <- 4
  }
  if(seeds$DATE[i] == "2021-02-03"){
    seeds$DATE.ORD[i] <- 5
  }
}

seeds$DATE.ORD <- factor(seeds$DATE.ORD, order = TRUE,
                         levels = c(1,2,3,4,5))

# Plot distribution
plotdist(seeds$DETECTIONS, histo = TRUE, demp = TRUE)
descdist(seeds$DETECTIONS, discrete=TRUE, boot=500) # NB

# Calculate 
m7 <- glmmTMB(DETECTIONS ~ TREATMENT + (1|BLOCK/DATE.ORD),
              data = seeds, family = nbinom2(link = "log"))
sim_m7 <- simulateResiduals(fittedModel = m7, n = 250)
plot(sim_m7)

car::Anova(m7)

seeds %>% 
  group_by(TREATMENT) %>% 
  summarize(DETECT = mean(DETECTIONS),
            SD = sd(DETECTIONS))

emtrends(m7, ~ TREATMENT.NUMB, var = "TREATMENT.NUMB", transform = "response")

emmeans(m7, ~TREATMENT,
        transform = "response"
)

## --------------- PRELIMINARY DATA --------------------------------------------


# prelim.seed.dat are from the baited vs un-baited feeders effect on seed count
prelim <- read.csv("data/initial_feeder.csv",
                                header = TRUE)

# Plot distribution
plotdist(prelim$rawseeds, histo = TRUE, demp = TRUE)
descdist(prelim$rawseeds, discrete=TRUE, boot=500) # NB/poisson

m8 <- glmmTMB(rawseeds ~ Treatment + (1|time/pair),
              data = prelim, family = nbinom2(link = "log"))

m8 <- glmer(rawseeds ~ Treatment + (1|pair),
            data = prelim, family = "poisson")

sim_m8 <- simulateResiduals(fittedModel = m8, n = 250)
plot(sim_m8)

car::Anova(m8)

emmeans(m8, ~Treatment,
        transform = "response"
)
