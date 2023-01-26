## --------------- HEADER ------------------------------------------------------
## Script name: 3a_Experiment-2-bird-rich.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliaton: University of Florida
## Date Created: 2022-10-15
## Date Last Modified: 2022-10-15
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script analyzes the camera trap data for birds
## in experiment 2.

## --------------- Set up workspace --------------------------------------------
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

rm(list=ls())

bird.dat.og <- read.csv(file = "data/bird_data2.0.csv",
                        header = TRUE)
head(bird.dat.og)

## --------------- Drop observations before acclimization ----------------------

bird.dat <- dplyr::select(bird.dat.og, -c("camera.card", "X", "notes", "photo.number"))
bird.dat$date <- as.Date(with(bird.dat, paste(year, month, day, sep = "-")), "%Y-%m-%d")

bird.week <- bird.dat %>% 
  group_by(week = week(date)) %>% 
  summarise(value = sum(presence))

bird.week
# looks like a 2 week acclimitizaton period
# observatons really picked up in January (observations are
# essentially increasing with time throughout our sampling period)

bird.dat.drop <- bird.dat %>% 
  mutate(week = week(date)) %>% 
  filter(week != 48 & week != 47)

# remove NA's
bird.clean <- bird.dat.drop %>% 
  drop_na(species)

bird.clean$treatment <- as.factor(bird.clean$treatment)

## --------------- Summarize the data ------------------------------------------

# count of individuals by species at each treatment and site
bird.sps.count.by.treat <- bird.clean %>% 
  group_by(site, treatment, species) %>% 
  summarise(count = n())

# total count of individuals by treatment and site
bird.obs <- bird.clean %>% 
  group_by(site, treatment) %>% 
  summarise(count = n())

bird.obs <- as.data.frame(bird.obs) %>% 
  dplyr::select(treatment,site,count)

missing <- data.frame(site = "six", 
                      treatment = "4",
                      count = 0)

bird.obs <- rbind(bird.obs,missing)

## --------------- Create the richness dataframe -------------------------------

sites <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,
           8,8,8,8,9,9,9,9,10,10,10,10)
treatment <- c(0,4,8,12,0,4,8,12,0,4,8,12,0,4,8,12,
               0,4,8,12,0,4,8,12,0,4,8,12,0,4,8,12,
               0,4,8,12,0,4,8,12)

bird.rich <- as.data.frame(cbind(sites,treatment))
bird.rich$richness <- c(2,2,3,3,2,2,2,3,3,3,4,4,3,3,3,3,4,4,5,5,
                        2,0,2,1,1,2,2,1,4,4,4,4,2,2,2,2,1,1,1,1)
bird.rich$treatment <- as.factor(bird.rich$treatment)
bird.rich <- as.data.frame(bird.rich) %>% 
  dplyr::select(treatment,sites,richness)

## --------------- Conduct the model -------------------------------------------

plotdist(bird.rich$richness, histo = TRUE, demp = TRUE)
descdist(bird.rich$richness, discrete=TRUE, boot=500) # poisson

# Convert to factor
bird.rich$treatment <- as_factor(bird.rich$treatment)
bird.rich$sites <- as_factor(bird.rich$sites)

bird.rich.lm <- lmer(richness ~ treatment + (1|sites),
                        data = bird.rich) # Normal fufilled assumptions best
summary(bird.rich.lm)
Anova(bird.rich.lm)

## --------------- Check assumptions -------------------------------------------

sim.bird.rich.lm <- simulateResiduals(bird.rich.lm)
plot(sim.bird.rich.lm) # increasing 'err'

# check some assumptions
plot(richness ~ treatment, data = bird.rich)
plot(residuals(bird.rich.lm) ~ bird.rich$treatment)
abline(a = 0, b = 0, col = "blue", lwd = 2)
hist(residuals(bird.rich.lm))

# plot residuals from prediction
resfit <- resid(bird.rich.lm)
hist(resfit)
plot(bird.rich$treatment, resfit, 
     ylab = "Residuals", 
     xlab = "Treatment") 
abline(0, 0)  

# look at predictions 
preds.lm <- predict(bird.rich.lm)
par(mfrow = c(1, 2))
plot(richness ~ treatment, 
     data = bird.rich)
plot(preds.lm ~ bird.rich$treatment) 

# look at predictions from model
# neg binom
predict(bird.rich.lm)
bird.rich$pred = predict(bird.rich.lm)
bird.rich$predicted = predict(bird.rich.lm)    # save the predicted values
bird.rich$residuals = residuals(bird.rich.lm)  # save the residual values

# quick look at the actual, predicted, and residual values
pred_df <- bird.rich %>% 
  dplyr::select(richness, predicted, residuals)

plot(fitted(bird.rich.lm) ~ bird.rich$richness)
abline(0, 1, col = "blue", lwd = 2)

## --------------- Export model validation -------------------------------------
E1 <- resid(bird.rich.lm)
F1 <- fitted(bird.rich.lm)

png("Figures/Experiment-2-bird-rich-model-validation.png", 
    width = 800, height = 800, units = "px")
par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x = F1, y = E1, xlab = "Fitted values", ylab = "Normalized residuals",
     cex.lab = 2)
abline(h = 0, lty = 2)
boxplot(E1 ~ treatment, data = bird.rich, xlab = "Treatment", ylab = "Normalized residuals",
        cex.lab = 2)
abline(h = 0, lty = 2)
boxplot(E1~sites, data = bird.rich, xlab = "Block", ylab = "Normalized residuals",
        cex.lab = 2)
abline(h = 0, lty = 2)
dev.off()

## --------------- Calculate and export means ----------------------------------

Anova(bird.rich.lm)

emmeans(bird.rich.lm, pairwise ~treatment, type = 'response',
        adjust = "none")

means <- as.data.frame(emmeans(bird.rich.lm, pairwise ~treatment, type = 'response')[1])
write.csv(means, "Model-output/experiment-2-bird-rich-emmeans.csv", row.names = FALSE)


# 

# Calculate r2
r.squaredGLMM(bird.rich.lm)

