## --------------- HEADER ------------------------------------------------------
## Script name: 2a_Experiment-2-bird-mag.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliaton: University of Florida
## Date Created: 2022-10-20
## Date Last Modified: 2022-10-20
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script analyses the bird detections data
## from experiment 2

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

## --------------- Acclimitization period --------------------------------------

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

## --------------- Model total detections --------------------------------------

plotdist(bird.obs$count, histo = TRUE, demp = TRUE)
descdist(bird.obs$count, discrete=TRUE, boot=500) # NB or poisson

birds.obs.fit.pois <- glmer(count ~ treatment + (1 | site),
                            data = bird.obs, family = poisson(link = "log"))
sim.birds.obs.fit.pois <- simulateResiduals(birds.obs.fit.pois)
plot(sim.birds.obs.fit.pois) # had to increase err??
Anova(birds.obs.fit.pois)
summary(birds.obs.fit.pois)

birds.obs.fit.nb <- glmer.nb(count ~ treatment + (1 | site),
                            data = bird.obs)
sim.birds.obs.fit.nb <- simulateResiduals(birds.obs.fit.nb)
plot(sim.birds.obs.fit.nb) # had to increase err??
Anova(birds.obs.fit.nb)
summary(birds.obs.fit.nb)

library('performance')
check_zeroinflation(birds.obs.fit.nb)
check_zeroinflation(birds.obs.fit.pois)

birds.obs.fit.zinb1 <- glmmTMB(count~treatment + (1|site),
                      family=nbinom1(), 
                      zi=~treatment, 
                      data=bird.obs)
summary(birds.obs.fit.zinb1)

birds.obs.fit.zinb2 <- glmmTMB(count~treatment + (1|site),
                               family=nbinom2(), 
                               zi=~treatment, 
                               data=bird.obs)
summary(birds.obs.fit.zinb2)
Anova(birds.obs.fit.zinb2)

birds.obs.fit.zinb2.sim <- simulateResiduals(birds.obs.fit.zinb2)
plot(birds.obs.fit.zinb2.sim) # not converging

check_zeroinflation(birds.obs.fit.zinb1)
check_zeroinflation(birds.obs.fit.zinb2)

# Selecting sim.birds.obs.fit.nb

## --------------- Check assumptions -------------------------------------------

# check some assumptions
plot(count ~ treatment, data = bird.obs)
plot(residuals(birds.obs.fit.nb, type = 'pearson') ~ bird.obs$treatment)
abline(a = 0, b = 0, col = "blue", lwd = 2)
hist(residuals(birds.obs.fit.nb, type = 'pearson'))

# plot residuals from prediction
resfit <- resid(birds.obs.fit.nb, type = 'pearson')
hist(resfit)
plot(bird.obs$treatment, resfit, 
     ylab = "Residuals", 
     xlab = "Treatment") 
abline(0, 0)  

# look at predictions 
preds.lm <- predict(birds.obs.fit.nb)
par(mfrow = c(1, 2))
plot(count ~ treatment, 
     data = bird.obs)
plot(preds.lm ~ bird.obs$treatment) 

# look at predictions from model
# neg binom
predict(birds.obs.fit.nb)
bird.obs$pred = predict(birds.obs.fit.nb)
bird.obs$predicted = predict(birds.obs.fit.nb)    # save the predicted values
bird.obs$residuals = residuals(birds.obs.fit.nb, type = 'pearson')  # save the residual values

# quick look at the actual, predicted, and residual values
pred_df <- bird.obs %>% 
  dplyr::select(count, predicted, residuals)

plot(fitted(birds.obs.fit.nb) ~ bird.obs$count)
abline(0, 1, col = "blue", lwd = 2)

## --------------- Export model validation -------------------------------------
E1 <- resid(birds.obs.fit.nb, type = 'pearson')
F1 <- fitted(birds.obs.fit.nb)

png("Figures/Experiment-2-bird-detect-model-validation.png", 
    width = 800, height = 800, units = "px")
par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x = F1, y = E1, xlab = "Fitted values", ylab = "Normalized residuals",
     cex.lab = 2)
abline(h = 0, lty = 2)
boxplot(E1 ~ treatment, data = bird.obs, xlab = "Treatment", ylab = "Normalized residuals",
        cex.lab = 2)
abline(h = 0, lty = 2)
boxplot(E1~site, data = bird.obs, xlab = "Block", ylab = "Normalized residuals",
        cex.lab = 2)
abline(h = 0, lty = 2)
dev.off()

## --------------- Calculate and export means ----------------------------------

Anova(birds.obs.fit.nb)

emmeans(birds.obs.fit.nb, pairwise ~treatment, type = 'response',
        adjust = "none")

means <- as.data.frame(emmeans(birds.obs.fit.nb, pairwise ~treatment, type = 'response')[1])
write.csv(means, "Model-output/experiment-2-bird-detect-emmeans.csv", row.names = FALSE)


# 

# Calculate r2
r.squaredGLMM(birds.obs.fit.nb)
