## --------------- HEADER ------------------------------------------------------
## Script name: 1a_Experiment-1-model.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2022-10-15
## Date Last Modified: 2023-01-26
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script analyzes the seed collection data for 
## experiment 1

## Overview: below you will find the model construction process, which ends with
## the final model that is used in the manuscript.
## 1] We start by exploring the distribution of the raw data. 
## 2] Based on this, we run a poisson model. This poisson model was overdispersed. 
## 3] We conducted a negative binomial model, which was heteroscedastic
## 4] We tried a conservative model without the "extreme" values associated with
##     with baited feeders at block 4 during the second sample. The model was
##     still heteroscedastic.

## --------------- Set up workspace --------------------------------------------
library(lattice)
library(fitdistrplus)
library(lme4)
library(tidyverse)
library(car)
library(DHARMa)
library(MuMIn)
library(emmeans)
library(lme4)
library(performance)

# Clear the deck
rm(list = ls())

# Bring in the data
d <- read.csv("Data/Experiment-1.csv")

# Convert treatment to factor
d$Treatment <- as_factor(d$Treatment)

## --------------- Check distribution of data ----------------------------------

# Visualize the data 
boxplot(rawseeds ~ Treatment,col=c("white","lightgray"), outline = TRUE, d)
dotplot(d$rawseeds) # huge outlier present

plotdist(d$rawseeds, histo = TRUE, demp = TRUE)
descdist(d$rawseeds, discrete=TRUE, boot=500)

# Make sure pair and time is a factor
d$pair <- as_factor(d$pair)
d$time <- as_factor(d$time)

# Check feeder treatment level
d |>
  filter(Treatment == 'Baited') |>
  dplyr::select(rawseeds) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # Negative binomial

d |>
  filter(Treatment == 'Control') |>
  dplyr::select(rawseeds) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500) # Negative binomial

# Time periods
d |>
  filter(time == '1') |>
  dplyr::select(rawseeds) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500) # poisson

d |>
  filter(time == '2') |>
  dplyr::select(rawseeds) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500) # poisson

## --------------- Full model --------------------------------------------------

seeds.fit.pois <- glmer(rawseeds ~ Treatment + time + (1 | pair),
                        data = d, family = poisson(link = "log"))
summary(seeds.fit.pois)
Anova(seeds.fit.pois)

# Check overdispersion
E1 <- resid(seeds.fit.pois, type = "pearson")
N <- nrow(d)
p <- length(fixef(seeds.fit.pois)) + 1
overdispersion <- sum(E1^2) / (N-p)
overdispersion 

# Overdispersion present, let's try negative binomial
seeds.fit.nb <- glmer.nb(rawseeds ~ Treatment + time + (1 | pair),
                         data = d)
summary(seeds.fit.nb)
Anova(seeds.fit.nb, type = 'II')

# Check overdispersion
E1 <- resid(seeds.fit.nb, type = "pearson")
N <- nrow(d)
p <- length(fixef(seeds.fit.nb)) + 1
overdispersion <- sum(E1^2) / (N-p)
overdispersion 

# check overall assumptions
dev.new()
check_model(seeds.fit.nb)

sim.seeds.fit.nb <- simulateResiduals(seeds.fit.nb)
dev.new()
plot(sim.seeds.fit.nb)

# Check values and residuals by treatment
plot(rawseeds ~ as_factor(Treatment), data = d)
plot(residuals(seeds.fit.nb) ~ as_factor(d$Treatment))
abline(a = 0, b = 0, col = "blue", lwd = 2)
hist(residuals(seeds.fit.nb))

# Plot residuals
resfit <- resid(seeds.fit.nb)
hist(resfit)
plot(d$rawseeds, resfit, 
     ylab = "Residuals", 
     xlab = "Detections") 
abline(0, 0)     

# Check predictions 
preds.nb <- predict(seeds.fit.nb)
par(mfrow = c(1, 2))
plot(rawseeds ~ as_factor(Treatment), 
     data = d)
plot(preds.nb ~ as_factor(d$Treatment))

predict(seeds.fit.nb)
d$pred = exp(predict(seeds.fit.nb))
d$predicted = predict(seeds.fit.nb)    # save the predicted values
d$residuals = residuals(seeds.fit.nb)  # save the residual values

# Quick look at the actual, predicted, and residual values
pred_df <- d %>% 
  dplyr::select(rawseeds, predicted, residuals)
pred_df$predicted = exp(pred_df$predicted)
pred_df$residuals = exp(pred_df$residuals)

# Plot fitted and residual values
plot(fitted(seeds.fit.nb) ~ d$rawseeds)
abline(0, 1, col = "blue", lwd = 2)

# Plot predicted and residual values
plot(predict(seeds.fit.nb), residuals(seeds.fit.nb, type = 'working'))

# Plot cooks distance
plot(cooks.distance(seeds.fit.nb), type='h')

# Plot predicted vs. actual
plot(density(d$rawseeds), xlim=c(0, 160), ylim=c(0, .08), main='M_0 y_hat')
lines(density(predict(seeds.fit.nb, type='response')), col='red')

# Check independence of observations
scatter.smooth(d$time, d$residuals) # Greater variance in residuals second sample

## --------------- Export coefficient estimates --------------------------------

# EMMEANS
emmeans(seeds.fit.nb, pairwise ~Treatment, type = 'response')
confint(emmeans(seeds.fit.nb, pairwise ~Treatment, type = 'response'))

means <- as.data.frame(emmeans(seeds.fit.nb, pairwise ~Treatment, 
                               type = 'response')[1], level = 0.95)
confint <- as.data.frame(confint(emmeans(seeds.fit.nb, 
                    pairwise ~Treatment, type = 'response')))

# Calculate r2
r.squaredGLMM(seeds.fit.nb)
# Marginal - conditional = 
# 0.7406428-0.5717837 = 0.1688591
# 0.8332853-0.6433047 = 0.1899806
# 0.7958567-0.6144094 = 0.1814473

## --------------- Export model validation -------------------------------------
E1 <- resid(seeds.fit.nb, type = 'pearson')
F1 <- fitted(seeds.fit.nb, type = 'response')

# png("Figures/Experiment-1-model-validation-full-data.png", 
    # width = 800, height = 800, units = "px")
par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x = F1, y = E1, xlab = "Fitted values", ylab = "Pearson residuals",
     cex.lab = 2)
abline(h = 0, lty = 2)
plot(d$time, y = E1, xlab = "Sampling period", ylab = "Pearson residuals",
     cex.lab = 2)
abline(h = 0, lty = 2)
boxplot(E1 ~ Treatment, data = d, xlab = "Treatment", ylab = "Pearson residuals",
        cex.lab = 2)
abline(h = 0, lty = 2)
boxplot(E1~pair, data = d, xlab = "Pair", ylab = "Pearson residuals",
        cex.lab = 2)
abline(h = 0, lty = 2)
dev.off()

# Below, we run the models again without the high value at
# block 4 during the second sampling period

## --------------- Conservative check distribution of data ----------------------------------

d.cons <- d[-c(17,18),]

# Visualize the data ##
boxplot(rawseeds ~ Treatment,col=c("white","lightgray"), outline = TRUE, d.cons)
dotplot(d.cons$rawseeds) # outlier in control

plotdist(d.cons$rawseeds, histo = TRUE, demp = TRUE)
descdist(d.cons$rawseeds, discrete=TRUE, boot=500) 

# Make sure pair and time is a factor
d.cons$pair <- as_factor(d.cons$pair)
d.cons$time <- as_factor(d.cons$time)

# Feeder treatment distributions
d.cons |>
  filter(Treatment == 'Baited') |>
  dplyr::select(rawseeds) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # Poisson

d.cons |>
  filter(Treatment == 'Control') |>
  dplyr::select(rawseeds) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500) # Poisson

# Sampling period distributions
d.cons |>
  filter(time == '1') |>
  dplyr::select(rawseeds) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # poisson

d.cons |>
  filter(time == '2') |>
  dplyr::select(rawseeds) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # normal?

# Distribution of conservative data (without the outlier) is a bit different than
# than the full data. Let's see if the negative binomial GLM can still produce a
# valid model for consistency.

## --------------- Conservative negative binomial GLM --------------------------

cons.seeds.fit.nb <- glmer.nb(rawseeds ~ Treatment + time + (1 | pair),
                         data = d.cons)
summary(cons.seeds.fit.nb)
Anova(cons.seeds.fit.nb, type = 'II')

# Check overall assumptions
dev.new()
check_model(cons.seeds.fit.nb) # good
sim.seeds.fit.nb <- simulateResiduals(cons.seeds.fit.nb)
plot(sim.seeds.fit.nb) # good

# Check overdispersion
E1 <- resid(cons.seeds.fit.nb, type = "pearson")
N <- nrow(d)
p <- length(fixef(cons.seeds.fit.nb)) + 1
overdispersion <- sum(E1^2) / (N-p)
overdispersion # Better

# Check values and residuals by treatment
plot(rawseeds ~ as_factor(Treatment), data = d.cons)
plot(residuals(cons.seeds.fit.nb) ~ as_factor(d.cons$Treatment))
abline(a = 0, b = 0, col = "blue", lwd = 2)
hist(residuals(cons.seeds.fit.nb))

# Plot residuals
resfit <- resid(cons.seeds.fit.nb)
hist(resfit)
plot(d.cons$rawseeds, resfit, 
     ylab = "Residuals", 
     xlab = "Detections") 
abline(0, 0)     

# Plot predictions 
preds.nb <- predict(cons.seeds.fit.nb)
par(mfrow = c(1, 2))
plot(rawseeds ~ as_factor(Treatment), 
     data = d.cons)
plot(preds.nb ~ as_factor(d.cons$Treatment))

predict(cons.seeds.fit.nb)
d.cons$pred = exp(predict(cons.seeds.fit.nb))
d.cons$predicted = predict(cons.seeds.fit.nb)    # save the predicted values
d.cons$residuals = residuals(cons.seeds.fit.nb)  # save the residual values

# Quick look at the actual, predicted, and residual values
pred_df <- d.cons %>% 
  dplyr::select(rawseeds, predicted, residuals)
pred_df$predicted = exp(pred_df$predicted)
pred_df$residuals = exp(pred_df$residuals)

plot(fitted(cons.seeds.fit.nb) ~ d.cons$rawseeds)
abline(0, 1, col = "blue", lwd = 2)

# Plot predictions against residuals
plot(predict(cons.seeds.fit.nb), residuals(cons.seeds.fit.nb, type = 'working'))

# Plot cooks distance
plot(cooks.distance(cons.seeds.fit.nb), type='h')

# Compare predicted vs. actual
plot(density(d.cons$rawseeds), xlim=c(0, 160), ylim=c(0, .08), main='M_0 y_hat')
lines(density(predict(cons.seeds.fit.nb, type='response')), col='red')

# Check independence of observations
scatter.smooth(d.cons$time, d.cons$residuals)

# Compare sampling events
boxplot(E1~time, data = d.cons, xlab = "Sample", ylab = "Pearson residuals",
        cex.lab = 2)
abline(h = 0, lty = 2)


## --------------- Conservative Export coefficient estimates --------------------------------

# EMMEANS
emmeans(cons.seeds.fit.nb, pairwise ~Treatment, type = 'response')
confint(emmeans(cons.seeds.fit.nb, pairwise ~Treatment, type = 'response'))

cons.means <- as.data.frame(emmeans(cons.seeds.fit.nb, pairwise ~Treatment, 
                                    type = 'response')[1], level = 0.95)

# Calculate r2
r.squaredGLMM(cons.seeds.fit.nb)
# 0.7341450-0.5880577 = 0.1460873
# 0.7790336-0.6240140 = 0.1550196
# 0.6705627-0.5371277 = 0.133435

# Export means
write.csv(cons.means, "Model-output/Experiment-1-emmeans.csv", row.names = FALSE)

## --------------- Conservative Export model validation -------------------------------------
E1 <- resid(cons.seeds.fit.nb, type = 'pearson')
F1 <- fitted(cons.seeds.fit.nb, type = 'response')

png("Figures/Experiment-1-model-validation-conserv.png", 
    width = 800, height = 800, units = "px")
par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x = F1, y = E1, xlab = "Fitted values", ylab = "Pearson residuals",
     cex.lab = 2)
abline(h = 0, lty = 2)
plot(d.cons$time, y = E1, xlab = "Sampling period", ylab = "Pearson residuals",
     cex.lab = 2)
abline(h = 0, lty = 2)
boxplot(E1 ~ Treatment, data = d.cons, xlab = "Treatment", ylab = "Pearson residuals",
        cex.lab = 2)
abline(h = 0, lty = 2)
boxplot(E1~pair, data = d.cons, xlab = "Pair", ylab = "Pearson residuals",
        cex.lab = 2)
abline(h = 0, lty = 2)
dev.off()
