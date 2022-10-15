## --------------- HEADER ------------------------------------------------------
## Script name: 1_Experiment-1-model.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliaton: University of Florida
## Date Created: 2022-10-15
## Date Last Modified: 2022-10-15
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script analyzes the seed collection data for 
## experiment 1

## --------------- Set up workspace --------------------------------------------

library(fitdistrplus)
library(lme4)
library(tidyverse)
library(car)
library(DHARMa)
library(lattice)
library(MuMIn)
library(emmeans)

# Clear the deck
rm(list = ls())

d <- read.csv("Data/initial_feeder.csv")

## --------------- Check distribution of data ----------------------------------

# Visualize the data ####
boxplot(rawseeds ~ Treatment,col=c("white","lightgray"), outline = TRUE, d)
dotplot(d$rawseeds) # huge outlier present

# Remove outlier
# d <- d[-c(17,18),]

plotdist(d$rawseeds, histo = TRUE, demp = TRUE)
descdist(d$rawseeds, discrete=TRUE, boot=500) # Poisson

## --------------- Poisson GLM -------------------------------------------------

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

# Overdispersion present

## --------------- Negative binomial GLM ---------------------------------------

seeds.fit.nb <- glmer.nb(rawseeds ~ Treatment + time + (1 | pair),
                         data = d)

summary(seeds.fit.nb)
Anova(seeds.fit.nb, type = 'II')

# Check overdispersion
E1 <- resid(seeds.fit.nb, type = "pearson")
N <- nrow(d)
p <- length(fixef(seeds.fit.nb)) + 1
overdispersion <- sum(E1^2) / (N-p)
overdispersion # Better

# check some assumptions
plot(rawseeds ~ Treatment, data = d)
plot(residuals(seeds.fit.nb) ~ d$Treatment)
abline(a = 0, b = 0, col = "blue", lwd = 2)
hist(residuals(seeds.fit.nb))

# plot residuals from prediction
resfit <- resid(seeds.fit.nb)
hist(resfit)
plot(d$rawseeds, resfit, 
     ylab = "Residuals", 
     xlab = "Detections") 
abline(0, 0)     

# look at predictions 
preds.nb <- predict(seeds.fit.nb)
par(mfrow = c(1, 2))
plot(rawseeds ~ Treatment, 
     data = d)
plot(preds.nb ~ d$Treatment)

# look at predictions from model
# neg binom
predict(seeds.fit.nb)
d$pred = exp(predict(seeds.fit.nb))
d$predicted = predict(seeds.fit.nb)    # save the predicted values
d$residuals = residuals(seeds.fit.nb)  # save the residual values

# quick look at the actual, predicted, and residual values
pred_df <- d %>% 
  dplyr::select(rawseeds, predicted, residuals)
pred_df$predicted = exp(pred_df$predicted)
pred_df$residuals = exp(pred_df$residuals)

plot(fitted(seeds.fit.nb) ~ d$rawseeds)
abline(0, 1, col = "blue", lwd = 2)

# Using dharma
sim.seeds.fit.nb <- simulateResiduals(seeds.fit.nb)
plot(sim.seeds.fit.nb)

## --------------- Bootstrap parameters ----------------------------------------

# look at our coefs 
coefs <- summary(seeds.fit.nb)$coef
coefs_est <- exp(coefs[, "Estimate"])
uprs <- exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"])
lwrs <- exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"])
(uprs - 1) * 100       # upr CI %'s
(coefs_est - 1) * 100  # coefficient estimates CI %'s
(lwrs - 1) * 100       # lwr CI %'s

par(mfrow = c(1, 1))
plot(exp(coefs[1:2,"Estimate"]), ylim = range(lwrs[1:2], uprs[1:2]))
segments(1:7, lwrs[1:2], 1:7, uprs[1:2])

coefs.fix <- fixef(seeds.fit.nb)
exp(coefs.fix[1])
exp(coefs.fix[2])

coefs_func <- function(.) {
  beta <- unname(fixef(.))
  treat <- exp(beta[1])             # mean count for control 
  control <-  exp(beta[1] + beta[2])   # mean count for treat 1 is this much greater than control
  c(control_mu = control, treat = treat)
  
}

start <- Sys.time()
rand <- bootMer(x = seeds.fit.nb, FUN = coefs_func, nsim = 1000)$t
Sys.time() - start

summ <- apply(rand, 2, function(x) c(mean = mean(x, na.rm = TRUE),
                                     quantile(x, c(0.025, 0.975), na.rm = TRUE)))

# EMMEANS
emmeans(seeds.fit.nb, pairwise ~Treatment, type = 'response')

means <- as.data.frame(emmeans(seeds.fit.nb, pairwise ~Treatment, type = 'response')[1])


# 

# Calculate r2
r.squaredGLMM(seeds.fit.nb)

## --------------- Export coefficient estimates --------------------------------

write.csv(summ, "Model-output/initial-feeder-boot-coeff.csv", row.names = FALSE)
write.csv(means, "Model-output/initial-feeder-emmeans.csv", row.names = FALSE)

## --------------- Export model validation -------------------------------------
E1 <- resid(seeds.fit.nb, type = 'pearson')
F1 <- fitted(seeds.fit.nb, type = 'response')

png("Figures/Experiment-1-model-validation.png", 
    width = 800, height = 800, units = "px")
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
