## --------------- HEADER ------------------------------------------------------
## Script name: 2b_Experiment-2-bird-count.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2022-10-20
## Date Last Modified: 2023-02-12
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script analyses the bird detections data
## from experiment 2

## In this script we model the distribution of the data and attempt
## models, beginning with the most complicated (block as a random factor)
## and interaction with time as a continuous variable). If this model
## performs poorly with reasonable distributions, we try less complicated
## models. Eventually we settle on a negative binomial model including block
## block (i.e., site) as a random factor and weeks since acclimatization
## as a factor with no interaction.

## After producing the model we check assumptions and export means and
## confidence intervals to use in the figure.

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
library(mvabund)
library(performance)

rm(list=ls())

bird.dat.og <- read.csv(file = "data/Experiment-2-birds-raw.csv",
                        header = TRUE)
head(bird.dat.og)

## --------------- Acclimatization period --------------------------------------

bird.dat <- dplyr::select(bird.dat.og, -c("camera.card", "X", "notes", "photo.number"))
bird.dat$date <- as.Date(with(bird.dat, paste(year, month, day, sep = "-")), "%Y-%m-%d")

bird.week <- bird.dat %>% 
  group_by(week = week(date)) %>% 
  summarise(value = sum(presence))

bird.week
# looks like a 2 week acclimatization period (Weeks 47 & 48)

bird.dat.drop <- bird.dat %>% 
  mutate(week = week(date)) %>% 
  filter(week != 48 & week != 47)

# Calculate days of coverage
range(bird.dat.drop$date)
as_date("2021-01-31")-as_date("2020-12-02")
# Time difference of 60 days

# remove NA's
bird.clean <- bird.dat.drop %>% 
  drop_na(species)

## --------------- Prepare data ------------------------------------------------

# Summarize data by week
bird.mod <- bird.clean |>
  group_by(site, week, treatment) |>
  summarise(count = n())

# Make sure week remains a factor
bird.mod$week <- as_factor(bird.mod$week)

# Some traps may not have had detections during a subset of weeks. We need
# to pivot wider to introduce true zeroes (i.e., )
bird.mod <- bird.mod |>
  pivot_wider(names_from = week, values_from = count)

# Replace NA with zero
bird.mod[is.na(bird.mod)] <- 0

# Bring it back to long form
bird.mod <- bird.mod |>
  pivot_longer(cols = 3:12, names_to = 'week', values_to = 'count')

# Check the data.
# We should have:
# 4 traps x 10 blocks = 40 traps
# 40 traps x 10 weeks = 400

# We are still missing true zeroes for site six, treatment 4
missing <- tibble(site = rep('six', 10),
                  week = as_factor(c('1','2','3','4','5','49','50','51','52','53')),
                  treatment = as_factor(rep('4', 10)),
                  count = rep(0,10))

bird.mod$treatment <- as_factor(bird.mod$treatment)
bird.mod$week <- as_factor(bird.mod$week)
bird.mod <- rbind(bird.mod, missing)

# Weeks are listed as week of year. It would be more useful for us to look at
# week since experiment started:
bird.mod$week.order <- NA
for(i in 1:nrow(bird.mod)){
  if(bird.mod$week[i] == 49){
    bird.mod$week.order[i] <- 1
  }
  if(bird.mod$week[i] == 50){
    bird.mod$week.order[i] <- 2
  }
  if(bird.mod$week[i] == 51){
    bird.mod$week.order[i] <- 3
  }
  if(bird.mod$week[i] == 52){
    bird.mod$week.order[i] <- 4
  }
  if(bird.mod$week[i] == 53){
    bird.mod$week.order[i] <- 5
  }
  if(bird.mod$week[i] == 1){
    bird.mod$week.order[i] <- 6
  }
  if(bird.mod$week[i] == 2){
    bird.mod$week.order[i] <- 7
  }
  if(bird.mod$week[i] == 3){
    bird.mod$week.order[i] <- 8
  }
  if(bird.mod$week[i] == 4){
    bird.mod$week.order[i] <- 9
  }
  if(bird.mod$week[i] == 5){
    bird.mod$week.order[i] <- 10
  }
}

# Reorder
bird.mod <- bird.mod |> dplyr::select(site, week.order, treatment, count)

# Make sure its a factor
bird.mod$week.order <- as_factor(bird.mod$week.order)
bird.mod$week.order <- factor(bird.mod$week.order, 
                      levels = c('1','2','3','4','5','6','7','8','9','10'))
bird.mod$treatment <- factor(bird.mod$treatment, 
                            levels = c("0", "4", "8", "12"))

# Export the data
write.csv(bird.mod, 'Data/Experiment-2-bird-observations.csv',
          row.names = FALSE)

## --------------- Model total detections --------------------------------------

# Check distribution
dev.new()
plotdist(bird.mod$count, histo = TRUE, demp = TRUE)
dev.new()
descdist(bird.mod$count, discrete=TRUE, boot=500) # negative binomial

bird.mod |>
  filter(treatment == '0') |>
  ungroup() |>
  dplyr::select(count) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # Negative binomial

bird.mod |>
  filter(treatment == '4') |>
  ungroup() |>
  dplyr::select(count) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # Negative binomial

bird.mod |>
  filter(treatment == '8') |>
  ungroup() |>
  dplyr::select(count) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # Negative binomial / poisson

bird.mod |>
  filter(treatment == '12') |>
  ungroup() |>
  dplyr::select(count) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # Negative binomial / poisson

# Negative binomial with week as continuous with interaction
birds.obs.fit.nb.1 <- glmer.nb(count ~ treatment * as.numeric(week.order) + (1 | site),
                      control = glmerControl(optimizer ="bobyqa"),
                      data = bird.mod) # Fails to converge with Nelder + bobyqa

# Poisson with week as continuous with interaction
birds.obs.fit.pois.1 <- glmer(count ~ treatment * as.numeric(week.order) + (1 | site),
                              control = glmerControl(optimizer ="Nelder_Mead"),
                               data = bird.mod, family = 'poisson')
sim.birds.obs.fit.pois.1 <- simulateResiduals(birds.obs.fit.pois.1) # Fails quantiles
plot(sim.birds.obs.fit.pois.1) # fails KS and quantiles
summary(birds.obs.fit.pois.1) # 8142.4

# Negative binomial with week as continuous with interaction no random
birds.obs.fit.nb.2 <- glm.nb(count ~ treatment * as.numeric(week.order) + site,
                               data = bird.mod)
sim.birds.obs.fit.nb.2 <- simulateResiduals(birds.obs.fit.nb.2) # Fails quantiles
dev.new()
plot(sim.birds.obs.fit.nb.2)
dev.new()
check_model(birds.obs.fit.nb.2) # potential collinearity, quantiles deeply skewed
summary(birds.obs.fit.nb.2) # AIC 1996.3

birds.obs.fit.pois.2 <- glm(count ~ treatment * as.numeric(week.order) + site,
                              data = bird.mod, family = 'poisson') 
sim.birds.obs.fit.pois.2 <- simulateResiduals(birds.obs.fit.pois.2)
plot(sim.birds.obs.fit.pois.2) # fails everything
summary(birds.obs.fit.pois.2) # AIC 8084

# Negative binomial with week as factor with no interaction
birds.obs.fit.nb.3 <- glmer.nb(count ~ treatment + week.order + (1 | site),
                               data = bird.mod)
sim.birds.obs.fit.nb.3 <- simulateResiduals(birds.obs.fit.nb.3) 
dev.new()
plot(sim.birds.obs.fit.nb.3) # Quantile deviations detected but ns
Anova(birds.obs.fit.nb.3)
summary(birds.obs.fit.nb.3) # AIC 1995.3
dev.new()
check_model(birds.obs.fit.nb.3)

# Clear out the model attempts
rm(birds.obs.fit.nb.1, birds.obs.fit.nb.2, sim.birds.obs.fit.nb.2,
   birds.obs.fit.nb.3, sim.birds.obs.fit.nb.3, birds.obs.fit.pois.1,
   sim.birds.obs.fit.pois.1, birds.obs.fit.pois.2, sim.birds.obs.fit.pois.2)

# Choosing model 3 (time as a factor, no interaction)
# Negative binomial with week as factor with no interaction
# This model did not have quantile issues, the observed residual variance better
# followed the predicted variance, and the residuals were normal. In this case,
# losing the interaction term was worth it.
birds.obs.fit.nb <- glmer.nb(count ~ treatment + week.order + (1 | site),
                               data = bird.mod)
Anova(birds.obs.fit.nb)
summary(birds.obs.fit.nb) 

sim.birds.obs.fit.nb <- simulateResiduals(birds.obs.fit.nb) 
dev.new()
plot(sim.birds.obs.fit.nb) # Quantile deviations detected but ns

dev.new()
check_model(birds.obs.fit.nb)

## --------------- Check assumptions -------------------------------------------

# Dharma package
testDispersion(birds.obs.fit.nb)
plotResiduals(birds.obs.fit.nb, form = bird.mod$site)
plotResiduals(birds.obs.fit.nb, form = bird.mod$treatment)
plotResiduals(birds.obs.fit.nb, form = bird.mod$week.order)
testOutliers(birds.obs.fit.nb)
testQuantiles(birds.obs.fit.nb)
testCategorical(birds.obs.fit.nb, catPred = bird.mod$treatment)
testCategorical(birds.obs.fit.nb, catPred = bird.mod$site)
testCategorical(birds.obs.fit.nb, catPred = bird.mod$week.order)


testUniformity(birds.obs.fit.nb, alternative = c('two.sided'))
testUniformity(birds.obs.fit.nb, alternative = c('less'))
testUniformity(birds.obs.fit.nb, alternative = c('greater'))


# Check differences among treatment
plot(count ~ treatment, data = bird.mod)
plot(residuals(birds.obs.fit.nb, type = 'pearson') ~ bird.mod$treatment)
abline(a = 0, b = 0, col = "blue", lwd = 2)
hist(residuals(birds.obs.fit.nb, type = 'pearson'))

# Plot residuals
resfit <- resid(birds.obs.fit.nb, type = 'pearson')
hist(resfit)
plot(bird.mod$treatment, resfit, 
     ylab = "Residuals", 
     xlab = "Treatment") 
abline(0, 0)  

plot(density(resid(birds.obs.fit.nb, type='pearson')))
plot(density(resid(birds.obs.fit.nb, type='deviance')))

# Look at predictions 
preds.lm <- predict(birds.obs.fit.nb)
par(mfrow = c(1, 2))
plot(count ~ treatment, 
     data = bird.mod)
plot(exp(preds.lm) ~ bird.mod$treatment) 

predict(birds.obs.fit.nb)
bird.mod$pred = predict(birds.obs.fit.nb)
bird.mod$predicted = predict(birds.obs.fit.nb)    # save the predicted values
bird.mod$residuals = residuals(birds.obs.fit.nb, type = 'pearson')  # save the residual values

# Look at the actual, predicted, and residual values
pred_df <- bird.mod %>% 
  dplyr::select(count, predicted, residuals)

plot(fitted(birds.obs.fit.nb) ~ bird.mod$count)
abline(0, 1, col = "blue", lwd = 2)  

plot(predict(birds.obs.fit.nb), residuals(birds.obs.fit.nb, type = 'working'))

# Compare predicted vs. actual
plot(density(bird.mod$count), xlim=c(0, 160), ylim=c(0, .08), main='M_0 y_hat')
lines(density(predict(birds.obs.fit.nb, type='response')), col='red')

# Check independence of observations
scatter.smooth(bird.mod$week.order, bird.mod$residuals)

# Check Cook's distance
library('influence.ME')
infl <- influence(birds.obs.fit.nb, obs = TRUE)
plot(infl, which = "cook", cutoff = 4/200)
# lots of influential values

infl <- influence(birds.obs.fit.nb, group = 'site')
plot(infl, which = "cook", cutoff = 4/10)
# 4 influential blocks

plot(cooks.distance(birds.obs.fit.nb), type='h')
# Handful of values over 1. Nothing dramatic
## --------------- Export model validation -------------------------------------
E1 <- resid(birds.obs.fit.nb, type = 'pearson')
F1 <- fitted(birds.obs.fit.nb)

# png("Figures/Experiment-2-bird-detect-model-validation.png", 
    # width = 800, height = 800, units = "px")
par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x = F1, y = E1, xlab = "Fitted values", ylab = "Normalized residuals",
     cex.lab = 2)
abline(h = 0, lty = 2)
boxplot(E1 ~ treatment, data = bird.mod, xlab = "Treatment", ylab = "Normalized residuals",
        cex.lab = 2)
abline(h = 0, lty = 2)
boxplot(E1~site, data = bird.mod, xlab = "Block", ylab = "Normalized residuals",
        cex.lab = 2)
abline(h = 0, lty = 2)
dev.off()

## --------------- Calculate and export means ----------------------------------

# Holm correction
emmeans(birds.obs.fit.nb, pairwise ~treatment, 
        type = 'response', adjust = 'holm')
confint(pairs(emmeans(birds.obs.fit.nb, ~treatment)), 
        type = "response", adjust = 'holm')

# No correction
emmeans(birds.obs.fit.nb, pairwise ~treatment, 
        type = 'response', adjust = 'none')
confint(pairs(emmeans(birds.obs.fit.nb, ~treatment)), 
        type = "response", adjust = 'none')

means <- as.data.frame(emmeans(birds.obs.fit.nb, pairwise ~treatment, 
                               type = 'response', adjust = 'none')[1])
write.csv(means, "Model-output/Experiment-2-weekly-bird-counts-emmeans.csv", row.names = FALSE)

# Calculate r2
r.squaredGLMM(birds.obs.fit.nb)
# Marginal - conditional = random
# 0.4696552-0.1584618 = 0.3111934
# 0.7868823-0.2654943 = 0.521388
# 0.6756451-0.2279628 = 0.4476823
