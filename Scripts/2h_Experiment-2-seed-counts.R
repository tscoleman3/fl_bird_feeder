## --------------- HEADER ------------------------------------------------------
## Script name: 2h_Experiment-2-seed-counts.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2022-10-21
## Date Last Modified: 2023-02-11
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script analyses the seed detection data by
## sampling date in experiment 2.

## In this script we model the distribution of the data and attempt
## models, beginning with the most complicated (block as a random factor)
## and interaction with time as a continuous variable). We settled on a negative 
## binomial model including block (i.e., site) as a random factor and the
## interaction term. Then, we run the model again with the filtered data.

## After producing the main model we check assumptions and export means and
## confidence intervals to use in the figure. For the filtered data, we export 
## the emmeans tables as pictures for the supplementary materials. We also export a ggplot object
## from emmip() for a supplementary figure describing interactive effects.

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
library(performance)

# Clear the decks
rm(list=ls())

# Bring in the data
seeds <- read.csv(file = "data/Experiment-2-clean-seeds.csv",
                  header = TRUE, stringsAsFactors = FALSE)

# Summarize the data
seeds.lg <- seeds |>
  pivot_longer(5:30, names_to = "SPECIES", values_to = "COUNT") |>
  group_by(DATE, BLOCK, TREATMENT) |>
  summarize(COUNT = sum(COUNT))

## --------------- Convert date to ordinal variable ----------------------------

# Add a numeric value for sampling period
seeds.lg$DATE.ORD <- NA
for(i in 1:nrow(seeds.lg)){
  if(seeds.lg$DATE[i] == "2020-11-27"){
    seeds.lg$DATE.ORD[i] <- 1
  }
  if(seeds.lg$DATE[i] == "2020-12-04"){
    seeds.lg$DATE.ORD[i] <- 2
  }
  if(seeds.lg$DATE[i] == "2020-12-11"){
    seeds.lg$DATE.ORD[i] <- 3
  }
  if(seeds.lg$DATE[i] == "2021-01-02"){
    seeds.lg$DATE.ORD[i] <- 4
  }
  if(seeds.lg$DATE[i] == "2021-02-03"){
    seeds.lg$DATE.ORD[i] <- 5
  }
} 

## --------------- Calculate days between samples ------------------------------
seeds.lg$DATE <- as_date(seeds.lg$DATE)

seeds.lg$SAMP.LENGTH <- NA
seeds.lg$EXP.DAYS <- NA

for (i in 1:nrow(seeds.lg)){
  if (seeds.lg$DATE[i] == '2020-11-27'){
    seeds.lg$SAMP.LENGTH[i] <- difftime('2020-11-27', '2020-11-20',  unit = c("days"))
    seeds.lg$EXP.DAYS[i] <- difftime('2020-11-27', '2020-11-20',  unit = c("days"))
  } 
  if (seeds.lg$DATE[i] == '2020-12-04'){
    seeds.lg$SAMP.LENGTH[i] <- difftime('2020-12-04', '2020-11-27',  unit = c("days"))
    seeds.lg$EXP.DAYS[i] <- difftime('2020-12-04', '2020-11-20',  unit = c("days"))
  } 
  if (seeds.lg$DATE[i] == '2020-12-11'){
    seeds.lg$SAMP.LENGTH[i] <- difftime('2020-12-11', '2020-12-04',  unit = c("days"))
    seeds.lg$EXP.DAYS[i] <- difftime('2020-12-11', '2020-11-20',  unit = c("days"))
  } 
  if (seeds.lg$DATE[i] == '2021-01-02'){
    seeds.lg$SAMP.LENGTH[i] <- difftime('2021-01-02', '2020-12-11',  unit = c("days"))
    seeds.lg$EXP.DAYS[i] <- difftime('2021-01-02', '2020-11-20',  unit = c("days"))
  } 
  if (seeds.lg$DATE[i] == '2021-02-03'){
    seeds.lg$SAMP.LENGTH[i] <- difftime('2021-02-03', '2021-01-02',  unit = c("days"))
    seeds.lg$EXP.DAYS[i] <- difftime('2021-02-03', '2020-11-20',  unit = c("days"))
    
  } 
}

## --------------- Check distribution ------------------------------------------

# Check the overall distribution
descdist(seeds.lg$COUNT, discrete = TRUE) # negative binomial

# Check the distribution for the treatment
seeds.lg |>
  filter(TREATMENT == 'Control') |>
  ungroup() |>
  dplyr::select(COUNT) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # negative binomial

seeds.lg |>
  filter(TREATMENT == 'Low') |>
  ungroup() |>
  dplyr::select(COUNT) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # Poisson

seeds.lg |>
  filter(TREATMENT == 'Medium') |>
  ungroup() |>
  dplyr::select(COUNT) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # negative binomial

seeds.lg |>
  filter(TREATMENT == 'High') |>
  ungroup() |>
  dplyr::select(COUNT) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # Poisson

## --------------- Full detections model ---------------------------------------

# Make sure everything is a factor
seeds.lg$TREATMENT <- as_factor(seeds.lg$TREATMENT)
seeds.lg$TREATMENT<- factor(seeds.lg$TREATMENT , 
                  levels = c("Control", "Low", "Medium", "High"))

seeds.lg$BLOCK <- as_factor(seeds.lg$BLOCK)
class(seeds.lg$EXP.DAYS)

# Run model
seed.detect.mod.nb.1 <- glmer.nb(COUNT ~  TREATMENT * EXP.DAYS + (1|BLOCK),
                            control = glmerControl(optimizer ="Nelder_Mead"),
                            data = seeds.lg)
Anova(seed.detect.mod.nb.1)
sim.seed.detect.mod.nb.1 <- simulateResiduals(seed.detect.mod.nb.1)
plot(sim.seed.detect.mod.nb.1) # KS test failed
dev.new()
check_model(seed.detect.mod.nb.1) # residuals skewed and residual
# variance doesn't follow predicted

seed.detect.mod.pois.1 <- glmer(COUNT ~  TREATMENT * EXP.DAYS + (1|BLOCK),
                            control = glmerControl(optimizer ="Nelder_Mead"),
                            data = seeds.lg, family = 'poisson')
Anova(seed.detect.mod.pois.1)
sim.seed.detect.mod.pois.1 <- simulateResiduals(seed.detect.mod.pois.1)
plot(sim.seed.detect.mod.pois.1) # KS and quantile failed
dev.new()
check_model(seed.detect.mod.pois.1) # residuals skewed and residual
# variance doesn't follow predicted

seed.detect.mod.nb.2 <- glm.nb(COUNT ~  TREATMENT * EXP.DAYS + BLOCK,
                               data = seeds.lg)
Anova(seed.detect.mod.nb.2) 
sim.seed.detect.mod.nb.2 <- simulateResiduals(seed.detect.mod.nb.2) # even worse
plot(sim.seed.detect.mod.nb.2) # good
dev.new()
check_model(seed.detect.mod.nb.2) # residuals worse

seed.detect.mod.pois.2 <- glm(COUNT ~  TREATMENT * EXP.DAYS + 1|BLOCK,
                                data = seeds.lg, family = 'poisson')
# no converge

seed.detect.mod.nb.3 <- glm.nb(COUNT ~  TREATMENT + as_factor(DATE.ORD) + BLOCK,
                               data = seeds.lg)
Anova(seed.detect.mod.nb.3) 
sim.seed.detect.mod.nb.3 <- simulateResiduals(seed.detect.mod.nb.3)
plot(sim.seed.detect.mod.nb.3) # good
dev.new()
check_model(seed.detect.mod.nb.3) # residuals skewed and residual
# variance doesn't follow predicted

seed.detect.mod.pois.3 <- glm(COUNT ~  TREATMENT + as_factor(DATE.ORD) + BLOCK,
                               data = seeds.lg, family = 'poisson')
Anova(seed.detect.mod.pois.3) 
sim.seed.detect.mod.pois.3 <- simulateResiduals(seed.detect.mod.pois.3)
plot(sim.seed.detect.mod.pois.3) # Fails all over
dev.new()
check_model(seed.detect.mod.pois.3) # residuals skewed 

# Best model is the initial negative binomial model with the interaction
# term and random effect of block. We chose this model because the simpler
# model did not perform much better than the complicated model.

rm(seed.detect.mod.nb.1, sim.seed.detect.mod.nb.1,
   seed.detect.mod.pois.1, sim.seed.detect.mod.pois.1,
   seed.detect.mod.nb.2, sim.seed.detect.mod.nb.2,
   seed.detect.mod.nb.3, sim.seed.detect.mod.nb.3,
   seed.detect.mod.pois.3, sim.seed.detect.mod.pois.3
)

# Run model
seed.detect.mod <- glmer.nb(COUNT ~  TREATMENT * EXP.DAYS + (1|BLOCK),
                                 control = glmerControl(optimizer ="Nelder_Mead"),
                                 data = seeds.lg)
Anova(seed.detect.mod)
sim.seed.detect.mod <- simulateResiduals(seed.detect.mod)
plot(sim.seed.detect.mod) # KS test failed

dev.new()
check_model(seed.detect.mod) # residuals skewed and residual
# variance doesn't follow predicted, collinearity

vif(seed.detect.mod)
# TREATMENT 1.702429
# EXP.DAYS (sampling period) 2.506455
# Interaction 2.076792

## --------------- Full model check assumptions --------------------------------

# Check overdispersion
E1 <- resid(seed.detect.mod, type = "pearson")
N <- nrow(seeds.lg)
p <- length(fixef(seed.detect.mod)) + 1
overdispersion <- sum(E1^2) / (N-p)
overdispersion # looks good

# Dharma package
testDispersion(seed.detect.mod)
dev.new()
plotResiduals(sim.seed.detect.mod, form = seeds.lg$BLOCK)
testOutliers(sim.seed.detect.mod)
testQuantiles(sim.seed.detect.mod)
testCategorical(sim.seed.detect.mod, catPred = seeds.lg$TREATMENT)
testCategorical(sim.seed.detect.mod, catPred = seeds.lg$BLOCK)

testUniformity(sim.seed.detect.mod, alternative = c('two.sided'))
testUniformity(sim.seed.detect.mod, alternative = c('less'))
testUniformity(sim.seed.detect.mod, alternative = c('greater'))


# Compare observed and residuals
plot(COUNT ~ TREATMENT, data = seeds.lg)
plot(residuals(seed.detect.mod, type = 'pearson') ~ seeds.lg$TREATMENT)
abline(a = 0, b = 0, col = "blue", lwd = 2)

# Plot residuals
resfit <- resid(seed.detect.mod, type = 'pearson')
hist(resfit)
plot(seeds.lg$TREATMENT, resfit, 
     ylab = "Residuals", 
     xlab = "Treatment") 
abline(0, 0)  

plot(density(resid(seed.detect.mod, type='pearson')))
plot(density(resid(seed.detect.mod, type='deviance')))

# Look at predictions 
preds.lm <- predict(seed.detect.mod)
par(mfrow = c(1, 2))
plot(COUNT ~ TREATMENT, 
     data = seeds.lg)
plot(exp(preds.lm) ~ seeds.lg$TREATMENT) 

predict(seed.detect.mod)
seeds.lg$pred = predict(seed.detect.mod)
seeds.lg$predicted = predict(seed.detect.mod) # save the predicted values
seeds.lg$residuals = residuals(seed.detect.mod, type = 'pearson')  # save the residual values

# Quick look at the actual, predicted, and residual values
pred_df <- seeds.lg %>% 
  dplyr::select(COUNT, predicted, residuals)

plot(fitted(seed.detect.mod) ~ seeds.lg$COUNT)
abline(0, 1, col = "blue", lwd = 2)  

# Residuals
plot(predict(seed.detect.mod), residuals(seed.detect.mod, type = 'working'))
# expected pattern

# Check Cook's distance
library('influence.ME')
infl <- influence(seed.detect.mod, obs = TRUE)
plot(infl, which = "cook", cutoff = 4/200)
# 14 values influential using 4/n. None above 1.

infl <- influence(seed.detect.mod, group = 'BLOCK')
plot(infl, which = "cook", cutoff = 4/10)
# 4/10 = 0.4

# Compare observed and predicted
plot(density(seeds.lg$COUNT), xlim=c(0, 160), ylim=c(0, .08), main='M_0 y_hat')
lines(density(predict(seed.detect.mod, type='response')), col='red')

# Check independence of observations
scatter.smooth(seeds.lg$EXP.DAYS, seeds.lg$residuals)

## --------------- Full export model validation --------------------------------
E1 <- resid(seed.detect.mod, type = 'pearson')
F1 <- fitted(seed.detect.mod)

# png("Figures/Experiment-2-seed-detect-model-validation.png", 
# width = 800, height = 800, units = "px")
par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x = F1, y = E1, xlab = "Fitted values", ylab = "Normalized residuals",
     cex.lab = 2)
abline(h = 0, lty = 2)
boxplot(E1 ~ TREATMENT, data = seeds.lg, xlab = "Treatment", ylab = "Normalized residuals",
        cex.lab = 2)
abline(h = 0, lty = 2)
boxplot(E1~BLOCK, data = seeds.lg, xlab = "Block", ylab = "Normalized residuals",
        cex.lab = 2)
abline(h = 0, lty = 2)
dev.off()

## --------------- Full calculate and export means -----------------------------

# Means Holm correction
emmeans(seed.detect.mod, pairwise ~TREATMENT, 
        type = 'response', adjust = 'holm')
confint(emmeans(seed.detect.mod, pairwise ~TREATMENT), 
        type = "response", adjust = 'holm') 

# Means no correction
emmeans(seed.detect.mod, revpairwise ~TREATMENT, 
        type = 'response', adjust = 'none')
confint(emmeans(seed.detect.mod, revpairwise~TREATMENT), 
        type = "response", adjust = 'none') 

means <- as.data.frame(emmeans(seed.detect.mod, pairwise ~TREATMENT, 
                               type = 'response', adjust = 'none')[1])
write.csv(means, "Model-output/Experiment-2-periodic-seed-counts-emmeans.csv", row.names = FALSE)

# Calculate interactive effects between treatment and week
emtrends(seed.detect.mod, pairwise~TREATMENT, var = "EXP.DAYS",
         adjust = "holm", type = 'response') # holm correction
confint(emtrends(seed.detect.mod, pairwise~TREATMENT, var = "EXP.DAYS",
                 adjust = "holm", type = 'response'))

emtrends(seed.detect.mod, revpairwise~TREATMENT, var = "EXP.DAYS",
         adjust = "none", type = 'response') # no correction
confint(emtrends(seed.detect.mod, revpairwise~TREATMENT, var = "EXP.DAYS",
                 adjust = "none", type = 'response'))

# Create a table for the trends
trends <- tibble(Treatment = c('Control', 'Low', 'Medium', 'High'),
                 Estimate = c(-0.0202, 0.0231, 0.0290, 0.0437),
                 LCL = c(-0.05173, -0.00156, 0.00620, 0.01954),
                 UCL = c(0.0113, 0.0478, 0.0519, 0.0679))
write.csv(trends, "Model-output/Experiment-2-periodic-seed-counts-emtrends.csv", row.names = FALSE)

# Calculate r2
r.squaredGLMM(seed.detect.mod)
# 0.21182160-0.11457391 = 0.09724769
# 0.40887481-0.22115963 = 0.1877152
# 0.05554567-0.03004455 = 0.02550112

# Fix order of factor levels for model
# seeds.lg$TREATMENT <- factor(seeds.lg$TREATMENT, 
#                              levels=c('Control', 'Low', 'Medium', 'High'),
#                              ordered = TRUE)
# # Run model
# seed.detect.mod <- glmer.nb(COUNT ~  TREATMENT * EXP.DAYS + (1|BLOCK),
#                             control = glmerControl(optimizer ="Nelder_Mead"),
#                             data = seeds.lg)
# Throwing up convergence issues after fixing factor order but still producing
# outcome

# Visualize the interactive effects

# emmip assigned treatment to a column called 'tvar'
colnames(seeds.lg)[3] <- 'tvar'

seed.count.interact <- emmip(seed.detect.mod, TREATMENT ~ EXP.DAYS, cov.reduce = range,
                             type = 'response')+
  scale_color_manual(values = c("darkgray",  "#00A9FF", "#00BF7D",  "#FF61CC"))+
  geom_jitter(data = seeds.lg, aes(x = EXP.DAYS, y = COUNT, color = tvar),
              shape = 20, height = 0.1, width = 2, alpha = 0.4, size = 3)+
  ylab('Seed count')+
  xlab('Days since start of experiment')+
  scale_y_continuous(
                     limits = c(-0.12,10),
                     oob = scales::squish)+
  scale_x_continuous(
                     limits = c(5,80),
                     breaks = seq(0,80,10))+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "none",
        axis.text = element_text(face="bold"),
        panel.grid.minor = element_blank())
saveRDS(seed.count.interact, file = "Model-output/periodic-mean-seed-count-interact.RDS")
# ggsave('Figures/Experiment-2-periodic-seed-count-interaction.png')

# Compare linear trends to actual values
ggplot(seeds.lg, aes(x = EXP.DAYS, y = COUNT, color = TREATMENT))+
  geom_jitter(size = 0.5, )+
  scale_color_manual(values = c("darkgray",  "#00A9FF", "#00BF7D",  "#FF61CC"))+
  ylab('Seed Count')+
  xlab('Days since start')+
  scale_x_continuous(expand = c(0,0),
                     limits = c(0,80),
                     breaks = seq(0,75,10))+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "right")+
  theme(axis.text = element_text(face="bold"))

## --------------- Filtered model ----------------------------------------------

# Clear the decks
rm(list = ls())

# Bring in the data
seeds <- read.csv(file = "data/Experiment-2-clean-seeds-filt.csv",
                  header = TRUE, stringsAsFactors = FALSE)

# Summarize
seeds.lg <- seeds |>
  pivot_longer(5:25, names_to = "SPECIES", values_to = "COUNT") |>
  group_by(DATE, BLOCK, TREATMENT) |>
  summarize(COUNT = sum(COUNT))

# Add the dates
seeds.lg$DATE.ORD <- NA
for(i in 1:nrow(seeds.lg)){
  if(seeds.lg$DATE[i] == "2020-11-27"){
    seeds.lg$DATE.ORD[i] <- 1
  }
  if(seeds.lg$DATE[i] == "2020-12-04"){
    seeds.lg$DATE.ORD[i] <- 2
  }
  if(seeds.lg$DATE[i] == "2020-12-11"){
    seeds.lg$DATE.ORD[i] <- 3
  }
  if(seeds.lg$DATE[i] == "2021-01-02"){
    seeds.lg$DATE.ORD[i] <- 4
  }
  if(seeds.lg$DATE[i] == "2021-02-03"){
    seeds.lg$DATE.ORD[i] <- 5
  }
} 

seeds.lg$DATE <- as_date(seeds.lg$DATE)
seeds.lg$SAMP.LENGTH <- NA
seeds.lg$EXP.DAYS <- NA

for (i in 1:nrow(seeds.lg)){
  if (seeds.lg$DATE[i] == '2020-11-27'){
    seeds.lg$SAMP.LENGTH[i] <- difftime('2020-11-27', '2020-11-20',  unit = c("days"))
    seeds.lg$EXP.DAYS[i] <- difftime('2020-11-27', '2020-11-20',  unit = c("days"))
  } 
  if (seeds.lg$DATE[i] == '2020-12-04'){
    seeds.lg$SAMP.LENGTH[i] <- difftime('2020-12-04', '2020-11-27',  unit = c("days"))
    seeds.lg$EXP.DAYS[i] <- difftime('2020-12-04', '2020-11-20',  unit = c("days"))
  } 
  if (seeds.lg$DATE[i] == '2020-12-11'){
    seeds.lg$SAMP.LENGTH[i] <- difftime('2020-12-11', '2020-12-04',  unit = c("days"))
    seeds.lg$EXP.DAYS[i] <- difftime('2020-12-11', '2020-11-20',  unit = c("days"))
  } 
  if (seeds.lg$DATE[i] == '2021-01-02'){
    seeds.lg$SAMP.LENGTH[i] <- difftime('2021-01-02', '2020-12-11',  unit = c("days"))
    seeds.lg$EXP.DAYS[i] <- difftime('2021-01-02', '2020-11-20',  unit = c("days"))
  } 
  if (seeds.lg$DATE[i] == '2021-02-03'){
    seeds.lg$SAMP.LENGTH[i] <- difftime('2021-02-03', '2021-01-02',  unit = c("days"))
    seeds.lg$EXP.DAYS[i] <- difftime('2021-02-03', '2020-11-20',  unit = c("days"))
    
  } 
}

# Make sure everything is a factor
seeds.lg$TREATMENT <- as_factor(seeds.lg$TREATMENT)
seeds.lg$BLOCK <- as_factor(seeds.lg$BLOCK)

seed.filt.detect.mod <- glmer.nb(COUNT ~  TREATMENT * EXP.DAYS + (1|BLOCK),
                                 data = seeds.lg) # No converge
Anova(seed.filt.detect.mod)

# Quick look at the model
sim.seed.detect.mod <- simulateResiduals(seed.filt.detect.mod)
plot(sim.seed.detect.mod)

dev.new()
check_model(seed.filt.detect.mod) 
# Fat tails / kurtosis

# Check model significance 
Anova(seed.filt.detect.mod, type = 2)

## --------------- Filtered Check assumptions ----------------------------------

# Compare treatments
plot(COUNT ~ TREATMENT, data = seeds.lg)
plot(residuals(seed.filt.detect.mod, type = 'pearson') ~ seeds.lg.filt$TREATMENT)
abline(a = 0, b = 0, col = "blue", lwd = 2)
hist(residuals(seed.filt.detect.mod, type = 'pearson'))

# Plot residuals
resfit <- resid(seed.filt.detect.mod, type = 'pearson')
hist(resfit)
plot(seeds.lg.filt$TREATMENT, resfit, 
     ylab = "Residuals", 
     xlab = "Treatment") 
abline(0, 0)  

plot(density(resid(seed.filt.detect.mod, type='pearson')))
plot(density(resid(seed.filt.detect.mod, type='deviance')))

# Look at predictions 
preds.lm <- predict(seed.filt.detect.mod)
par(mfrow = c(1, 2))
plot(COUNT ~ TREATMENT, 
     data = seeds.lg.filt)
plot(exp(preds.lm) ~ seeds.lg.filt$TREATMENT) 

predict(seed.filt.detect.mod)
seeds.lg.filt$pred = predict(seed.filt.detect.mod)
seeds.lg.filt$predicted = predict(seed.filt.detect.mod)    # save the predicted values
seeds.lg.filt$residuals = residuals(seed.filt.detect.mod, type = 'pearson')  # save the residual values

# Quick look at the actual, predicted, and residual values
pred_df <- seeds.lg.filt %>% 
  dplyr::select(COUNT, predicted, residuals)

# Plot fitted against observed
plot(fitted(seed.filt.detect.mod) ~ seeds.lg.filt$COUNT)
abline(0, 1, col = "blue", lwd = 2)  

# Plot predictions against residuals
plot(predict(seed.filt.detect.mod), residuals(seed.filt.detect.mod, type = 'working'))

# Check Cook's distances
plot(cooks.distance(seed.filt.detect.mod), type='h')

# Compare predicted vs. actual
plot(density(seeds.lg.filt$COUNT), xlim=c(0, 160), ylim=c(0, .08), main='M_0 y_hat')
lines(density(predict(seed.filt.detect.mod, type='response')), col='red')

# Check independence of observations
scatter.smooth(seeds.lg.filt, seeds.lg.filt$residuals)

## --------------- Filtered Export model validation ----------------------------
E1 <- resid(seed.filt.detect.mod, type = 'pearson')
F1 <- fitted(seed.filt.detect.mod)

# png("Figures/Experiment-2-seed-count-filt-model-validation.png", 
# width = 800, height = 800, units = "px")
par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x = F1, y = E1, xlab = "Fitted values", ylab = "Normalized residuals",
     cex.lab = 2)
abline(h = 0, lty = 2)
boxplot(E1 ~ TREATMENT, data = seeds.lg.filt, xlab = "Treatment", ylab = "Normalized residuals",
        cex.lab = 2)
abline(h = 0, lty = 2)
boxplot(E1~BLOCK, data = seeds.lg.filt, xlab = "Block", ylab = "Normalized residuals",
        cex.lab = 2)
abline(h = 0, lty = 2)
dev.off()

## --------------- Filtered Calculate and export means -------------------------
library(gridExtra)
# Calculate means averaged over entire sampling periods and export the outpu
# for supplementary materials

# Filtered data emmeans with Holm correction
filt.means.holm.em <- emmeans(seed.filt.detect.mod, pairwise ~TREATMENT, 
                              type = 'response', adjust = 'holm') # back-transform the values, holm correction
filt.means.holm <- as.data.frame(filt.means.holm.em$emmeans)
# Save the file in model output
png("Model-output/Experiment-2-periodic-seed-counts-filtered-holm.png", 
    height = 50*nrow(filt.means.holm), width = 200*ncol(filt.means.holm))
grid.table(filt.means.holm)
dev.off()

# Filtered data emmeans confidence intervals with Holm correction
filt.means.holm.conf <- as.data.frame(confint(filt.means.holm.em$contrasts))
# Make a column for p values
pval <- tibble(p.value = c(1.0000,0.7227,0.4368,0.7227,0.6552,1.0000))
filt.means.holm.conf <- cbind(filt.means.holm.conf,pval)
png("Model-output/Experiment-2-periodic-seed-counts-filtered-holm-confint.png", 
    height = 50*nrow(filt.means.holm.conf), width = 200*ncol(filt.means.holm.conf))
grid.table(filt.means.holm.conf)
dev.off()

# Filtered data emmeans without Holm correction
filt.means.no.correct.em <- emmeans(seed.filt.detect.mod, pairwise ~TREATMENT, 
                                    type = 'response', adjust = 'none') # back-transform the values, holm correction
filt.means.no.correct <- as.data.frame(filt.means.no.correct.em$emmeans)
png("Model-output/Experiment-2-periodic-seed-counts-filtered-no-correct.png", 
    height = 50*nrow(filt.means.no.correct), width = 200*ncol(filt.means.no.correct))
grid.table(filt.means.no.correct)
dev.off()

# Filtered data emmeans confidence intervals with Holm correction
filt.means.no.correct.conf <- as.data.frame(confint(filt.means.no.correct.em$contrasts))
pval <- tibble(p.value = c(0.9053,0.1807,0.0728,0.2360,0.1310,0.6291))
filt.means.no.correct.conf <- cbind(filt.means.no.correct.conf,pval)
png("Model-output/Experiment-2-periodic-seed-counts-filtered-no-correct-confint.png", 
    height = 50*nrow(filt.means.no.correct.conf), width = 200*ncol(filt.means.no.correct.conf))
grid.table(filt.means.no.correct.conf)
dev.off()

# Calculate interactive effects between treatment and week correction
filt.trends.holm.em <- emtrends(seed.filt.detect.mod, pairwise~TREATMENT, var = "EXP.DAYS",
                                adjust = "holm", type = 'response')
filt.trends.holm <- as.data.frame(filt.trends.holm.em$emtrends)
png("Model-output/Experiment-2-periodic-seed-counts-trends-filtered-holm.png", 
    height = 50*nrow(filt.trends.holm), width = 200*ncol(filt.trends.holm))
grid.table(filt.trends.holm)
dev.off()

# Confidence intervals correction
filt.trends.holm.confint <- as.data.frame(confint(filt.trends.holm.em$contrasts))
pval <- tibble(p.value = c(0.0353,0.1301,0.1856,0.5552,0.3622,0.5862))
filt.trends.holm.confint <- cbind(filt.trends.holm.confint,pval)
png("Model-output/Experiment-2-periodic-seed-counts-trends-filtered-holm-confint.png", 
    height = 50*nrow(filt.trends.holm.confint), width = 200*ncol(filt.trends.holm.confint))
grid.table(filt.trends.holm.confint)
dev.off()

# Calculate interactive effects between treatment and week no correction
filt.trends.no.correct.em <- emtrends(seed.filt.detect.mod, pairwise~TREATMENT, var = "EXP.DAYS",
                                      adjust = "none", type = 'response') # No adjustments
filt.trends.no.correct <- as.data.frame(filt.trends.no.correct.em$emtrends)
png("Model-output/Experiment-2-periodic-seed-counts-trends-filtered-no-correct.png", 
    height = 50*nrow(filt.trends.no.correct), width = 200*ncol(filt.trends.no.correct))
grid.table(filt.trends.no.correct)
dev.off()

# Confidence no correction
filt.trends.no.correct.confint <- as.data.frame(confint(filt.trends.no.correct.em$contrasts))
pval <- tibble(p.value = c(0.0059,0.0260,0.0464,0.2776,0.1207,0.5862))
filt.trends.no.correct.confint <- cbind(filt.trends.no.correct.confint,pval)
png("Model-output/Experiment-2-periodic-seed-counts-trends-filtered-no-correct-confint.png", 
    height = 50*nrow(filt.trends.no.correct.confint), width = 200*ncol(filt.trends.no.correct.confint))
grid.table(filt.trends.no.correct.confint)
dev.off()
