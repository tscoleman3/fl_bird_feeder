## --------------- HEADER ------------------------------------------------------
## Script name: 2k_Experiment-2-seed-rich.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2022-10-15
## Date Last Modified: 2023-02-11
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script analyzes the seed collection data for 
## experiment 2

## In this script we model both total and weekly richness. In both cases
## we begin by preparing data and calculating richness. Then we assess the 
## distribution of the data before attempting models, beginning with the 
## most complicated. For the total richness model this there is temporal variable. 

## For the total species richness model we settled on a poisson distribution
## including block as a random effect.

## For the weekly richness model, we settled on a negative binomial model
## with the continuous temporal interaction and block as a random effect.

## After producing the main model we check assumptions and export means and
## confidence intervals to use in the figure. We also export a ggplot object
## from emmip() for a supplementary figure describing interactive effects.

## Then we run the models again with the filtered data. For the filtered data 
## models, we export the emmeans tables as pictures for the supplementary materials. 


## --------------- Set up workspace --------------------------------------------
library(fitdistrplus)
library(lme4)
library(tidyverse)
library(car)
library(DHARMa)
library(lattice)
library(MuMIn)
library(emmeans)
library(vegan)
library(lubridate)
library(forcats)
library(gridExtra)

# Clear the deck
rm(list = ls())

# Bring in the data
seeds <- read.csv(file = "data/Experiment-2-clean-seeds.csv",
                  header = TRUE, stringsAsFactors = FALSE)

## --------------- Calculate total richness ------------------------------------

# Pivot the matrix into long format and summarize across sampling periods
seeds.lg <- seeds |>
  pivot_longer(cols = 5:30, names_to = 'SPECIES', values_to = 'COUNT') |>
  group_by(BLOCK, TRAP, TREATMENT, SPECIES) |>
  dplyr::summarize(COUNT = sum(COUNT))

# Make it wide again
seeds.wd <- seeds.lg |>
  pivot_wider(names_from = SPECIES, values_from = COUNT)

# Separate species matrix from predictors
comm.mat <- seeds.wd[,4:29]
predictors <- seeds.wd[,1:3]

# Convert to binary matrix
comm.mat.pa <- decostand(comm.mat, "pa")

# Recombine with predictors
seed.rich <- cbind(predictors[,1:3], comm.mat.pa)

# Calculate richness
seed.rich$RICH <- NA
for(i in 1:nrow(seed.rich)){
  seed.rich$RICH[i] <- sum(seed.rich[i,4:29])
}  

write.csv(seed.rich, "Data/Experiment-2-seed-richness.csv",
          row.names = FALSE)

## --------------- Model total richness ----------------------------------------

seed.rich$TREATMENT <- as_factor(seed.rich$TREATMENT)

# Visualize the data ##
boxplot(RICH ~ TREATMENT,col=c("white","lightgray",
                                      "darkgray", "black"), outline = TRUE, seed.rich)
                                      
# Plot distribution
plotdist(seed.rich$RICH, histo = TRUE, demp = TRUE)
descdist(seed.rich$RICH, discrete=TRUE, boot=500) # negative binomial

# Visualize data by treatment
seed.rich |>
  filter(TREATMENT == 'Control') |>
  ungroup() |>
  dplyr::select(RICH) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500) 
# Error in plot.window(...) : need finite 'xlim' values

seed.rich |>
  filter(TREATMENT == 'Control') |>
  ungroup() |>
  dplyr::select(RICH) |>
  as_vector() |>
  hist()
# Not normal

seed.rich |>
  filter(TREATMENT == 'Low') |>
  ungroup() |>
  dplyr::select(RICH) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # Normal

seed.rich |>
  filter(TREATMENT == 'Medium') |>
  ungroup() |>
  dplyr::select(RICH) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # Poisson                                      

seed.rich |>
  filter(TREATMENT == 'High') |>
  ungroup() |>
  dplyr::select(RICH) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # Normal                                      

# Feeder treatments are distributed differently. The model
# might have issues

seed.rich$BLOCK <- as_factor(seed.rich$BLOCK)

seeds.tot.rich.nb.1 <- glmer.nb(RICH ~ TREATMENT + (1|BLOCK), 
                                data = seed.rich) 
# Having trouble calculating dispersion parameter

seeds.tot.rich.pois.1 <- glmer(RICH ~ TREATMENT + (1|BLOCK), data = seed.rich,
                               family = 'poisson') 
sim.seeds.tot.rich.pois.1 <- simulateResiduals(seeds.tot.rich.pois.1)
dev.new()
plot(sim.seeds.tot.rich.pois.1) # quantile test fails
summary(seeds.tot.rich.pois.1) # 130.2

seeds.tot.rich.lm.1 <- lmer(RICH ~ TREATMENT + (1|BLOCK), data = seed.rich) 
sim.seeds.tot.rich.lm.1 <- simulateResiduals(seeds.tot.rich.lm.1)
plot(sim.seeds.tot.rich.lm.1)
AIC(seeds.tot.rich.lm.1) # 148.3612
seeds.tot.rich.nb.2 <- glm.nb(RICH ~ TREATMENT + BLOCK, data = seed.rich) 
# Having trouble calculating dispersion parameter

seeds.tot.rich.pois.2 <- glm(RICH ~ TREATMENT + BLOCK, data = seed.rich,
                             family = 'poisson')
sim.seeds.tot.rich.pois.2 <- simulateResiduals(seeds.tot.rich.pois.2)
plot(sim.seeds.tot.rich.pois.2) # problems converging, quantile test fails
summary(seeds.tot.rich.pois.2) # 130.12

seeds.tot.rich.nb.2 <- glm.nb(RICH ~ TREATMENT + BLOCK, data = seed.rich)
# No converging

# Remove old models
rm(seeds.tot.rich.pois.1, sim.seeds.tot.rich.pois.1, seeds.tot.rich.pois.2,
   sim.seeds.tot.rich.pois.2, seeds.tot.rich.nb.1, seeds.tot.rich.nb.2,
   seeds.tot.rich.lm.1, seeds.tot.rich.nb.2)

# Choose model
# Having trouble estimating theta in the negative binomial models and not
# getting convergence with poisson model. The best model is the poisson
# with the random effect for block
seeds.tot.rich.pois <- glmer(RICH ~ TREATMENT + (1|BLOCK), data = seed.rich,
                             family = 'poisson') 
Anova(seeds.tot.rich.pois)

## --------------- Check total model assumptions -------------------------------

# Check overall assumptions
sim.seeds.tot.rich.pois <- simulateResiduals(seeds.tot.rich.pois)
dev.new()
# Outer newtons do not converge fully
plot(sim.seeds.tot.rich.pois) # hit stop or wait

# Test quantiles
testQuantiles(seeds.tot.rich.pois) # hit stop or wait


# Dharma package
testDispersion(seeds.tot.rich.pois)
dev.new()
plotResiduals(seeds.tot.rich.pois, form = seed.rich$BLOCK)
plotResiduals(seeds.tot.rich.pois, form = seed.rich$TREATMENT)
testOutliers(seeds.tot.rich.pois)
testCategorical(seeds.tot.rich.pois, catPred = seed.rich$TREATMENT)
testCategorical(seeds.tot.rich.pois, catPred = seed.rich$BLOCK)

testUniformity(seeds.tot.rich.pois, alternative = c('two.sided'))
testUniformity(seeds.tot.rich.pois, alternative = c('less'))
testUniformity(seeds.tot.rich.pois, alternative = c('greater'))

# Overall assumptions
dev.new()
check_model(seeds.tot.rich.pois)

# Check overdispersion
E1 <- resid(seeds.tot.rich.pois, type = "pearson")
N <- nrow(seed.rich)
p <- length(coef(seeds.tot.rich.pois)) + 1
overdispersion <- sum(E1^2) / (N-p)
overdispersion # good

# Explore residuals and observed values by treatment
plot(RICH ~ TREATMENT, data = seed.rich)
plot(residuals(seeds.tot.rich.pois) ~ seed.rich$TREATMENT)
abline(a = 0, b = 0, col = "blue", lwd = 2)
hist(residuals(seeds.tot.rich.pois))

# Plot residuals 
resfit <- resid(seeds.tot.rich.pois)
hist(resfit)
plot(seed.rich$TREATMENT, resfit, 
     ylab = "Residuals", 
     xlab = "Detections") 
abline(0, 0)  

# Look at predictions 
preds.pois <- predict(seeds.tot.rich.pois)
par(mfrow = c(1, 2))
plot(RICH ~ TREATMENT, 
     data = seed.rich)
plot(preds.pois ~ seed.rich$TREATMENT) # funnel is problem

predict(seeds.tot.rich.pois)
seed.rich$pred = exp(predict(seeds.tot.rich.pois))
seed.rich$predicted = predict(seeds.tot.rich.pois)    # save the predicted values
seed.rich$residuals = residuals(seeds.tot.rich.pois)  # save the residual values

# Quick look at the actual, predicted, and residual values
pred_df <- seed.rich %>% 
  dplyr::select(RICH, predicted, residuals)
pred_df$predicted = exp(pred_df$predicted)
pred_df$residuals = exp(pred_df$residuals)

# Plot fitted by observed
plot(fitted(seeds.tot.rich.pois) ~ seed.rich$RICH)
abline(0, 1, col = "blue", lwd = 2)

# Plot prediced by residuals
plot(predict(seeds.tot.rich.pois), residuals(seeds.tot.rich.pois, type = 'working'))

# Check Cook's distance
plot(cooks.distance(seeds.tot.rich.pois), type='h')

# Compare predicted vs. actual
plot(density(seed.rich$RICH), xlim=c(0, 160), ylim=c(0, .08), main='M_0 y_hat')
lines(density(predict(seeds.tot.rich.pois, type='response')), col='red')

# Check independence of observations
scatter.smooth(seed.rich$BLOCK, seed.rich$residuals)

# Check Cook's distance
library('influence.ME')
infl <- influence(seeds.tot.rich.pois, obs = TRUE)
plot(infl, which = "cook", cutoff = 4/200)
# lots of influential values but few much higher

infl <- influence(seeds.tot.rich.pois, group = 'BLOCK')
plot(infl, which = "cook", cutoff = 4/10)
# 4/10 = 0.4

## --------------- Export total model outputs ----------------------------------

# EMMEANS
emmeans(seeds.tot.rich.pois, revpairwise ~TREATMENT, type = 'response',
        adjust = 'holm') # With correction
confint(emmeans(seeds.tot.rich.pois, pairwise ~TREATMENT, type = 'response',
                adjust = 'holm'))

emmeans(seeds.tot.rich.pois, revpairwise ~TREATMENT, type = 'response',
        adjust = 'none') # No correction
confint(emmeans(seeds.tot.rich.pois, revpairwise ~TREATMENT, type = 'response',
                adjust = 'none'))

# Calculate explained variance
r.squaredGLMM(seeds.tot.rich.pois)
# 0.2381983-0.10170547 = 0.1364928
# 0.2929250-0.12507259 = 0.1678524
# 0.1789071-0.07638944 = 0.1025177

# Export the means
means <- as.data.frame(emmeans(seeds.tot.rich.pois, pairwise ~TREATMENT, 
                               type = 'response', adjust = 'none')[1])
means2 <- as.data.frame(emmeans(seeds.tot.rich.pois, pairwise ~TREATMENT, 
                               type = 'response', adjust = 'holm')[1])
write.csv(means, "Model-output/Experiment-2-total-seed-rich-emmeans.csv", row.names = FALSE)

## --------------- Calculate weekly richness -----------------------------------

# Clear the deck
rm(list = ls())

# Bring in the data
seeds <- read.csv(file = "data/Experiment-2-clean-seeds.csv",
                  header = TRUE, stringsAsFactors = FALSE)

# Separate species matrix from predictors
comm.mat <- seeds[,5:30]
predictors <- seeds[,1:4]

# Convert to binary and recombine with predictors
comm.mat.pa <- decostand(comm.mat, "pa")
seed.rich <- cbind(predictors[,1:4], comm.mat.pa)

# Calculate richness
seed.rich$RICH <- NA
for(i in 1:nrow(seed.rich)){
  seed.rich$RICH[i] <- sum(seed.rich[i,5:30])
}  

## --------------- Convert predictors to ordinal -------------------------------

# Convert date to a dummy variable
seed.rich$DATE.ORD <- NA
for(i in 1:nrow(seed.rich)){
  if(seed.rich$DATE[i] == "2020-11-27"){
    seed.rich$DATE.ORD[i] <- 1
  }
  if(seed.rich$DATE[i] == "2020-12-04"){
    seed.rich$DATE.ORD[i] <- 2
  }
  if(seed.rich$DATE[i] == "2020-12-11"){
    seed.rich$DATE.ORD[i] <- 3
  }
  if(seed.rich$DATE[i] == "2021-01-02"){
    seed.rich$DATE.ORD[i] <- 4
  }
  if(seed.rich$DATE[i] == "2021-02-03"){
    seed.rich$DATE.ORD[i] <- 5
  }
} 

## --------------- Calculate days between samples ------------------------------
seed.rich$DATE <- as_date(seed.rich$DATE)

seed.rich$SAMP.LENGTH <- NA
seed.rich$EXP.DAYS <- NA

for (i in 1:nrow(seed.rich)){
  if (seed.rich$DATE[i] == '2020-11-27'){
    seed.rich$SAMP.LENGTH[i] <- difftime('2020-11-27', '2020-11-20',  unit = c("days"))
    seed.rich$EXP.DAYS[i] <- difftime('2020-11-27', '2020-11-20',  unit = c("days"))
  } 
  if (seed.rich$DATE[i] == '2020-12-04'){
    seed.rich$SAMP.LENGTH[i] <- difftime('2020-12-04', '2020-11-27',  unit = c("days"))
    seed.rich$EXP.DAYS[i] <- difftime('2020-12-04', '2020-11-20',  unit = c("days"))
  } 
  if (seed.rich$DATE[i] == '2020-12-11'){
    seed.rich$SAMP.LENGTH[i] <- difftime('2020-12-11', '2020-12-04',  unit = c("days"))
    seed.rich$EXP.DAYS[i] <- difftime('2020-12-11', '2020-11-20',  unit = c("days"))
  } 
  if (seed.rich$DATE[i] == '2021-01-02'){
    seed.rich$SAMP.LENGTH[i] <- difftime('2021-01-02', '2020-12-11',  unit = c("days"))
    seed.rich$EXP.DAYS[i] <- difftime('2021-01-02', '2020-11-20',  unit = c("days"))
  } 
  if (seed.rich$DATE[i] == '2021-02-03'){
    seed.rich$SAMP.LENGTH[i] <- difftime('2021-02-03', '2021-01-02',  unit = c("days"))
    seed.rich$EXP.DAYS[i] <- difftime('2021-02-03', '2020-11-20',  unit = c("days"))
    
  } 
}

## --------------- Weekly richness model ---------------------------------------

# Make sure the factors are factors
seed.rich$BLOCK <- forcats::as_factor(seed.rich$BLOCK)
seed.rich$TRAP <- forcats::as_factor(seed.rich$TRAP)

# Make sure the factors are in the correct order
seed.rich$TREATMENT <- factor(seed.rich$TREATMENT, order = TRUE,
                              levels = c('Control', 'Low', 'Medium', 'High'))

# Visualize the data ##
boxplot(RICH ~ TREATMENT,col=c("white","lightgray", 
                                      "darkgray", "black"), outline = TRUE, seed.rich)
                                      
# Plot distribution
plotdist(seed.rich$RICH, histo = TRUE, demp = TRUE)
descdist(seed.rich$RICH, discrete=TRUE, boot=500) # Negative binomial

# Check distribution by treatment
seed.rich |>
  filter(TREATMENT == 'Control') |>
  ungroup() |>
  dplyr::select(RICH) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500) # Poisson

seed.rich |>
  filter(TREATMENT == 'Low') |>
  ungroup() |>
  dplyr::select(RICH) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # Poisson

seed.rich |>
  filter(TREATMENT == 'Medium') |>
  ungroup() |>
  dplyr::select(RICH) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # Negative binomial                                      

seed.rich |>
  filter(TREATMENT == 'High') |>
  ungroup() |>
  dplyr::select(RICH) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # Poisson                                      

# Treatments follow different distributions, which may produce
# problems in the analysis

seeds.rich.nb <- glmer.nb(RICH ~ TREATMENT * EXP.DAYS + (1|BLOCK),
                            control = glmerControl(optimizer ="bobyqa"),
                            data = seed.rich) 
sim.seeds.rich.nb <- simulateResiduals(seeds.rich.nb)
plot(sim.seeds.rich.nb) # Looks good
Anova(seeds.rich.nb)

dev.new()
check_model(seeds.rich.nb) # fat tails, observed residual variance not following
# predicted residual variance at high values

# generalized variance inflation factors comparable
vif(seeds.rich.nb) 
# Treatment 1.952123, EXP.DAYS 1.240154, TREATMENT:EXP.DAYS, 2.028202,

r.squaredGLMM(seeds.rich.nb)
# 0.16258068-0.12329856 = 0.03928212
# 0.30207280-0.22908713 = 0.07298567
# 0.05187482-0.03934103 = 0.01253379
## --------------- Check model assumptions -------------------------------------

# Simulated residuals
sim.seeds.rich.nb <- simulateResiduals(seeds.rich.nb)
plot(sim.seeds.rich.nb)

testDispersion(seeds.rich.nb)
dev.new()
plotResiduals(sim.seeds.rich.nb, form = seed.rich$EXP.DAYS)
plotResiduals(sim.seeds.rich.nb, form = seed.rich$BLOCK)
testOutliers(sim.seeds.rich.nb)
testQuantiles(sim.seeds.rich.nb)
testCategorical(sim.seeds.rich.nb, catPred = seed.rich$TREATMENT)
testUniformity(sim.seeds.rich.nb, alternative = c('two.sided'))
testUniformity(sim.seeds.rich.nb, alternative = c('less'))
testUniformity(sim.seeds.rich.nb, alternative = c('greater'))

# Performance package
dev.new()
check_model(seeds.rich.nb)

# Check overdispersion
E1 <- resid(seeds.rich.nb, type = "pearson")
N <- nrow(seed.rich)
p <- length(coef(seeds.rich.nb)) + 1
overdispersion <- sum(E1^2) / (N-p)
overdispersion 

# check some assumptions
plot(RICH ~ TREATMENT, data = seed.rich)
plot(residuals(seeds.rich.nb) ~ seed.rich$TREATMENT)
abline(a = 0, b = 0, col = "blue", lwd = 2)
hist(residuals(seeds.rich.nb))

# plot residuals from prediction
resfit <- resid(seeds.rich.nb)
hist(resfit)
plot(seed.rich$TREATMENT, resfit, 
     ylab = "Residuals", 
     xlab = "Detections") 
abline(0, 0)  

# look at predictions 
preds.nb <- predict(seeds.rich.nb)
par(mfrow = c(1, 2))
plot(RICH ~ TREATMENT, 
     data = seed.rich)
plot(preds.nb ~ seed.rich$TREATMENT) # funnel is problem

# look at predictions from model
# neg binom
predict(seeds.rich.nb)
seed.rich$pred = exp(predict(seeds.rich.nb))
seed.rich$predicted = predict(seeds.rich.nb)    # save the predicted values
seed.rich$residuals = residuals(seeds.rich.nb)  # save the residual values

# quick look at the actual, predicted, and residual values
pred_df <- seed.rich %>% 
  dplyr::select(RICH, predicted, residuals)
pred_df$predicted = exp(pred_df$predicted)
pred_df$residuals = exp(pred_df$residuals)

plot(fitted(seeds.rich.nb) ~ seed.rich$RICH)
abline(0, 1, col = "blue", lwd = 2)

# Check Cook's distance
library('influence.ME')
infl <- influence(seeds.rich.nb, obs = TRUE)
plot(infl, which = "cook", cutoff = 4/200)
# lots of influential values
infl <- influence(seeds.rich.nb, group = 'BLOCK')
plot(infl, which = "cook", cutoff =  4/10)
# 3 Blocks over cut-off

## --------------- Export model validation -------------------------------------
E1 <- resid(seeds.rich.nb, type = 'pearson')
F1 <- fitted(seeds.rich.nb, type = 'response')

# png("Figures/Experiment-2-seed-rich-model-validation.png", 
# width = 800, height = 800, units = "px")

par(mfrow = c(2,2), mar = c(5,5,2,2))
plot(x = F1, y = E1, xlab = "Fitted values", ylab = "Pearson residuals",
     cex.lab = 2)
abline(h = 0, lty = 2)
plot(seeds$DATE.ORD, y = E1, xlab = "Sampling period", ylab = "Pearson residuals",
     cex.lab = 2)
abline(h = 0, lty = 2)
boxplot(E1 ~ TREATMENT, data = seeds, xlab = "Treatment", ylab = "Pearson residuals",
        cex.lab = 2)
abline(h = 0, lty = 2)
boxplot(E1~BLOCK, data = seeds, xlab = "Block", ylab = "Pearson residuals",
        cex.lab = 2)
abline(h = 0, lty = 2)
dev.off()

## --------------- Calculate and export means ----------------------------------

# Mean weekly seed richness
emmeans(seeds.rich.nb, pairwise ~TREATMENT, type = 'response', 
        adjust = "holm")
confint(pairs(emmeans(seeds.rich.nb, pairwise ~TREATMENT, type = 'response', 
        adjust = "holm")))

emmeans(seeds.rich.nb, pairwise ~TREATMENT, type = 'response', 
        adjust = "none") # no correction
confint(emmeans(seeds.rich.nb, pairwise ~TREATMENT, type = 'response', 
        adjust = "none")) # Pairs will not allow no adjustment?

means <- as.data.frame(emmeans(seeds.rich.nb, pairwise ~TREATMENT, type = 'response',
                               adjust = "none", at=list(week.order=c(1,2,3,4,5)))[1])
write.csv(means, "Model-output/Experiment-2-periodic-seed-rich-emmeans.csv", row.names = FALSE)


# Interactive effect between treatment and sampling period
emtrends(seeds.rich.nb, pairwise ~ TREATMENT, 
         var = 'EXP.DAYS', type = 'response', adjust = 'holm')
confint(emtrends(seeds.rich.nb, pairwise ~ TREATMENT, 
                 var = 'EXP.DAYS', type = 'response', adjust = 'holm'))
emtrends(seeds.rich.nb, pairwise ~ TREATMENT, 
         var = 'EXP.DAYS', type = 'response', adjust = 'none')
confint(emtrends(seeds.rich.nb, pairwise ~ TREATMENT, 
                 var = 'EXP.DAYS', type = 'response', adjust = 'none'))

trends <- tibble(Treatment = c('Control', 'Low', 'Medium', 'High'),
                 Estimate = c(-0.0222, 0.0180, 0.0304, 0.0339),
                 LCL = c(-0.05488, -0.00219, 0.01289, 0.01590),
                 UCL = c(0.0104, 0.0382, 0.0480, 0.0520))

write.csv(trends, "Model-output/Experiment-2-periodic-seed-rich-emtrends.csv", row.names = FALSE)

# Export the interaction for supplemental information

# emmip assigned treatment to a column called 'tvar'
colnames(seed.rich)[4] <- 'tvar'
  
seed.rich.interact <- emmip(seeds.rich.nb, TREATMENT ~ EXP.DAYS, cov.reduce = range,
                            type = 'response')+
  scale_color_manual(values = c("darkgray",  "#00A9FF", "#00BF7D",  "#FF61CC"))+
  geom_jitter(data = seed.rich, aes(x = EXP.DAYS, y = RICH, color = tvar),
             shape = 20, height = 0.1, width = 2, alpha = 0.4, size = 3)+
  ylab('Seed richness')+
  scale_y_continuous(limits = c(-0.11,4.1))+
  scale_x_continuous(limits = c(5,80),
                     breaks = seq(0,80,10))+
  xlab('Days since start of experiment')+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "none",
        axis.text = element_text(face="bold"),
        panel.grid.minor = element_blank())
saveRDS(seed.rich.interact, file = "Model-output/periodic-mean-seed-rich-interact.RDS")

# ggsave("Figures/Experiment-2-periodic-seed-rich-interaction.png", width = 5, height = 7, units = "in")

# Compare linear trends to actual values
ggplot(seed.rich, aes(x = EXP.DAYS, y = RICH, color = TREATMENT))+
  geom_jitter(size = 0.5, )+
  scale_color_manual(values = c("darkgray",  "#00A9FF", "#00BF7D",  "#FF61CC"))+
  ylab('Seed Rich')+
  xlab('Days Since Start')+
  scale_x_continuous(expand = c(0,0),
                     limits = c(0,80),
                     breaks = seq(0,75,10))+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "right")+
  theme(axis.text = element_text(face="bold"))



## --------------- FILTERED TOTAL MODEL ----------------------------------------

rm(list = ls())

# Bring in the data
seeds <- read.csv(file = "data/Experiment-2-clean-seeds-filt.csv",
                  header = TRUE, stringsAsFactors = FALSE)

comm.mat <- seeds[,5:25]
predictors <- seeds[,1:4]

comm.mat.pa <- decostand(comm.mat, "pa")
seeds <- cbind(predictors[,1:4], comm.mat.pa)

# Calculate richness
seeds.lg <- seeds |>
  pivot_longer(cols = 5:25, names_to = 'SPECIES', values_to = 'COUNT') |>
  group_by(BLOCK, TRAP, TREATMENT, SPECIES) |>
  dplyr::summarize(COUNT = sum(COUNT))

seeds.wd <- seeds.lg |>
  pivot_wider(names_from = SPECIES, values_from = COUNT)

# Separate species matrix from predictors
comm.mat <- seeds.wd[,4:24]
predictors <- seeds.wd[,1:3]

comm.mat.pa <- decostand(comm.mat, "pa")
seed.rich <- cbind(predictors[,1:3], comm.mat.pa)

seed.rich$RICH <- NA
for(i in 1:nrow(seed.rich)){
  seed.rich$RICH[i] <- sum(seed.rich[i,4:24])
}  

# Model
seeds.tot.rich.pois <- glmer(RICH ~ TREATMENT + (1|BLOCK), data = seed.rich,
                             family = 'poisson') 
Anova(seeds.tot.rich.pois) # Anova

# Exporting tables for supplementary materials

# Filtered data emmeans with Holm correction
filt.means.holm.emm <- emmeans(seeds.tot.rich.pois, pairwise ~TREATMENT, type = 'response',
                               adjust = 'holm') # Holm correction
filt.means.holm <- as.data.frame(filt.means.holm.emm$emmeans)
# Save the file in model output
png("Model-output/Experiment-2-total-seed-rich-filtered-holm.png", 
    height = 50*nrow(filt.means.holm), width = 200*ncol(filt.means.holm))
grid.table(filt.means.holm)
dev.off()

# Filtered data emmeans confidence intervals with Holm correction
filt.means.holm.confint <- as.data.frame(confint(filt.means.holm.emm$contrasts))
pval <- tibble(p.value = c(0.1963, 0.8268, 0.1628, 0.8268, 0.8474, 0.8268))
filt.means.holm.confint <- cbind(filt.means.holm.confint, pval)
png("Model-output/Experiment-2-total-seed-rich-filtered-holm-confint.png", 
    height = 50*nrow(filt.means.holm.confint), width = 200*ncol(filt.means.holm.confint))
grid.table(filt.means.holm.confint)
dev.off()

# Filtered data emmeans without Holm correction
filt.means.no.correct.emm <- emmeans(seeds.tot.rich.pois, pairwise ~TREATMENT, type = 'response',
                                     adjust = 'none') # Holm correction
filt.means.no.correct <- as.data.frame(filt.means.no.correct.emm$emmeans)
# Save the file in model output
png("Model-output/Experiment-2-total-seed-rich-filtered-no-correct.png", 
    height = 50*nrow(filt.means.no.correct), width = 200*ncol(filt.means.no.correct))
grid.table(filt.means.no.correct)
dev.off()

# Filtered data emmeans confidence intervals without Holm correction
filt.means.no.correct.confint <- as.data.frame(confint(filt.means.no.correct.emm$contrasts))
pval <- tibble(p.value = c(0.0393, 0.2577, 0.0271, 0.2799, 0.8474, 0.2067))
filt.means.no.correct.confint <- cbind(filt.means.no.correct.confint, pval)
png("Model-output/Experiment-2-total-seed-rich-filtered-no-correct-confint.png", 
    height = 50*nrow(filt.means.no.correct.confint), width = 200*ncol(filt.means.no.correct.confint))
grid.table(filt.means.no.correct.confint)
dev.off()

emmeans(seeds.tot.rich.pois, pairwise ~TREATMENT, type = 'response',
        adjust = 'none') # No correction
confint(emmeans(seeds.tot.rich.pois, pairwise ~TREATMENT, type = 'response',
                adjust = 'none')) # No correction


## --------------- FILTERED WEEKLY MODEL ---------------------------------------

rm(list = ls())

# Bring in the data
seeds <- read.csv(file = "data/Experiment-2-clean-seeds-filt.csv",
                  header = TRUE, stringsAsFactors = FALSE)

# Separate species matrix from predictors
comm.mat <- seeds[,5:25]
predictors <- seeds[,1:4]

comm.mat.pa <- decostand(comm.mat, "pa")
seed.rich <- cbind(predictors[,1:4], comm.mat.pa)

seed.rich$RICH <- NA
for(i in 1:nrow(seed.rich)){
  seed.rich$RICH[i] <- sum(seed.rich[i,5:25])
}  

# Convert date to a dummy variable
seed.rich$DATE.ORD <- NA
for(i in 1:nrow(seed.rich)){
  if(seed.rich$DATE[i] == "2020-11-27"){
    seed.rich$DATE.ORD[i] <- 1
  }
  if(seed.rich$DATE[i] == "2020-12-04"){
    seed.rich$DATE.ORD[i] <- 2
  }
  if(seed.rich$DATE[i] == "2020-12-11"){
    seed.rich$DATE.ORD[i] <- 3
  }
  if(seed.rich$DATE[i] == "2021-01-02"){
    seed.rich$DATE.ORD[i] <- 4
  }
  if(seed.rich$DATE[i] == "2021-02-03"){
    seed.rich$DATE.ORD[i] <- 5
  }
} 

seed.rich$DATE <- as_date(seed.rich$DATE)
seed.rich$SAMP.LENGTH <- NA
seed.rich$EXP.DAYS <- NA

for (i in 1:nrow(seed.rich)){
  if (seed.rich$DATE[i] == '2020-11-27'){
    seed.rich$SAMP.LENGTH[i] <- difftime('2020-11-27', '2020-11-20',  unit = c("days"))
    seed.rich$EXP.DAYS[i] <- difftime('2020-11-27', '2020-11-20',  unit = c("days"))
  } 
  if (seed.rich$DATE[i] == '2020-12-04'){
    seed.rich$SAMP.LENGTH[i] <- difftime('2020-12-04', '2020-11-27',  unit = c("days"))
    seed.rich$EXP.DAYS[i] <- difftime('2020-12-04', '2020-11-20',  unit = c("days"))
  } 
  if (seed.rich$DATE[i] == '2020-12-11'){
    seed.rich$SAMP.LENGTH[i] <- difftime('2020-12-11', '2020-12-04',  unit = c("days"))
    seed.rich$EXP.DAYS[i] <- difftime('2020-12-11', '2020-11-20',  unit = c("days"))
  } 
  if (seed.rich$DATE[i] == '2021-01-02'){
    seed.rich$SAMP.LENGTH[i] <- difftime('2021-01-02', '2020-12-11',  unit = c("days"))
    seed.rich$EXP.DAYS[i] <- difftime('2021-01-02', '2020-11-20',  unit = c("days"))
  } 
  if (seed.rich$DATE[i] == '2021-02-03'){
    seed.rich$SAMP.LENGTH[i] <- difftime('2021-02-03', '2021-01-02',  unit = c("days"))
    seed.rich$EXP.DAYS[i] <- difftime('2021-02-03', '2020-11-20',  unit = c("days"))
    
  } 
}

seed.rich$BLOCK <- as_factor(seed.rich$BLOCK)

seeds.rich.nb <- glmer.nb(RICH ~ TREATMENT * EXP.DAYS + (1|BLOCK),
                          control = glmerControl(optimizer ="bobyqa"),
                          data = seed.rich) 
Anova(seeds.rich.nb) # iteration/alteration limit reached

# Mean weekly seed richness
emmeans(seeds.rich.nb, pairwise ~TREATMENT, type = 'response', 
        adjust = "holm")
emmeans(seeds.rich.nb, pairwise ~TREATMENT, type = 'response', 
        adjust = "none") # no correction

# Filtered data emmeans with Holm correction
filt.means.holm.em <- emmeans(seeds.rich.nb, pairwise ~TREATMENT, 
                              type = 'response', adjust = 'holm')# averaged over dates
filt.means.holm <- as.data.frame(filt.means.holm.em$emmeans)
# Save the file in model output
png("Model-output/Experiment-2-periodic-seed-rich-filtered-holm.png", 
    height = 50*nrow(filt.means.holm), width = 200*ncol(filt.means.holm))
grid.table(filt.means.holm)
dev.off()

# Filtered data emmeans confidence intervals with Holm correction
filt.means.holm.conf <- filt.means.holm <- as.data.frame(confint(filt.means.holm.em$contrasts))
# Make a column for p values
pval <- tibble(p.value = c(1.0000,1.0000,1.0000,1.0000,1.0000,1.0000))
filt.means.holm.conf <- cbind(filt.means.holm.conf,pval)
png("Model-output/Experiment-2-periodic-seed-rich-filtered-holm-confint.png", 
    height = 50*nrow(filt.means.holm.conf), width = 200*ncol(filt.means.holm.conf))
grid.table(filt.means.holm.conf)
dev.off()

# Filtered data emmeans without Holm correction
filt.means.no.correct.em <- emmeans(seeds.rich.nb, pairwise ~TREATMENT, 
                                    type = 'response', adjust = 'none') # back-transform the values, holm correction
filt.means.no.correct <- as.data.frame(filt.means.no.correct.em$emmeans)
png("Model-output/Experiment-2-periodic-seed-rich-filtered-no-correct.png", 
    height = 50*nrow(filt.means.no.correct), width = 200*ncol(filt.means.no.correct))
grid.table(filt.means.no.correct)
dev.off()

# Filtered data emmeans confidence intervals with Holm correction
filt.means.no.correct.conf <- as.data.frame(confint(filt.means.no.correct.em$contrasts))
pval <- tibble(p.value = c(0.7189,0.4181,0.2574,0.3088,0.2059,0.7002))
filt.means.no.correct.conf <- cbind(filt.means.no.correct.conf,pval)
png("Model-output/Experiment-2-periodic-seed-rich-filtered-no-correct-confint.png", 
    height = 50*nrow(filt.means.no.correct.conf), width = 200*ncol(filt.means.no.correct.conf))
grid.table(filt.means.no.correct.conf)
dev.off()

# Calculate interactive effects between treatment and week correction
filt.trends.holm.em <- emtrends(seeds.rich.nb, pairwise~TREATMENT, var = "EXP.DAYS",
                                adjust = "holm", type = 'response')
filt.trends.holm <- as.data.frame(filt.trends.holm.em$emtrends)
png("Model-output/Experiment-2-periodic-seed-rich-trends-filtered-holm.png", 
    height = 50*nrow(filt.trends.holm), width = 200*ncol(filt.trends.holm))
grid.table(filt.trends.holm)
dev.off()

# Confidence intervals correction
filt.trends.holm.confint <- as.data.frame(confint(filt.trends.holm.em$contrasts))
pval <- tibble(p.value = c(0.0297,0.2564,0.1212,0.2564,0.3980,0.3980))
filt.trends.holm.confint <- cbind(filt.trends.holm.confint,pval)
png("Model-output/Experiment-2-periodic-seed-rich-trends-filtered-holm-confint.png", 
    height = 50*nrow(filt.trends.holm.confint), width = 200*ncol(filt.trends.holm.confint))
grid.table(filt.trends.holm.confint)
dev.off()

# Calculate interactive effects between treatment and week no correction
filt.trends.no.correct.em <- emtrends(seeds.rich.nb, pairwise~TREATMENT, var = "EXP.DAYS",
                                      adjust = "none", type = 'response') # No adjustments
filt.trends.no.correct <- as.data.frame(filt.trends.no.correct.em$emtrends)
png("Model-output/Experiment-2-periodic-seed-rich-trends-filtered-no-correct.png", 
    height = 50*nrow(filt.trends.no.correct), width = 200*ncol(filt.trends.no.correct))
grid.table(filt.trends.no.correct)
dev.off()

# Confidence no correction
filt.trends.no.correct.confint <- as.data.frame(confint(filt.trends.no.correct.em$contrasts))
pval <- tibble(p.value = c(0.0049,0.0841,0.0242,0.0641,0.1990,0.3765))
filt.trends.no.correct.confint <- cbind(filt.trends.no.correct.confint,pval)
png("Model-output/Experiment-2-periodic-seed-rich-trends-filtered-no-correct-confint.png", 
    height = 50*nrow(filt.trends.no.correct.confint), width = 200*ncol(filt.trends.no.correct.confint))
grid.table(filt.trends.no.correct.confint)
dev.off()
