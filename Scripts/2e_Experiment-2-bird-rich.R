## --------------- HEADER ------------------------------------------------------
## Script name: 2e_Experiment-2-bird-rich.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2022-10-15
## Date Last Modified: 2023-02-22
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script analyzes the camera trap data for birds
## in experiment 2.

## In this script we model both total and weekly richness. In both cases
## we begin by preparing data and calculating richness. Then we assess the 
## distribution of the data before attempting models, beginning with the 
## most complicated. For the total richness model this there is temporal variable. 
## If this complicated model performs poorly with reasonable distributions, 
## we try less complicated models. Eventually we settle on a mixed linear model
## with a random blocking effect for the overall richness model and a poisson
## model with week.since.acclimatization(continuous)*treatment interaction and 
## block as a fixed effect.

## After producing the model we check assumptions and export means and
## confidence intervals to use in the figure. We also export a ggplot object
## from emmip() for a supplementary figure describing interactive effects.

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
library(pscl)
library(performance)

rm(list=ls())

bird.dat.og <- read.csv(file = "data/Experiment-2-birds-raw.csv",
                        header = TRUE)
head(bird.dat.og)

## --------------- Drop observations before acclimatization --------------------

bird.dat <- dplyr::select(bird.dat.og, -c("camera.card", "X", "notes", "photo.number"))
bird.dat$date <- as.Date(with(bird.dat, paste(year, month, day, sep = "-")), "%Y-%m-%d")

bird.week <- bird.dat %>% 
  group_by(week = week(date)) %>% 
  summarise(value = sum(presence))

bird.week
# looks like a 2 week acclimatization period
# observations really picked up in January (observations are
# essentially increasing with time throughout our sampling period)

bird.dat.drop <- bird.dat %>% 
  mutate(week = week(date)) %>% 
  filter(week != 48 & week != 47)

# remove NA's
bird.clean <- bird.dat.drop %>% 
  drop_na(species)

bird.clean$treatment <- as.factor(bird.clean$treatment)

## --------------- Prepare total richness --------------------------------------

# count of individuals by species at each treatment and site
bird.sps.count.by.treat <- bird.clean %>% 
  group_by(site, treatment, species) %>% 
  summarise(count = n())

bird.wd <- bird.sps.count.by.treat |>
  pivot_wider(names_from = 'species', values_from = 'count')

# Missing true zeroes
bird.wd[is.na(bird.wd)] <- 0

missing <- data.frame(site = "six", 
                      treatment = "4",
                      chipping_sparrow_spizella_passerina = 0,
                      eastern_phoebe_sayornis_phoebe = 0,
                      mourning_dove_zenaida_macroura = 0,
                      northern_cardinal_cardinalis_cardinalis = 0,
                      gray_catbird_dumetella_carolinensis = 0,
                      red_shoulder_hawk_buteo_lineatus = 0,
                      pine_warbler_setophaga_pinus = 0,
                      sedge_wren_cistothorus_stellaris = 0,
                      eastern_screech_owl_megascops_asio = 0,
                      barred_owl_strix_varia = 0)

bird.wd <- rbind(bird.wd, missing)

## --------------- Create total richness dataframe -----------------------------

spec <- bird.wd[,3:12]
pred <- bird.wd[,1:2]

# Convert to binary
spec <- as.matrix((spec>0)+0)

# Bring back together
bird.wd <- cbind(pred,spec)

bird.rich <- bird.wd |>
  mutate(richness = chipping_sparrow_spizella_passerina+
           eastern_phoebe_sayornis_phoebe+
           mourning_dove_zenaida_macroura+
           northern_cardinal_cardinalis_cardinalis+
           gray_catbird_dumetella_carolinensis+
           red_shoulder_hawk_buteo_lineatus+
           pine_warbler_setophaga_pinus+
           sedge_wren_cistothorus_stellaris+
           eastern_screech_owl_megascops_asio+
           barred_owl_strix_varia)

bird.rich$site <- as_factor(bird.rich$site)

# Save the richness dataframe
write.csv(bird.rich, 'Data/Experiment-2-bird-richness.csv',
          row.names = FALSE)
## --------------- Model total richness ----------------------------------------

# Visualize the data 
boxplot(richness ~ as_factor(treatment), outline = TRUE, data = bird.rich)

plotdist(bird.rich$richness, histo = TRUE, demp = TRUE)
descdist(bird.rich$richness, discrete=TRUE, boot=500)

bird.rich |>
  filter(treatment == '0') |>
  ungroup() |>
  dplyr::select(richness) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # Normal

bird.rich |>
  filter(treatment == '4') |>
  ungroup() |>
  dplyr::select(richness) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # Normal/Poisson

bird.rich |>
  filter(treatment == '8') |>
  ungroup() |>
  dplyr::select(richness) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # Poisson

bird.rich |>
  filter(treatment == '12') |>
  ungroup() |>
  dplyr::select(richness) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # Normal/Poisson

bird.rich$treatment <- factor(bird.rich$treatment, 
                              levels = c("0", "4", "8", "12"))

bird.rich.pois <- glmer(richness ~ treatment + (1|site),
                     data = bird.rich,
                     family = 'poisson') #
summary(bird.rich.pois)
Anova(bird.rich.pois)
dev.new()
check_model(bird.rich.pois) # Not performing well


bird.rich.lm <- lmer(richness ~ treatment + (1|site),
                        data = bird.rich)
summary(bird.rich.lm)
Anova(bird.rich.lm)
dev.new()
check_model(bird.rich.lm) # Much better

## --------------- Check total richness assumptions ----------------------------

# Dharma package
testDispersion(bird.rich.lm)
dev.new()
plotResiduals(bird.rich.lm, form = bird.rich$site)
plotResiduals(bird.rich.lm, form = bird.rich$treatment)
testOutliers(bird.rich.lm)
testCategorical(bird.rich.lm, catPred = bird.rich$treatment)
testCategorical(bird.rich.lm, catPred = bird.rich$site)

testUniformity(bird.rich.lm, alternative = c('two.sided'))
testUniformity(bird.rich.lm, alternative = c('less'))
testUniformity(bird.rich.lm, alternative = c('greater'))


# Check overdispersion
E1 <- resid(bird.rich.lm, type = "pearson")
N <- nrow(bird.rich)
p <- length(fixef(bird.rich.lm)) + 1
overdispersion <- sum(E1^2) / (N-p)
overdispersion 

# check other assumptions
dev.new()
check_model(bird.rich.lm)

sim.bird.rich.lm <- simulateResiduals(bird.rich.lm)
plot(sim.bird.rich.lm)

# plot residuals
plot(richness ~ as_factor(treatment), data = bird.rich)
plot(residuals(bird.rich.lm) ~ as_factor(bird.rich$treatment))
abline(a = 0, b = 0, col = "blue", lwd = 2)
hist(residuals(bird.rich.lm))

# plot residuals from prediction
resfit <- resid(bird.rich.lm)
hist(resfit)
plot(bird.rich$richness, resfit, 
     ylab = "Residuals", 
     xlab = "Detections") 
abline(0, 0)     

# look at predictions 
preds.lm <- predict(bird.rich.lm)
par(mfrow = c(1, 2))
plot(richness ~ as_factor(treatment), 
     data = bird.rich)
plot(preds.lm ~ as_factor(bird.rich$treatment))

# look at predictions from model
# neg binom
predict(bird.rich.lm)
bird.rich$pred = exp(predict(bird.rich.lm))
bird.rich$predicted = predict(bird.rich.lm)    # save the predicted values
bird.rich$residuals = residuals(bird.rich.lm)  # save the residual values

# quick look at the actual, predicted, and residual values
pred_df <- bird.rich %>% 
  dplyr::select(richness, predicted, residuals)
pred_df$predicted = exp(pred_df$predicted)
pred_df$residuals = exp(pred_df$residuals)

plot(fitted(bird.rich.lm) ~ bird.rich$richness)
abline(0, 1, col = "blue", lwd = 2)

plot(predict(bird.rich.lm), residuals(bird.rich.lm, type = 'working'))
plot(cooks.distance(bird.rich.lm), type='h') # All less than 1

# compare predicted vs. actual
plot(density(bird.rich$richness), xlim=c(0, 160), ylim=c(0, .08), main='M_0 y_hat')
lines(density(predict(bird.rich.lm, type='response')), col='red')

# independence of observations
scatter.smooth(bird.rich$site, bird.rich$residuals)

## --------------- Export coefficient estimates --------------------------------

# EMMEANS
emmeans(bird.rich.lm, revpairwise ~treatment, type = 'response',
        adjust = 'holm')
confint(emmeans(bird.rich.lm, revpairwise ~treatment, type = 'response',
        adjust = 'holm'))

emmeans(bird.rich.lm, pairwise ~treatment, type = 'response',
        adjust = 'none')
confint(emmeans(bird.rich.lm, pairwise ~treatment, type = 'response',
                adjust = 'none'))

means <- as.data.frame(emmeans(bird.rich.lm, pairwise ~treatment, 
                               type = 'response', adjust = 'none')[1])
write.csv(means, "Model-output/Experiment-2-total-bird-rich-emmeans.csv", row.names = FALSE)


# Calculate r2
r.squaredGLMM(bird.rich.lm)
# 0.8574247-0.02106278 = 0.8363619

## --------------- Prepare weekly richness dataframe ---------------------------

# count of individuals by species at each treatment and site
bird.sps.count.by.treat <- bird.clean %>% 
  group_by(site, week, treatment, species) %>% 
  summarise(count = n())

# Fix week 
bird.sps.count.by.treat$week.order <- NA
for(i in 1:nrow(bird.sps.count.by.treat)){
  if(bird.sps.count.by.treat$week[i] == 49){
    bird.sps.count.by.treat$week.order[i] <- 1
  }
  if(bird.sps.count.by.treat$week[i] == 50){
    bird.sps.count.by.treat$week.order[i] <- 2
  }
  if(bird.sps.count.by.treat$week[i] == 51){
    bird.sps.count.by.treat$week.order[i] <- 3
  }
  if(bird.sps.count.by.treat$week[i] == 52){
    bird.sps.count.by.treat$week.order[i] <- 4
  }
  if(bird.sps.count.by.treat$week[i] == 53){
    bird.sps.count.by.treat$week.order[i] <- 5
  }
  if(bird.sps.count.by.treat$week[i] == 1){
    bird.sps.count.by.treat$week.order[i] <- 6
  }
  if(bird.sps.count.by.treat$week[i] == 2){
    bird.sps.count.by.treat$week.order[i] <- 7
  }
  if(bird.sps.count.by.treat$week[i] == 3){
    bird.sps.count.by.treat$week.order[i] <- 8
  }
  if(bird.sps.count.by.treat$week[i] == 4){
    bird.sps.count.by.treat$week.order[i] <- 9
  }
  if(bird.sps.count.by.treat$week[i] == 5){
    bird.sps.count.by.treat$week.order[i] <- 10
  }
}

# Reorder
bird.sps.count.by.treat <- bird.sps.count.by.treat |> 
  ungroup() |>
  dplyr::select(site, week.order, treatment, species, count)

# Make sure its a factor
bird.sps.count.by.treat$week.order <- as_factor(bird.sps.count.by.treat$week.order)
bird.sps.count.by.treat$week.order <- factor(bird.sps.count.by.treat$week.order, 
                                             levels = c('1','2','3','4','5','6','7','8','9','10'))
# Pivot wider
bird.wd <- bird.sps.count.by.treat |>
  pivot_wider(names_from = 'species', values_from = 'count')

# Missing true zeroes
bird.wd[is.na(bird.wd)] <- 0

## --------------- Create the weekly richness dataframe ------------------------

spec <- bird.wd[,4:13]
pred <- bird.wd[,1:3]

# Convert to binary
spec <- as.matrix((spec>0)+0)

# Bring back together
bird.wd <- cbind(pred,spec)

bird.rich <- bird.wd |>
  mutate(richness = chipping_sparrow_spizella_passerina+
           eastern_phoebe_sayornis_phoebe+
           mourning_dove_zenaida_macroura+
           northern_cardinal_cardinalis_cardinalis+
           gray_catbird_dumetella_carolinensis+
           red_shoulder_hawk_buteo_lineatus+
           pine_warbler_setophaga_pinus+
           sedge_wren_cistothorus_stellaris+
           eastern_screech_owl_megascops_asio+
           barred_owl_strix_varia)

bird.rich$site <- as_factor(bird.rich$site)

bird.rich <- bird.rich[,c(1,2,3,14)]

# Fix missing data
bird.rich <- complete(bird.rich, site, week.order, treatment)

bird.rich[is.na(bird.rich)] <- 0
## --------------- Model weekly richness ---------------------------------------

dev.new()
descdist(bird.rich$richness, discrete = TRUE) # Poisson

bird.rich |>
  filter(treatment == '0') |>
  ungroup() |>
  dplyr::select(richness) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # Poisson

bird.rich |>
  filter(treatment == '4') |>
  ungroup() |>
  dplyr::select(richness) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # Normal

bird.rich |>
  filter(treatment == '8') |>
  ungroup() |>
  dplyr::select(richness) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # Poisson

bird.rich |>
  filter(treatment == '12') |>
  ungroup() |>
  dplyr::select(richness) |>
  as_vector() |>
  descdist(discrete=TRUE, boot=500)  # Poisson

# Treatment distributions are different, which may cause issues with the
# models moving forward

# Poisson with numeric week and interaction
bird.rich$week.order <- as.numeric(bird.rich$week.order)
bird.rich.pois <- glmer(richness ~ treatment * week.order + (1|site), 
                          control = glmerControl(optimizer ="bobyqa"),
                          family = poisson(link = "log"), data = bird.rich)
sim.bird.rich.pois <- simulateResiduals(bird.rich.pois)
plot(sim.bird.rich.pois)
testResiduals(sim.bird.rich.pois)
dev.new()
check_model(bird.rich.pois)
Anova(bird.rich.pois)

## --------------- Check weekly richness assumptions ---------------------------

# Dharma package
testDispersion(bird.rich.pois)
dev.new()
plotResiduals(bird.rich.pois, form = bird.rich$site)
plotResiduals(bird.rich.pois, form = bird.rich$treatment)
testOutliers(bird.rich.pois)
testCategorical(bird.rich.pois, catPred = bird.rich$treatment)
testCategorical(bird.rich.pois, catPred = bird.rich$site)

testUniformity(bird.rich.pois, alternative = c('two.sided'))
testUniformity(bird.rich.pois, alternative = c('less'))
testUniformity(bird.rich.pois, alternative = c('greater'))

# Check overdispersion
E1 <- resid(bird.rich.pois, type = "pearson")
N <- nrow(bird.rich)
p <- length(coef(bird.rich.pois)) + 1
overdispersion <- sum(E1^2) / (N-p)
overdispersion 

# Model dispersion
testDispersion(sim.bird.rich.pois, alternative = c('less'))
testDispersion(sim.bird.rich.pois, alternative = c('greater'))
testDispersion(sim.bird.rich.pois, alternative = c('two.sided'))

# check other assumptions
dev.new()
check_model(bird.rich.pois) # Fat tails collinearity

# Dharma package
sim.bird.rich.pois <- simulateResiduals(bird.rich.pois)
dev.new()
plot(sim.bird.rich.pois) # Significant deviation and adjusted quantile 
testDispersion(sim.bird.rich.pois)
dev.new()
plotResiduals(sim.bird.rich.pois, form = bird.rich$site)
plotResiduals(sim.bird.rich.pois, form = bird.rich$week.order)
testOutliers(sim.bird.rich.pois)
testQuantiles(sim.bird.rich.pois)
testCategorical(sim.bird.rich.pois, catPred = bird.rich$treatment)
testUniformity(sim.bird.rich.pois, alternative = c('two.sided'))
testUniformity(sim.bird.rich.pois, alternative = c('less'))
testUniformity(sim.bird.rich.pois, alternative = c('greater'))

# plot residuals
plot(richness ~ as_factor(treatment), data = bird.rich)
plot(residuals(bird.rich.pois) ~ as_factor(bird.rich$treatment))
abline(a = 0, b = 0, col = "blue", lwd = 2)
hist(residuals(bird.rich.pois))

# plot residuals from prediction
resfit <- resid(bird.rich.pois)
hist(resfit)
plot(bird.rich$richness, resfit, 
     ylab = "Residuals", 
     xlab = "Detections") 
abline(0, 0)     

# look at predictions 
preds.lm <- predict(bird.rich.pois)
par(mfrow = c(1, 2))
plot(richness ~ as_factor(treatment), 
     data = bird.rich)
plot(preds.lm ~ as_factor(bird.rich$treatment))

# look at predictions from model
# neg binom
predict(bird.rich.pois)
bird.rich$pred = exp(predict(bird.rich.pois))
bird.rich$predicted = predict(bird.rich.pois)    # save the predicted values
bird.rich$residuals = residuals(bird.rich.pois)  # save the residual values

# quick look at the actual, predicted, and residual values
pred_df <- bird.rich %>% 
  dplyr::select(richness, predicted, residuals)
pred_df$predicted = exp(pred_df$predicted)
pred_df$residuals = exp(pred_df$residuals)

plot(fitted(bird.rich.pois) ~ bird.rich$richness)
abline(0, 1, col = "blue", lwd = 2)

plot(predict(bird.rich.pois), residuals(bird.rich.pois, type = 'working'))
plot(cooks.distance(bird.rich.pois), type='h')

# compare predicted vs. actual
plot(density(bird.rich$richness), xlim=c(0, 160), ylim=c(0, .08), main='M_0 y_hat')
lines(density(predict(bird.rich.pois, type='response')), col='red')

# independence of observations
scatter.smooth(bird.rich$site, bird.rich$residuals)

# Check Cook's distance
library('influence.ME')
infl <- influence(bird.rich.pois, obs = TRUE)
plot(infl, which = "cook", cutoff = 4/200)
# lots of influential values

infl <- influence(bird.rich.pois, group = 'site')
plot(infl, which = "cook", cutoff = 4/10)
# 4 influential blocks

## --------------- Export coefficient estimates --------------------------------

# Calculate r2
r.squaredGLMM(bird.rich.pois)
# 0.4032949-0.04921284 = 0.3540821
# 0.4921386-0.06005416 = 0.4320844
# 0.2933441-0.03579588 = 0.2575482

# Calculate variance inflation factors
vif(bird.rich.pois)
# treatment 2.355904
# week since acclimatization 2.344432
# interaction 2.651140

# Calculate means averaged over entire sampling periods
emmeans(bird.rich.pois, pairwise ~treatment, type = 'response',
        adjust = "holm") # Correction
confint(pairs(emmeans(bird.rich.pois, ~treatment)), 
        type = "response", adjust = 'holm') # Confidence intervals
emmeans(bird.rich.pois, pairwise ~treatment, type = 'response',
        adjust = "none") # No correction
confint(pairs(emmeans(bird.rich.pois, ~treatment)), 
        type = "response", adjust = 'none') # Confidence intervals

means <- as.data.frame(emmeans(bird.rich.pois, pairwise ~treatment, type = 'response',
                               adjust = "none")[1])
write.csv(means, "Model-output/Experiment-2-weekly-bird-rich-emmeans.csv", row.names = FALSE)

# Calculate interactive effects between treatment and week
emtrends(bird.rich.pois, pairwise~treatment, var = "week.order",
         adjust = "holm", type = 'response')
confint(pairs(emtrends(bird.rich.pois, pairwise~treatment, var = "week.order",
         adjust = "holm", type = 'response')))
emtrends(bird.rich.pois, pairwise~treatment, var = "week.order",
         adjust = "none", type = 'response') # No adjustments
confint(emtrends(bird.rich.pois, pairwise~treatment, var = "week.order",
         adjust = "none", type = 'response')) # Removed pairs for no correction

trends <- tibble(Treatment = c('Control', 'Low', 'Medium', 'High'),
                 Estimate = c(0.0132, 0.0903, 0.0709, 0.0981),
                 LCL = c(-0.06855, 0.02032, 0.00621, 0.03120),
                 UCL = c(0.0949, 0.1602, 0.1355, 0.1650))

write.csv(trends, "Model-output/Experiment-2-weekly-bird-rich-emtrends.csv", row.names = FALSE)

# Export the interaction graphs for supplemental information
cols <- c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC")

# emmip assigned treatment to a column called 'tvar'
colnames(bird.rich)[3] <- 'tvar'

bird.rich.interact <- emmip(bird.rich.pois, treatment ~ week.order, cov.reduce = range,
      type = 'response')+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10))+
  scale_y_continuous(breaks = seq(0, 5, 1))+
  scale_color_manual(values = c("darkgray",  "#00A9FF", "#00BF7D",  "#FF61CC"))+
  geom_jitter(data = bird.rich, aes(x = week.order, y = richness, color = tvar),
              shape = 20, height = 0.25, width = 0.25, alpha = 0.3, size = 2.5)+
  ylab('Bird Richness')+
  xlab('Weeks since acclimitization period')+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "right",
        axis.text = element_text(face="bold"),
        panel.grid.minor = element_blank())

saveRDS(bird.rich.interact, file = "Model-output/weekly-mean-bird-rich-interact.RDS")
# ggsave('Figures/Experiment-2-weekly-bird-rich-interaction.png')

ggplot(bird.rich, aes(x = week.order, y = richness, color = treatment))+
  geom_jitter(size = 0.5, )+
  scale_color_manual(values = c("darkgray",  "#00A9FF", "#00BF7D",  "#FF61CC"))+
  ylab('Bird Richness')+
  xlab('Weeks since acclimitization period')+
  emmip(bird.rich.pois, treatment ~ week.order, cov.reduce = range,
        type = 'response')+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "right")+
  theme(axis.text = element_text(face="bold"))


