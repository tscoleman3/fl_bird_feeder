## --------------- HEADER ------------------------------------------------------
## Script name: 4_Experiment-2-seed-rich.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliaton: University of Florida
## Date Created: 2022-10-15
## Date Last Modified: 2022-10-15
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script analyzes the seed collection data for 
## experiment 2

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

# Clear the deck
rm(list = ls())

# Bring in the data
seeds <- read.csv(file = "data/clean_seed_traps.csv",
                  header = TRUE, stringsAsFactors = FALSE)


comm.mat <- seeds[,5:31]
predictors <- seeds[,1:4]

comm.mat.pa <- decostand(comm.mat, "pa")
seeds <- cbind(predictors[,1:4], comm.mat.pa)

## --------------- Calculate richness ------------------------------------------

seeds$RICH <- NA
for(i in 1:nrow(seeds)){
  seeds$RICH[i] <- sum(seeds[i,5:31])
}  

## --------------- Calculate days between samples ------------------------------

seeds$DATE <- as_date(seeds$DATE)

seeds$SAMP.LENGTH <- NA
seeds$EXP.DAYS <- NA

for (i in 1:nrow(seeds)){
  if (seeds$DATE[i] == '2020-11-27'){
      seeds$SAMP.LENGTH[i] <- difftime('2020-11-27', '2020-11-20',  unit = c("days"))
      seeds$EXP.DAYS[i] <- difftime('2020-11-27', '2020-11-20',  unit = c("days"))
  } 
  if (seeds$DATE[i] == '2020-12-04'){
      seeds$SAMP.LENGTH[i] <- difftime('2020-12-04', '2020-11-27',  unit = c("days"))
      seeds$EXP.DAYS[i] <- difftime('2020-12-04', '2020-11-20',  unit = c("days"))
  } 
  if (seeds$DATE[i] == '2020-12-11'){
      seeds$SAMP.LENGTH[i] <- difftime('2020-12-11', '2020-12-04',  unit = c("days"))
      seeds$EXP.DAYS[i] <- difftime('2020-12-11', '2020-11-20',  unit = c("days"))
  } 
  if (seeds$DATE[i] == '2021-01-02'){
      seeds$SAMP.LENGTH[i] <- difftime('2021-01-02', '2020-12-11',  unit = c("days"))
      seeds$EXP.DAYS[i] <- difftime('2021-01-02', '2020-11-20',  unit = c("days"))
  } 
  if (seeds$DATE[i] == '2021-02-03'){
    seeds$SAMP.LENGTH[i] <- difftime('2021-02-03', '2021-01-02',  unit = c("days"))
    seeds$EXP.DAYS[i] <- difftime('2021-02-03', '2020-11-20',  unit = c("days"))
    
  } 
}

## --------------- Convert predictors to ordinal -------------------------------

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

library('forcats')
seeds$BLOCK <- forcats::as_factor(seeds$BLOCK)
seeds$TRAP <- forcats::as_factor(seeds$TRAP)

seeds$TREATMENT <- factor(seeds$TREATMENT, order = TRUE,
                             levels = c('Control', 'Low', 'Medium', 'High'))

## --------------- Richness model ----------------------------------------------

seeds <- seeds %>%
  dplyr::select(DATE.ORD, SAMP.LENGTH, EXP.DAYS, 
                BLOCK, TRAP, TREATMENT, TREATMENT.NUMB, RICH)

# Visualize the data ####
boxplot(RICH ~ TREATMENT,col=c("white","lightgray", 
                                    "darkgray", "black"), outline = TRUE, seeds)
dotplot(seeds$RICH) 

# Plot distribution
plotdist(seeds$RICH, histo = TRUE, demp = TRUE)
descdist(seeds$RICH, discrete=TRUE, boot=500) # NB/poisson

seeds.rich.nb <- glmer.nb(RICH ~ TREATMENT * DATE.ORD + (1|BLOCK), data = seeds) 

library('optimx')

# seeds.rich.nb <- glmer(RICH ~ TREATMENT * EXP.DAYS + (1|BLOCK),
                          # family = 'poisson', data = seeds) 

# Predictions are funnel shaped when using treatment as an ordinal variable
# Model fail assumptions using poisson with factored treatment

sim.seeds.rich.nb <- simulateResiduals(seeds.rich.nb)
plot(sim.seeds.rich.nb)

summary(seeds.rich.nb)
car::Anova(seeds.rich.nb)
r.squaredGLMM(seeds.rich.nb)

## --------------- Check model assumptions -------------------------------------

# Check overdispersion
E1 <- resid(seeds.rich.nb, type = "pearson")
N <- nrow(seeds)
p <- length(fixef(seeds.rich.nb)) + 1
overdispersion <- sum(E1^2) / (N-p)
overdispersion 

# check some assumptions
plot(RICH ~ TREATMENT, data = seeds)
plot(residuals(seeds.rich.nb) ~ seeds$TREATMENT)
abline(a = 0, b = 0, col = "blue", lwd = 2)
hist(residuals(seeds.rich.nb))

# plot residuals from prediction
resfit <- resid(seeds.rich.nb)
hist(resfit)
plot(seeds$TREATMENT, resfit, 
     ylab = "Residuals", 
     xlab = "Detections") 
abline(0, 0)  

# look at predictions 
preds.nb <- predict(seeds.rich.nb)
par(mfrow = c(1, 2))
plot(RICH ~ TREATMENT, 
     data = seeds)
plot(preds.nb ~ seeds$TREATMENT) # funnel is problem

# look at predictions from model
# neg binom
predict(seeds.rich.nb)
seeds$pred = exp(predict(seeds.rich.nb))
seeds$predicted = predict(seeds.rich.nb)    # save the predicted values
seeds$residuals = residuals(seeds.rich.nb)  # save the residual values

# quick look at the actual, predicted, and residual values
pred_df <- seeds %>% 
  dplyr::select(RICH, predicted, residuals)
pred_df$predicted = exp(pred_df$predicted)
pred_df$residuals = exp(pred_df$residuals)

plot(fitted(seeds.rich.nb) ~ seeds$RICH)
abline(0, 1, col = "blue", lwd = 2)

library('performance')
check_zeroinflation(seeds.rich.nb) 

## --------------- Export model validation -------------------------------------
E1 <- resid(seeds.rich.nb, type = 'pearson')
F1 <- fitted(seeds.rich.nb, type = 'response')

png("Figures/Experiment-2-seed-rich-model-validation.png", 
    width = 800, height = 800, units = "px")

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

emmeans(seeds.rich.nb, pairwise ~TREATMENT, type = 'response',
        adjust = "none")

emtrends(seeds.rich.nb, pairwise ~ TREATMENT, var = "DATE.ORD", type = 'response')
# emtrends(seeds.rich.nb, pairwise ~ TREATMENT, var = "EXP.DAYS")

cols <- c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC")
emmip(seeds.rich.nb, TREATMENT ~ DATE.ORD, 
      cov.reduce = range, CIs = TRUE, type = 'response',
      style = 'factor')+
  ylab(expression(paste("Fitted mean seed richness", " (m"^{2}, ")")))+
  xlab('')+
  scale_x_discrete(labels = c("Sample 1", "Sample 5"))+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "none")+
  theme(axis.text = element_text(face="bold"))+
  scale_color_manual(values=cols)

ggsave("Figures/Experiment-2-seed-rich.png", width = 5, height = 7, units = "in")

