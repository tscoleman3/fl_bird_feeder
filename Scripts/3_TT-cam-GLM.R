## --------------- HEADER ------------------------------------------------------
## Script name: 3a_TT-detection-model.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliaton: University of Florida
## Date Created: 2022-08-11
## Date Last Modified: 2023-07-26
## Copyright (c) David S. Mason, 2023
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script conducts a GLM on annual avian detections
## for the Tall Timbers experiment

## --------------- SET—UP WORKSPACE --------------------------------------------
library(tidyverse)
library(lubridate)
library(tidylog)
library(styler)
library(fitdistrplus)
library(glmmTMB)
library(DHARMa)
library(performance)
library(emmeans)
library(car)
library(MuMIn) # r2
library(lme4) # VarCorr

# Wipe the decks
rm(list = ls())

# Bring in the data
TT.cam <- read.csv("Birds/Clean-data/1_TT-annual.csv")
coverage <- read.csv("Birds/Clean-data/Coverage.csv")
range(TT.cam$Detections)
head(TT.cam)

## --------------- PREPARE THE DATA --------------------------------------------

# Concatenate the treatment and trap information into a single ID for the
# paired analysis
TT.cam$ID <- paste(TT.cam$Site, TT.cam$Treatment, TT.cam$Trap)

# Standardize terminology for burn year
colnames(TT.cam)[2] <- 'Burn.year'

# Create a value for unit (paired clusters of three seed traps)
TT.cam$Unit <- c(rep('A', 6), rep('B', 6), rep('C', 6), rep('D', 6), rep('E', 6),
								rep('F', 6), rep('G', 6), rep('H', 6), rep('I', 6), rep('J', 6),
								rep('K', 6), rep('L', 6))

# Make sure factors are factors
TT.cam$Site <- as_factor(TT.cam$Site)
TT.cam$Burn.Year <- as_factor(TT.cam$Burn.year)
TT.cam$ID <- as_factor(TT.cam$ID)
TT.cam$Trap <- as_factor(TT.cam$Trap)
TT.cam$Unit <- as_factor(TT.cam$Unit)
TT.cam$Burn <- as_factor(TT.cam$Burn)

# Merge coverage
colnames(coverage)[2] <- 'Burn.year'
TT.cam <- merge(TT.cam, coverage)
rm(coverage)

TT.cam <- TT.cam |>
	mutate(Coverage = Potential.days-Total.missing.days)

## --------------- DETERMINE DISTRIBUTION --------------------------------------

descdist(TT.cam$Detections, boot = 1000, discrete = TRUE) # negative binomial

# Explore mean-variance relationship to assess which negative binomial
# distribution to use
mean.var <- TT.cam |>
	group_by(Site, Unit, Burn) |>
	dplyr::summarize(Mean.det = mean(Detections),
						Var.det = var(Detections))
q1 <- qplot(Mean.det,Var.det,data=mean.var)

print(q1+
	## linear (quasi-Poisson/NB1) fit
	geom_smooth(method="lm",formula=y~x-1)+
	## smooth (loess)
	geom_smooth(colour="red")+
	## semi-quadratic (NB2/LNP)
	geom_smooth(method="lm",formula=y~I(x^2)+offset(x)-1,colour="purple")+
	## Poisson (v=m)
	geom_abline(a=0,b=1,lty=2))
# blue = nb1
# purple = nb2

# Hard to tell which value fits better

## --------------- CONDUCT INITIAL MODEL ---------------------------------------

# Run both versions of negative binomial and compare outputs
m1 <- glmmTMB(Detections ~ Burn + (1|Site/Unit/ID) + offset(log(Coverage)), data = TT.cam,
						 family = 'nbinom1') # variance = µ * phi
m2 <- glmmTMB(Detections ~ Burn + (1|Site/Unit/ID) + offset(log(Coverage)), data = TT.cam,
						family = 'nbinom2') # variance = µ(1+µ/k)
check_model(m1)
check_model(m2)
# Residual variance looks better with variance = µ(1+µ/k) (nbinom2). Without it
# there are diverging trends in variance at the highest values.

summary(m2)
Anova(m2)
r.squaredGLMM(m2) # fixed effects explain 3-4% of variance, random = 52-63%
# throws up a warning.
check_overdispersion(m2) # No overdispersion detected
check_singularity(m2) # False

# Compare to null
null <- glmmTMB(Detections ~ 1 + (1|Site/Unit/ID) + offset(log(Coverage)), data = TT.cam,
						family = 'nbinom2')
anova(null, m2) # 0.023

sim.m2 <- simulateResiduals(m2)
plot(sim.m2) # ns
# plotResiduals(sim.m2, TT.cam$Coverage) # ns
plotResiduals(sim.m2, TT.cam$Burn) # ns
# plotResiduals(sim.m2, TT.cam$Burn.year) # ns
plotResiduals(sim.m2, TT.cam$Site) # ns
testDispersion(sim.m2) # ns
testQuantiles(sim.m2)
testUniformity(sim.m2, alternative = c('two.sided'))

# Plot residuals
resfit <- resid(m2, type = 'pearson')
plot(TT.cam$Burn, resfit, 
     ylab = "Residuals", 
     xlab = "Treatment") 
abline(0, 0) 

# delete
plot(TT.cam$Burn.year, resfit, 
     ylab = "Residuals", 
     xlab = "Treatment") 
abline(0, 0)  

plot(TT.cam$Site, resfit, 
     ylab = "Residuals", 
     xlab = "Treatment") 
abline(0, 0)  

plot(TT.cam$ID, resfit, 
     ylab = "Residuals", 
     xlab = "Treatment") 
abline(0, 0)  

plot(density(resid(m2, type='pearson')))

# look at the predictions
predict(m2)
TT.cam$pred = predict(m2)
TT.cam$residuals = residuals(m2, type = 'pearson')  # save the residual values

pred_df <- TT.cam %>% 
  dplyr::select(Detections, pred, residuals)
plot(fitted(m2) ~ TT.cam$Detections)
abline(0, 1, col = "blue", lwd = 2)  
plot(predict(m2), residuals(m2, type = 'working'))

# Compare predicted vs. actual
plot(density(TT.cam$Detections), xlim=c(0, 160), ylim=c(0, .08), main='M_0 y_hat')
lines(density(predict(m2, type='response')), col='red')

# Check independence of observations
scatter.smooth(TT.cam$Site, pred_df$residuals)				

# Check for influential observations
source(system.file("other_methods","influence_mixed.R", package="glmmTMB"))
m2_influence <- influence_mixed(m2, groups="Site")
4/6 # 0.6666667
car::infIndexPlot(m2_influence)

m2_influence <- influence_mixed(m2, groups="Unit")
4/12 # 0.3333333
car::infIndexPlot(m2_influence)

m2_influence <- influence_mixed(m2, groups="ID")
4/36 # 0.1111111
car::infIndexPlot(m2_influence)

m2_influence <- influence_mixed(m2, groups=".case")
4/72 # 0.05555556
car::infIndexPlot(m2_influence)

## --------------- MODEL OUTPUTS -----------------------------------------------

if (requireNamespace("car") && getRversion() >= "3.6.0") {
car::Anova(m2) ## default type II
}

# Model output
Anova(m2)
model.output <- tibble('Fixed effect' = c('Burn'),
											 Chisq = c(5.4332),
											 df = c(1),
											 pval = c('0.020'))
write.csv(model.output, "Birds/Results/Detection-mod-output.csv",
					row.names = FALSE)

total.sd <- sqrt(0.34286^2 + 0.22907^2 + 0.67002^2)
detection.means <- as.data.frame(emmeans(m2, ~Burn, type = "response",
												bias.adjust = TRUE, sigma = total.sd))
write.csv(detection.means, "Birds/Results/Detection-means.csv",
					row.names = FALSE)

# Means comp
emmeans(m2, revpairwise ~ Burn, type = "response",
												bias.adjust = TRUE, sigma = total.sd)
emmean.output <- tibble(Comparison = c('Contrast'),
										 Ratio = c(1.97),
										 se = c(0.345),
										 df = (66),
										 null = (1.31),
										 t.ratio = c(2.331),
										 pval = c('0.023'))
write.csv(emmean.output, "Birds/Results/Detection-means-test.csv",
					row.names = FALSE)

## --------------- CONDUCT MODEL EXPLORING BURN YEAR ---------------------------

rm(list=ls()[! ls() %in% c("TT.cam")])

# Run both versions of negative binomial and compare outputs
m1 <- glmmTMB(Detections ~ Burn * Burn.year + (1|Site/Unit/ID) + offset(log(Coverage)), data = TT.cam,
						 family = 'nbinom1') # variance = µ * phi
m2 <- glmmTMB(Detections ~ Burn * Burn.year + (1|Site/Unit/ID) + offset(log(Coverage)), data = TT.cam,
						 family = 'nbinom2') # variance = µ(1+µ/k)
check_model(m1)
check_model(m2)
# Residual variance looks better with variance = µ(1+µ/k) (nbinom2). Without it
# there are diverging trends in variance at the highest values.

summary(m2)
Anova(m2)
r.squaredGLMM(m2) # fixed effects explain 3-4% of variance, random = 52-63%
# throws up a warning.
check_overdispersion(m2) # No overdispersion detected
check_singularity(m2) # False

# Compare to null
null <- glmmTMB(Detections ~ 1 + (1|Site/Unit/ID) + offset(log(Coverage)), data = TT.cam,
						family = 'nbinom2')
anova(null, m2)

sim.m2 <- simulateResiduals(m2)
plot(sim.m2) # ns
# plotResiduals(sim.m2, TT.cam$Coverage) # ns
plotResiduals(sim.m2, TT.cam$Burn) # ns
# plotResiduals(sim.m2, TT.cam$Burn.year) # ns
plotResiduals(sim.m2, TT.cam$Site) # ns
testDispersion(sim.m2) # ns
testQuantiles(sim.m2)
testUniformity(sim.m2, alternative = c('two.sided'))

# Plot residuals
resfit <- resid(m2, type = 'pearson')
plot(TT.cam$Burn, resfit, 
     ylab = "Residuals", 
     xlab = "Treatment") 
abline(0, 0) 

plot(TT.cam$Site, resfit, 
     ylab = "Residuals", 
     xlab = "Treatment") 
abline(0, 0)  

plot(TT.cam$ID, resfit, 
     ylab = "Residuals", 
     xlab = "Treatment") 
abline(0, 0)  

plot(density(resid(m2, type='pearson')))

# look at the predictions
predict(m2)
TT.cam$pred = predict(m2)
TT.cam$residuals = residuals(m2, type = 'pearson')  # save the residual values

pred_df <- TT.cam %>% 
  dplyr::select(Detections, pred, residuals)
plot(fitted(m2) ~ TT.cam$Detections)
abline(0, 1, col = "blue", lwd = 2)  
plot(predict(m2), residuals(m2, type = 'working'))

# Compare predicted vs. actual
plot(density(TT.cam$Detections), xlim=c(0, 160), ylim=c(0, .08), main='M_0 y_hat')
lines(density(predict(m2, type='response')), col='red')

# Check independence of observations
scatter.smooth(TT.cam$Site, pred_df$residuals)				

# Check for influential observations
source(system.file("other_methods","influence_mixed.R", package="glmmTMB"))
m2_influence <- influence_mixed(m2, groups="Site")
4/6 # 0.6666667
car::infIndexPlot(m2_influence)

m2_influence <- influence_mixed(m2, groups="Unit")
4/12 # 0.3333333
car::infIndexPlot(m2_influence)

m2_influence <- influence_mixed(m2, groups="ID")
4/36 # 0.1111111
car::infIndexPlot(m2_influence)

m2_influence <- influence_mixed(m2, groups=".case")
4/72 # 0.05555556
car::infIndexPlot(m2_influence)

## --------------- MODEL OUTPUTS -----------------------------------------------

if (requireNamespace("car") && getRversion() >= "3.6.0") {
car::Anova(m2, type = 2) ## default type II
}

# Model output
Anova(m2)
model.output <- tibble('Fixed effect' = c('Burn', 'Burn year', 'Interaction'),
											 Chisq = c(6.5067, 0.7048, 12.5432),
											 df = c(1, 1, 1),
											 pval = c('0.011', '0.401', '<0.001'))
write.csv(model.output, "Birds/Results/Detection-interaction-mod-output.csv",
					row.names = FALSE)

# Means values
total.sd <- sqrt(0.41486^2 + 0.28185^2 + 0.65937^2)
detection.means <- as.data.frame(emmeans(m1, ~ Burn*Burn.year, type = "response",
																bias.adjust = TRUE, sigma = total.sd))
write.csv(detection.means, "Birds/Results/Detection-interaction-means.csv",
					row.names = FALSE)

# Means comp

emmeans(m2, revpairwise ~ Burn|Burn.year, type = "response",
				bias.adjust = FALSE, sigma = total.sd)
emmean.output <- tibble(Comparison = c('Year 1 contrast', 'Year 2 contrast'),
										 Ratio = c(1.13, 3.29),
										 se = c(0.246, 0.689),
										 df = c(64, 64),
										 null = c(1.34, 1.34),
										 t.ratio = c(-0.780, 4.294),
										 pval = c('0.438', '<0.001'))

write.csv(emmean.output, "Birds/Results/Detection-interaction-means-test.csv",
					row.names = FALSE)

