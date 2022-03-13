##### JIMMY HOLDGRAFER ###########
##### SEED DISPERSAL PROJECT #####

# H01: as food diversity increases, bird and seed richness increase
# H02: as food diversity increases, bird and seed observations increase


##### PACKAGES #####
####################
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


###### DATA ########
####################
bird.dat.og <- read.csv(file = "data/bird_data2.0.csv",
                   header = TRUE)
head(bird.dat.og)

seed.dat.og <- read.csv(file = "data/seed_traps2.0.csv",
                   header = TRUE)
head(seed.dat.og)


### ACCLIM PER #####
####################

# Birds
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

# Seeds
seed.dat <- pivot_longer(seed.dat.og, cols = 5:35, names_to = "SPECIES",
												 values_to = "SEEDS")

seed.dat$DATE <- as.Date(seed.dat$DATE)

seed.dat <- seed.dat %>% 
	mutate(week = week(DATE)) %>% 
	filter(week != 48 & week != 47)


## SUMMARIZE DATA ##
####################
#### SUMMARIZE #####
####################
# Birds
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

# Seeds
seed.obs <- seed.dat %>% 
  dplyr::group_by(TREATMENT, BLOCK) %>% 
  dplyr::summarise(sum(SEEDS)) 
seed.obs <- dplyr::rename(seed.obs, SEEDS = "sum(SEEDS)")

# total count of species by treatment and site
seed.species.count <- seed.dat %>% 
  group_by(TREATMENT, BLOCK, SPECIES) %>% 
  summarise(sum(SEEDS)) 

seed.species.count <- dplyr::rename(seed.species.count, SEEDS = "sum(SEEDS)")

seed.species.count$TREATMENT <- factor(seed.species.count$TREATMENT, 
																		levels=c('Control', 'Low', 'Medium', 'High'))

seed.species.count <- seed.species.count %>% 
	filter(SEEDS > 0)

seed.rich <- as.data.frame(cbind(sites,treatment))
seed.rich$richness <- c(1,2,5,2,0,0,0,2,1,0,0,5,0,3,4,3,0,2,2,3,
												1,1,2,2,1,1,1,1,0,1,1,3,0,0,4,4,2,0,0,0)
seed.rich <- as.data.frame(seed.rich) %>% 
	dplyr::select(treatment,sites,richness)

####################
### DISTRIBUTION ###
####################
plotdist(bird.obs$count, histo = TRUE, demp = TRUE)
descdist(bird.obs$count, discrete=TRUE, boot=500) # NB or poisson

plotdist(bird.rich$richness, histo = TRUE, demp = TRUE)
descdist(bird.rich$richness, discrete=TRUE, boot=500) # poisson

plotdist(seed.obs$SEEDS, histo = TRUE, demp = TRUE)
descdist(seed.obs$SEEDS, discrete=TRUE, boot=500) # poisson

plotdist(seed.rich$richness, histo = TRUE, demp = TRUE)
descdist(seed.rich$richness, discrete=TRUE, boot=500) # poisson


####################
## BIRD OBS MODELS #
####################
# Bird observations
# we want to look at the effect of treatment on bird count
# while accounting for site as a random effect
# neg binom fit
birds.obs.fit.nb <- glmer.nb(count ~ treatment + (1 | site),
                             data = bird.obs)
# poisson fit
birds.obs.fit.pois <- glmer(count ~ treatment + (1 | site),
                            data = bird.obs, family = poisson(link = "log"))
glmmPQL(count ~ treatment, random = ~1|site,
        data = bird.obs, family = poisson(link = "log"))
# birds.obs.fit.pois <- glm(count ~ treatment * site,
#                             data = bird.obs)

# check some assumptions
plot(count ~ as.numeric(treatment), 
     data = bird.obs)
plot(residuals(birds.obs.fit.nb) ~ as.numeric(bird.obs$treatment))
abline(a = 0, b = 0, col = "blue", lwd = 2)
plot(residuals(birds.obs.fit.pois) ~ as.numeric(bird.obs$treatment))
abline(a = 0, b = 0, col = "blue", lwd = 2)
hist(residuals(birds.obs.fit.nb))
hist(residuals(birds.obs.fit.pois))
# plot residuals from prediction to test for assumptions
resfit <- resid(birds.obs.fit.nb)
resfit <- resid(birds.obs.fit.pois)
hist(resfit)
plot(bird.obs$count, resfit, 
     ylab = "Residuals", 
     xlab = "Count") 
abline(0, 0)        
#' we meet the following assumptions based on our residual plots:
#' - normality
#' - linearity
#' - homoscedasticity
#' - no autocorrelation 
#' perfect

# look at predictions 
preds.nb <- predict(birds.obs.fit.nb)
preds.pois <- predict(birds.obs.fit.pois)
par(mfrow = c(2, 2))
plot(count ~ as.numeric(treatment), 
     data = bird.obs)
plot(preds.nb ~ as.numeric(bird.obs$treatment))
plot(count ~ as.numeric(treatment), 
     data = bird.obs)
plot(preds.pois ~ as.numeric(bird.obs$treatment))

# look at predictions from model
# neg binom
predict(birds.obs.fit.nb)
bird.obs$pred = exp(predict(birds.obs.fit.nb))
bird.obs$predicted = predict(birds.obs.fit.nb)    # save the predicted values
bird.obs$residuals = residuals(birds.obs.fit.nb)  # save the residual values
# poisson
predict(birds.obs.fit.pois)
bird.obs$pred = exp(predict(birds.obs.fit.pois))
bird.obs$predicted = predict(birds.obs.fit.pois)    # save the predicted values
bird.obs$residuals = residuals(birds.obs.fit.pois)  # save the residual values
# quick look at the actual, predicted, and residual values
pred_df <- bird.obs %>% 
  dplyr::select(count, predicted, residuals)
pred_df$predicted = exp(pred_df$predicted)
pred_df$residuals = exp(pred_df$residuals)
par(mfrow = c(1, 2))
plot(fitted(birds.obs.fit.pois) ~ bird.obs$count)
abline(0, 1, col = "blue", lwd = 2)
plot(fitted(birds.obs.fit.nb) ~ bird.obs$count)
abline(0, 1, col = "blue", lwd = 2)

# look at our coefs 
summary(birds.obs.fit.pois)
anova(birds.obs.fit.pois)
coefs <- summary(birds.obs.fit.pois)$coef
coefs_est <- exp(coefs[, "Estimate"])
uprs <- exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"])
lwrs <- exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"])
(uprs - 1) * 100       # upr CI %'s
(coefs_est - 1) * 100  # coefficient estimates CI %'s
(lwrs - 1) * 100       # lwr CI %'s

par(mfrow = c(1, 1))
plot(exp(coefs[1:4,"Estimate"]), ylim = range(lwrs[1:4], uprs[1:4]))
segments(1:7, lwrs[1:4], 1:7, uprs[1:4])

coefs.fix <- fixef(birds.obs.fit.pois)
exp(coefs.fix[1])
exp(coefs.fix[1]) * exp(coefs.fix[2])
exp(coefs.fix[1]) * exp(coefs.fix[3])
exp(coefs.fix[1]) * exp(coefs.fix[4])

coefs_func <- function(.) {
  beta <- unname(fixef(.))
  
  control <- exp(beta[1])             # mean count for control 
  treat1 <-  exp(beta[1] + beta[2])   # mean count for treat 1 is this much greater than control
  treat2 <-  exp(beta[1] + beta[3])   # mean count for treat 2 is this much greater than control
  treat3 <-  exp(beta[1] + beta[4])   # mean count for treat 3 is this much greater than control
  
  c(control_mu = control, treat1 = treat1, treat2 = treat2, treat3 = treat3)
  
}

start <- Sys.time()
rand <- bootMer(x = birds.obs.fit.pois, FUN = coefs_func, nsim = 1000)$t
Sys.time() - start

summ <- apply(rand, 2, function(x) c(mean = mean(x, na.rm = TRUE),
                                     quantile(x, c(0.025, 0.975), na.rm = TRUE)))
summ
####################
# BIRD RICH MODELS #
####################
bird.rich.lm <- lmer(richness~treatment + (1|sites), data = bird.rich)
anova(bird.rich.lm)
qqnorm(resid(bird.rich.lm))

Betas <- fixef(bird.rich.lm)
SE <- sqrt(diag(vcov(bird.rich.lm)))
pval <- 2 * pnorm(-abs(Betas / SE))
Output <- cbind(Betas, SE, pval)
print(Output, digits = 3)

Resid <- resid(bird.rich.lm)
Fitted <- fitted(bird.rich.lm)

# Checking values given by fitted
Betas <- fixef(bird.rich.lm)
X <- model.matrix(bird.rich.lm)
FitManual <- X %*% Betas
RE <- ranef(bird.rich.lm)$sites$'(Intercept)'
ALLRE <- RE[as.numeric(bird.rich$sites)]
FitManual + ALLRE - fitted(bird.rich.lm)

# Checking assumptions
# Zurr says these residuals should be normalized but the objects
# called in Sigmas <- as.numeric(summary(bird.rich.lm@REmat[,4]))
# do not exist

plot(richness ~ as.numeric(treatment), 
     data = bird.rich)
plot(residuals(bird.rich.lm) ~ as.numeric(bird.rich$treatment))
abline(a = 0, b = 0, col = "blue", lwd = 2)
hist(residuals(bird.rich.lm))

resfit <- resid(bird.rich.lm)
hist(resfit)
plot(bird.rich$rich, resfit, 
     ylab = "Residuals", 
     xlab = "Richness") 
abline(0, 0) # better than poisson model but not perfect

# LRT
bird.rich.lm <- lmer(richness~treatment + (1|sites), 
                     data = bird.rich, REML = TRUE)
bird.rich.lm.null <- update(bird.rich.lm, .~. - treatment)
anova(bird.rich.lm,bird.rich.lm.null)

####################
# SEED OBS MODELS #
####################
# neg binom fit
seeds.obs.fit.nb <- glmer.nb(SEEDS ~ TREATMENT + (1 | BLOCK),
                             data = seed.obs)
# poisson fit
seed.obs.fit.pois <- glmer(SEEDS ~ TREATMENT + (1 | BLOCK),
                            data = seed.obs, family = poisson(link = "log"))
glmmPQL(SEEDS ~ TREATMENT, random = ~1|BLOCK,
        data = seed.obs, family = poisson(link = "log"))

# check some assumptions
plot(SEEDS ~ as.numeric(TREATMENT), 
     data = seed.obs)
plot(residuals(seeds.obs.fit.nb) ~ as.numeric(seed.obs$TREATMENT))
abline(a = 0, b = 0, col = "blue", lwd = 2)
plot(residuals(seed.obs.fit.pois) ~ as.numeric(seed.obs$TREATMENT))
abline(a = 0, b = 0, col = "blue", lwd = 2)
hist(residuals(seeds.obs.fit.nb))
hist(residuals(seed.obs.fit.pois))

# plot residuals from prediction to test for assumptions
resfit <- resid(seeds.obs.fit.nb)
resfit <- resid(seed.obs.fit.pois)
hist(resfit)
plot(seed.obs$SEEDS, resfit, 
     ylab = "Residuals", 
     xlab = "Count") 
abline(0, 0)        

#' we do not meet the following assumptions based on our residual plots:
#' - normality N
#' - linearity Y
#' - homoscedasticity N
#' - no autocorrelation N

# Trying a zero-inflation model
library("performance")
check_zeroinflation(seeds.obs.fit.nb) # error
check_zeroinflation(seed.obs.fit.pois) # probably zero-inflation

# Zero-inflation fit
library(glmmTMB)
seed.obs.fit.pois.ZI <- glmmTMB(SEEDS ~ TREATMENT + (1 | BLOCK),
                                family = poisson,
                                zi=~TREATMENT,
                                data=seed.obs)
  
summary(seed.obs.fit.pois.ZI)

# Looking at simulated residuals
sim.seed.obs.fit.pois.ZI <- simulateResiduals(fittedModel = seed.obs.fit.pois.ZI,
                                              n = 250)
plot(sim.seed.obs.fit.pois.ZI)

# check some assumptions
plot(SEEDS ~ as.numeric(TREATMENT), 
     data = seed.obs)
plot(residuals(seed.obs.fit.pois.ZI) ~ as.numeric(seed.obs$TREATMENT))
abline(a = 0, b = 0, col = "blue", lwd = 2)
hist(residuals(seed.obs.fit.pois.ZI))

# plot residuals from prediction to test for assumptions
resfit <- resid(seed.obs.fit.pois.ZI)
hist(resfit)
plot(seed.obs$SEEDS, resfit, 
     ylab = "Residuals", 
     xlab = "Count") 
abline(0, 0)        

#' we do not meet the following assumptions based on our residual plots:
#' - normality N but better
#' - linearity Y
#' - homoscedasticity N but better
#' - no autocorrelation N but better

# look at predictions 
preds.pois <- predict(seed.obs.fit.pois.ZI)
par(mfrow = c(1, 2))
plot(SEEDS ~ as.numeric(TREATMENT), 
     data = seed.obs)
plot(preds.pois ~ as.numeric(seed.obs$TREATMENT))

# look at predictions from model
# poisson
predict(seed.obs.fit.pois.ZI)
seed.obs$pred = exp(predict(seed.obs.fit.pois.ZI))
seed.obs$predicted = predict(seed.obs.fit.pois.ZI)    # save the predicted values
seed.obs$residuals = residuals(seed.obs.fit.pois.ZI)  # save the residual values

# quick look at the actual, predicted, and residual values
pred_df <- seed.obs %>% 
  dplyr::select(SEEDS, predicted, residuals)
pred_df$predicted = exp(pred_df$predicted)
pred_df$residuals = exp(pred_df$residuals)
plot(fitted(seed.obs.fit.pois.ZI) ~ seed.obs$SEEDS)
abline(0, 1, col = "blue", lwd = 2)

# look at our coefs 
summary(seed.obs.fit.pois.ZI)
# anova(seed.obs.fit.pois.ZI) not appropriate for ZI
Anova(seed.obs.fit.pois.ZI, type = "III") # Not the right Anova?

# RUNNING INTO PROBLEMS HERE #

coefs <- summary(seed.obs.fit.pois.ZI)$coef
coefs_est <- exp(coefs[, "Estimate"])
uprs <- exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"])
lwrs <- exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"])
(uprs - 1) * 100       # upr CI %'s
(coefs_est - 1) * 100  # coefficient estimates CI %'s
(lwrs - 1) * 100       # lwr CI %'s

####################
# SEED RICH MODELS #
####################
# neg binom fit
seeds.rich.fit.nb <- glmer.nb(richness ~ treatment + (1 | sites),
                             data = seed.rich) # iteration limit reached
# poisson fit
seeds.rich.fit.pois <- glmer(richness ~ treatment + (1 | sites),
                            data = seed.rich, family = poisson(link = "log"))
glmmPQL(richness ~ treatment, random = ~1|sites,
        data = seed.rich, family = poisson(link = "log"))

# check some assumptions
plot(richness ~ as.numeric(treatment), 
     data = seed.rich)
plot(residuals(seeds.rich.fit.nb) ~ as.numeric(seed.rich$treatment))
abline(a = 0, b = 0, col = "blue", lwd = 2)
plot(residuals(seeds.rich.fit.pois) ~ as.numeric(seed.rich$treatment))
abline(a = 0, b = 0, col = "blue", lwd = 2)
hist(residuals(seeds.rich.fit.nb))
hist(residuals(seeds.rich.fit.pois))
# plot residuals from prediction to test for assumptions
resfit <- resid(seeds.rich.fit.nb)
resfit <- resid(seeds.rich.fit.pois)
hist(resfit)
plot(seed.rich$richness, resfit, 
     ylab = "Residuals", 
     xlab = "Richness") 
abline(0, 0)        
#' we do not meet the following assumptions based on our residual plots:
#' - normality N
#' - linearity Y
#' - homoscedasticity N
#' - no autocorrelation N

# Trying a zero-inflation model
library("performance")
check_zeroinflation(seeds.rich.fit.pois) # probably zero-inflation

library(glmmTMB)
## Poisson with Zero Inflation (note no overdispersion). 
# Is it enough to account for zero inflation without correcting for overdispersion? 
seeds.rich.fit.pois.ZI <- glmmTMB(richness ~ treatment + (1 | sites),
                                family = poisson,
                                zi=~treatment,
                                data=seed.rich)

summary(seeds.rich.fit.pois.ZI)

sim.seeds.rich.fit.pois.ZI <- simulateResiduals(fittedModel = seeds.rich.fit.pois.ZI,
                                              n = 250)
plot(sim.seeds.rich.fit.pois.ZI) # Quantile deviations detected
# we are modeling a unimodal relationship as a linear process

# check some assumptions
plot(richness ~ as.numeric(treatment), 
     data = seed.rich)
plot(residuals(seeds.rich.fit.pois.ZI) ~ as.numeric(seed.rich$treatment))
abline(a = 0, b = 0, col = "blue", lwd = 2)
hist(residuals(seeds.rich.fit.pois.ZI))
# plot residuals from prediction to test for assumptions
resfit <- resid(seeds.rich.fit.pois.ZI)
hist(resfit)
plot(seed.rich$richness, resfit, 
     ylab = "Residuals", 
     xlab = "Count") 
abline(0, 0)        
#' we do not meet the following assumptions based on our residual plots:
#' - normality Y 
#' - linearity Y
#' - homoscedasticity Y 
#' - no autocorrelation N but better

# look at predictions 
preds.pois <- predict(seeds.rich.fit.pois.ZI)
par(mfrow = c(1, 2))
plot(richness ~ as.numeric(treatment), 
     data = seed.rich)
plot(preds.pois ~ as.numeric(seed.rich$treatment))

# look at predictions from model
# poisson
predict(seeds.rich.fit.pois.ZI)
seed.rich$pred = exp(predict(seeds.rich.fit.pois.ZI))
seed.rich$predicted = predict(seeds.rich.fit.pois.ZI)    # save the predicted values
seed.rich$residuals = residuals(seeds.rich.fit.pois.ZI)  # save the residual values

# quick look at the actual, predicted, and residual values
pred_df <- seed.rich %>% 
  dplyr::select(richness, predicted, residuals)
pred_df$predicted = exp(pred_df$predicted)
pred_df$residuals = exp(pred_df$residuals)
plot(fitted(seeds.rich.fit.pois.ZI) ~ seed.rich$richness)
abline(0, 1, col = "blue", lwd = 2)

# look at our coefs 
summary(seeds.rich.fit.pois.ZI)
# anova(seed.obs.fit.pois.ZI) not appropriate for ZI
Anova(seeds.rich.fit.pois.ZI, type = "III") # NEED TO CHECK IF APPROPRIATE

# RUNNING INTO PROBLEMS HERE #

coefs <- summary(seed.obs.fit.pois.ZI)$coef
coefs_est <- exp(coefs[, "Estimate"])
uprs <- exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"])
lwrs <- exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"])
(uprs - 1) * 100       # upr CI %'s
(coefs_est - 1) * 100  # coefficient estimates CI %'s
(lwrs - 1) * 100       # lwr CI %'s

