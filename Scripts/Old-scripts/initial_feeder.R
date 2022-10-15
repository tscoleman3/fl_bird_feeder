# 24 March 2019
# Amazing R, User David Mason
# Feeders as an Evolutionary Trap for Plants
# Set-up Workspace ####
d <- read.csv("data/initial_feeder.csv")

library(fitdistrplus)
library(lme4)
library(tidyverse)
library(car)
library(DHARMa)




# Visualize the data ####
boxplot(Seeds ~ Treatment,col=c("white","lightgray"), outline = TRUE, d)
boxplot(rawseeds ~ Treatment,col=c("white","lightgray"), outline = FALSE,d)

### DISTRIBUTION ####
plotdist(d$rawseeds, histo = TRUE, demp = TRUE)
descdist(d$rawseeds, discrete=TRUE, boot=500) # NB or poisson

# neg binom fit
seeds.fit.nb <- glmer.nb(rawseeds ~ Treatment + (1 | time/pair),
                             data = d)

summary(seeds.fit.nb)
Anova(seeds.fit.nb)
# poisson fit
seeds.fit.pois <- glmer(rawseeds ~ Treatment + time + (1 | pair),
                            data = d, family = poisson(link = "log"))
summary(seeds.fit.pois)
Anova(seeds.fit.pois)

E1 <- resid(seeds.fit.pois, type = "pearson")
N <- nrow(d)
p <- length(fixef(seeds.fit.pois)) + 1
overdispersion <- sum(E1^2) / (N-p)

d <- d[-c(17,18),]

plotdist(d$rawseeds, histo = TRUE, demp = TRUE)
descdist(d$rawseeds, discrete=TRUE, boot=500) # NB or poisson

seeds.fit.pois <- glmer(rawseeds ~ Treatment + time + (1 | pair),
                        data = d, family = poisson(link = "log"))
summary(seeds.fit.pois)
Anova(seeds.fit.pois)

E1 <- resid(seeds.fit.pois, type = "pearson")
N <- nrow(d)
p <- length(fixef(seeds.fit.pois)) + 1
overdispersion <- sum(E1^2) / (N-p)
overdispersion

# check some assumptions
plot(rawseeds ~ Treatment, data = d)
plot(residuals(seeds.fit.nb) ~ d$Treatment)
abline(a = 0, b = 0, col = "blue", lwd = 2)
plot(residuals(seeds.fit.pois) ~ d$Treatment)
abline(a = 0, b = 0, col = "blue", lwd = 2)
hist(residuals(seeds.fit.nb))
hist(residuals(seeds.fit.pois))

# plot residuals from prediction to test for assumptions
resfit <- resid(seeds.fit.nb)
resfit <- resid(seeds.fit.pois)
hist(resfit)
plot(d$rawseeds, resfit, 
     ylab = "Residuals", 
     xlab = "Detections") 
abline(0, 0)     

# look at predictions 
preds.nb <- predict(seeds.fit.nb)
preds.pois <- predict(seeds.fit.pois)
par(mfrow = c(2, 2))
plot(rawseeds ~ Treatment, 
     data = d)
plot(preds.nb ~ d$Treatment)
plot(rawseeds ~ Treatment, 
     data = d)
plot(preds.pois ~ d$Treatment)

# look at predictions from model
# neg binom
predict(seeds.fit.nb)
d$pred = exp(predict(seeds.fit.nb))
d$predicted = predict(seeds.fit.nb)    # save the predicted values
d$residuals = residuals(seeds.fit.nb)  # save the residual values
# poisson
predict(seeds.fit.pois)
d$pred = exp(predict(seeds.fit.pois))
d$predicted = predict(seeds.fit.pois)    # save the predicted values
d$residuals = residuals(seeds.fit.pois)  # save the residual values
# quick look at the actual, predicted, and residual values
pred_df <- d %>% 
  dplyr::select(rawseeds, predicted, residuals)
pred_df$predicted = exp(pred_df$predicted)
pred_df$residuals = exp(pred_df$residuals)
par(mfrow = c(1, 2))
plot(fitted(seeds.fit.nb) ~ d$rawseeds)
abline(0, 1, col = "blue", lwd = 2)
plot(fitted(seeds.fit.pois) ~ d$rawseeds)
abline(0, 1, col = "blue", lwd = 2)

# look at our coefs 
summary(seeds.fit.pois)
anova(seeds.fit.pois)
coefs <- summary(seeds.fit.pois)$coef
coefs_est <- exp(coefs[, "Estimate"])
uprs <- exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"])
lwrs <- exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"])
(uprs - 1) * 100       # upr CI %'s
(coefs_est - 1) * 100  # coefficient estimates CI %'s
(lwrs - 1) * 100       # lwr CI %'s

par(mfrow = c(1, 1))
plot(exp(coefs[1:2,"Estimate"]), ylim = range(lwrs[1:2], uprs[1:2]))
segments(1:7, lwrs[1:2], 1:7, uprs[1:2])

coefs.fix <- fixef(seeds.fit.pois)
exp(coefs.fix[1])
exp(coefs.fix[2])


coefs_func <- function(.) {
  beta <- unname(fixef(.))
  control <- exp(beta[1])             # mean count for control 
  treat <-  exp(beta[1] - beta[2])   # mean count for treat 1 is this much greater than control
  c(control_mu = control, treat = treat)
  
}

start <- Sys.time()
rand <- bootMer(x = seeds.fit.pois, FUN = coefs_func, nsim = 1000)$t
Sys.time() - start

summ <- apply(rand, 2, function(x) c(mean = mean(x, na.rm = TRUE),
                                     quantile(x, c(0.025, 0.975), na.rm = TRUE)))
summ # these are backwards at the moment

write.csv(summ, "data/initial-feeder-mod-coeff.csv", row.names = FALSE)

library(MuMIn)
r.squaredGLMM(seeds.fit.pois)

# OLD #######
# Build a model 
mod <- aov(Seeds ~ Treatment + pair + time, data=raw)
mod2 <- glm(Seeds ~ Treatment + pair + time, data=raw)
mod3 <- glm(seeds ~ trap, family = poisson, data=d)
mod4 <- glm(seeds ~ trap + pair, family = poisson, data=d)
mod5 <- glm(seeds ~ pair, family = poisson, data=d)
mod6 <- glm(seeds ~ time, family = poisson, data=d)
mod7 <- glm(seeds ~ 1, family = poisson, data=d)
install.packages("fitdistrplus")
hist(raw$Seeds)
descdist(raw$Seeds, discrete=FALSE, boot=500)
Cand.models.group<-list("null1"=mod2,"null2"=mod3, "null3"=mod4,"null4"=mod5,
                        "null5"=mod6,"null6"= mod7)


#Model selection
library("AICcmodavg")

# Model selection table based on AIC #
aictab(Cand.models.group, modnames = NULL,
       second.ord = TRUE, nobs = NULL, sort = TRUE)

shapiro.test(mod$residuals)
shapiro.test(mod2$residuals)
plot(mod)
plot(mod2)
hist(mod$residuals)
hist(mod2$residuals)
bartlett.test(mod$residuals, g = raw$Treatment)
shapiro.test(raw$Seeds)

summary(aov(Seeds~Treatment+pair, d = raw))

# Nice GG plot #
install.packages("ggplot2")
install.packages("Rmisc")
install.packages("Hmisc")
install.packages("wesanderson")
library(wesanderson)

d1 <- d[which(d$time== "1"),]
d2 <- d[which(d$time== "2"),]

ggplot(data_avgs, aes(Condition, Avg)) +
  geom_col(fill = "royalblue1") + 
  geom_errorbar(aes(ymin = L_CI, ymax = U_CI), width = .25)

ggplot(d1, aes(x = trap, y = seeds, fill = trap)) +
  geom_boxplot() + ggtitle("Box plot")
ggplot(d2, aes(x = trap, y = seeds, fill = trap)) +
  geom_boxplot() + ggtitle("Box plot")

d$time <- factor(d$time)
d$pair <- factor(d$pair)
d1$time <- factor(d1$time)
d1$pair <- factor(d1$pair)
d2$time <- factor(d2$time)
d2$pair <- factor(d2$pair)
raw$pair <- factor(raw$pair)
raw$time <- factor(raw$time)

p <- ggplot(raw, aes(x = Treatment, y = Seeds, fill = Treatment)) +
  geom_violin(trim = FALSE) + scale_fill_manual(values = c("darkseagreen2","peru")) + xlab("") + ylab("Seeds (log transformed)")

p2<- p + ggtitle("Seeds caught in traps") + theme_classic() + theme(legend.position ="top") +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    theme(legend.title = element_blank())+
    theme(legend.position = c(0.8, 0.9),
            legend.direction = "vertical")+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text = element_text(size=10))+
    theme(plot.title = element_text(size = 16))+
    theme(axis.title.y = element_text(size = 14))


p2 + geom_boxplot(width=0.1, fill="white")
p2 + geom_dotplot(binaxis='y', stackdir='center', dotsize=1, binwidth = 0.2, fill = "grey")

ggplot(raw, aes(x = time, y = Seeds, fill = Treatment)) +
  geom_violin(trim = TRUE) + scale_fill_manual(values = c("darkseagreen2","peru")) + xlab("Sampling Event") + ylab("Seeds (log transformed)")+
  geom_boxplot(width=0.2, fill="white")

p3 <- p + ggtitle("Seeds caught in traps") + theme_classic() + theme(legend.position ="top") +
  theme(legend.position = c(0.9, 0.9),
        legend.direction = "vertical")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(size=10))+
  theme(plot.title = element_text(size = 16))+
  theme(axis.title.y = element_text(size = 14))

p3 <- p3 + scale_y_continuous(name="Seeds (log transformed)", limits=c(0, 8))

ggplot(raw, aes(x=time, y=Seeds, fill=Treatment)) +
  geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=1, binwidth = 0.2, fill = "grey")

data()
data(package = .packages(all.available = TRUE))

# Survey data Main #
surv <- read.csv("survey_main.csv")
surv <- surv[1:394,1:4]
city <- surv[which(surv$Location=="City"),]
rural <- surv[which(surv$Location=="Rural"),]

p <- ggplot(surv, aes(Proximity, Location)) + geom_jitter(aes(color = Surface),width = 0.12, size = 1)+
      xlab("Distance to natural area (km)") + ylab("Location")+scale_color_manual(values=c('violetred', 'peru'))

p + theme_classic() + ggtitle("Feeder Location") +
  guides(color = guide_legend(reverse = TRUE))+
  theme(legend.position = "bottom",
        legend.direction = "horizontal")+
  theme(axis.text = element_text(size=10))+
  theme(plot.title = element_text(size = 16))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.y = element_text(size = 14))
# Waffle ##
w <- read.csv("waffle.csv")
w <- w[1:373,1:5]

var <- w$Managed  # the categorical data 
nrows <- 10
df <- expand.grid(y = 1:nrows, x = 1:nrows)
categ_table <- round(table(var) * ((nrows*nrows)/(length(var))))
df$Managed <- factor(rep(names(categ_table), categ_table))  

p <- ggplot(df, aes(x = x, y = y, fill = Managed)) + 
  geom_tile(color = "black", size = 0.5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), trans = 'reverse') +
  scale_fill_manual(values = c("darkseagreen2","lightblue"))+
  labs(title="Soil Surface Beneath Feeders", subtitle = "73% Managed, 27% Unmanaged")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(axis.text = element_text(size=10))+
  theme(plot.title = element_text(size = 16))

p + guides(fill = guide_legend(reverse=TRUE))
  

# Bar ##
bar <- read.csv("bars_easy2.csv")
bar$SEASONS <- as.factor(bar$SEASONS)
g <- ggplot(bar, aes(SEASON))
g + geom_bar(aes(fill=SEASONS), width = 0.5) + theme_classic() + coord_flip() + 
  ylab("Number of responses") + theme(axis.text = element_text(size=10)) +
  theme(plot.title = element_text(size = 16)) + ggtitle("Seasons when respondents stock feeders") + theme(plot.title = element_text(hjust = 0.5)) 
