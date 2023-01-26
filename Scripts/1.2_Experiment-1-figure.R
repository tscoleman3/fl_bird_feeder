## --------------- HEADER ------------------------------------------------------
## Script name: 2_Experiment-1-figure.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliaton: University of Florida
## Date Created: 2022-10-15
## Date Last Modified: 2022-10-15
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script analyzes the seed collection data for 
## experiment 1

## --------------- SET UP WORKPLACE --------------------------------------------

library(tidyverse)

# Bootstrapped coefficient
boot <- read.csv("data/initial-feeder-boot-coeff.csv")
means <- read.csv("data/initial-feeder-emmeans.csv")

## --------------- EMMEANS -----------------------------------------------------

colnames(means)[1] <- "Treatment"
colnames(means)[2] <- "Mean"
colnames(means)[5] <- "LCL"
colnames(means)[6] <- "UCL"

ggplot(means, aes(x = Treatment, y = Mean/14))+
  geom_errorbar(aes(ymin = LCL/14, ymax = UCL/14), width = 0, size = 1.2)+
  geom_point(aes(color = Treatment, fill = Treatment),
             shape = 21, size = 3, color = "Black")+
  scale_fill_manual(values = c("Black", "White"))+
  scale_y_continuous(breaks = c(0,5,10,15,20,25))+
  xlab("")+
  ylab(expression(paste("Mean daily seed arrival", " (0.5 m"^{2}, ")")))+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "none")+
  theme(axis.text = element_text(face="bold"))

ggsave("Figures/Experiment-1-emmeans.png", width = 5, height = 7, units = "in")


## --------------- BOOTSTRAPPED COEFFICIENTS -----------------------------------

rownames(boot)[1] = "Mean"
rownames(boot)[2] = "LCL"
rownames(boot)[3] = "UCL"

boot <- t(boot)

boot2 <- data.frame("Treatment" = c("Control", "Baited"))

boot <- cbind(boot, boot2)

ggplot(boot, aes(x = Treatment, y = Mean/14))+
  geom_errorbar(aes(ymin = LCL/14, ymax = UCL/14), width = 0, size = 1.2)+
  geom_point(aes(color = Treatment, fill = Treatment),
             shape = 21, size = 3, color = "Black")+
  scale_fill_manual(values = c("Black", "White"))+
  xlab("")+
  ylab(expression(paste("Mean daily seed arrival", " (0.5 m"^{2}, ")")))+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "none")+
  theme(axis.text = element_text(face="bold"))

ggsave("Figures/Experiment-1-boot.png", width = 5, height = 7, units = "in")
