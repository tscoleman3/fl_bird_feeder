## --------------- HEADER ------------------------------------------------------
## Script name: 1b_Experiment-1-figure.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2022-10-15
## Date Last Modified: 2023-01-26
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script analyzes the seed collection data for 
## experiment 1

## --------------- SET UP WORKPLACE --------------------------------------------

library(tidyverse)

rm(list=ls())

# Bring in the data
d <- read.csv("Data/Experiment-1.csv") # raw data
means <- read.csv("Model-output/experiment-1-emmeans.csv") # means

## --------------- EMMEANS -----------------------------------------------------

# Rename the columns
colnames(means)[1] <- "Treatment"
colnames(means)[2] <- "Mean"
colnames(means)[5] <- "LCL"
colnames(means)[6] <- "UCL"

# Combine the data
comb <- merge(d, means)
comb$time <- as_factor(comb$time)

ggplot(comb, aes(x = Treatment, y = Mean))+
  geom_jitter(aes(x = Treatment, y = rawseeds, shape = time), color = 'black', size = 3,
              fill = 'gray', alpha = 0.5, width = 0.3)+
  scale_shape_manual(values = c(21,23))+
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 1.2)+
  geom_point(aes(color = Treatment, fill = Treatment),
             shape = 21, size = 5, color = "Black")+
  scale_fill_manual(values = c("Black", "White"))+
  scale_y_continuous(limits = c(0,2000))+
  xlab("")+
  ylab(expression(paste("Mean total seed arrival", " (0.5 m"^{2}, ")")))+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "none")+
  theme(axis.text = element_text(face="bold"))


ggsave("Figures/Experiment-1-emmeans.png", width = 5, height = 7, units = "in")


# OLD
ggplot(means, aes(x = Treatment, y = Mean))+
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 1.2)+
  geom_point(aes(color = Treatment, fill = Treatment),
             shape = 21, size = 3, color = "Black")+
  scale_fill_manual(values = c("Black", "White"))+
  scale_y_continuous(limits = c(0,275),
                     breaks = c(0,25,50,75,100,125,150,175,200,
                                225,250,275))+
  xlab("")+
  ylab(expression(paste("Mean total seed arrival", " (0.5 m"^{2}, ")")))+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "none")+
  theme(axis.text = element_text(face="bold"))