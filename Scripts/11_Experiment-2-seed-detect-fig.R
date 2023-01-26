## --------------- HEADER ------------------------------------------------------
## Script name: 4b_Experiment-2-seed-detect-fig.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliaton: University of Florida
## Date Created: 2022-10-21
## Date Last Modified: 2022-10-21
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script analyses the seed detection data by
## sampling date in experiment 2.

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

# Clear the decks
rm(list=ls())

# Bring in the data
seeds <- read.csv(file = "data/seed_traps2.0.csv",
                  header = TRUE, stringsAsFactors = FALSE)


# Make the figure with time series
d <- seeds %>% 
  dplyr::select(-NOTES) %>% 
  pivot_longer(5:35, names_to = "SPECIES", values_to = "DETECTIONS") %>% 
  group_by(DATE, BLOCK, TREATMENT) %>% 
  summarize(DETECT = sum(DETECTIONS)) %>% 
  group_by(DATE, TREATMENT) %>%
  summarize(mean = mean(DETECT),
            sd = sd(DETECT),
            n = n(),
            se = sd/sqrt(n),
            margin = qt(0.975,df=n-1)*sd/sqrt(n),
            LCL = mean - margin,
            UCL = mean + margin)

d$TREATMENT <- as_factor(d$TREATMENT)
d$TREATMENT <- factor(d$TREATMENT, 
                      levels=c('Control', 'Low', 'Medium', 'High'))

d$DATE <- as_date(d$DATE)

cols <- c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC")


ggplot(d, aes(x = TREATMENT, y = mean, col = TREATMENT))+
  geom_point(position=position_dodge(width=3.5))+
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 1.2,
                position=position_dodge(width=3.5))+
  scale_color_manual(values = cols)+
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 8))+
  xlab("Sampling period")+
  ylab('Mean total seed detections')+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "none")+
  theme(axis.text = element_text(face="bold"))  

ggsave("Figures/Experiment-2-seed-detect.png", width = 5, height = 7, units = "in")

# Make the figure with no time series
d <- seeds %>% 
  dplyr::select(-NOTES) %>% 
  pivot_longer(5:35, names_to = "SPECIES", values_to = "DETECTIONS") %>% 
  group_by(BLOCK, TREATMENT) %>% 
  summarize(DETECT = sum(DETECTIONS)) %>% 
  group_by(TREATMENT) %>%
  summarize(mean = mean(DETECT),
            sd = sd(DETECT),
            n = n(),
            se = sd/sqrt(n),
            margin = qt(0.975,df=n-1)*sd/sqrt(n),
            LCL = mean - margin,
            UCL = mean + margin)

d$TREATMENT <- as_factor(d$TREATMENT)
d$TREATMENT <- factor(d$TREATMENT, 
                      levels=c('Control', 'Low', 'Medium', 'High'),
                      ordered = TRUE)
# R won't change the order of the factor, despite it being a factor
# and despite the code clearly asking it to follow that order.

# SO, I am putting these colors in....the wrong order!
cols <- c("darkgray", "#FF61CC", "#00A9FF", "#00BF7D")

ggplot(d, aes(x = TREATMENT, y = mean))+
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 1.2,
                color = cols)+
  scale_color_manual(values = c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC"))+
  geom_point(aes(color = TREATMENT, fill = TREATMENT),
             shape = 21, size = 6, color = "Black", stroke = 1.7)+
  scale_fill_manual(values = c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC"))+
  xlab("Feeder resource richness")+
  ylab('Mean total seed rain detections')+
  scale_y_continuous(limits = c(0,12),
                     breaks = c(0,2,4,6,8,10,12),
                     labels = c(0,2,4,6,8,10,12))+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "none")+
  theme(axis.text = element_text(face="bold"))

ggsave("Figures/Experiment-2-seed-detect-totals.png", width = 5, height = 7, units = "in")
