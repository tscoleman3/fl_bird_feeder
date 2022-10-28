## --------------- HEADER ------------------------------------------------------
## Script name: 6b_Experiment-2-seed-rich.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliaton: University of Florida
## Date Created: 2022-10-15
## Date Last Modified: 2022-10-15
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script generates a total seed richness figure

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

# Clear the deck
rm(list = ls())

seed.dat.og <- read.csv(file = "data/seed_traps2.0.csv",
                        header = TRUE)
head(seed.dat.og)

## --------------- Prepare the data --------------------------------------------
sites <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,
           8,8,8,8,9,9,9,9,10,10,10,10)
treatment <- c(0,4,8,12,0,4,8,12,0,4,8,12,0,4,8,12,
               0,4,8,12,0,4,8,12,0,4,8,12,0,4,8,12,
               0,4,8,12,0,4,8,12)
seed.rich <- as.data.frame(cbind(sites,treatment))
seed.rich$richness <- c(1,2,5,2,0,0,0,2,1,0,0,5,0,3,4,3,0,2,2,3,
                        1,1,2,2,1,1,1,1,0,1,1,3,0,0,4,4,2,0,0,0)
seed.rich <- as.data.frame(seed.rich) %>% 
  dplyr::select(treatment,sites,richness)

seed.rich$treatment <- as.factor(seed.rich$treatment)
seed.rich$treatment <- recode_factor(seed.rich$treatment ,
        `0` = "Control", `4` = "Low", `8` = "Medium", `12` = "High")

# Make the figure with no time series
d <- seed.rich %>% 
  group_by(treatment) %>%
  summarize(mean = mean(richness),
            sd = sd(richness),
            n = n(),
            se = sd/sqrt(n),
            margin = qt(0.975,df=n-1)*sd/sqrt(n),
            LCL = mean - margin,
            UCL = mean + margin)


## --------------- Visualize total richness ------------------------------------

ggplot(d, aes(x = treatment, y = mean))+
  geom_errorbar(aes(ymin = LCL, ymax = UCL, 
                    color = treatment), width = 0, size = 1.2)+
  scale_color_manual(values = c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC"))+
  geom_point(aes(color = treatment, fill = treatment),
             shape = 21, size = 6, color = "Black", stroke = 1.7)+
  scale_fill_manual(values = c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC"))+
  xlab("Feeder resource richness")+
  ylab('Mean total seed rain richness')+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "none")+
  theme(axis.text = element_text(face="bold"))

ggsave("Figures/Experiment-2-seed-rich-totals.png", width = 5, height = 7, units = "in")
