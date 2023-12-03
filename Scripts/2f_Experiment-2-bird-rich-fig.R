## --------------- HEADER ------------------------------------------------------
## Script name: 2f_Experiment-2-bird-rich-fig.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2022-10-15
## Date Last Modified: 2023-02-18
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script generates a figure for the bird richness
## data in experiment 2

## We prepare and then plot the data. Then, we save the ggplot object
## so that we can produce a combined plot with the other response 
## variables.

library(tidyverse)

# Clear the decks
rm(list=ls())

# Bring in emmeans and raw data
means <- read.csv("Model-output/Experiment-2-total-bird-rich-emmeans.csv")
raw <- read.csv("Data/Experiment-2-bird-richness.csv") |>
  dplyr::select(site, treatment, richness)

## --------------- Clean the means dataframe -----------------------------------

# Fix column names for means df
colnames(means)[1] <- "Treatment.num"
colnames(means)[2] <- "Mean"
colnames(means)[3] <- "se"
colnames(means)[4] <- "df"
colnames(means)[5] <- "LCL"
colnames(means)[6] <- "UCL"

# sort data for means df
means <- means[order(means$Treatment.num),]

# add treatment description for means df
treat <- tibble(Treatment = c("Control", "Low",
                                  "Medium", "High"))
means <- cbind(means,treat)
means$Treatment <- as_factor(means$Treatment)

## --------------- Clean and add the raw dataframe -----------------------------

# Fix column names for raw df
colnames(raw)[2] <- "Treatment.num"

# Combine the dataframes
d <- merge(raw, means)

# Remove other dataframes
rm(means, raw, treat)

## --------------- Create the new figure ---------------------------------------

cols <- c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC")

dev.new()
bird.rich <- ggplot(d, aes(x = Treatment, y = richness, fill = Treatment))+
  geom_jitter(shape = 21, width = 0.2, height = 0.2, size = 3, alpha = 0.5)+
  scale_fill_manual(values = c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC"))+
  geom_point(aes(x = Treatment, y = Mean, size = 5))+
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 1)+  
  scale_y_continuous(limits = c(0.5,5), breaks = c(0,1,2,3,4,5),
                     oob = scales::squish)+
  xlab("Feeder resource richness")+
  ylab('Mean total bird richness')+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "none",
        axis.text = element_text(face="bold"),
        panel.grid = element_blank())

saveRDS(bird.rich, file = "Model-output/total-mean-bird-rich.RDS")

## --------------- Create the old figure ---------------------------------------

cols <- c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC")

dev.new()
bird.rich <- ggplot(d, aes(x = Treatment, y = Mean))+
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 1.2,
                color = cols)+
  geom_point(aes(color = Treatment, fill = Treatment),
             shape = 21, size = 6, color = "Black", stroke = 1.7)+
  scale_fill_manual(values = c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC"))+
  scale_y_continuous(limits = c(0,4), breaks = c(0,1,2,3,4))+
  xlab("Feeder resource richness")+
  ylab('Mean total bird richness')+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "none")+
  theme(axis.text = element_text(face="bold"))

# saveRDS(bird.rich, file = "Model-output/total-mean-bird-rich.RDS")
# ggsave("Figures/Experiment-2-bird-total-rich.png", width = 5, height = 7, units = "in")
