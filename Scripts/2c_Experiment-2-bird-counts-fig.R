## --------------- HEADER ------------------------------------------------------
## Script name: 2c_Experiment-2-bird-counts-fig.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2022-10-20
## Date Last Modified: 2023-01-28
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script creates a figure for the bird 
## detection data from experiment 2

## We prepare and then plot the data. Then, we save the ggplot object
## so that we can produce a combined plot with the other response 
## variables. We multiply the mean and CI for weekly counts by the total weeks
## since acclimatization period to generate an overall total.

library(tidyverse)

# Clear the decks
rm(list=ls())

# Import the data
means <- read.csv("Model-output/Experiment-2-weekly-bird-counts-emmeans.csv")
raw <- read.csv("Data/Experiment-2-bird-observations.csv")

## --------------- Clean the means dataframe -----------------------------------

# Drop some variables we do not need
means <- means[,c(1,2,5,6)]

# Give the columns meaningful or clear names
colnames(means)[1] <- "Treatment.num"
colnames(means)[2] <- "Mean"
colnames(means)[3] <- "LCL"
colnames(means)[4] <- "UCL"

# Add a verbal descriptor column for the treatment
treat <- tibble(Treatment = c("Control", "Low",
                                  "Medium", "High"))
means <- cbind(means,treat)

# Make sure the factor is a factor
means$Treatment <- as_factor(means$Treatment)

## --------------- Clean and add the raw dataframe -----------------------------

# Fix column names for raw df
colnames(raw)[3] <- "Treatment.num"

# Create totals
summed <- raw |>
  group_by(site, Treatment.num) |>
  summarize(Count = sum(count))

# Combine the dataframes
d <- merge(summed, means)

# Remove other dataframes
rm(means, raw, treat, summed)

## --------------- Create the new figure ---------------------------------------

cols <- c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC")

dev.new()
bird.count <- ggplot(d, aes(x = Treatment, y = log(Count+1), fill = Treatment))+
  geom_jitter(shape = 21, width = 0.2, height = 0.2, size = 3, alpha = 0.5)+
  scale_fill_manual(values = c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC"))+
  geom_point(aes(x = Treatment, y = log(Mean*10), size = 5))+
  geom_errorbar(aes(ymin = log(LCL*10), ymax = log(UCL*10)), width = 0, size = 1)+
  scale_y_continuous(limits = c(0,6.6),
                     breaks = c(0, 2.302585, 4.60517, 6.907755),
                     labels = c("0", "10", "100", "1000"),
                     oob = scales::squish)+
  xlab("Feeder resource richness")+
  ylab('Mean total bird count')+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "none",
        axis.text = element_text(face="bold"),
        panel.grid = element_blank())

saveRDS(bird.count, file = "Model-output/total-mean-bird-count.RDS")

## --------------- Create the old figure ---------------------------------------

# Create a vector of colors for the figure
cols <- c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC")

bird.count <- ggplot(d, aes(x = Treatment, y = Mean*10))+
  geom_errorbar(aes(ymin = LCL*10, ymax = UCL*10), width = 0, size = 1.2,
                color = cols)+
  geom_point(aes(color = Treatment, fill = Treatment),
             shape = 21, size = 6, color = "Black", stroke = 1.7)+
  scale_fill_manual(values = c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC"))+
  xlab("Feeder resource richness")+
  ylab('Mean total bird count')+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "none")+
  theme(axis.text = element_text(face="bold"))
# saveRDS(bird.count, file = "Model-output/total-mean-bird-count.RDS")
# ggsave("Figures/Experiment-2-bird-counts-total.png", width = 5, height = 7, units = "in")


