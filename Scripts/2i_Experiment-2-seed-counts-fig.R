## --------------- HEADER ------------------------------------------------------
## Script name: 2i_Experiment-2-seed-counts-fig.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2022-10-21
## Date Last Modified: 2023-02-18
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script analyses the seed detection data by
## sampling date in experiment 2.

## We prepare and then plot the data. Then, we save the ggplot object
## so that we can produce a combined plot with the other response 
## variables. We multiply the mean and CI per sampling by the total
## number of sampling periods (5) to generate an overall total.

# Load packages
library(tidyverse)

# Clear the decks
rm(list=ls())

# Bring in the data
means <- read.csv(file = "Model-output/Experiment-2-periodic-seed-counts-emmeans.csv",
                  header = TRUE, stringsAsFactors = FALSE)
raw <- read.csv(file = "data/Experiment-2-clean-seeds.csv",
                  header = TRUE, stringsAsFactors = FALSE) |>
  pivot_longer(5:30, names_to = "SPECIES", values_to = "COUNT") |>
  group_by(DATE, BLOCK, TREATMENT) |>
  summarize(COUNT = sum(COUNT))

## --------------- Clean the means dataframe -----------------------------------

means <- means[, c(1,2,5,6)]

colnames(means)[1] <- "Treatment"
colnames(means)[2] <- "Mean"
colnames(means)[3] <- "LCL"
colnames(means)[4] <- "UCL"

means$Treatment <- as_factor(means$Treatment)
means$Treatment <- factor(means$Treatment, 
                      levels=c('Control', 'Low', 'Medium', 'High'))

## --------------- Clean and add the raw dataframe -----------------------------

# Fix column names for raw df
colnames(raw)[3] <- "Treatment"

# Create totals
summed <- raw |>
  group_by(BLOCK, Treatment) |>
  summarize(Count = sum(COUNT))

# Combine the dataframes
d <- merge(summed, means)

# Remove other dataframes
rm(means, raw, summed)

# Make sure treatment is a factor
d$Treatment <- as_factor(d$Treatment)

# sort data for means df
d <- d[order(d$Treatment.num),]

d$Treatment <- fct_relevel(d$Treatment, c('Control', 'Low', 'Medium', 'High'))

## --------------- Create the new figure ---------------------------------------

# Set the colors
cols <- c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC")

seed.count <- ggplot(d, aes(x = Treatment, y = Count, fill = Treatment))+
  geom_jitter(shape = 21, width = 0.2, height = 0.2, size = 3, alpha = 0.5)+
  scale_fill_manual(values = c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC"))+
  geom_point(aes(x = Treatment, y = Mean*5, size = 5))+
  geom_errorbar(aes(ymin = LCL*5, ymax = UCL*5), width = 0, size = 1)+
  scale_y_continuous(limits = c(0,20), oob = scales::squish)+
  xlab("Feeder resource richness")+
  ylab('Mean total seed count')+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "none",
        axis.text = element_text(face="bold"),
        panel.grid = element_blank())
saveRDS(seed.count, file = "Model-output/total-mean-seed-count.RDS")

## --------------- Create the old figure ---------------------------------------

# Set the colors
cols <- c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC")

seed.count <- ggplot(d, aes(x = Treatment, y = Mean*5))+
  geom_errorbar(aes(ymin = LCL*5, ymax = UCL*5), width = 0, size = 1.2,
                color = cols)+
  geom_point(aes(color = Treatment, fill = Treatment),
             shape = 21, size = 6, color = "Black", stroke = 1.7)+
  scale_fill_manual(values = c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC"))+
  scale_y_continuous(limits=c(0,6))+
  xlab("Feeder resource richness")+
  ylab('Mean total seed count')+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "none")+
  theme(axis.text = element_text(face="bold"))
# saveRDS(seed.count, file = "Model-output/total-mean-seed-count.RDS")
# ggsave("Figures/Experiment-2-seed-count-totals.png", width = 5, height = 7, units = "in")
