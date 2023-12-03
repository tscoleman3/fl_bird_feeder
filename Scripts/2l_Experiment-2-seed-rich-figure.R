## --------------- HEADER ------------------------------------------------------
## Script name: 2l_Experiment-2-seed-rich-figure.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2022-10-15
## Date Last Modified: 2022-02-18
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script generates a total seed richness figure

## We prepare and then plot the data. Then, we save the ggplot object
## so that we can produce a combined plot with the other response 
## variables.

library(tidyverse)

# Clear the deck
rm(list = ls())

# Import the data
means <- read.csv('Model-output/Experiment-2-total-seed-rich-emmeans.csv')
raw <- read.csv('Data/Experiment-2-seed-richness.csv') |>
  dplyr::select(BLOCK, TREATMENT, RICH)

## --------------- Clean the means dataframe -----------------------------------

# Give columns meaningful or clear names
colnames(means)[1] <- "Treatment"
colnames(means)[2] <- "Mean"
colnames(means)[3] <- "se"
colnames(means)[4] <- "df"
colnames(means)[5] <- "LCL"
colnames(means)[6] <- "UCL"

# Make sure the treatment factor is correct
means$Treatment <- as_factor(means$Treatment)
means$Treatment <- factor(means$Treatment,
                      levels = c("Control", "Low", "Medium", "High"))

## --------------- Clean and add the raw dataframe -----------------------------

# Fix column names for raw df
colnames(raw)[2] <- "Treatment"

# Combine the dataframes
d <- merge(raw, means)

# Remove other dataframes
rm(means, raw)

# Convert to factor
d$Treatment <- as_factor(d$Treatment)

# Reorder factor level
d$Treatment <- fct_relevel(d$Treatment, c('Control', 'Low', 'Medium', 'High'))

## --------------- Create the new figure ---------------------------------------

# Create a vector of colors for each treatment in the figure
cols <- c("#FF61CC","darkgray","#00A9FF","#00BF7D")

seed.rich <- ggplot(d, aes(x = Treatment, y = RICH, fill = Treatment))+
  geom_jitter(shape = 21, width = 0.2, height = 0.2, size = 3, alpha = 0.5)+
  scale_fill_manual(values = c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC"))+
  geom_point(aes(x = Treatment, y = Mean, size = 5))+
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 1)+  
  xlab("Feeder resource richness")+
  ylab('Mean total seed richness')+
  scale_y_continuous(limits = c(0,7), breaks = c(0,1,2,3,4,5,6,7),
                     oob = scales::squish)+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "none",
        axis.text = element_text(face="bold"),
        panel.grid = element_blank())

saveRDS(seed.rich, file = "Model-output/total-mean-seed-rich.RDS")

## --------------- Create the old figure ---------------------------------------

# Create a vector of colors for each treatment in the figure
cols <- c("#FF61CC","darkgray","#00A9FF","#00BF7D")

seed.rich <- ggplot(d, aes(x = Treatment, y = Mean))+
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 1.2,
                color = cols)+
  geom_point(aes(color = Treatment, fill = Treatment),
             shape = 21, size = 6, color = "Black", stroke = 1.7)+
  scale_fill_manual(values = c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC"))+
  scale_y_continuous(limits=c(0,3))+
  xlab("Feeder resource richness")+
  ylab('Mean total seed richness')+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "none")+
  theme(axis.text = element_text(face="bold"))

# saveRDS(seed.rich, file = "Model-output/total-mean-seed-rich.RDS")
# ggsave("Figures/Experiment-2-total-seed-rich.png", 
#         width = 5, height = 7, units = "in")


