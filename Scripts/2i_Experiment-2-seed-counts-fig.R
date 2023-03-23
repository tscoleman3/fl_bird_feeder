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
d <- read.csv(file = "Model-output/Experiment-2-periodic-seed-counts-emmeans.csv",
                  header = TRUE, stringsAsFactors = FALSE)

d <- d[, c(1,2,5,6)]

colnames(d)[1] <- "Treatment"
colnames(d)[2] <- "Mean"
colnames(d)[3] <- "LCL"
colnames(d)[4] <- "UCL"

d$Treatment <- as_factor(d$Treatment)
d$Treatment <- factor(d$Treatment, 
                      levels=c('Control', 'Low', 'Medium', 'High'))

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
saveRDS(seed.count, file = "Model-output/total-mean-seed-count.RDS")
# ggsave("Figures/Experiment-2-seed-count-totals.png", width = 5, height = 7, units = "in")
