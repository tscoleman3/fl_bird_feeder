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
d <- read.csv('Model-output/Experiment-2-total-seed-rich-emmeans.csv')

# Give columns meaningful or clear names
colnames(d)[1] <- "Treatment"
colnames(d)[2] <- "Mean"
colnames(d)[3] <- "se"
colnames(d)[4] <- "df"
colnames(d)[5] <- "LCL"
colnames(d)[6] <- "UCL"

# Make sure the treatment factor is correct
d$Treatment <- as_factor(d$Treatment)
d$Treatment <- factor(d$Treatment,
                      levels = c("Control", "Low", "Medium", "High"))

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

saveRDS(seed.rich, file = "Model-output/total-mean-seed-rich.RDS")
# ggsave("Figures/Experiment-2-total-seed-rich.png", width = 5, height = 7, units = "in")


