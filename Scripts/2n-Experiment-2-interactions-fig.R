## --------------- HEADER ------------------------------------------------------
## Script name: 2n_Experiment-2-interactions-figure.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2022-02-12
## Date Last Modified: 2022-02-18
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script generates a combined 3-panel figure from
## the bird and seed count and richness interactions

## We bring in the ggplot objects for the temporal variable*treatment interaction 
## effects for bird rich, seed counts, and seed rich. We also create
## figures depicting the coefficient value and confidence intervals for this
## each interaction. Then, we combine these six objects using patchwork, and 
## we export the image.

library(patchwork)
library(tidyverse)

# Clear the decks
rm(list=ls())

# Bring in the data
seed.count <- read.csv('Model-output/Experiment-2-periodic-seed-counts-emtrends.csv')
seed.rich <- read.csv('Model-output/Experiment-2-periodic-seed-rich-emtrends.csv')
bird.rich <- read.csv('Model-output/Experiment-2-weekly-bird-rich-emtrends.csv')

# Create the coefficient estimate figures

# Bird rich
bird.rich$Treatment <- factor(bird.rich$Treatment, 
                              levels = c("Control", "Low", "Medium", "High"))
cols <- c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC")

bird.rich.coef <- ggplot(bird.rich, aes(x = Treatment, y = Estimate))+
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 1.2,
                color = cols)+
  geom_point(aes(color = Treatment, fill = Treatment),
             shape = 21, size = 6, color = "Black", stroke = 1.7)+
  scale_fill_manual(values = c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC"))+  
  scale_y_continuous(limits = c(-0.1, 0.2))+
  ylab('Parameter Estimate')+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "none")+
  theme(axis.text = element_text(face="bold"))

# Seed counts
seed.count$Treatment <- as_factor(seed.count$Treatment)
cols <- c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC")

seed.count.coef <- ggplot(seed.count, aes(x = Treatment, y = Estimate))+
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 1.2,
                color = cols)+
  geom_point(aes(color = Treatment, fill = Treatment),
             shape = 21, size = 6, color = "Black", stroke = 1.7)+
  scale_fill_manual(values = c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC"))+  
  scale_y_continuous(limits = c(-0.075, 0.075))+
  ylab('Parameter Estimate')+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "none")+
  theme(axis.text = element_text(face="bold"))

# Seed rich
seed.rich$Treatment <- factor(seed.rich$Treatment, 
                    levels = c("Control", "Low", "Medium", "High"))

seed.rich.coef <- ggplot(seed.rich, aes(x = Treatment, y = Estimate))+
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 1.2,
                color = cols)+
  geom_point(aes(color = Treatment, fill = Treatment),
             shape = 21, size = 6, color = "Black", stroke = 1.7)+
  scale_fill_manual(values = c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC"))+  
  ylab('Parameter Estimate')+
  scale_y_continuous(limits = c(-0.06, 0.06))+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "none")+
  theme(axis.text = element_text(face="bold"))

# Bring in the interaction plots
bird.rich.interact <- readRDS(file = "Model-output/weekly-mean-bird-rich-interact.RDS")
seed.rich.interact <- readRDS(file = "Model-output/periodic-mean-seed-rich-interact.RDS")
seed.count.interact <- readRDS(file = "Model-output/periodic-mean-seed-count-interact.RDS")

dev.new()
((bird.rich.interact+theme(plot.margin = margin(0.5,0.5,0.5,0.5,'cm'))+theme(legend.position = "none"))+(bird.rich.coef+labs(x=NULL)+theme(axis.text.x=element_blank())))/
  ((seed.count.interact+labs(x=NULL)+theme(plot.margin = margin(0.5,0.5,0.5,0.5,'cm'))) + (seed.count.coef+labs(x=NULL)+theme(axis.text.x=element_blank())))/
    ((seed.rich.interact+theme(plot.margin = margin(0.5,0.5,0.5,0.5,'cm'))) + seed.rich.coef)+plot_annotation(tag_levels = c('A', '1'))
ggsave("Figures/Experiment-2-comb-interaction.png", width = 12, height = 11, units = "in")


