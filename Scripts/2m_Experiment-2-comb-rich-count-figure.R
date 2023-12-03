## --------------- HEADER ------------------------------------------------------
## Script name: 2m_Experiment-2-combined-rich-count-figure.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2022-02-12
## Date Last Modified: 2022-02-12
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script generates a combined 4-panel figure from
## the bird and seed count and richness figures.

## We bring in the ggplot objects for bird counts, bird rich, seed counts,
## seed rich and combine them using patchwork. Then we export the image.

library(tidyverse)
library(patchwork)

# Clear the decks
rm(list = ls())

# Bring in the figures
bird.count <- readRDS(file = "Model-output/total-mean-bird-count.RDS")
bird.rich <- readRDS(file = "Model-output/total-mean-bird-rich.RDS")
seed.count <- readRDS(file = "Model-output/total-mean-seed-count.RDS")
seed.rich <- readRDS(file = "Model-output/total-mean-seed-rich.RDS")

(bird.count + labs(x = NULL) + theme(axis.text.x = element_blank()) +
  bird.rich + labs(x = NULL) + theme(axis.text.x = element_blank())) / (seed.count + seed.rich) +
  plot_annotation(tag_levels = c("A", "1"))
ggsave("Figures/Experiment-2-comb-rich-count.png", width = 9, height = 11, units = "in")
