## --------------- HEADER ------------------------------------------------------
## Script name: 2b_Experiment-2-bird-detect-fig.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliaton: University of Florida
## Date Created: 2022-10-20
## Date Last Modified: 2022-10-20
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script creates a figure for the bird 
## detection data from experiment 2

rm(list=ls())

d <- read.csv("Model-output/experiment-2-bird-detect-emmeans.csv")

colnames(d)[1] <- "Treatment.num"
colnames(d)[2] <- "Mean"
colnames(d)[3] <- "se"
colnames(d)[4] <- "df"
colnames(d)[5] <- "LCL"
colnames(d)[6] <- "UCL"

treat <- data_frame(Treatment = c("Control", "Low",
                                  "Medium", "High"))

d <- cbind(d,treat)

d$Treatment <- as_factor(d$Treatment)

cols <- c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC")

ggplot(d, aes(x = Treatment, y = Mean))+
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, size = 1.2,
                color = cols)+
  geom_point(aes(color = Treatment, fill = Treatment),
             shape = 21, size = 6, color = "Black", stroke = 1.7)+
  scale_fill_manual(values = c("darkgray", "#00A9FF", "#00BF7D",  "#FF61CC"))+
  xlab("Feeder resource richness")+
  ylab('Mean total bird detections')+
  theme_bw()+
  theme(text = element_text(size = 20),
        legend.position = "none")+
  theme(axis.text = element_text(face="bold"))

ggsave("Figures/Experiment-2-bird-detect.png", width = 5, height = 7, units = "in")
