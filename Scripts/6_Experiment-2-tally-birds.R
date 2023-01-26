## --------------- HEADER ------------------------------------------------------
## Script name: 6_tally-birds.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2023-01-23
## Date Last Modified: 2023-01-23
## Copyright (c) David S. Mason, 2023
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script summarizes bird data from experiment 2

# Clear the decks
rm(list=ls())

# Bring in the data
bird.dat.og <- read.csv(file = "data/bird_data.csv",
                        header = TRUE)
head(bird.dat.og)

# Prep the data
bird.dat <- dplyr::select(bird.dat.og, -c("camera.card", "X", "notes", "photo.number"))
bird.dat$date <- as.Date(with(bird.dat, paste(year, month, day, sep = "-")), "%Y-%m-%d")

# Get the total counts
bird.dat |>
  mutate(week = week(date)) |>
  filter(week != 48 & week != 47) |>
  na.omit(species) |>
  summarize(count = n())

# 4863 total counts

# Create a data frame with summarized value for each species
spec.sum <- bird.dat |>
  mutate(week = week(date)) |>
  filter(week != 48 & week != 47) |>
  na.omit(species) |>
  group_by(species) |>
  summarize(count = n(),
            percent = (count/4863) * 100)

# Check total
sum(spec.sum$percent)
