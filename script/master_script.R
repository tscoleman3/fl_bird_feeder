##########################################
##### FLORIDA SEED DISPERSAL PROJECT #####
##########################################
##### JIMMY HOLDHOLDGRAFER #####
######### DAVID MASON ##########
##### TYLER STEVEN COLEMAN #####





# ----------------------- #
##### PACKAGES NEEDED #####
# ----------------------- #
library(dplyr)          # data wrangling
library(tidyverse)      # data wrangling
library(lme4)           # glmm's
library(fitdistrplus)   # distribution test





# ---------------- #
##### RAW DATA #####
# ---------------- #

# bird.dat are the presence of species in all camera traps ever recorded
raw.bird.dat <- read.csv("data/bird_data2.0.csv",
                     header = TRUE)
head(raw.bird.dat)

# seed.dat are the seeds found in the traps throughout the study
raw.seed.dat <- read.csv("data/seed_traps2.0.csv",
                     header = TRUE)
head(raw.seed.dat)

# prelim.seed.dat are from the baited vs un-baited feeders effect on seed count
raw.prelim.seed.dat <- read.csv("data/feeder_all.csv",
                            header = TRUE)
head(raw.prelim.seed.dat)





# --------------------- #
##### DATA CLEANING #####
# --------------------- #

### bird presence/observation data ###

# remove unwanted rows 
bird.dat <- dplyr::select(raw.bird.dat, -c("camera.card", "X", "notes", "photo.number"))

# add date column
bird.dat$date <- as.Date(with(bird.dat, paste(year, month, day, sep = "-")), "%Y-%m-%d")
head(bird.dat)

# 
bird
