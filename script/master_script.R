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
library(lubridate)      # date time stuff





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





# ----------------------------------- #
##### DATA CLEANING & SUMMARIZING #####
# ----------------------------------- #

##### bird presence/observation data #####
head(raw.bird.dat)

# remove unwanted rows 
bird.dat <- dplyr::select(raw.bird.dat, -c("camera.card", "X", "notes", "photo.number"))

# add date column
bird.dat$date <- as.Date(with(bird.dat, paste(year, month, day, sep = "-")), "%Y-%m-%d")
head(bird.dat)

#' here we want to make a summary table that groups 
#' observations of birds on the camera traps by week
#' to see if there was a clear acclimation period
#' from when the feeders were originally put out on 
#' 20 November 2020 (i.e., week 47)
bird.week <- bird.dat %>% 
  group_by(week = week(date)) %>% 
  summarise(value = sum(presence))
bird.week
# there seems to be a 2 week acclimation period 

# remove those first two weeks of acclimation
bird.dat.clean <- bird.dat %>% 
  mutate(week = week(date)) %>% 
  filter(week != 48 & week != 47)
unique(bird.week$week)        # weeks 47 and 48 exist 
unique(bird.dat.clean$week)   # weeks 47 and 48 should be gone

# remove NA's 
bird.dat.clean <- bird.dat.clean %>% 
  drop_na(species)

# make treatment a factor 
str(bird.dat.clean)
bird.dat.clean$treatment <- as.factor(bird.dat.clean$treatment)
str(bird.dat.clean)

#' we want to summarize the count of birds present 
#' by species at each treatment and site
bird.sps.count.by.treat <- bird.dat.clean %>% 
  group_by(site, treatment, species) %>% 
  summarise(count = n())
bird.sps.count.by.treat

#' we want to summarize the count of birds 
#' by treatment and site to give us the total
#' observations for analysis
bird.obs <- bird.dat.clean %>% 
  group_by(site, treatment) %>% 
  summarise(count = n())
bird.obs <- as.data.frame(bird.obs)   # make the tibble table a df

#' make sure we have observations for all sites;
#' we should have a count of 4 (4 treatments) for
#' each site
bird.obs %>% 
  group_by(site) %>% 
  summarise(count = n())
# we do not have 4 for each

# add in the missing 0 (0's ARE IMPORTANT)
missing <- data.frame(site = "six", 
                      treatment = "4",
                      count = 0)
bird.obs <- rbind(bird.obs, missing)

# now check it again
bird.obs %>% 
  group_by(site) %>% 
  summarise(count = n())
# 4 for each site

# add in richness 
bird <- bird.obs %>% arrange(site)

# david originally had this:
# sites <- c(1,1,1,1,
#            2,2,2,2,
#            3,3,3,3,
#            4,4,4,4,
#            5,5,5,5,
#            6,6,6,6,
#            7,7,7,7,
#            8,8,8,8,
#            9,9,9,9,
#            10,10,10,10)
# treatment <- c(0,4,8,12,
#                0,4,8,12,
#                0,4,8,12,
#                0,4,8,12,
#                0,4,8,12,
#                0,4,8,12,
#                0,4,8,12,
#                0,4,8,12,
#                0,4,8,12,
#                0,4,8,12)
# bird.rich$richness <- c(2,2,3,3,   # 1
#                         2,2,2,3,   # 2
#                         3,3,4,4,   # 3
#                         3,3,3,3,   # 4
#                         4,4,5,5,   # 5
#                         2,0,2,1,   # 6
#                         1,2,2,1,   # 7
#                         4,4,4,4,   # 8
#                         2,2,2,2,   # 9
#                         1,1,1,1)   # 10

# now let's add richness to our bird data
unique(bird$site)   # this order; alphabetical 
bird$rich <- c(4, 4, 4, 4,   # 8
               4, 4, 5, 5,   # 5
               3, 3, 3, 3,   # 4
               2, 2, 2, 2,   # 9
               2, 2, 3, 3,   # 1
               1, 2, 2, 1,   # 7
               2, 0, 2, 1,   # 6
               1, 1, 1, 1,   # 10
               3, 3, 4, 4,   # 3
               2, 2, 2, 3)   # 2
bird





##### seed data #####
head(raw.seed.dat)

# remove unwanted rows 
seed.dat <- dplyr::select(raw.seed.dat, -c("NOTES"))
head(seed.dat)

# change data layout
seed.dat <- pivot_longer(seed.dat, cols = 5:35,
                         names_to = "SPECIES",
                         values_to = "SEEDS")
head(seed.dat)

# change date format
seed.dat$DATE <- as.Date(seed.dat$DATE)

# remove the acclimation weeks 47 and 48
seed.dat.clean <- seed.dat %>% 
  mutate(week = week(DATE)) %>% 
  filter(week != 47 & week != 48)
