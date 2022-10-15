# October 31st 2020—David Mason
# Code for selecting random assortments of seeds for each site
# Set up ####
library(dplyr)

species <- c('wheat', 'rye', 'brown top millet', 'white millet', 'oats', 
						 'corn', 'sunflower', 'barely', 'safflower', 'peanuts', 'nyjer',
						 'grubs')

set.seed(34)
# Make function to generate a site ####
sel_spec <- function(species){
		# four traps per site, up to 12 species
		site_feeders <- matrix(nrow = 4, ncol = 12)
		# name the columns
		colnames(site_feeders) <- paste0((seq(1,12,1)), c("st", "nd", "rd", rep("th", 9))) #
		# select four species for the low richness treatment
		site_feeders[2, 1:4] <- sample(species, 4)
		# select eight species for the low richness treatment
		site_feeders[3, 1:8] <- sample(species, 8)
		# select twelve species for the low richness treatment
		site_feeders[4, 1:12] <- sample(species, 12)
		return(site_feeders)
		}
# For loop to create 10 sites ####
# Create empty list to store vectors
all_feeders <- list() # create a list to store all sites
n <- c(1:10) # create a sequence for the loop
for(i in n){
	site_name <- paste('Site', i, sep = '') # name each site
	assign(site_name, sel_spec(species)) # run the species selecting function
	all_feeders[[i]]<-get(site_name) # get all the sites in the list
	}

# Convert the list to  dataframes ####
# This is where my code gets a REALLY clunky. I was struggling to get the list
# of sites into a readable format. I ended up just pulling them one at a time.

Site_1 <- as.data.frame(all_feeders[[1]]) # pull and convert
Site_1$Site <- as.character(c("1")) # add value for site column
Site_1$Feeder <- c("Empty", "Low", "Medium", "High") # add values for feeder
Site_1 <- select(Site_1, Site, Feeder, everything()) # reorder columns

Site_2 <- as.data.frame(all_feeders[[2]])
Site_2$Site <- as.character(c("2"))
Site_2$Feeder <- c("Empty", "Low", "Medium", "High")
Site_2 <- select(Site_2, Site, Feeder, everything())

Site_3 <- as.data.frame(all_feeders[[3]])
Site_3$Site <- as.character(c("3"))
Site_3$Feeder <- c("Empty", "Low", "Medium", "High")
Site_3 <- select(Site_3, Site, Feeder, everything())

Site_4 <- as.data.frame(all_feeders[[4]])
Site_4$Site <- as.character(c("4"))
Site_4$Feeder <- c("Empty", "Low", "Medium", "High")
Site_4 <- select(Site_4, Site, Feeder, everything())

Site_5 <- as.data.frame(all_feeders[[5]])
Site_5$Site <- as.character(c("5"))
Site_5$Feeder <- c("Empty", "Low", "Medium", "High")
Site_5 <- select(Site_5, Site, Feeder, everything())

Site_6 <- as.data.frame(all_feeders[[6]])
Site_6$Site <- as.character(c("6"))
Site_6$Feeder <- c("Empty", "Low", "Medium", "High")
Site_6 <- select(Site_6, Site, Feeder, everything())

Site_7 <- as.data.frame(all_feeders[[7]])
Site_7$Site <- as.character(c("7"))
Site_7$Feeder <- c("Empty", "Low", "Medium", "High")
Site_7 <- select(Site_7, Site, Feeder, everything())

Site_8 <- as.data.frame(all_feeders[[8]])
Site_8$Site <- as.character(c("8"))
Site_8$Feeder <- c("Empty", "Low", "Medium", "High")
Site_8 <- select(Site_8, Site, Feeder, everything())

Site_9 <- as.data.frame(all_feeders[[9]])
Site_9$Site <- as.character(c("9"))
Site_9$Feeder <- c("Empty", "Low", "Medium", "High")
Site_9 <- select(Site_9, Site, Feeder, everything())

Site_10 <- as.data.frame(all_feeders[[10]])
Site_10$Site <- as.character(c("10"))
Site_10$Feeder <- c("Empty", "Low", "Medium", "High")
Site_10 <- select(Site_10, Site, Feeder, everything())

all_feeder_df <- rbind(Site_1, Site_2, Site_3, Site_4,
											 Site_5, Site_6, Site_7, Site_8,
											 Site_9, Site_10)

# Create a table ####
write.csv(all_feeder_df, "selected_species.csv")
