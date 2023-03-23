## --------------- HEADER ------------------------------------------------------
## Script name: 2a_Experiment-2-treatment-randomization.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2020-10-31
## Date Last Modified: 2023-02-10
## Copyright (c) David S. Mason, 2023
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script assigns randomized treatments for the 
## experiment 2.

# Clear the decks
rm(list=ls())

# Create a vector for feeder resources
species <- c('wheat', 'rye', 'brown top millet', 'white millet', 'oats', 
						 'corn', 'sunflower', 'barely', 'safflower', 'peanuts', 'nyjer',
						 'grubs')

# Set seed for randomization
set.seed(34)

# Make function to generate a block ####
sel_spec <- function(species){
		# four traps per block, up to 12 species
		block_feeders <- matrix(nrow = 4, ncol = 12)
		# name the columns
		colnames(block_feeders) <- paste0((seq(1,12,1)), c("st", "nd", "rd", rep("th", 9))) #
		# select four species for the low richness treatment
		block_feeders[2, 1:4] <- sample(species, 4)
		# select eight species for the low richness treatment
		block_feeders[3, 1:8] <- sample(species, 8)
		# select twelve species for the low richness treatment
		block_feeders[4, 1:12] <- sample(species, 12)
		return(block_feeders)
}

# For loop to create 10 blocks ####
# Create empty list to store vectors
all_feeders <- list() # create a list to store all blocks
n <- c(1:10) # create a sequence for the loop
for(i in n){
	block_name <- paste('Block', i, sep = '') # name each block
	assign(block_name, sel_spec(species)) # run the species selecting function
	all_feeders[[i]]<-get(block_name) # get all the blocks in the list
	}

# Convert the list to  dataframes ####
Block_1 <- as.data.frame(all_feeders[[1]]) # pull and convert
Block_1$Block <- as.character(c("1")) # add value for block column
Block_1$Feeder <- c("Empty", "Low", "Medium", "High") # add values for feeder
Block_1 <- dplyr::select(Block_1, Block, Feeder, everything()) # reorder columns

Block_2 <- as.data.frame(all_feeders[[2]])
Block_2$Block <- as.character(c("2"))
Block_2$Feeder <- c("Empty", "Low", "Medium", "High")
Block_2 <- dplyr::select(Block_2, Block, Feeder, everything())

Block_3 <- as.data.frame(all_feeders[[3]])
Block_3$Block <- as.character(c("3"))
Block_3$Feeder <- c("Empty", "Low", "Medium", "High")
Block_3 <- dplyr::select(Block_3, Block, Feeder, everything())

Block_4 <- as.data.frame(all_feeders[[4]])
Block_4$Block <- as.character(c("4"))
Block_4$Feeder <- c("Empty", "Low", "Medium", "High")
Block_4 <- dplyr::select(Block_4, Block, Feeder, everything())

Block_5 <- as.data.frame(all_feeders[[5]])
Block_5$Block <- as.character(c("5"))
Block_5$Feeder <- c("Empty", "Low", "Medium", "High")
Block_5 <- dplyr::select(Block_5, Block, Feeder, everything())

Block_6 <- as.data.frame(all_feeders[[6]])
Block_6$Block <- as.character(c("6"))
Block_6$Feeder <- c("Empty", "Low", "Medium", "High")
Block_6 <- dplyr::select(Block_6, Block, Feeder, everything())

Block_7 <- as.data.frame(all_feeders[[7]])
Block_7$Block <- as.character(c("7"))
Block_7$Feeder <- c("Empty", "Low", "Medium", "High")
Block_7 <- dplyr::select(Block_7, Block, Feeder, everything())

Block_8 <- as.data.frame(all_feeders[[8]])
Block_8$Block <- as.character(c("8"))
Block_8$Feeder <- c("Empty", "Low", "Medium", "High")
Block_8 <- dplyr::select(Block_8, Block, Feeder, everything())

Block_9 <- as.data.frame(all_feeders[[9]])
Block_9$Block <- as.character(c("9"))
Block_9$Feeder <- c("Empty", "Low", "Medium", "High")
Block_9 <- dplyr::select(Block_9, Block, Feeder, everything())

Block_10 <- as.data.frame(all_feeders[[10]])
Block_10$Block <- as.character(c("10"))
Block_10$Feeder <- c("Empty", "Low", "Medium", "High")
Block_10 <- dplyr::select(Block_10, Block, Feeder, everything())

# Combine dataframe
all_feeder_df <- rbind(Block_1, Block_2, Block_3, Block_4,
											 Block_5, Block_6, Block_7, Block_8,
											 Block_9, Block_10)

# Write to csv ###
# write.csv(all_feeder_df, "selected_species.csv")
