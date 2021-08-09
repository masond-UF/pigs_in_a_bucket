############ 7 MAY 2021—DAVID MASON—GARBAGE CAN PIG DECOMP ANALYSIS ############
############ LOAD PACKAGES AND BRING IN DATA ###################################  
d.tmp <- read.csv("Data/garbage_can.csv")
head(d)
summary(d)

library(tidyverse)
############ CALCULATE CHANGE IN BIOMASS #######################################
d <- d.tmp %>% 
	mutate(LOST_BIO = ABS(PIG_FINAL-PIG_START))
					