## --------------- Header ------------------------------------------------------
## Script name: 2c_future-decomp-biomass-ES.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2021-08-16
## Date Last Modified: 2023-03-11
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script analyzes the future decomposition data

## --------------- Set-up workspace --------------------------------------------
library(ggplot2)
require(tidyverse)
require(MASS)

# Clear the decks
rm(list=ls())

# Bring in the data
d <- read.csv("Clean-data/2a_future-decomp-stage.csv")

d$Pig.biomass <- round(d$Pig.biomass/2.205)

## --------------- Prepare the data --------------------------------------------

# Fix the exclusion labels
d$Fence <- c("Open", "Open", "Open",
						"Fenced","Fenced", "Fenced")
d$Fence <- as.factor(d$Fence)
is.factor(d$Fence)

## --------------- Visualize the model -----------------------------------------

# Adjust the weight values to create space for visualization
d[1,2] <- 80
d[4,2] <- 30

p1 <- ggplot(d,aes(y=Stage,x=Pig.biomass,color=Fence,fill=Fence))+
	geom_point(size = 6.5,shape=21,stroke=2)+
	scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
	scale_color_manual(values=c("black","black"))+
	theme_classic()+
	theme(legend.title = element_blank())+
	scale_x_continuous(breaks = c(0, 200, 400, 600, 800))+
	scale_y_continuous(breaks = c(1,2,3,4,5),
										 limits = c(1,5))+
	xlab("Initial biomass exposure (kg)")+
	ylab("Decomposition stage")+
	theme_classic()+
	theme(axis.title = element_text(face="bold"))+
	theme(axis.text = element_text(size = 20),
				axis.title = element_text(size = 25),
				legend.position = 'none')

# Output the ggplot object
saveRDS(p1, file = "Output/2a_future-decomp-stage-fig.RDS")

