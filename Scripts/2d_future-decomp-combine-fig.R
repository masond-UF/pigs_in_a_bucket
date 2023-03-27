## --------------- Header ------------------------------------------------------
## Script name: 2d_future-decomp-combine-fig.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2021-08-16
## Date Last Modified: 2023-03-11
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script analyzes the future decomposition data

## --------------- Set-up workspace --------------------------------------------
library(tidyverse)
library(patchwork)

# Clear the decks
rm(list=ls())

# Bring in the figures
decomp.stage <- readRDS(file = "Output/2a_future-decomp-stage-fig.RDS")
biomass <- readRDS(file = "Output/2b_future-biomass-fig.RDS")
biomass.ES <- readRDS(file = "Output/2c_future-biomass-ES-fig.RDS")

# annotate_figure(figure,bottom = textGrob("Carrion biomass (kg)", gp = gpar(cex = 2)))

dev.new()
comb <- (decomp.stage+ labs(x = NULL))+biomass+(biomass.ES+labs(x = NULL))

ggsave('Figures/2_decomp-rate.png',width = 20, height = 10, 
			 units = 'in', dpi = 300)
