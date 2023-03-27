## --------------- Header ------------------------------------------------------
## Script name: 2c_future-decomp-biomass-ES.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2021-05-07
## Date Last Modified: 2023-03-11
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script analyzes the future decomposition data

## --------------- Set-up workspace --------------------------------------------
library(tidyverse)
library(fitdistrplus)
library(performance)
library(devtools)
library(emmeans)

# Clear the deck
rm(list=ls())

# Bring in the data
d <- read.csv("Clean-data/2c_future-decomp-biomass-ES.csv")
# Effect size is log ratio 
# log(biomass change in fenced/biomass change in open)

## --------------- Run the models ----------------------------------------------

mod <- lm(effect.magnitude~Treatment, d = d)
summary(mod)
test(emtrends(mod, ~0, var = 'Treatment'), null = 0)

mod.ln <- lm(effect.magnitude~log(Treatment), d = d)
summary(mod.ln)
test(emtrends(mod.ln, ~0, var = 'Treatment'), null = 0)

anova(mod, mod.ln) # 0.0026427

## --------------- Effect size calculations ------------------------------------

(0.0430672/100)
# 1% increase in biomass increases effect size 0.000430672
# doubling biomass increases effect size by 0.04
# At 55 kg effect size = 0.17, at 110 kg = 0.21
# At 250 kg effect size = 0.23, at 500 kg = 0.27
# At 750 kg effect size = 0.28, at 1000 kg = 0.32
## --------------- Plot the log model ------------------------------------------

# New
dev.new()
p1 <- ggplot(d,aes(y=effect.magnitude,x=Treatment))+
	geom_point(size = 6.5,shape=21,stroke=2)+
	geom_smooth(method="lm",se=FALSE, color = "black",
							formula = y ~ log(x))+
	scale_y_continuous(name="Effect magnitude",
										 breaks = c(0.1,0.15,0.2,0.25,0.3,0.35),
										 limits = c(0.1,0.38))+
	theme_classic()+
	theme(legend.title = element_blank())+
	scale_x_continuous(breaks = c(0, 250, 500, 750, 1000,
																1250, 1500))+
	xlab("Initial biomass exposure (kg)")+
	theme_classic()+
	theme(axis.title = element_text(face="bold"))+
	theme(axis.text = element_text(size = 20),
				axis.title = element_text(size = 25),
				legend.position = 'none')+
	annotate('text', x = 180, y = 0.38, 
					 label = "p = 0.076")

# Output the ggplot object
saveRDS(p1, file = "Output/2c_future-biomass-fig-ES.RDS")

# ggsave('effect_size.jpg', width = 6, height = 6, dpi = 300)
