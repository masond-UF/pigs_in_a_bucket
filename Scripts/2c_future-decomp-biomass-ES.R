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

# mod <- lm(effect.magnitude~Treatment, d = d)
# summary(mod)
# test(emtrends(mod, ~0, var = 'Treatment'), null = 0)
# 
# mod.ln <- lm(effect.magnitude~log(Treatment), d = d)
# summary(mod.ln)
# test(emtrends(mod.ln, ~0, var = 'Treatment'), null = 0)
# 
# anova(mod, mod.ln) # 0.0026427

## --------------- Effect size calculations ------------------------------------

min.data <- d %>% filter(Treatment == min(Treatment)) # 25 kg
max.data <- d %>% filter(Treatment == max(Treatment)) # 725 kg

total.effect.change <- max.data$effect.magnitude - min.data$effect.magnitude
total.treatment.change.kg <- max.data$Treatment - min.data$Treatment # 725 - 25 = 700

# Average increase in effect size per 1 kg increase in Treatment
overall.rate.of.change <- total.effect.change / total.treatment.change.kg

print(paste("Overall rate of change (per 1 kg):", round(overall.rate.of.change, 7)))

# Calculate the "doubling" effect (e.g., doubling from 25 kg)
# Change in biomass for this doubling = 25 kg
increase.for.doubling.from.start <- overall.rate.of.change * 25
print(paste("Estimated increase when doubling from 25 kg:", round(increase.for.doubling.from.start, 4)))

# Report the rate of change per kg
print(paste("Increase in effect size per 1 kg change:", round(overall.rate.of.change, 7)))

# Calculate the estimated increase for a 100 kg change
increase.per.100kg <- overall.rate.of.change * 100

# Print the result
print(paste("Estimated increase in effect size per 100 kg change:", round(increase.per.100kg, 4)))

## --------------- Plot the log model ------------------------------------------

# New
p1 <- ggplot(d, aes(y = effect.magnitude, x = Treatment)) +
  geom_point(size = 6.5, shape = 21, stroke = 2) +
  scale_y_continuous(name = "Effect magnitude",
                     breaks = c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35),
                     limits = c(0.1, 0.38)) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  scale_x_continuous(breaks = c(0, 200, 400, 600, 800)) + 
  xlab("Initial biomass exposure (kg)") +
  theme_classic() +
  theme(axis.title = element_text(face = "bold")) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        legend.position = 'none')

# Output the ggplot object
saveRDS(p1, file = "Output/2c_future-biomass-fig-ES.RDS")

