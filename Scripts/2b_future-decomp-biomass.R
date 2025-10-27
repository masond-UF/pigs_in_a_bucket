## --------------- Header ------------------------------------------------------
## Script name: 2b_future-decomp-biomass.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2021-05-07
## Date Last Modified: 2023-03-10
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script analyzes the future decomposition data

## --------------- Set-up workspace --------------------------------------------
library(tidyverse)
library(fitdistrplus)
library(performance)
library(car)
library(emmeans)
library(patchwork)

# Clear the deck
rm(list=ls())

# Load the data
d <- read.csv("Clean-data/2b_future-decomp-biomass.csv")
head(d)
summary(d)
is.factor(d$EXPOSURE)

## --------------- Prepare the data --------------------------------------------

# Convert to factor
d$FENCE <- as_factor(d$FENCE)

d$EXPOSURE.KG <- round(d$EXPOSURE / 2.205)

# Fix the labels
d$FENCE <- c("Open","Open","Open",
						 "Fenced","Fenced","Fenced")

# Calculate the change in biomass
d <- d %>% 
	mutate(LOST_BIO = abs(PIG_FINAL-PIG_START))

## --------------- Explore the data --------------------------------------------

library(performance)
descdist(d$LOST_BIO, discrete = FALSE) # beta or log-normal
descdist(log(d$LOST_BIO), discrete = FALSE) # closer to beta or log-normal
descdist(log10(d$LOST_BIO), discrete = FALSE) # closer to beta or log-normal

## --------------- Create the model --------------------------------------------

mod <- lm(LOST_BIO~FENCE*EXPOSURE.KG, d)
anova(lm(LOST_BIO~FENCE*EXPOSURE.KG, d))

hist(mod$residuals)
plot(mod$residuals)
summary(mod)

check_model(mod)

# Means
emmeans(mod, pairwise~FENCE*EXPOSURE)
confint(emmeans(mod, pairwise~FENCE))

# Trends
emtrends(mod, pairwise~FENCE, var = 'EXPOSURE')
confint(emtrends(mod, pairwise~FENCE, var = 'EXPOSURE'))

# Coefficients
coef <- tibble(Fence = c('Fenced', 'Open'),
							 Value = c(0.00361, 0.00215),
							 LCL = c(0.00294, 0.00147),
							 UCL = c(0.00428, 0.00282))

# The model fails assumptions.

## --------------- Descriptive statistics ---------------------------------------

# Means
d |>
  group_by(FENCE) |>
  summarize(
    mean_lost_bio = mean(LOST_BIO),
    .groups = 'drop'
  )

slope <- d |>
  # Select the columns we need
  dplyr::select(FENCE, EXPOSURE.KG, LOST_BIO) |>
  
  # Group by fence type
  group_by(FENCE) |>
  
  # Find the min/max values for exposure and the corresponding bio loss
  summarize(
    low_exp_bio = LOST_BIO[EXPOSURE.KG == min(EXPOSURE.KG)],
    high_exp_bio = LOST_BIO[EXPOSURE.KG == max(EXPOSURE.KG)],
    low_exp_kg = min(EXPOSURE.KG),
    high_exp_kg = max(EXPOSURE.KG)
  ) |>
  
  # Now, calculate the slopes and the "per 100kg" value
  mutate(
    # Rise (change in biomass lost)
    total_bio_change = high_exp_bio - low_exp_bio,
    
    # Run (change in kg of exposure)
    total_exp_change_kg = high_exp_kg - low_exp_kg,
    
    # Descriptive slope (kg lost / kg exposure)
    descriptive_slope = total_bio_change / total_exp_change_kg,
    
    # Calculate the "per 100kg" value
    loss_per_100kg = descriptive_slope * 100
  )

## --------------- Visualize the model -----------------------------------------

# Main figure
p1 <- ggplot(data = d, aes(x = EXPOSURE.KG, y = LOST_BIO, color = FENCE, fill = FENCE)) +
  geom_line(aes(group = FENCE), linewidth = 1.75) + 
  geom_point(size = 8, shape = 21, color = 'black', stroke = 2) +
  scale_color_manual(values = c("Fenced" = "#E7B800", "Open" = "#00AFBB")) +
  scale_fill_manual(values = c("Fenced" = "#E7B800", "Open" = "#00AFBB")) +
  scale_x_continuous(breaks = c(0, 200, 400, 600, 800)) + 
  ylab("Biomass loss (kg)") + 
  xlab("Initial biomass exposure (kg)") + 
  theme_classic() +
  theme(axis.title = element_text(face = "bold")) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        legend.position = 'none')

# Output the ggplot object
saveRDS(p1, file = "Output/2b_future-biomass-fig.RDS")


