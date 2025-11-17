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

# what does this mean???
# experiment starts on July 7, 2016.
# pigs moved on July 22nd, 2016.

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

library(performance)
check_model(mod)

# Means
library(emmeans)
emmeans(mod, pairwise~FENCE*EXPOSURE.KG)
confint(emmeans(mod, pairwise~FENCE))

# Trends
emtrends(mod, pairwise~FENCE, var = 'EXPOSURE.KG')
confint(emtrends(mod, pairwise~FENCE, var = 'EXPOSURE.KG'))

# Coefficients
coef <- tibble(Fence = c('Fenced', 'Open'),
							 Value = c(0.00796, 0.00473),
							 LCL = c(0.00647, 0.00324),
							 UCL = c(0.00944, 0.00621))

## --------------- Descriptive statistics ---------------------------------------

# Means
d |>
  group_by(FENCE) |>
  dplyr::summarize(
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

# Create new data frame to predict from
pred.dat <- expand.grid(
  EXPOSURE.KG = c(25, 181, 726),
  FENCE = c("Open", "Fenced")
)

pred.dat <- pred.dat[order(pred.dat$FENCE, pred.dat$EXPOSURE.KG), ]

# Predictions from model
preds <- predict(mod, 
								 newdata = pred.dat, 
								 se = T)

# Combine predictions to new data frame for plotting
pred.dat <- cbind(pred.dat, fit = preds$fit)
pred.dat <- cbind(pred.dat, se.fit = preds$se.fit)

# Calculate 95% CI for predictions from predicted standard errors
pred.dat$LCL <- pred.dat$fit - (1.96*pred.dat$se.fit) # Correct?
pred.dat$UCL <- pred.dat$fit + (1.96*pred.dat$se.fit)

# Main figure
p1 <- ggplot(data=pred.dat, aes(x = EXPOSURE.KG, y = fit, color = FENCE))+
	geom_line(linewidth = 1.75)+
	scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
	scale_color_manual(values = c('black', 'black'))+
	geom_point(data = d, aes(x = EXPOSURE.KG, y = LOST_BIO, fill = FENCE), 
						 size = 8, shape = 21, color = 'black', stroke = 2)+
	scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
	scale_x_continuous(breaks = c(0, 200, 400, 800))+
	ylab("Biomass loss (kg)")+
	xlab("Initial biomass exposure (kg)")+ 
	theme_classic()+
	theme(axis.title = element_text(face="bold"))+
	theme(axis.text = element_text(size = 20),
				axis.title = element_text(size = 25),
				legend.position = 'none')

# Coefficient comparisons
p2 <- ggplot(coef, aes(x=Fence,y=Value,color=Fence,fill=Fence))+
	geom_errorbar(aes(ymin=LCL,ymax=UCL),
								position = position_dodge(width = 0.5),color='black', width=0.2)+
	geom_point(size=6.5,shape=21,stroke=2,
						 position=position_dodge(width = 0.5),color='black')+
	scale_y_continuous(limits = c(0.003,0.01))+
	scale_color_manual(values=c("#00AFBB", "#E7B800"))+
	scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
	theme_classic()+
	theme(legend.position = 'none')+
	xlab("")+
	ylab('Slope estimate')+
	theme(axis.title = element_text(face="bold"))+
	annotate('text', x = 2, y = 0.01, 
					 label = "p = 0.022")+
	theme(plot.title = element_text(hjust = 0.5))+
	theme(axis.text = element_text(size = 10),
				axis.title = element_text(size = 15))+
	theme(plot.title = element_text(hjust = 0.5),
				axis.text.x = element_blank(),
				axis.ticks.x = element_blank())

# Combine the figures
comb <- p1 +
  inset_element(
    p2,
    left   = 0.07,
    bottom = 0.65,
    right  = 0.3,
    top    = 0.98,
    align_to = "full"
  )

# Output the ggplot object
saveRDS(p1, file = "Output/2b_future-biomass-fig.RDS")


