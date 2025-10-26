## --------------- Header ------------------------------------------------------
## Script name: 1_fly-survey.R
## Author: David S. Mason, UF D.E.E.R. Lab
## Department: Wildlife Ecology and Conservation
## Affiliation: University of Florida
## Date Created: 2021-06-18
## Date Last Modified: 2023-03-10
## Copyright (c) David S. Mason, 2022
## Contact: masond@ufl.edu, @EcoGraffito
## Purpose of script: This script analyzes and visualizes the fly survey data
## from experiment 1.

## --------------- Set-up workspace --------------------------------------------
# Data munging and visualization
library(tidyverse)
library(lubridate)
# Modeling
library(fitdistrplus)
library(car)
library(emmeans)
library(performance)
library(DHARMa)
# Combing figures
library(patchwork)
# Exporting tables
library(knitr) 
library(htmlTable)
library(kableExtra)
library(magick)
library(broom)
library(broom.mixed)
library(tidymodels)
library(lme4)
options(scipen = 999)

# Clean the decks
rm(list=ls())

# Bring in the data
d <- read.csv("Clean-data/1_fly-survey.csv", stringsAsFactors=T)
head(d)
summary(d)

# The MME at John Starr started July 5, 2016. We exposed fresh pigs to those 
# plots for 4hrs on July 22, then they were moved to the common garden location.
## --------------- Prepare the data --------------------------------------------

# Fix the biomass
d <- d |>
	mutate(Biomass = round(Biomass/2.205, 2))

# Create a Date column
d$Date <- paste(d$Year, d$Month, d$Day, sep="-") %>% ymd() %>% as.Date()

# Pivot data into long format
d.lg <- d |> 
	dplyr::select(-NumObs, -AVGrate, -SUMrate) |>
	pivot_longer(9:16,
							 names_to = c('Obs', 'Metric'),
							 names_sep = "_") |>
	dplyr::filter(Metric != 'tot')

colnames(d.lg)[12] ="Flies.sec"

# Flies per hour
d.lg <- d.lg |>
	dplyr::select(-Metric) |>
	mutate(Flies.hr = Flies.sec*60)
d.lg$Flies.hr <- round(d.lg$Flies.hr)

# Days since first observation of flies
d.lg <- d.lg |>
	mutate(Days.since.start = as.numeric(ymd(Date) - ymd('2016-07-05')))

d.lg$Fence <- as_factor(d.lg$Fence)
d.lg$Biomass <- as.numeric(d.lg$Biomass)
d.lg$Days.since.start <- as.numeric(d.lg$Days.since.start)

d.lg <- d.lg |>
	na.omit(Flies.hr)

## --------------- Data exploration --------------------------------------------

# Flies caught per hour by date for each biomass treatment
ggplot(d.lg, aes(y=Flies.hr, x=Date, color=Fence))+
	geom_point()+
	facet_wrap(~Biomass)


## --------------- Build the model -----------------------------------------------

descdist(d.lg$Flies.hr, discrete = TRUE)

# Explore mean-variance relationship to assess which negative binomial
# distribution to use for the raw seeds
mean.var <- d.lg |>
	group_by(Biomass, Fence) |>
	dplyr::summarize(Mean.det = mean(Flies.hr),
						Var.det = var(Flies.hr))
q1 <- qplot(Mean.det,Var.det,data=mean.var)

print(q1+
	## linear (quasi-Poisson/NB1) fit
	geom_smooth(method="lm",formula=y~x-1)+
	## smooth (loess)
	geom_smooth(colour="red")+
	## semi-quadratic (NB2/LNP)
	geom_smooth(method="lm",formula=y~I(x^2)+offset(x)-1,colour="purple")+
	## Poisson (v=m)
	geom_abline(intercept=0,slope=1,lty=2))
# blue = nb1
# purple = nb2

rm(mean.var, q1)

library(glmmTMB)
m1 <- glmmTMB(Flies.hr ~ Biomass*Fence*Days.since.start, data = d.lg, family = 'nbinom2')
ggplot(d.lg, aes(x = Days.since.start, y = residuals(m1))) +
  geom_point() +
  geom_smooth(method = "loess")


d.lg$Days_squared <- d.lg$Days.since.start^2
d.lg$c.Days.since.start <- scale(d.lg$Days.since.start, center = TRUE, scale = FALSE)
d.lg$c.Days_squared <- d.lg$c.Days.since.start^2

m2 <- glmmTMB(Flies.hr ~ Biomass*Fence + c.Days.since.start*Fence + c.Days_squared, family = 'nbinom2',
							data = d.lg)
ggplot(d.lg, aes(x = Days.since.start, y = residuals(m2))) +
  geom_point() +
  geom_smooth(method = "loess")


# Compare models
AIC(m1, m2)  # m2 should be much better

plot(d.lg$Days.since.start, residuals(m2))
lines(lowess(d.lg$Days.since.start, residuals(m2)), col="red")

Anova(m2, type = 3)

ggplot(d.lg,aes(x=Biomass,y=Flies.hr,color=Fence)) +
	geom_point() +
	stat_smooth(method="glm")

## --------------- Check the model ---------------------------------------------

# Overall checks
dev.new()
png(filename="Output/1_Fly-surveys/Check-model.png")
check_model(m2)
dev.off()

# VIF
check_collinearity(m2)

# Simulated residuals
sim.m2 <- simulateResiduals(m2)
plot(sim.m2)

# Check overdispersion
E1 <- resid(m2, type = "pearson")
N <- nrow(d.lg)
p <- length(coef(m2)) + 1
overdispersion <- sum(E1^2) / (N-p)
overdispersion # looks good

# Dharma package
testDispersion(sim.m2)
dev.new()
plotResiduals(sim.m2, d.lg$Days.since.start, quantreg = TRUE)
plotResiduals(sim.m2, form = d.lg$Biomass)
testOutliers(sim.m2)
testQuantiles(sim.m2)
testCategorical(sim.m2, catPred = d.lg$Fence)
testUniformity(sim.m2, alternative = c('two.sided'))

# Compare observed and residuals
plot(Flies.hr ~ Fence, data = d.lg)
plot(residuals(m1, type = 'pearson') ~ d.lg$Fence)
abline(a = 0, b = 0, col = "blue", lwd = 2)

plot(residuals(m1, type = 'pearson') ~ d.lg$Biomass)
abline(a = 0, b = 0, col = "blue", lwd = 2)

plot(residuals(m1, type = 'pearson') ~ d.lg$Days.since.start)
abline(a = 0, b = 0, col = "blue", lwd = 2)

# Visualize residuals
plot(density(resid(m1, type='pearson')))
plot(density(resid(m1, type='deviance')))

# Look at predictions 
preds.lm <- predict(m1)
par(mfrow = c(1, 2))
plot(Flies.hr ~ Fence, 
		 data = d.lg)
plot(exp(preds.lm) ~ d.lg$Fence) 

par(mfrow = c(1, 2))
plot(Flies.hr ~ Biomass, 
		 data = d.lg)
plot(exp(preds.lm) ~ d.lg$Biomass) 

par(mfrow = c(1, 2))
plot(Flies.hr ~ Days.since.start, 
		 data = d.lg)
plot(exp(preds.lm) ~ d.lg$Days.since.start) 

plot(predict(m1), residuals(m1, type = 'working'))

# Check cooks
# dev.new()
# png(filename="Output/1_Fly-surveys/Cooks-distance.png")
# plot(cooks.distance(m1), main = 'Fly survey model Cook distance')
# abline(h = 4/nrow(d.lg), lty = 2, col = "steelblue")
# dev.off()
# 5 influential values

# Compare observed and predicted
plot(density(d.lg$Flies.hr), xlim=c(0, 160), ylim=c(0, .08), main='M_0 y_hat')
lines(density(predict(m1, type='response')), col='red')

# Check independence of observations
d.lg$residuals = residuals(m1, type = 'pearson')  # save the residual values
scatter.smooth(d.lg$Days.since.start, d.lg$residuals)
scatter.smooth(d.lg$Biomass, d.lg$residuals)

check_singularity(m2)

## --------------- Model results -----------------------------------------------

# Explained variance
logLik_full <- logLik(m2)
null <- glmmTMB(Flies.hr ~ 1, family = 'nbinom2', data = d.lg)
logLik_null <- logLik(null)
1 - (as.numeric(logLik_full) / as.numeric(logLik_null))

# Effect sizes
mean_day <- mean(d.lg$Days.since.start) # mean
target_day <- 19 # target day
c_day_19 <- target_day - mean_day
c_day_19_sq <- c_day_19^2

emtrends(m2, pairwise ~ Fence, var = "Biomass", 
         at = list(c.Days.since.start = c_day_19, c.Days_squared = c_day_19_sq),
         regrid = "response")

# Days since
coefs <- fixef(m2)$cond
b_days <- coefs["c.Days.since.start"]
b_days_sq <- coefs["c.Days_squared"]

peak_day_centered <- -b_days / (2 * b_days_sq)
peak_day_centered + mean_day # 19 day

# 500 kg
0.1191 * 500 # Fence
0.0593 * 500 # Open

# Rough visualizations
emmip(m2,  ~ Fence~Biomass, mult.name = "Fence", cov.reduce = FALSE)

## --------------- Visualize predict -------------------------------------------

# Mean
mean_day <- mean(d.lg$Days.since.start)

# Target day
target_day <- 19

# Calculate the centered values FOR Day 19
c_day_19 <- target_day - mean_day
c_day_19_sq <- c_day_19^2

# Create pred
pred.dat.midpoint <- expand.grid(
  Biomass = seq(min(d.lg$Biomass), max(d.lg$Biomass), length.out = 100),
  Fence = c('F', 'O'),
  c.Days.since.start = c_day_19,    
  c.Days_squared = c_day_19_sq    
)

preds <- predict(m2, 
                 newdata = pred.dat.midpoint, 
                 se.fit = TRUE, 
                 type = "link")

# Combine predictions
pred.dat.midpoint$fit <- preds$fit
pred.dat.midpoint$se.fit <- preds$se.fit

pred.dat.midpoint$LCL <- pred.dat.midpoint$fit - (1.96 * pred.dat.midpoint$se.fit)
pred.dat.midpoint$UCL <- pred.dat.midpoint$fit + (1.96 * pred.dat.midpoint$se.fit)

# Back-transform
pred.dat.midpoint$fit.resp <- exp(pred.dat.midpoint$fit)
pred.dat.midpoint$LCL.resp <- exp(pred.dat.midpoint$LCL)
pred.dat.midpoint$UCL.resp <- exp(pred.dat.midpoint$UCL)

## --------------- Visualize Fence*Biomass -------------------------------------

biomass.int <- ggplot(data = pred.dat.midpoint, aes(x = Biomass, y = fit.resp, color = Fence)) +
    geom_ribbon(aes(ymin = LCL.resp, ymax = UCL.resp, fill = Fence), 
                alpha = 0.4, color = NA) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = c("#00AFBB", "#E7B800")) +
    geom_jitter(data = d.lg, aes(x = Biomass, y = Flies.hr, 
                                group = Fence, fill = Fence),
                height = 0, width = 20, size = 2, stroke = 0.75, pch = 21,
                color = 'black') +
    scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
    scale_y_continuous(limits = c(0, 250)) +
	  scale_x_continuous(limits = c(0, 800)) +
    theme_classic() +
    theme(legend.position = 'none') +
    ylab("Flies caught per hour") +
    xlab('Carrion biomass (kg)') +
    theme(axis.title = element_text(face = "bold")) +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 25))+
		theme(aspect.ratio = 1.2)
	
## --------------- Visualize Fence*Time ----------------------------------------

mean_day <- mean(d.lg$Days.since.start)
mean_biomass <- mean(d.lg$Biomass)

pred.dat.midbiomass <- expand.grid(
  Days.since.start = seq(8, 30, length.out = 100),
  Fence = c('F', 'O'),
  Biomass = mean_biomass
)

pred.dat.midbiomass$c.Days.since.start <- pred.dat.midbiomass$Days.since.start - mean_day
pred.dat.midbiomass$c.Days_squared <- pred.dat.midbiomass$c.Days.since.start^2

preds <- predict(m2, 
                 newdata = pred.dat.midbiomass, 
                 se.fit = TRUE, 
                 type = "link")

# Combine predictions
pred.dat.midbiomass$fit <- preds$fit
pred.dat.midbiomass$se.fit <- preds$se.fit

# Calculate 95% CI on the LINK (log) scale
pred.dat.midbiomass$LCL <- pred.dat.midbiomass$fit - (1.96 * pred.dat.midbiomass$se.fit)
pred.dat.midbiomass$UCL <- pred.dat.midbiomass$fit + (1.96 * pred.dat.midbiomass$se.fit)

# Back-transform to the RESPONSE (fly count) scale
pred.dat.midbiomass$fit.resp <- exp(pred.dat.midbiomass$fit)
pred.dat.midbiomass$LCL.resp <- exp(pred.dat.midbiomass$LCL)
pred.dat.midbiomass$UCL.resp <- exp(pred.dat.midbiomass$UCL)

days <- ggplot(data = pred.dat.midbiomass, aes(x = Days.since.start, y = fit.resp, color = Fence)) +
  geom_ribbon(aes(ymin = LCL.resp, ymax = UCL.resp, fill = Fence), 
              alpha = 0.4, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  geom_jitter(data = d.lg, aes(x = Days.since.start, y = Flies.hr, 
                              group = Fence, fill = Fence),
              size = 2, stroke = 0.75, pch = 21, height = 0, width = 0.5,
              color = 'black') +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  scale_y_continuous(limits = c(0, 250)) +
  scale_x_continuous(limits = c(5, 31),
                     breaks = c(5, 10, 15, 20, 25, 30)) +
  theme_classic() +
  theme(legend.position = 'none') +
  ylab("") +
  xlab('Days since deployment') +
  theme(axis.title = element_text(face = "bold")) +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 30))
ggsave('Figures/2_fly-surveys-time-ALL-SI.png',width = 7, height = 11, units = 'in', dpi = 300)

## --------------- Visualize Fence*Time at different Biomass levels ----------------

mean_day <- mean(d.lg$Days.since.start)

biomass_levels <- c(24.94, 58.96, 181.41, 362.81, 725.62)
biomass_labels <- c("25 kg", "60 kg", "180 kg", "360 kg", "725 kg")

pred.dat.faceted <- expand.grid(
  Days.since.start = seq(8, 30, length.out = 100), # Smooth time axis
  Fence = c('F', 'O'),
  Biomass = biomass_levels  # 5 user-defined levels for faceting
)
																				 
pred.dat.faceted$c.Days.since.start <- pred.dat.faceted$Days.since.start - mean_day
pred.dat.faceted$c.Days_squared <- pred.dat.faceted$c.Days.since.start^2

pred.dat.faceted$Biomass_label <- factor(pred.dat.faceted$Biomass,
                                         levels = biomass_levels,
                                         labels = biomass_labels)

preds <- predict(m2, 
                 newdata = pred.dat.faceted, 
                 se.fit = TRUE, 
                 type = "link")

pred.dat.faceted$fit <- preds$fit
pred.dat.faceted$se.fit <- preds$se.fit
pred.dat.faceted$LCL <- pred.dat.faceted$fit - (1.96 * pred.dat.faceted$se.fit)
pred.dat.faceted$UCL <- pred.dat.faceted$fit + (1.96 * pred.dat.faceted$se.fit)

pred.dat.faceted$fit.resp <- exp(pred.dat.faceted$fit)
pred.dat.faceted$LCL.resp <- exp(pred.dat.faceted$LCL)
pred.dat.faceted$UCL.resp <- exp(pred.dat.faceted$UCL)

days.int <- ggplot(data = pred.dat.faceted, aes(x = Days.since.start, y = fit.resp, color = Fence)) +
  geom_ribbon(aes(ymin = LCL.resp, ymax = UCL.resp, fill = Fence), 
              alpha = 0.4, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  geom_jitter(data = d.lg, aes(x = Days.since.start, y = Flies.hr, 
                              group = Fence, fill = Fence),
              size = 2, stroke = 0.75, pch = 21, height = 0, width = 0.5,
              color = 'black', alpha = 0.3) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  scale_y_continuous(limits = c(0, 250)) +
  scale_x_continuous(limits = c(5, 31),
                     breaks = c(5, 10, 15, 20, 25, 30)) +
  theme_classic() +
  theme(legend.position = 'none') +
  ylab("") +
  xlab('Days since deployment') +
  theme(axis.title = element_text(face = "bold")) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25)) +
  facet_wrap(~Biomass_label, ncol = 5) + 
  theme(
    strip.text.x = element_text(size = 14, face = "bold"), 
    strip.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    plot.background = element_rect(fill = "transparent", colour = NA)
  )
ggsave('Figures/2_fly-surveys-time-SI.png',width = 7, height = 11, units = 'in', dpi = 300)

## --------------- Visualize coefficients raw ----------------------------------

# Biomass 
emtrends(m2, pairwise ~ Fence, var = "Biomass", 
         at = list(c.Days.since.start = c_day_19, c.Days_squared = c_day_19_sq),
         regrid = "response")

biomass <- tibble(
  Group = c("Fenced", "Open", "All"),
  Metric = c("Biomass Slope (at Day 19)", "Biomass Slope (at Day 19)", "Peak Fly Abundance"),
  Value = c(0.1191, 0.0593, 19.41),
  LCL = c(0.0782, 0.0329, NA),
  UCL = c(0.1601, 0.0857, NA),
  Units = c("Flies/kg", "Flies/kg", "Days")
)

# Days
coefs <- fixef(m2)$cond
b_days <- coefs["c.Days.since.start"]
b_days_sq <- coefs["c.Days_squared"]

peak_day_centered <- -b_days / (2 * b_days_sq)
mean_day <- mean(d.lg$Days.since.start)
peak_day_unscaled <- peak_day_centered + mean_day

## --------------- Visualize coefficients back-transformed ---------------------

emtrends(m1, pairwise~Fence*Days.since.start, var = 'Biomass',
				 regrid = "response")
confint(emtrends(m1, pairwise~Fence*Days.since.start, var = 'Biomass',
								 regrid = "response")
 )

emtrends(m1, pairwise~Biomass*Fence, var = 'Days.since.start')

# Biomass coefficients are back-transformed 
coef <- tibble(Fence = c('Fenced', 'Open'),
							 Trend = c('Biomass', 'Biomass'),
							 Value = c(0.0319, 0.0163),
							 LCL = c(0.02165, 0.00917),
							 UCL = c(0.0422, 0.0235))

biomass.coef <- ggplot(coef, aes(x=Fence,y=Value,color=Fence,fill=Fence))+
	geom_errorbar(aes(ymin=LCL,ymax=UCL),
								position = position_dodge(width = 0.5),color='black', width=0.2)+
	geom_point(size=4.5,shape=21,stroke=2,position=position_dodge(width = 0.5),color='black')+
	scale_color_manual(values=c("#00AFBB", "#E7B800"))+
	scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
	theme_classic()+
	theme(legend.position = 'none')+
	xlab("")+
	ylab('Slope estimate')+
	theme(axis.title = element_text(face="bold"))+
	# annotate('text', x = 2.1, y = 0.04, 
	# 				 label = "p = 0.0236")+
	theme(plot.title = element_text(hjust = 0.5))+
	theme(axis.text = element_text(size = 15),
				axis.title = element_text(size = 17))+
	theme(plot.title = element_text(hjust = 0.5),
				axis.text.x = element_blank(),
				axis.ticks.x = element_blank())

## --------------- Combine figures ---------------------------------------------

library(patchwork)
# Raw biomass
biomass.int+inset_element(biomass.coef, 0.2, 0.6, 0.6, 1, align_to = 'full')
ggsave('Figures/1_fly-surveys.png',width = 10, height = 8, units = 'in', dpi = 300)

# Raw days
dev.new()
days/days.int

# Add a new label to the 'average' prediction data
pred_data_avg <- pred.dat.midbiomass %>% 
  mutate(Biomass_label = "Averaged")

# Combine the 5-level data and the average data
plot_data_combined <- rbind(pred.dat.faceted, pred_data_avg)

# Copy the original raw data and give it the new label
d.lg_avg <- d.lg %>% 
  mutate(Biomass_label = "Averaged")

# Now, we need to add the correct labels to your *original* raw data
# (This assumes your 'biomass_labels' object is still in your environment)
d.lg_faceted <- d.lg %>% 
  mutate(Biomass_label = factor(Biomass,
                                levels = biomass_levels,
                                labels = biomass_labels))

# Combine the raw data
d.lg_combined <- rbind(d.lg_faceted, d.lg_avg)

all_labels <- c("Averaged", "25 kg", "60 kg", "180 kg", "360 kg", "725 kg")

# 2. Re-order the 'Biomass_label' column in BOTH dataframes
plot_data_combined$Biomass_label <- factor(plot_data_combined$Biomass_label, levels = all_labels)
d.lg_combined$Biomass_label      <- factor(d.lg_combined$Biomass_label,      levels = all_labels)

ggplot(data = plot_data_combined, aes(x = Days.since.start, y = fit.resp, color = Fence)) +
  geom_ribbon(aes(ymin = LCL.resp, ymax = UCL.resp, fill = Fence), 
              alpha = 0.4, color = NA) +
  geom_line(linewidth = 1) +
 
  # Use the combined raw data
  geom_jitter(data = d.lg_combined, aes(x = Days.since.start, y = Flies.hr, 
                                        group = Fence, fill = Fence),
              size = 2, stroke = 0.75, pch = 21, height = 0, width = 0.5,
              color = 'black', alpha = 0.3) +
 
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  scale_y_continuous(limits = c(0, 250)) +
  scale_x_continuous(limits = c(5, 31), breaks = c(5, 10, 15, 20, 25, 30)) +
  theme_classic() +
  ylab("Flies caught per hour") + 
  xlab('Days since deployment') +
 
  # Tell facet_wrap to make 2 columns (for a 3x2 grid)
  facet_wrap(~Biomass_label, ncol = 2) + 
 
  # --- Updated Theme Section ---
  theme(
    legend.position = 'none',
    aspect.ratio = 1, # Added from your other plot
    
    # Text sizes from your other plot
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 25, face = "bold"),
    
    # Facet theme elements
    strip.text.x = element_text(size = 20, face = "bold"), # Matched to axis.text
    strip.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

ggsave('Figures/2_fly-surveys-time-SI.png',width = 10, height = 8, units = 'in', dpi = 300)


## --------------- Export model information ------------------------------------

library(broom)       # For tidy() and glance()
library(htmlTable)   # For htmlTable()
library(kableExtra)  # For save_kable()
library(dplyr)       # For filtering
library(broom.mixed) # For glmmmTMB

# Parameter estimates
m2.sum <- tidy(m2) %>%
  filter(component == "cond") # <-- This is the main change needed

m2.sum$estimate <- as.numeric(round(m2.sum$estimate, 2))
m2.sum$std.error <- as.numeric(round(m2.sum$std.error, 2))
m2.sum$statistic <- as.numeric(round(m2.sum$statistic, 2))
m2.sum$p.value <- as.numeric(round(m2.sum$p.value, 3))

# Select and rename columns to your preference
m2.sum <- m2.sum %>%
  dplyr::select(term, estimate, std.error, statistic, p.value)

colnames(m2.sum) <- c('Term', 'Estimate', 'SE', 'Z value', 'p value')

# Save the table
htmlTable(m2.sum, align = 'l') %>%
  save_kable(file = 'Output/1_Fly-surveys/Model-summary-m2.png')

# Other model information
library(webshot2)
m2.glance <- glance(m2) %>% round(2)

htmlTable(m2.glance) %>%
  save_kable(file = 'Output/1_Fly-surveys/Model-glance-m2.png')
