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

# Days since start of experiment [15 days would be 2016-06-28]
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

m1 <- glm.nb(Flies.hr ~ Biomass*Fence*Days.since.start, data = d.lg)
Anova(m1, type = 2)

ggplot(d.lg,aes(x=Biomass,y=Flies.hr,color=Fence)) +
	geom_point() +
	stat_smooth(method="glm")
## --------------- Check the model ---------------------------------------------

# Overall checks
dev.new()
png(filename="Output/1_Fly-surveys/Check-model.png")
check_model(m1)
dev.off()

# VIF
m1.VIF <- tidy(vif(m1))
colnames(m1.VIF)[1] <- 'Term'
colnames(m1.VIF)[2] <- 'VIF'
m1.VIF$VIF <- round(m1.VIF$VIF, digits = 1)
htmlTable(m1.VIF) %>%
	save_kable(file = 'Output/1_Fly-surveys/VIF.png')

sim.m1 <- simulateResiduals(m1)
plot(sim.m1)

# Check overdispersion
E1 <- resid(m1, type = "pearson")
N <- nrow(d.lg)
p <- length(coef(m1)) + 1
overdispersion <- sum(E1^2) / (N-p)
overdispersion # looks good

# Dharma package
testDispersion(sim.m1)
dev.new()
plotResiduals(sim.m1, d.lg$Days.since.start, quantreg = TRUE)
plotResiduals(sim.m1, form = d.lg$Biomass)
testOutliers(sim.m1)
testQuantiles(sim.m1)
testCategorical(sim.m1, catPred = d.lg$Fence)
testUniformity(sim.m1, alternative = c('two.sided'))

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
dev.new()
png(filename="Output/1_Fly-surveys/Cooks-distance.png")
plot(cooks.distance(m1), main = 'Fly survey model Cook distance')
abline(h = 4/nrow(d.lg), lty = 2, col = "steelblue")
dev.off()
# 5 influential values

# Compare observed and predicted
plot(density(d.lg$Flies.hr), xlim=c(0, 160), ylim=c(0, .08), main='M_0 y_hat')
lines(density(predict(m1, type='response')), col='red')

# Check independence of observations
d.lg$residuals = residuals(m1, type = 'pearson')  # save the residual values
scatter.smooth(d.lg$Days.since.start, d.lg$residuals)
scatter.smooth(d.lg$Biomass, d.lg$residuals)

## --------------- Model results -----------------------------------------------

1 - (m1$deviance/m1$null.deviance) # 34.8%

emmeans(m1, pairwise~Fence, type ='response')
confint(emmeans(m1, pairwise~Fence, type ='response'))

emtrends(m1, pairwise~Fence|Days.since.start, var = 'Biomass',
				 regrid = c("response"))
confint(emtrends(m1, pairwise~Fence|Days.since.start, var = 'Biomass',
				 regrid = c("response")))

emtrends(m1, pairwise~Fence|Biomass, var = 'Days.since.start',
				 regrid = c("response"))
confint(emtrends(m1, pairwise~Fence|Biomass, var = 'Days.since.start',
				 regrid = c("response")))

# Rough visualizations
emmip(m1,  ~ Fence~Biomass, mult.name = "Fence", cov.reduce = FALSE)

ggplot(data=d.lg, aes(x = Biomass, y = log(flies.hr), group = Fence, color = Fence))+
	geom_jitter()+
	geom_abline(aes(slope=0.001551,intercept=log(2.9568697985)), color = 'red')+
	geom_abline(aes(slope=0.000983,intercept=log(2.9568697985-0.1531859253)), 
							color = 'blue')+
	theme_bw() 

## --------------- Visualize predict -------------------------------------------

# Create new data frame to predict from
pred.dat <- data.frame(Days.since.start = rep(seq(8,30,1),3092),
											 Biomass = rep(seq(55,1600,1), 46),
											 Fence = c(rep('F', 35558),rep('O',35558)))

# Predictions from model
preds <- predict(m1, 
								 newdata = pred.dat, 
								 se = T)

# Combine predictions to new data frame for plotting
pred.dat <- cbind(pred.dat, fit = preds$fit)
pred.dat <- cbind(pred.dat, se.fit = preds$se.fit)

# Calculate 95% CI for predictions from predicted standard errors
pred.dat$LCL <- pred.dat$fit - (1.96*pred.dat$se.fit) # Correct?
pred.dat$UCL <- pred.dat$fit + (1.96*pred.dat$se.fit)

## --------------- Visualize Fence*Biomass -------------------------------------

# Average by date
pred.dat.avg.date <- pred.dat |> 
	group_by(Biomass, Fence)|>
	summarize(fit = mean(fit),
						LCL = mean(LCL),
						UCL = mean(UCL))

dev.new()
biomass.int <- ggplot(data=pred.dat.avg.date, aes(x = Biomass, y = exp(fit), color = Fence))+
	geom_ribbon(aes(ymin = exp(LCL), ymax = exp(UCL), fill = Fence), 
							alpha = 0.4, color = NA)+
	geom_line()+
	scale_color_manual(values=c("#00AFBB", "#E7B800"))+
	geom_jitter(data=d.lg, aes(x = Biomass, y = Flies.hr, 
														 group = Fence, fill = Fence),
							height = 0, width = 80, size = 2, stroke = 0.75, pch = 21,
							color = 'black')+
	scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
	scale_y_continuous(limits = c(0,200))+
	theme_classic()+
	theme(legend.position = 'none')+
	ylab("Flies caught per hour")+
	xlab('Carrion biomass')+
	theme(axis.title = element_text(face="bold"))+
	theme(axis.text = element_text(size = 20),
				axis.title = element_text(size = 25))
	
## --------------- Visualize Fence*Time ----------------------------------------

# Average by biomass
pred.dat.avg.biomass <- pred.dat |> 
	group_by(Days.since.start, Fence)|>
	summarize(fit = mean(fit),
						LCL = mean(LCL),
						UCL = mean(UCL))


# Days alone no interaction with weight
days <- ggplot(data=pred.dat.avg.biomass, aes(x = Days.since.start, y = exp(fit), color = Fence))+
	geom_ribbon(aes(ymin = exp(LCL), ymax = exp(UCL), fill = Fence), 
							alpha = 0.4, color = NA)+
	geom_line()+
	scale_color_manual(values=c("#00AFBB", "#E7B800"))+
	geom_jitter(data=d.lg, aes(x = Days.since.start, y = Flies.hr, 
														 group = Fence, fill = Fence),
							size = 2, stroke = 0.75, pch = 21, height = 0, width = 0.5,
							color = 'black')+
	scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
	scale_y_continuous(limits = c(0,200))+
	scale_x_continuous(limits = c(5,31),
										 breaks = c(5,10,15,20,25, 30))+
	theme_classic()+
	theme(legend.position = 'none')+
	ylab("")+
	xlab('Days since deployment')+
	theme(axis.title = element_text(face="bold"))+
	theme(axis.text = element_text(size = 25),
				axis.title = element_text(size = 30))
ggsave('Figures/2_fly-surveys-time-ALL-SI.png',width = 7, height = 11, units = 'in', dpi = 300)

dev.new()
days.int <- ggplot(data=pred.dat.avg.biomass, aes(x = Days.since.start, y = exp(fit), color = Fence))+
	geom_ribbon(aes(ymin = exp(LCL), ymax = exp(UCL), fill = Fence), 
							alpha = 0.4, color = NA)+
	geom_line()+
	scale_color_manual(values=c("#00AFBB", "#E7B800"))+
	geom_jitter(data=d.lg, aes(x = Days.since.start, y = Flies.hr, 
														 group = Fence, fill = Fence),
							size = 2, stroke = 0.75, pch = 21, height = 0, width = 0.5,
							color = 'black')+
	scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
	scale_y_continuous(limits = c(0,200))+
	scale_x_continuous(limits = c(5,31),
										 breaks = c(5,10,15,20,25, 30))+
	theme_classic()+
	theme(legend.position = 'none')+
	ylab("")+
	xlab('Days since deployment')+
	theme(axis.title = element_text(face="bold"))+
	theme(axis.text = element_text(size = 20),
				axis.title = element_text(size = 25))+
	facet_wrap(~Biomass, ncol = 2)+
	theme(
		strip.text.x = element_blank(),
		panel.border = element_rect(colour = "black", fill=NA, linewidth = 1),
		plot.background = element_rect(fill = "transparent", colour = NA)
	)

## --------------- Visualize coefficients raw --------------------------------------

emtrends(m1, pairwise~Fence*Days.since.start, var = 'Biomass',
				 type = "response")

emtrends(m1, pairwise~Biomass*Fence, var = 'Days.since.start')

# Coefficients are not back-transformed 
coef <- tibble(Fence = c('Fenced', 'Open', 'Fenced', 'Open'),
							 Trend = c('Biomass', 'Biomass', "Days", "Days"),
							 Value = c(0.001551, 0.000983, -0.0343, -0.0179),
							 LCL = c(0.001172, 0.000604, -0.0637, -0.0474),
							 UCL = c(0.00193, 0.00136, -0.00485, 0.01155))

# Biomass interpretation
(exp(0.001551)-1)*100
# For every 1 kg in biomass, flies increases by 0.16% in fenced plots
# For every 100 kg in biomass, flies increase by 16% in fenced plots
# For every 500 kg in biomass, flies increase by 78% in fenced plots
# At 0 kg flies = 10, at 500 kg flies = 18
# At 1000 kg flies = 40, at 1500 kg flies = 75

(exp(0.000983)-1)*100
# For every 1 kg in biomass, flies increases by 0.10% in fenced plots
# For every 100 kg in biomass, flies increase by 10% in fenced plots
# For every 500 kg in biomass, flies increase by 49% in fenced plots
# At 0 kg flies = 10, at 500 kg flies = 18
# At 1000 kg flies = 25, at 1500 kg flies = 75

# Days interpretation
(exp(-0.0343)-1)*100
# For every 1 day, flies decreased by 3% in fenced plots
# For every 5 days, flies decreased by 17% in fenced plots
# For every 10 days, flies decreased by 34% in fenced plots
# At 4 days flies = 45, at 20 days flies = 30

(exp(-0.0179)-1)*100
# For every 1 day, flies decreased by 2% in fenced plots
# For every 5 days, flies decreased by 9% in fenced plots
# For every 10 days, flies decreased by 18% in fenced plots
# At 4 days flies = 30, at 20 days flies = 25

biomass <- coef |> filter(Trend == 'Biomass')
days <- coef |> filter(Trend == 'Days')

# Biomass interaction
dev.new()
biomass.coef <- ggplot(biomass, aes(x=Fence,y=Value,color=Fence,fill=Fence))+
	geom_errorbar(aes(ymin=LCL,ymax=UCL),
								position = position_dodge(width = 0.5),color='black', width=0.2)+
	geom_point(size=6.5,shape=21,stroke=2,position=position_dodge(width = 0.5),color='black')+
	scale_y_continuous(limits = c(0.0005,0.002))+
	scale_color_manual(values=c("#00AFBB", "#E7B800"))+
	scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
	theme_classic()+
	theme(legend.position = 'none')+
	xlab("")+
	ylab('Estimate')+
	theme(axis.title = element_text(face="bold"))+
	annotate('text', x = 2.1, y = 0.002, 
					 label = "p = 0.038")+
	theme(plot.title = element_text(hjust = 0.5))+
	theme(axis.text = element_text(size = 10),
				axis.title = element_text(size = 15))+
	theme(plot.title = element_text(hjust = 0.5),
				axis.text.x = element_blank(),
				axis.ticks.x = element_blank())

# Time interaction
dev.new()
days.coef <- ggplot(days, aes(x=Fence,y=Value,color=Fence,fill=Fence))+
	geom_errorbar(aes(ymin=LCL,ymax=UCL),
								position = position_dodge(width = 0.5),color='black', width=0.2)+
	geom_point(size=6.5,shape=21,stroke=2,position=position_dodge(width = 0.5),color='black')+
	scale_y_continuous(limits = c(-0.07,0.03),
										 breaks = c(-0.06, -0.03, 0.00,
										 					 0.03))+
	scale_color_manual(values=c("#00AFBB", "#E7B800"))+
	scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
	theme_classic()+
	theme(legend.position = 'none')+
	xlab("")+
	ylab('Estimate')+
	theme(axis.title = element_text(face="bold"))+
	annotate('text', x = 0.7, y = 0.03, 
					 label = "p = 0.442")+
	theme(plot.title = element_text(hjust = 0.5))+
	theme(axis.text = element_text(size = 10),
				axis.title = element_text(size = 15),
				axis.text.x = element_blank(),
				axis.ticks.x = element_blank())

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

# Biomass interaction
dev.new()

biomass.coef <- ggplot(coef, aes(x=Fence,y=Value,color=Fence,fill=Fence))+
	geom_errorbar(aes(ymin=LCL,ymax=UCL),
								position = position_dodge(width = 0.5),color='black', width=0.2)+
	geom_point(size=6.5,shape=21,stroke=2,position=position_dodge(width = 0.5),color='black')+
	scale_color_manual(values=c("#00AFBB", "#E7B800"))+
	scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
	theme_classic()+
	theme(legend.position = 'none')+
	xlab("")+
	ylab('Estimate')+
	theme(axis.title = element_text(face="bold"))+
	annotate('text', x = 2.1, y = 0.04, 
					 label = "p = 0.015")+
	theme(plot.title = element_text(hjust = 0.5))+
	theme(axis.text = element_text(size = 12.5),
				axis.title = element_text(size = 15))+
	theme(plot.title = element_text(hjust = 0.5),
				axis.text.x = element_blank(),
				axis.ticks.x = element_blank())

# Time interaction

emtrends(m1, pairwise~Fence*Biomass, var = 'Days.since.start',
				 regrid = "response")
confint(emtrends(m1, pairwise~Fence*Biomass, var = 'Days.since.start',
								 regrid = "response")
)

# Day coefficients are back-transformed 
coef <- tibble(Fence = c('Fenced', 'Open'),
							 Trend = c('Biomass', 'Biomass'),
							 Value = c(-0.704, -0.297),
							 LCL = c(-1.33, -0.79),
							 UCL = c(-0.0834, 0.1950))

dev.new()
days.coef <- ggplot(coef, aes(x=Fence,y=Value,color=Fence,fill=Fence))+
	geom_errorbar(aes(ymin=LCL,ymax=UCL),
								position = position_dodge(width = 0.5),color='black', width=0.2)+
	geom_point(size=6.5,shape=21,stroke=2,position=position_dodge(width = 0.5),color='black')+
	scale_color_manual(values=c("#00AFBB", "#E7B800"))+
	scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
	scale_y_continuous(limits = c(-1.4,0.6))+
	theme_classic()+
	theme(legend.position = 'none')+
	xlab("")+
	ylab('Estimate')+
	theme(axis.title = element_text(face="bold"))+
	annotate('text', x = 0.7, y = 0.5, 
					 label = "p = 0.442")+
	theme(plot.title = element_text(hjust = 0.5))+
	theme(axis.text = element_text(size = 10),
				axis.title = element_text(size = 15),
				axis.text.x = element_blank(),
				axis.ticks.x = element_blank())

## --------------- Combine figures ---------------------------------------------

# Raw biomass
dev.new()
biomass.int+inset_element(biomass.coef, 0.2, 0.6, 0.6, 1, align_to = 'full')
ggsave('Figures/1_fly-surveys.png',width = 12, height = 10, units = 'in', dpi = 300)

# New
dev.new()
biomass.int+inset_element(biomass.coef, 0.2, 0.6, 0.6, 1, align_to = 'full')
# ggsave('Figures/1_fly-surveys.png', width = 6, height = 12, units = 'in', dpi = 300)

# Raw days
dev.new()
days.int+inset_element(days.coef, 0.5445, 0.04, 0.95, 0.355, align_to = 'full')
ggsave('Figures/2_fly-surveys-time-SI.png',width = 12, height = 10, units = 'in', dpi = 300)

# New
dev.new()
days.int+inset_element(days.coef, 0.5445, 0.04, 0.95, 0.355, align_to = 'full')
# ggsave('Figures/2_fly-surveys-time-SI.png',width = 12, height = 10, units = 'in', dpi = 300)



## --------------- Export model information ------------------------------------

# Parameter estimates
m1.sum <- tidy(m1)
m1.sum$estimate <- as.numeric(round(m1.sum$estimate,2))
m1.sum$std.error <- as.numeric(round(m1.sum$std.error,2))
m1.sum$statistic <- as.numeric(round(m1.sum$statistic,2))
m1.sum$p.value <- as.numeric(round(m1.sum$p.value,3))
colnames(m1.sum) <- c('Term', 'Estimate', 'SE', 'Z value', 'p value')
htmlTable(m1.sum, align = 'l') %>%
	save_kable(file = 'Output/1_Fly-surveys/Model-summary.png')

# Other model information
m1.glance <- glance(m1) |> round(2)
htmlTable(m1.glance) %>%
	save_kable(file = 'Output/1_Fly-surveys/Model-glance.png')