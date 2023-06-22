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

# Fix the labels
d$FENCE <- c("Open","Open","Open",
						 "Fenced","Fenced","Fenced")

# Calculate the change in biomass
d <- d %>% 
	mutate(LOST_BIO = abs(PIG_FINAL-PIG_START))

## --------------- Explore the data --------------------------------------------

descdist(d$LOST_BIO, discrete = FALSE) # beta or log-normal
descdist(log(d$LOST_BIO), discrete = FALSE) # closer to beta or log-normal
descdist(log10(d$LOST_BIO), discrete = FALSE) # closer to beta or log-normal

## --------------- Create the model --------------------------------------------

mod <- lm(LOST_BIO~FENCE*EXPOSURE, d)
Anova(lm(LOST_BIO~FENCE*EXPOSURE, d))
hist(mod$residuals)
plot(mod$residuals)
summary(mod)

dev.new()
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

## --------------- Visualize the model -----------------------------------------

# Create new data frame to predict from
pred.dat <- data.frame(EXPOSURE = c(55,400,1600,55,400,1600),
											 FENCE = c('Open','Open','Open','Fenced','Fenced','Fenced'))

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

# Old
ggplot(d,aes(y=LOST_BIO,x=log(EXPOSURE),color=FENCE,fill=FENCE))+
	geom_point(size = 6.5,shape=21)+
	scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
	scale_color_manual(values=c("black","black"))+
	stat_smooth(method="lm",se=FALSE,show.legend = FALSE)+
	scale_x_continuous(name="Carrion exposure (kg)",
										 breaks=c(4.007333,5.991465,7.377759),
										 labels=c("55","400","1600"))+
	ylab("Biomass reduction (kg)")+
	theme_classic()+
	theme(legend.title = element_blank())+
	theme(legend.position = c(0.17, 0.92))+
	theme(text = element_text(size = 25))+
	theme(axis.title.x = element_blank())+
	theme(axis.title.y = element_text(face="bold", vjust=0.7))+
		 guides(fill=guide_legend(
                 keyheight=0.4,
                 default.unit="inch"))

# ggsave('garbage_can.jpg', width = 6, height = 6, dpi = 300)

# Main figure
p1 <- ggplot(data=pred.dat, aes(x = EXPOSURE, y = fit, color = FENCE))+
	geom_line(linewidth = 1.75)+
	# geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = FENCE), 
							# alpha = 0.4, color = NA)+
	scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
	scale_color_manual(values = c('black', 'black'))+
	geom_point(data = d, aes(x = EXPOSURE, y = LOST_BIO, fill = FENCE), 
						 size = 8, shape = 21, color = 'black', stroke = 2)+
	scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
	scale_x_continuous(breaks = c(0, 250, 500, 750, 1000,
																1250, 1500))+
	ylab("Biomass loss after X days (kg)")+
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
	scale_y_continuous(limits = c(0.001,0.005))+
	scale_color_manual(values=c("#00AFBB", "#E7B800"))+
	scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
	theme_classic()+
	theme(legend.position = 'none')+
	xlab("")+
	ylab('Slope estimate')+
	theme(axis.title = element_text(face="bold"))+
	annotate('text', x = 1.3, y = 0.005, 
					 label = "p = 0.022")+
	theme(plot.title = element_text(hjust = 0.5))+
	theme(axis.text = element_text(size = 10),
				axis.title = element_text(size = 15))+
	theme(plot.title = element_text(hjust = 0.5),
				axis.text.x = element_blank(),
				axis.ticks.x = element_blank())

# Combine the figures
dev.new()
comb <- p1+inset_element(p2, 0.14, 0.64, 0.45, 0.97, align_to = 'full')

# Output the ggplot object
saveRDS(comb, file = "Output/2b_future-biomass-fig.RDS")


