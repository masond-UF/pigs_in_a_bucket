############ 7 MAY 2021—DAVID MASON—GARBAGE CAN PIG DECOMP ANALYSIS ############
############ LOAD PACKAGES AND BRING IN DATA ###################################  
d.tmp <- read.csv("Data/garbage_can.csv")
head(d)
summary(d)

library(tidyverse)
library(RVAideMemoire)
library(simpleboot)
library(sm)
library(lmPerm)
library(multcomp)
############ CALCULATE CHANGE IN BIOMASS #######################################
d <- d.tmp %>% 
	mutate(LOST_BIO = abs(PIG_FINAL-PIG_START))

############ MODEL CHANGE IN BIOMASS ###########################################

# Basic model
mod <- lm(LOST_BIO~FENCE*EXPOSURE, d)
anova(lm(LOST_BIO~FENCE*EXPOSURE, d))

# Bootstrapped model
lm.mod <- lm.boot(mod, R = 1000, rows = FALSE) 
summary(lm.mod)


set.seed(600)

# Two-way anova 
perm.anova(LOST_BIO~FENCE+EXPOSURE, data = d, nperm = 1000)

## Ancova
mod3 <- lmp(LOST_BIO~FENCE*EXPOSURE, data=d)
summary(mod3) 
summary(mod3)$r.squared
Anova(mod3)

# no intercept?
# Fence increased biomass lost by 1.19 kg
# 1 kg increase in biomass exposure increased biomass lost by 0.003 kg
# That slope is 0.001 higher in Fence (0.004 kg) 

sem <- function(x, na.rm=FALSE){
		out <- sd(x, na.rm = na.rm)/sqrt(length(x))
		return(out)
}

d %>% group_by(FENCE) %>% 
	summarise(mean = mean(LOST_BIO),
						sd = mean(LOST_BIO),
						se = sem(LOST_BIO))

d %>% group_by(EXPOSURE) %>% 
	summarise(mean = mean(LOST_BIO),
						sd = mean(LOST_BIO),
						se = sem(LOST_BIO))

# Just to demonstrate that aovp is the same as lmp/anova
mod4 <- aovp(LOST_BIO~FENCE*EXPOSURE, data=d)
summary(mod4)

############ VISUALIZE MODEL ###################################################
levels(d$FENCE) <- c("Fenced","Open")

ggplot(d,aes(y=LOST_BIO,x=log(EXPOSURE),color=FENCE,fill=FENCE))+
	geom_point(size = 6.5,shape=21,stroke=2)+
	scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
	scale_color_manual(values=c("black","black"))+
	stat_smooth(method="lm",se=FALSE,show.legend = FALSE)+
	scale_x_continuous(name="Biomass treatment (kg)",
										 breaks=c(4.007333,5.991465,7.377759),
										 labels=c("55","400","1600"))+
	ylab("Lost biomass (kg)")+
	theme_classic()+
	theme(legend.title = element_blank())+
	theme(legend.position = c(0.17, 0.92))+
	theme(text = element_text(size = 25))+
	theme(axis.title.x = element_text(face="bold"))+
	theme(axis.title.y = element_text(face="bold", vjust=0.7))+
		 guides(fill=guide_legend(
                 keyheight=0.4,
                 default.unit="inch"))

ggsave('garbage_can.jpg', width = 6, height = 6, dpi = 300)

