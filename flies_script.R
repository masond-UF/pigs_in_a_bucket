############ 18 JUNE 2021—DAVID MASON—GARBAGE FLY ABUNDANCE ANALYSIS ###########
############ LOAD PACKAGES AND BRING IN DATA ################################### 
d <- read.csv("Data/flies.csv")
head(d)
summary(d)

library(tidyverse)
library(lme4)
library(lubridate)
library(MuMIn)

d$Date <- paste(d$Year, d$Month, d$Day, sep="-") %>% ymd() %>% as.Date()
############ EXPLORE DATA VISUALLY ############################################# 
ggplot(d, aes(y=AVGrate, x=Date, color=Fence))+
	geom_point()+
	facet_wrap(~Biomass)
############ BUILD A MODEL ##################################################### 
mod1 <- lm(AVGrate ~ Date * Biomass * Fence, data = d)
plot(mod1)
hist(mod1$residuals, breaks=10)
anova(mod1)

library(car)
mod2 <- lmer(AVGrate ~ Biomass * Fence + (1|Date), data = d)
Anova(mod2)
############ VISUALIZE THE DATE ################################################
library(emmeans)
emmeans(mod, pairwise~Fence)
emtrends(mod, pairwise~Fence, var = 'Biomass') 

d.NA <- d[!is.na(d$AVGrate),]

sem <- function(x, na.rm = FALSE) {
  out <-sd(x, na.rm = na.rm)/sqrt(length(x))
  return(out)
}

d.means <- d.NA %>% 
							group_by(Biomass,Fence) %>% 
							summarise(mean = mean(AVGrate),
								sd = sd(AVGrate),
								se = sem(AVGrate))


log(55) # 4.007333
log(130) # 4.867534
log(400) # 5.991465
log(800) # 6.684612 
log(1600) # 7.377759 
log(0.00001)

levels(d.means$Fence) <- c("Fenced","Open")

ggplot(d.means, aes(x=log(Biomass),y=mean,color=Fence,fill=Fence))+
	geom_errorbar(aes(ymin=mean-se,ymax=mean+se),
								position = position_dodge(width = 0.5),color='black', width=0.2)+
	geom_point(size=6.5,shape=21,stroke=2,position=position_dodge(width = 0.5),color='black')+
	scale_color_manual(values=c("#00AFBB", "#E7B800"))+
	scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
	scale_x_continuous(name="Biomass treatment (kg)",
										 breaks=c(4.007333,4.867534,5.991465, 6.684612, 7.377759),
										 labels=c("55", "130", "400", "800", "1600"))+
	scale_y_continuous(name="Mean flies per minute",
										 breaks=c(0.0,0.25,0.50,0.75,1.00,1.25,1.5))+
	theme_classic()+
	theme(legend.title = element_blank())+
	theme(legend.position = c(0.18, 0.92))+
	theme(text = element_text(size = 25))+
	theme(axis.title.x = element_text(face="bold"))+
	theme(axis.title.y = element_text(face="bold", vjust=0.7))+
	 guides(fill=guide_legend(
                 keyheight=0.4,
                 default.unit="inch"))

ggsave('flies.jpg', width = 6, height = 6, dpi = 300)
############ REPEATED MODEL DROPPING OUTLIER ###################################
d.tmp <- filter(d, !Date=='2016-07-13')

mod3 <- lm(AVGrate ~ Date * Biomass * Fence, data = d.tmp)
plot(mod3)
hist(mod3$residuals, breaks=10)
anova(mod3)

library(DHARMa)
mod4 <- lmer(AVGrate ~ Biomass * Fence + (1|Date), data = d.tmp)
mod4.sim <- simulateResiduals(fittedModel = mod4, n = 250)
plot(mod4.sim)

Anova(mod4)

## test for overdispersion (ie. more variance in the data than expected)
testDispersion(mod4.sim)

## test for zero-inflation (ie. more zeros in the data than expected)
testZeroInflation(mod4.sim)

## emmeans, correcting for bias (only needed for GLMMs) by calculating the total SD of the random effect variance components
modSD <- VarCorr(mod4) %>% as.data.frame() %>% summarize(totSD=sum(vcov)) %>% mutate(totSD=sqrt(totSD))

emmeans(mod4, pairwise~Fence, type="response", bias.adj=T, sigma=modSD$totSD)
# flies caught per minute is 0.204 (SE = 0.0805) higher in fence than control
# fence mean is 0.596 (LCL = 0.231, UCL = 0.961)
# open mean is 0.392 (LCL = 0.027, UCL = 0.758)
emtrends(mod4, pairwise~Fence, var = 'Biomass') 
# each 1 kg of biomass increased mean flys caught per minute by 0.000451
# (SE = 0.000139) for fence in comparison to control.
# fence per unit (kg) increase in flies is 0.000804 (LCL = 0.000603, UCL = 0.001004)
# open per unit (kg) increase in flies is 0.000353 (LCL = 0.000152, UCL = 0.000554)

r.squaredGLMM(mod4) # fixed efects = 0.557324, entire model = 0.759619
VarCorr(mod4) # random effects = 0.22769, unexplained variance = 0.24820
############ REPEATED VISUALIZE DROPPING OUTLIER ###################################
library(emmeans)

d.NA <- d.tmp[!is.na(d.tmp$AVGrate),]

sem <- function(x, na.rm = FALSE) {
  out <-sd(x, na.rm = na.rm)/sqrt(length(x))
  return(out)
}

d.means <- d.NA %>% 
							group_by(Biomass,Fence) %>% 
							summarise(mean = mean(AVGrate*60),
								sd = sd(AVGrate*60),
								se = sem(AVGrate*60))


log(55) # 4.007333
log(130) # 4.867534
log(400) # 5.991465
log(800) # 6.684612 
log(1600) # 7.377759 
log(0.00001)

levels(d.means$Fence) <- c("Fenced","Open")

ggplot(d.means, aes(x=log(Biomass),y=mean,color=Fence,fill=Fence))+
	geom_errorbar(aes(ymin=mean-se,ymax=mean+se),
								position = position_dodge(width = 0.5),color='black', width=0.2)+
	geom_point(size=6.5,shape=21,stroke=2,position=position_dodge(width = 0.5),color='black')+
	scale_color_manual(values=c("#00AFBB", "#E7B800"))+
	scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
	scale_x_continuous(name="Biomass treatment (kg)",
										 breaks=c(4.007333,4.867534,5.991465, 6.684612, 7.377759),
										 labels=c("55", "130", "400", "800", "1600"))+
	scale_y_continuous(name="Mean flies per hour",
										 breaks=c(0.0*60,0.25*60,0.50*60,0.75*60,1.00*60,1.25*60,
										 				 1.5*60,1.75*60,2*60),
										 limits = c(0,2*60))+
	theme_classic()+
	theme(legend.title = element_blank())+
	theme(legend.position = c(0.17, 0.92))+
	theme(text = element_text(size = 25))+
	theme(axis.title.x = element_text(face="bold"))+
	theme(axis.title.y = element_text(face="bold", vjust=0.7))+
	 guides(fill=guide_legend(
                 keyheight=0.4,
                 default.unit="inch"))

ggsave('flies_drop_first.jpg', width = 6, height = 6, dpi = 300)
