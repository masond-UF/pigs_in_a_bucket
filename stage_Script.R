# Decomposition stage

library(ggplot2)
library(lmPerm)
library(car)
library(lme4)

d <- read.csv("decomp_stage.csv")

d$Fence <- c("Open", "Open", "Open",
						"Fenced","Fenced", "Fenced")

d$Fence <- as.factor(d$Fence)
is.factor(d$Fence)

ggplot(d,aes(y=Stage,x=log(Pig.biomass),color=Fence,fill=Fence))+
	geom_point(size = 6.5,shape=21,stroke=2)+
	scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
	scale_color_manual(values=c("black","black"))+
	scale_x_continuous(name="Carrion biomass (kg)",
										 breaks=c(4.007333,5.991465,7.377759),
										 labels=c("55","400","1600"))+
	scale_y_continuous(name="Decomposition stage",
										 breaks=c(0,1,2,3,4,5),
										 labels=c("0","1","2","3","4","5"))+
	theme_classic()
	theme(legend.title = element_blank())+
	theme(legend.position = c(0.17, 0.92))+
	theme(text = element_text(size = 25))+
	theme(axis.title.x = element_text(face="bold"))+
	theme(axis.title.y = element_text(face="bold", vjust=0.7))+
		 guides(fill=guide_legend(
                 keyheight=0.4,
                 default.unit="inch"))

ggsave('decomp_stage.jpg', width = 6, height = 6, dpi = 300)

require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)

m <- polr(Stage ~ Fence + Pig.biomass, data = d, Hess=TRUE)
summary(m)
