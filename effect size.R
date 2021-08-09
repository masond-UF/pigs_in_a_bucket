# effect size

d <- read.csv("es.csv")

is.numeric(d$effect.magnitude)

ggplot(d,aes(y=effect.magnitude,x=log(Treatment)))+
	geom_point(size = 6.5,shape=21,stroke=2)+
	geom_smooth(method="lm",se=FALSE, color = "black")+
	scale_x_continuous(name="Carrion biomass (kg)",
										 breaks=c(4.007333,5.991465,7.377759),
										 labels=c("55","400","1600"))+
	scale_y_continuous(name="Effect magnitude of vertebrate scavengers",
										 breaks = c(0.1,0.15,0.2,0.25,0.3,0.35),
										 limits = c(0.1,0.35))+
	theme_classic()+
	theme(legend.title = element_blank())+
	theme(legend.position = c(0.17, 0.92))+
	theme(text = element_text(size = 25))+
	theme(axis.title.x = element_text(face="bold"))+
	theme(axis.title.y = element_text(face="bold", vjust=0.7, size = 18))+
		 guides(fill=guide_legend(
                 keyheight=0.4,
                 default.unit="inch"))

ggsave('effect_size.jpg', width = 6, height = 6, dpi = 300)
