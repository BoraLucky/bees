rm(list=ls())
library(ggplot2)
library(ggpubr)
library(dplyr)
args <- commandArgs(TRUE)
statistics_file <- args[1]
file_table <- read.table(statistics_file,header=T)


dxy_pi_plot <- 	file_table  %>% 
		ggplot(aes(avg_dxy,avg_pi.pop,color=avg_wc_fst)) +
		geom_point(size=0.1) +
		facet_wrap(~population,nrow=2)+
		scale_color_gradientn(colors = c("grey50", "#BB2020"),limits = c(0, 1)) +
		theme_classic() +  
		theme(strip.background = element_rect(color = "white"))+  
		theme(plot.title = element_text(hjust = 0.5),legend.position = "top",legend.key.width=unit(0.5,'cm'),legend.key.size=unit(0.1,'cm'))+ 
		geom_smooth(method = "lm",formula = y ~ x,se = FALSE,size=0.3)+ 
		stat_cor(method="spearman")+  
		xlab("Dxy")+                   
		ylab("Pi")+
		ylim(0, 0.04)+       
		xlim(0, 0.04)+          
		geom_segment(aes(x =0, y = 0, xend =  0.04, yend =  0.04), colour= "black",size=0.2)+  
		ggsave("Fig3B.pdf",width=11,height =4)
