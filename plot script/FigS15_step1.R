rm(list=ls())
library(ggplot2)
library(ggpubr)
library(dplyr)

args <- commandArgs(TRUE)
input_file <- args[1]
Dxy_SV_boxplot <- args[2]
Fst_SV_boxplot <- args[3]
Pop_Pi_boxplot <- args[4]


file_table <- read.table(input_file,header=T)
###########################################dxy###########################################
dxy_stat <- file_table  %>%  ggplot(aes(population,avg_dxy,fill=SV_heterogeneity))+
  geom_boxplot(outlier.shape = NA,width = 0.5)+   
  coord_cartesian(ylim=c(0, 0.03))+   
  theme_classic() +
  stat_compare_means(aes(group = SV_heterogeneity),method ="wilcox.test",label = "p.signif",label.y = 0.025)+
  scale_fill_manual(values = c("grey50", "#BB2020"))+ 
  ggsave(Dxy_SV_boxplot,width=7,height = 4)

##########################################Fst############################################
dxy_stat <- file_table  %>%  ggplot(aes(population,avg_wc_fst,fill=SV_heterogeneity))+
	geom_boxplot(outlier.shape = NA,width = 0.5)+   
        coord_cartesian(ylim=c(0, 1))+    
	theme_classic() +
	stat_compare_means(aes(group = SV_heterogeneity),method ="wilcox.test",label = "p.signif",label.y = 0.025)+
	scale_fill_manual(values = c("grey50", "#BB2020"))+ 
	ggsave(Fst_SV_boxplot,width=7,height = 4)

###########################################pi#########################################
pop_pi_stat <- file_table  %>%  ggplot(aes(population,avg_pi.pop,fill=SV_heterogeneity))+
  geom_boxplot(outlier.shape = NA,width = 0.5)+    
  coord_cartesian(ylim=c(0, 0.03))+
  theme_classic() +
  stat_compare_means(aes(group = SV_heterogeneity), method ="wilcox.test",label = "p.signif",label.y = 0.025)+
  scale_fill_manual(values = c("grey50", "#BB2020"))+ 
  ggsave(Pop_Pi_boxplot,width=7,height = 4)




