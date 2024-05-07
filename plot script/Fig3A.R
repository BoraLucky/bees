rm(list=ls())
library(ggplot2)
library(ggpubr)
library(dplyr)

args <- commandArgs(TRUE)
input_file <- args[1]
file_table <- read.table(input_file,header=T)
###########################################dxy###########################################
dxy_stat <- file_table  %>%  ggplot(aes(population,avg_dxy,fill=heterogeneity))+
  geom_boxplot(outlier.shape = NA,width = 0.5)+    
  coord_cartesian(ylim=c(0, 0.03))+    
  theme_classic() +
  stat_compare_means(aes(group = heterogeneity),method ="wilcox.test",label = "p.signif",label.y = 0.025)+
  scale_fill_manual(values = c("grey50", "#BB2020"))+               
  ggsave("Fig3A1.pdf",width=7,height = 4)

pop_pi_stat <- file_table  %>%  ggplot(aes(population,avg_pi.pop,fill=heterogeneity))+
  geom_boxplot(outlier.shape = NA,width = 0.5)+    
  coord_cartesian(ylim=c(0, 0.03))+
  theme_classic() +
  stat_compare_means(aes(group = heterogeneity), method ="wilcox.test",label = "p.signif",label.y = 0.025)+
  scale_fill_manual(values = c("grey50", "#BB2020"))+ 
  ggsave("Fig3A2.pdf",width=7,height = 4)


