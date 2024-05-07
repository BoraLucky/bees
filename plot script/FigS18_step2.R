
rm(list=ls())
library(ggplot2)
library(ggpubr)
library(dplyr)

args <- commandArgs(TRUE)
input_file <- args[1]
output_file <- args[2]
file_table <- read.table(input_file,header=T)

###########################################Rho##########################################
dxy_stat <- file_table  %>%  ggplot(aes(population,Rho,fill=selscan_heterogeneity))+
  geom_boxplot(outlier.shape = NA,width = 0.5)+    
  coord_cartesian(ylim=c(0, 500))+   
  theme_classic() +
  stat_compare_means(aes(group = selscan_heterogeneity),method ="wilcox.test",label = "p.signif",label.y = 0.025)+
  scale_fill_manual(values = c("grey50", "#BB2020"))+
  ggsave(output_file,width=7,height = 4)

