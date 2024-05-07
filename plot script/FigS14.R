rm(list=ls())
library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(Hmisc)
library(corrplot)

args <- commandArgs(TRUE)
background_outlier_file <- args[1]       
Pi_Dxy_Fst_Rho <- read.table(background_outlier,sep="\t", header=TRUE)
############################standardized nucleotide diversity#################################################
  
AllPop_dxy_pi  <- all_windows %>% 
                  na.omit()   %>%
                  dplyr::filter(grepl("outlier", heterogeneity)) %>%
                  mutate(pi_dxy=avg_pi.pop/avg_dxy) %>%
                  ggplot(aes(x=pi_dxy)) + 
                  geom_density(fill="#FF8C64", colour="grey50", size=.2)+
                  theme_classic()+xlab("pi/dxy")+ylab("density")+
                  scale_x_continuous(limits = c(0,2),breaks = c(0,1,2)) +
                  facet_grid(population ~ ., scales="free_y", space="free_y") 
ggsave(AllPop_dxy_pi,filename="Fig14.pdf",width = 4,height = 3)