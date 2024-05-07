library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(Hmisc)
library(corrplot)
library(UpSetR)

args <- commandArgs(TRUE)
xpnsl_heterogeneity_table <- args[1] 
Window10K_Site10.xpnsl <- read.table(xpnsl_heterogeneity_table,header=T,stringsAsFactors=F)  
bin_boundary.txt <- args[2] 
qbin <- read.table(bin_boundary.txt,header=T,stringsAsFactors=F)
block_boundary <-  qbin$bin_boundary

##############################################outlier distribution#####################################
block_boundary_frame <- data.frame(X_line=block_boundary[2:10]) 

point_plot <- Window10K_Site10.xpnsl  %>%  
  select(Pop,sites_in_win,frac_scores_gt_threshold,heterogeneity) %>%
  ggplot(aes(x = sites_in_win, y = frac_scores_gt_threshold, color = heterogeneity ))+
  geom_point(size = 1.5,stroke = 0)+
  xlab("Number of sites in 10K window")+
  ylab("Proportion of scores >2  in windows")+
  scale_color_manual(values = c("grey50", "#BB2020"))+
  theme_classic()+
  theme(strip.text.x = element_text(color = "black", face = "italic",angle=0,size = 7),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position ="none")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1.1))+
  geom_vline(data=block_boundary_frame,aes(xintercept=X_line),colour="black", linetype="dashed",size = 0.4)+
  facet_grid(Pop ~ ., scales="free_y", space="free_y")+
ggsave(filename="FigS17.pdf",width = 6,height = 4)

