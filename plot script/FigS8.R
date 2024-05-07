library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(Hmisc)
library(corrplot)
library(UpSetR)


args <- commandArgs(TRUE)
input_path <- args[1]

file_table <- read.table(input_path,header=T)

all_tabls <- file_table %>% 
  select(population,chromosome,window_pos_1,window_pos_2,heterogeneity) %>%
  spread(population, value = heterogeneity) %>% 
  filter_all(any_vars(. == "outlier")) %>%   
  select(-chromosome,-window_pos_1,-window_pos_2)  %>%
  na.omit()

all_tabls[all_tabls=="background"]<- 0 
all_tabls[all_tabls=="outlier"]<- 1
all_tabls <- as.data.frame(lapply(all_tabls,as.numeric)) 


pdf(file='pops_overlap_outlier.pdf',onefile=FALSE,height = 6,width =12)
upset(all_tabls,nsets = 6, nintersects =80,mb.ratio=c(0.7, 0.3),order.by="freq",  
      show.numbers="yes",number.angles=-30,main.bar.color="blue",
      matrix.color="blue", point.size=3,line.size=1,shade.alpha=1,         
      sets.x.label="outlier number",set_size.show="TRUE",sets.bar.color="blue",set_size.numbers_size=7)

dev.off()


