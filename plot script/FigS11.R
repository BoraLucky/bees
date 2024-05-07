library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(Hmisc)
library(corrplot)


args <- commandArgs(TRUE)
SC2N2mG_pramer_file <- args[1]       
SC2N2mG_pramer <- read.table(SC2N2mG_pramer_file,sep="\t", header=T)

cb_palette <- c("#A65628","#FF7F00","#377EB8","#F781BF","#E41A1C","#00BFFF")

SC2N2mG_prame_long  <-     SC2N2mG_pramer  %>% 
			   mutate(All_M12 = m12+me12)  %>%
			   mutate(All_M21 = m21+me21)  %>%
		       gather(-Pop,key = "statistic", value = "value") #long form

island_no_island_geneflow <- SC2N2mG_prame_long %>%
                      dplyr::filter(grepl("All_M12", statistic) | grepl("All_M21", statistic) | 
                      grepl("m12", statistic) | grepl("m21", statistic) |
                      grepl("me12", statistic) | grepl("me21", statistic)) %>%  
                      ggplot(aes(Pop,value,fill=statistic))+
                      geom_boxplot(outlier.shape = NA,width = 0.5)+    
                      theme_classic() +
                      scale_fill_brewer(palette = "Set2")+
                      ggsave("island_no_island_geneflow.pop.pdf",width=7,height = 3) 

