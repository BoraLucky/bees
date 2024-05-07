
library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(Hmisc)
library(corrplot)
args <- commandArgs(TRUE)
all_pop.file <- args[1] 

file_table <- read.table(all_pop.file,header = T)

##########################################Fst#######################################################
Fst_corr_params <- file_table %>% 
		   select(population,chromosome,window_pos_1,window_pos_2,avg_wc_fst) %>%
                   spread(population, value = avg_wc_fst) %>% 
                   select(-chromosome,-window_pos_1,-window_pos_2) %>% 
                   na.omit()  %>% 
                   as.matrix()%>% 
                   rcorr(type = 'spearman') 
      
pdf("Fst_window_correlated.pdf",width = 4,height = 4)
corrplot(Fst_corr_params$r, number.cex = 0.8, diag = FALSE, 
           tl.cex = 0.8,col = c("darkorange", "steelblue"),bg = "white",addCoef.col = "black")
dev.off
##########################################Dxy#######################################################
Dxy_corr_params <- file_table %>% 
		   select(population,chromosome,window_pos_1,window_pos_2,avg_dxy) %>%
                   spread(population, value = avg_dxy) %>% 
                   select(-chromosome,-window_pos_1,-window_pos_2) %>%
		   filter(heterogeneity=="outlier") %>%    
                   na.omit()  %>% 
                   as.matrix()%>%
                   rcorr(type = 'spearman') 
pdf("Dxy_window_correlated.pdf",width = 4,height = 4)
corrplot(Dxy_corr_params$r, number.cex = 0.8, diag = FALSE,
         tl.cex = 0.8,col = c("darkorange", "steelblue"),bg = "white",addCoef.col = "black")
dev.off

########################################mean_Pi###########################################################
Mean_Pi_corr_params <- file_table %>% 
		  mutate(mean_pi=(avg_pi.pop+avg_pi.central)/2) %>%
		  select(population,chromosome,window_pos_1,window_pos_2,mean_pi) %>%
		  spread(population, value = mean_pi) %>%
	          select(-chromosome,-window_pos_1,-window_pos_2) %>%
		  filter(heterogeneity=="outlier") %>%    
                  na.omit()  %>%
	          as.matrix()%>% 
		  rcorr(type = 'spearman')
pdf("Mean_pi_window_correlated.pdf",width = 4,height = 4) 
corrplot(Mean_Pi_corr_params$r, number.cex = 0.8, diag = FALSE,
	tl.cex = 0.8,col = c("darkorange", "steelblue"),bg = "white",addCoef.col = "black")
dev.off

