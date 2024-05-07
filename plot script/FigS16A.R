library(Hmisc)
library(gridExtra)
library(tidyr)
library(magrittr)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyverse)
library(ggpubr)
library(corrplot)

args <- commandArgs(TRUE)
pixy_df_file  <-  args[1] 
xpnsl_file    <-  args[2] 
pixy_df <- read.table(pixy_df_file, sep="\t", header=TRUE)
xpnsl_snp <- read.table(xpnsl_file, sep="\t", header=TRUE)



a <- aggregate(norm_xpnsl~population+chromosome+window_pos_1+window_pos_2,xpnsl_snp,mean)   
b <- aggregate(norm_xpnsl~population+chromosome+window_pos_1+window_pos_2,xpnsl_snp,median)  

combined <- a %>%
  	    inner_join(b,by=c("population","chromosome","window_pos_1","window_pos_2"),suffix = c('.mean','.median')) %>%
  	    dplyr::rename(mean_score = norm_xpnsl.mean, median_score = norm_xpnsl.median)  %>%
  	    arrange(population,chromosome, window_pos_1,window_pos_2) 

rm(a,b)
allpop_xpnsl_win_stat <- xpnsl_snp %>%
  		         select(-pos,-norm_xpnsl) %>%
			 distinct(.keep_all = TRUE) %>%  
			 inner_join(combined,by=c("population","chromosome","window_pos_1","window_pos_2")) %>%
			 mutate(window_pos_2=window_pos_2 - 1)  %>%
			 select(population,chromosome,window_pos_1,window_pos_2,sites_in_win,max_score,min_score,mean_score,median_score,xpnsl_heterogeneity)

pixy_xpnsl_stat <-  	allpop_xpnsl_win_stat %>%
  		    	inner_join(pixy_df,by=c("population","chromosome","window_pos_1","window_pos_2")) %>%
		    	mutate(π_ratio=avg_pi.central/avg_pi.pop)
write.table(pixy_xpnsl_stat,file="allpop_pixy_xpnsl_win_stat.tab",quote=F,row.names=F,sep="\t")

AllPop_Fst_xpnsl_corr  <- pixy_xpnsl_stat %>%
                          filter(mean_score > 0 ) %>%
  			  ggplot(aes(avg_wc_fst,mean_score)) + geom_point(size=0.1,colour="grey50") +  
			  facet_wrap(~ population, nrow = 2)+
			  geom_smooth(method = "lm",formula = y ~ x,se = FALSE,size=0.5,color="blue")+
			  stat_cor(method="spearman",size = 2)+
			  theme_classic() +
			  ylab("XP-nSL") +
			  xlab("Fst") +
			  theme(strip.background = element_rect(color = "white"))+ 
			  theme(text =  element_text(size = 13,family = "sans")) 
ggsave("FigS16A.pdf",width=7,height = 3)
