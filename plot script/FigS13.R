library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(Hmisc)
library(corrplot)



args <- commandArgs(TRUE)
Pi_Dxy_Fst_Rho_file <- args[1]       
Pi_Dxy_Fst_Rho <- read.table(Pi_Dxy_Fst_Rho_file,sep="\t", header=TRUE)


for (pop in c("AB","BM","QH","NE","HN","TW")){
  assign(paste("Central_",pop,sep=""), 
         Pi_Dxy_Fst_Rho %>% filter(population == pop) %>% select(-heterogeneity)) 
  
  assign(paste("Central_",pop,"_fst_q99",sep=""),
         quantile(get(paste("Central_",pop,sep=""))$avg_wc_fst,prob=0.99)[[1]])
  
  assign(paste(pop,"_NewOutlier",sep=""),
         get(paste("Central_",pop,sep="")) %>% 
           mutate(heterogeneity = ifelse(avg_wc_fst>=get(paste("Central_",pop,"_fst_q99",sep="")), "outlier", "background"))) #heterogeneity
  
}

Pi_Dxy_Fst_Rho_NewOutlier <- dplyr::bind_rows(AB_NewOutlier,
                                              BM_NewOutlier,
                                              QH_NewOutlier,
                                              NE_NewOutlier,
                                              HN_NewOutlier,
                                              TW_NewOutlier)

dxy_stat <- Pi_Dxy_Fst_Rho_NewOutlier  %>%  
  ggplot(aes(population,Rho.pop,fill=heterogeneity))+
  geom_boxplot(outlier.shape = NA,width = 0.5)+   
  coord_cartesian(ylim=c(20, 500))+    
  theme_classic() +
  stat_compare_means(aes(group = heterogeneity),method ="wilcox.test",label = "p.signif",label.y = 0.025)+
  scale_fill_manual(values = c("grey50", "#BB2020"))+         
  ggsave("FigS13A.pdf",width=7,height = 4)
  

dxy_stat <- Pi_Dxy_Fst_Rho_NewOutlier  %>%  
  mutate(Pi_Rho=Rho.pop/avg_pi.pop) %>% 
  ggplot(aes(population,Pi_Rho,fill=heterogeneity))+
  geom_boxplot(outlier.shape = NA,width = 0.5)+   
  coord_cartesian(ylim=c(20, 50000))+    
  theme_classic() +
  stat_compare_means(aes(group = heterogeneity),method ="wilcox.test",label = "p.signif",label.y = 0.025)+
  scale_fill_manual(values = c("grey50", "#BB2020"))+          
  ggsave("FigS13B.pdf",width=7,height = 4)



