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
library(patchwork) 

args <- commandArgs(TRUE)
pixy_df_file  <-  args[1] 
xpnsl_file    <-  args[2] 
pixy_df <- read.table(pixy_df_file, sep="\t", header=TRUE)
xpnsl_result <- read.table(xpnsl_file, sep="\t", header=TRUE)



names(xpnsl_result) <- c("population","chromosome","window_pos_1","window_pos_2","sites_in_win","max_score","min_score","heterogeneity")
Fst_sites_overlap  <-   xpnsl_result %>%
                        mutate(window_pos_2 = window_pos_2 - 1) %>%
                        inner_join(pixy_df,by=c("population","chromosome","window_pos_1","window_pos_2"),suffix = c(paste(".","xpnsl",sep=""),'.Fst')) %>%
                        select(-heterogeneity.Fst,-max_score,-min_score,-heterogeneity.xpnsl) %>%
                        mutate(π_ration = avg_pi.pop/avg_pi.central)
						
############################，Fst top5 #################################
for (pop in c("AB","BM","QH","NE","HN","TW")){
  assign(paste("Central_",pop,sep=""),
         Fst_sites_overlap %>% filter(population == pop))
  assign(paste("Central_",pop,"_fst_q95",sep=""),
         quantile(get(paste("Central_",pop,sep=""))$avg_wc_fst,prob=0.95)[[1]])
  assign(paste(pop,"_NewOutlier",sep=""),
         get(paste("Central_",pop,sep="")) %>%
         mutate(heterogeneity_top5 = ifelse(avg_wc_fst>=get(paste("Central_",pop,"_fst_q95",sep="")), "outlier", "background")))
}

Fst_Top5 <- dplyr::bind_rows(AB_NewOutlier,
                             BM_NewOutlier,
                             QH_NewOutlier,
                             NE_NewOutlier,
                             HN_NewOutlier,
                             TW_NewOutlier)
rm(list=objects(pattern='^Central.*'))  
Fst_xpnsl_conflicts <-  xpnsl_result  %>%
                        mutate(window_pos_2 = window_pos_2 - 1) %>%
                        inner_join(Fst_Top5,by=c("population","chromosome","window_pos_1","window_pos_2")) %>%
                        select(-sites_in_win.y) %>%
                        dplyr::rename(sites_in_win = sites_in_win.x,heterogeneity.xpnsl=heterogeneity,heterogeneity.Fst=heterogeneity_top5) %>%
                        mutate(Conflicts = ifelse(grepl("outlier", heterogeneity.Fst) & grepl("background", heterogeneity.xpnsl),"conflicts", "conformity"))  %>% 
                        select(-heterogeneity.xpnsl, -heterogeneity.Fst)

                         ################Random100#################
Fst_xpnsl_conflicts_Random100 <-    Fst_xpnsl_conflicts  %>%
                                    filter(Conflicts=="conflicts") %>%
                                    group_by(population) %>%  sample_n(size = 100) %>%  
                                    mutate(Conflicts="Random100_conflicts") %>%
                                    right_join(Fst_xpnsl_conflicts,
				     by=c("population","chromosome","window_pos_1","window_pos_2","sites_in_win",
                                    "max_score","min_score","avg_dxy","avg_wc_fst","avg_pi.pop","avg_pi.central","π_ration")) %>%
				     mutate(Conflicts = ifelse(grepl("Random100_conflicts", Conflicts.x) & grepl("conflicts", Conflicts.y),"Random100_conflicts", "conformity")) %>%
                                    select(-Conflicts.x,-Conflicts.y)  %>%
                                    arrange(population,chromosome, window_pos_1,window_pos_2)

conflicts_Random100_central_pi <- Fst_xpnsl_conflicts_Random100 %>%
                               ggplot(aes(population,avg_pi.central,fill=Conflicts))+
                               geom_boxplot(outlier.shape = NA,width = 0.5)+    
                               coord_cartesian(ylim=c(0, 0.03))+  
                               theme_classic() +
                               stat_compare_means(aes(group = Conflicts),method ="wilcox.test",label = "p.signif",label.y = 0.025)+
                               scale_fill_manual(values = c("grey50", "#BB2020"))

conflicts_Random100_density <- Fst_xpnsl_conflicts_Random100 %>%
　　　　　　　　　　　　　　　 ggplot(aes(population,sites_in_win,fill=Conflicts))+
　　　　　　　　　　　　　　　 geom_boxplot(outlier.shape = NA,width = 0.5)+    
　　　　　　　　　　　　　　　 coord_cartesian(ylim=c(10, 500))+  
　　　　　　　　　　　　　　　 theme_classic() +
　　　　　　　　　　　　　　　 stat_compare_means(aes(group = Conflicts),method ="wilcox.test",label = "p.signif",label.y = 0.025)+
　　　　　　　　　　　　　　　 scale_fill_manual(values = c("grey50", "#BB2020")) 
　

conflicts_Random100_combined <- conflicts_Random100_central_pi + conflicts_Random100_density
ggsave("Fig16B.pdf",width=14,height = 4)
