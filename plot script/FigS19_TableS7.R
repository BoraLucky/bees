library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(Hmisc)
library(corrplot)
library(UpSetR)
library(patchwork) 
library(gg.gap)  

args <- commandArgs(TRUE)
all_pop.file <- args[1]


Mean_Fst <- aggregate(avg_wc_fst~population,file_table,mean) %>% dplyr::rename(Fst_mean = avg_wc_fst) 
Sd_Fst <-   aggregate(avg_wc_fst~population,file_table,sd) %>% dplyr::rename(Fst_sd = avg_wc_fst)  
Fst_mean_Sd <- Mean_Fst %>% inner_join(Sd_Fst,by="population") 


for (pop in c("AB","BM","QH","NE","HN","TW")){
  assign(paste("MeanFst_",pop,sep=""),(Fst_mean_Sd %>% filter(grepl(pop, population)) %>% select("Fst_mean"))$Fst_mean)
  assign(paste("SdFst_",pop,sep=""),(Fst_mean_Sd %>% filter(grepl(pop, population)) %>% select("Fst_sd"))$Fst_sd)
  assign(paste(pop,"_table",sep=""),
         file_table %>% dplyr::filter(grepl(pop, population)) %>%
           select(-heterogeneity) %>%
           mutate(Zfst = (avg_wc_fst - get(paste("MeanFst_",pop,sep="")))/get(paste("SdFst_",pop,sep=""))) %>% 
           mutate(Zfst_heterogeneity = ifelse(Zfst >= 2, "outlier", "background")))  
}
all_Zfst <- dplyr::bind_rows(AB_table,BM_table,QH_table,NE_table,HN_table,TW_table) 


dxy_stat <- all_Zfst  %>%
  ggplot(aes(population,avg_dxy,fill=Zfst_heterogeneity))+
  geom_boxplot(outlier.shape = NA,width = 0.5)+    
  coord_cartesian(ylim=c(0, 0.03))+    
  theme_classic() +
  stat_compare_means(aes(group = Zfst_heterogeneity),method ="wilcox.test",label = "p.signif",label.y = 0.025)+
  scale_fill_manual(values = c("grey50", "#BB2020")) 
ggsave("Zfst_outlier_Dxy_boxplot.pdf",width=7,height = 4)


#####################################Zfst outlier##################################
Mean_Dxy <- aggregate(avg_dxy~population,all_Zfst,mean) %>% dplyr::rename(dxy_mean = avg_dxy) 
Sd_Dxy <- aggregate(avg_dxy~population,all_Zfst,sd) %>% dplyr::rename(dxy_sd = avg_dxy) 
Dxy_mean_Sd <- Mean_Dxy %>%
  inner_join(Sd_Dxy,by="population") %>% 
  mutate(below_2SD= dxy_mean-2*dxy_sd) %>% 
  mutate(above_2SD= dxy_mean+2*dxy_sd)

for (pop in c("AB","BM","QH","NE","HN","TW")){
  assign(paste("below_2SD_",pop,sep=""),(Dxy_mean_Sd %>% filter(grepl(pop, population)) %>% select("below_2SD"))$below_2SD)
  assign(paste("above_2SD_",pop,sep=""),(Dxy_mean_Sd %>% filter(grepl(pop, population)) %>% select("above_2SD"))$above_2SD)
  assign(paste(pop,"_clssify",sep=""), all_Zfst %>%
           filter(grepl(pop, population)) %>%
           mutate(sites = ifelse(grepl("background", Zfst_heterogeneity), "background",
                                 ifelse(avg_dxy >= get(paste("above_2SD_",pop,sep="")), "above_2sd",
                                        ifelse(avg_dxy <= get(paste("below_2SD_",pop,sep="")),"below_2sd","Middle_2sd")))))
}
all_pop_clssify <- dplyr::bind_rows(AB_clssify,BM_clssify,QH_clssify,NE_clssify,HN_clssify,TW_clssify)  

################################################Table S7_2SD##########################################################################
outlier_count <- all_Zfst %>% filter(Zfst_heterogeneity=="outlier") %>% group_by(population) %>% count()
Dxy_below_2sd_count <- all_pop_clssify %>% filter(sites=="below_2sd") %>% group_by(population) %>% count()
Dxy_up_2sd_count <- all_pop_clssify %>% filter(sites=="above_2sd") %>% group_by(population) %>% count()
Dxy_count_table <- outlier_count %>%
  inner_join(Dxy_below_2sd_count,by="population") %>%
  inner_join(Dxy_up_2sd_count,by="population") %>%
  rename(all_outlier_n=n.x,below_2sd_n=n.y,above_2sd_n=n) %>% 
  mutate(up_2sd_frac=above_2sd_n/all_outlier_n) %>%
  mutate(below_2sd_frac=below_2sd_n/all_outlier_n)
write.table(Dxy_count_table,file="Table_S7_Dxy_2SD.tab",quote=F,row.names=F,sep="\t")¡¡¡¡

#################################################Fig19######################################################################
plot_classify_his0 <- all_pop_clssify %>%
  ggplot(aes(avg_wc_fst,fill = Zfst_heterogeneity))+
  geom_histogram(bins = 50)+
  facet_wrap(~ population, nrow = 1)+
  theme_classic() +
  scale_fill_manual(values = c("grey50", "#BB2020")) + 
  ggsave("Fst_background_outlier.pdf",width=8,height = 2)
plot_classify_his1 <- all_pop_clssify %>%
  ggplot(aes(avg_dxy,fill = Zfst_heterogeneity))+
  geom_histogram(bins = 50)+
  facet_wrap(~ population, nrow = 1)+
  theme_classic() +
  scale_fill_manual(values = c("grey50", "#BB2020")) + 
  theme(strip.background = element_blank()) + 
  ggsave("Fst_background_outlier_dxy_2SD_2classify.pdf",width=8,height = 2)

fig <- plot_classify_his0 + plot_classify_his1 + plot_layout(ncol=1) +
  ggsave("Fig19AB.pdf",width=12,height = 5)  

#########################################################
plot_classify_his2 <- all_pop_clssify %>%
  filter(grepl("outlier", Zfst_heterogeneity)) %>%
  ggplot(aes(avg_dxy, fill = sites))+
  geom_histogram(bins = 50)+
  facet_wrap(~ population, nrow = 1)+
  theme_classic() +
  scale_fill_manual(values = c("#BB2020","blue","grey50")) +
  theme(strip.background = element_blank(),strip.text = element_text(size = 0)) + 
  ggsave("Fst_outlier_dxy_2SD_3classify.pdf",width=8,height = 2)

picture <- gg.gap(plot = plot_classify_his2,
                  segments = c(10, 30), 
                  tick_width = 10,        
                  rel_heights = c(0.4, 0, 0.4),  
                  ylim = c(0,180)) +        
  ggsave("Fig19C_truncation.pdf",width=8,height = 3) 
