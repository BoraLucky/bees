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
library(qpdf)


args <- commandArgs(TRUE)
avg_pi_tab <- args[1]  
pixy_tab <- args[2] 
FastEPRR_tab <- args[3]   
background_outlier_tab <- args[4] 

pop_avg_pi <- read.table(avg_pi_tab, sep="\t", header=TRUE)
Allpop_pair_pixy <- read.table(pixy_tab, sep="\t", header=TRUE)
FastEPRR_result <- read.table(FastEPRR_tab, sep="\t", header=TRUE)
pixy_df <- read.table(background_outlier_tab, sep="\t", header=TRUE)
names(Allpop_pair_pixy) <- c("pop1","pop2","chromosome","window_pos_1","window_pos_2","summary_statistic","summary_value") 

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(pop_pair1 = rownames(cormat)[row(cormat)[ut]],
             pop_pair2 = rownames(cormat)[col(cormat)[ut]],
             cor =(cormat)[ut],
             p = pmat[ut])
}   
 

             ##############################[corr(ρi,ρj), corr(πi,πj)]###########
pdf("a)between_populations_Rho_π_correlated.pdf",width = 6,height = 3)  
opar <- par(no.readonly = TRUE)  
par(mfrow=c(1,2))
 
Rho_corr_params <- FastEPRR_result  %>%
　　　　　　　　　　select(population,chromosome,window_pos_1,window_pos_2,Rho) %>%
　　　　　　　　　　spread(population, value = Rho) %>%
　　　　　　　　　　select(-chromosome,-window_pos_1,-window_pos_2,-ML) %>%
　　　　　　　　　　na.omit()  %>% 
　　　　　　　　　　as.matrix()%>% 
　　　　　　　　　  rcorr(type = 'spearman') 
corrplot(Rho_corr_params$r, number.cex = 0.8, diag = FALSE, title="Rho", mar=c(0, 0, 1, 0),
         tl.cex = 0.8,tl.col = "black",
         col = c("darkorange", "steelblue"),bg = "white",addCoef.col = "black")
 
π_corr_params <- pop_avg_pi  %>%
　　　　　　　　　　select(population,chromosome,window_pos_1,window_pos_2,avg.pi) %>%
　　　　　　　　　　spread(population, value = avg.pi) %>%
　　　　　　　　　　select(-chromosome,-window_pos_1,-window_pos_2,-ML) %>%
　　　　　　　　　　na.omit()  %>% 
　　　　　　　　　　as.matrix() %>% 
　　　　　　　　　　rcorr(type = 'spearman') 
corrplot(π_corr_params$r, number.cex = 0.8, diag = FALSE, title="pi", mar=c(0, 0, 1, 0),
         tl.cex = 0.8,tl.col = "black",
         col = c("darkorange", "steelblue"),bg = "white",addCoef.col = "black")
par(opar) 
dev.off()
Rho_corr <-  flattenCorrMatrix(Rho_corr_params$r, Rho_corr_params$P) %>% 
　　　　　　 dplyr::rename(pop1=pop_pair1,pop2=pop_pair2)
π_corr <- flattenCorrMatrix(π_corr_params$r, π_corr_params$P) %>% 
　　　　     dplyr::rename(pop1=pop_pair1,pop2=pop_pair2)
between_populations_params <- inner_join(Rho_corr,π_corr,by=c('pop1','pop2'),suffix = c(paste(".","Rho",sep=""),'.π')) %>% 
　　　　　　　　　　　　　　　	  arrange(pop1, pop2) %>% 　　　　　　　　　　　　　　　
　　                           mutate(pop_pair = paste(get("pop1"),get("pop2"),sep="_")) %>% 　　　　　　　　　　　　　　　　　　
                              select(pop_pair,cor.Rho,cor.π,p.Rho,p.π) %>% 　　　　　　　　　　　　　　　　　
　                 	          mutate(pop_pair = str_replace_all(pop_pair,c("Central_HN"="HN_Central","Central_NE"="NE_Central", 　　　　　　　　　　　　　　　　　　　　　　　　　          "Central_QH"="QH_Central","Central_TW"="TW_Central"))) 
rm(Rho_corr,π_corr,Rho_corr_params,π_corr_params) 

                #############################[Figure S12: corr(ρi,πi) ]##################################
within_Pop_π_Rho_corr <- inner_join(pop_avg_pi,FastEPRR_result,by=c('population','chromosome','window_pos_1','window_pos_2')) %>%
　　　　　　　　　　　　 filter(population!="ML") %>%
　　　　　　　　　　　　 select(-chromosome,-window_pos_1,-window_pos_2,-CIL,-CIR) %>% # 
　　　　　　　　　　　　 na.omit() %>%
　　　　　　　　　　　　 group_by(population) %>%
　　　　　　　　　　　　 dplyr::summarize("Rho~π"=cor(avg.pi,Rho,method = "spearman")) %>%
　　　　　　　　　　　　 arrange(population)   
π_Rho <-  inner_join(pop_avg_pi,FastEPRR_result,
                     by=c('population','chromosome','window_pos_1','window_pos_2')) %>%
  　　　　filter(population!="ML" & population!="Central")   %>%
  　　　  select(-CIL,-CIR)  
AllPop_pi_Rho_pic  <- π_Rho %>%
 　　　　　　　　　　 ggplot(aes(Rho,avg.pi)) + geom_point(size=0.1,colour="grey50") +  # Create basic ggplot
　　　　　　　　　　　facet_wrap(~ population, nrow = 2)+
　　　　　　　　　　　geom_smooth(method = "lm",formula = y ~ x,se = FALSE,size=0.5,color="blue")+
　　　　　　　　　　　stat_cor(method="spearman",size = 2)+
　　　　　　　　　　　theme_classic() +
　　　　　　　　　　　ylab("pi") +
　　　　　　　　　　　theme(strip.background = element_rect(color = "white"))+  
　　　　　　　　　　　theme(text =  element_text(size = 13,family = "sans"))  
ggsave("b)within_population_Rho_π_correlations.pdf",width=7,height = 3)
rm(π_Rho,AllPop_pi_Rho_pic)



            #################### [e.g. corr(πi/j,FstI=i,j),  (πi/j,DxyI=i,j)]######
π <- inner_join(pop_avg_pi,FastEPRR_result,by=c('population','chromosome','window_pos_1','window_pos_2')) %>%
     filter(population!="ML") %>% select(-Rho,-CIL,-CIR)  
                      ##############################Fst~π##############################
Fst_π <-  Allpop_pair_pixy %>%
          filter(summary_statistic == "avg_wc_fst")  %>%
          mutate(pop_pair=paste(get("pop1"),get("pop2"),sep="_"))  %>%
          select(pop_pair,chromosome,window_pos_1,window_pos_2,summary_value) %>%
          dplyr::rename(Fst=summary_value)%>%
          inner_join(π,by=c('chromosome','window_pos_1','window_pos_2'))  %>%  
          filter(population != "ML") %>%
          filter(
                  (grepl("AB", pop_pair) & grepl("AB", population))  |
                  (grepl("BM", pop_pair) & grepl("BM", population))  |
                  (grepl("QH", pop_pair) & grepl("QH", population))  |
                  (grepl("NE", pop_pair) & grepl("NE", population))  |
                  (grepl("HN", pop_pair) & grepl("HN", population))  |
                  (grepl("TW", pop_pair) & grepl("TW", population))  |
                  (grepl("Central", pop_pair) & grepl("Central", population))
                )     
for (target_pop in c("AB","BM","QH","HN","TW","NE","Central")){
    assign(paste(target_pop,"Fst_π",sep="_"), Fst_π %>%
          filter(population==target_pop) %>%
          na.omit()  %>%
          group_by(pop_pair) %>%
          dplyr::summarize(single_pop_π=cor(Fst,avg.pi,method = "spearman")) %>%
          `names<-`(c("pop_pair",paste(target_pop,"π~Fst",sep="."))))  
}  
 
Fst_samePair_π <- full_join(AB_Fst_π,BM_Fst_π,by="pop_pair") %>%
                  full_join(QH_Fst_π,by="pop_pair") %>%
                  full_join(NE_Fst_π,by="pop_pair") %>%
                  full_join(HN_Fst_π,by="pop_pair") %>%
                  full_join(TW_Fst_π,by="pop_pair") %>%
                  full_join(Central_Fst_π,by="pop_pair") %>%
                  arrange(pop_pair)
rm(AB_Fst_π,BM_Fst_π,QH_Fst_π,NE_Fst_π,HN_Fst_π,TW_Fst_π,Central_Fst_π)
 
                  #################################Dxy~π#################################
Dxy_π <-  Allpop_pair_pixy %>%
          filter(summary_statistic == "avg_dxy")  %>%
          mutate(pop_pair=paste(get("pop1"),get("pop2"),sep = "_"))  %>%
          select(pop_pair,chromosome,window_pos_1,window_pos_2,summary_value) %>%
          dplyr::rename(Dxy=summary_value)%>%
          inner_join(π,by=c('chromosome','window_pos_1','window_pos_2'))  %>%  
          filter(population != "ML") %>%
          filter(
                  (grepl("AB", pop_pair) & grepl("AB", population))  |
                  (grepl("BM", pop_pair) & grepl("BM", population))  |
                  (grepl("QH", pop_pair) & grepl("QH", population))  |
                  (grepl("NE", pop_pair) & grepl("NE", population))  |
                  (grepl("HN", pop_pair) & grepl("HN", population))  |
                  (grepl("TW", pop_pair) & grepl("TW", population))  |
                  (grepl("Central", pop_pair) & grepl("Central", population))
                )     
 
for (target_pop in c("AB","BM","QH","HN","TW","NE","Central")){
  assign(paste(target_pop,"Dxy_π",sep="_"), Dxy_π %>%
        filter(population==target_pop) %>%
        na.omit()  %>%
        group_by(pop_pair) %>%
        dplyr::summarize(single_pop_π=cor(Dxy,avg.pi,method = "spearman")) %>%
        `names<-`(c("pop_pair",paste(target_pop,"π~Dxy",sep=".")))) 
}  
 
Dxy_samePair_π <- full_join(AB_Dxy_π,BM_Dxy_π,by="pop_pair") %>%
                  full_join(QH_Dxy_π,by="pop_pair") %>%
                  full_join(HN_Dxy_π,by="pop_pair") %>%
                  full_join(TW_Dxy_π,by="pop_pair") %>%
                  full_join(NE_Dxy_π,by="pop_pair") %>%
                  full_join(Central_Dxy_π,by="pop_pair") %>%
                  arrange(pop_pair)
combined_SamPair_Fst_π_Dxy_corr <- full_join(Fst_samePair_π,Dxy_samePair_π,by="pop_pair") 
rm(AB_Dxy_π,BM_Dxy_π,QH_Dxy_π,NE_Dxy_π,HN_Dxy_π,TW_Dxy_π,Central_Dxy_π,Dxy_samePair_π,Fst_samePair_π,π)

##################Figure 2B:  [e.g. corr(πi/j,FstI=i,j), corr(ρi/j,FstI=i,j), (πi/j,DxyI=i,j), corr(ρi/j,DxyI=i,j))]##################
combineed_peripatric_Rho <- inner_join(pixy_df,FastEPRR_result,
                                       by=c('population','chromosome','window_pos_1','window_pos_2'))  %>% select(-CIL,-CIR)  
 
Central_Rho <- FastEPRR_result  %>% filter(population  %in% "Central") %>% select(-CIL,-CIR) 
 
all_pixy_Rho <- inner_join(combineed_peripatric_Rho,Central_Rho,
                           by=c('chromosome','window_pos_1','window_pos_2'),
                           suffix = c(paste(".","pop",sep=""),'.central'))  %>% 
 　　　　　　　　 select(-population.central)  %>%
 　　　　　　　　 select(population.pop,chromosome,window_pos_1,window_pos_2,
        　　　　　　　　 avg_dxy,avg_wc_fst,avg_pi.pop,avg_pi.central,
         　　　　　　　　Rho.pop,Rho.central,heterogeneity) %>% 
  　　　　　　　　dplyr::rename(population=population.pop)  
 
pdf("c)between_populations_Fst_dxy_π_ρ_correlated.pdf",width = 7,height = 4) 
opar <- par(no.readonly = TRUE)  
par(mfrow=c(2,3)) 
for (target_pop in c("AB","BM","QH","HN","TW","NE")){
  assign(paste(target_pop,"summary_corr_params",sep="_"), all_pixy_Rho %>%
           na.omit()  %>%
           select(population,avg_dxy,avg_wc_fst,avg_pi.pop,avg_pi.central,Rho.pop,Rho.central) %>%
           filter(population==target_pop) %>%
           select(-population)%>%
           as.matrix()%>% 
           rcorr(type = 'spearman')) 
  corrplot(get(paste(target_pop,"summary_corr_params",sep="_"))$r, number.cex = 0.8, diag = FALSE,
           title=paste("CT",target_pop,sep="_"), mar=c(0, 0, 1, 0),
           tl.cex = 0.8,tl.col = "black",
           col = c("darkorange", "steelblue"),bg = "white",addCoef.col = "black")
}  
par(opar)  
dev.off() 
rm(AB_summary_corr_params,BM_summary_corr_params,QH_summary_corr_params,NE_summary_corr_params,HN_summary_corr_params,TW_summary_corr_params,combineed_peripatric_Rho,Central_Rho,all_pixy_Rho)
 
               #################################FigureS9#############
pdf("d)population_pairs_Fst_dxy_landscapes_correlated.pdf",width = 6,height = 3) 
opar <- par(no.readonly = TRUE)  
par(mfrow=c(1,2)) 
Fst_corr_params <- pixy_df %>%
  　　　　select(population,chromosome,window_pos_1,window_pos_2,avg_wc_fst) %>%
 　　　　 spread(population, value = avg_wc_fst) %>%
  　　　　select(-chromosome,-window_pos_1,-window_pos_2) %>%
  　　　　na.omit()  %>%
  　　　　as.matrix()%>% 
  　　　　rcorr(type = 'spearman') 
corrplot(Fst_corr_params$r, number.cex = 0.8, diag = FALSE,
         title="CT-Peripatric pairs Fst",mar=c(0, 0, 1, 0),
         tl.cex = 0.8, tl.col = "black",
         col = c("darkorange", "steelblue"),bg = "white",addCoef.col = "black")
Dxy_corr_params <- pixy_df %>%
  　　　　select(population,chromosome,window_pos_1,window_pos_2,avg_dxy) %>%
  　　　　spread(population, value = avg_dxy) %>%
  　　　　select(-chromosome,-window_pos_1,-window_pos_2) %>%
  　　　　na.omit()  %>%
  　　　　as.matrix()%>% 
 　　　　 rcorr(type = 'spearman') 
corrplot(Dxy_corr_params$r, number.cex = 0.8, diag = FALSE,
         tl.cex = 0.8,tl.col = "black",
         title="CT-Peripatric pairs Dxy",mar=c(0, 0, 1, 0),
         col = c("darkorange", "steelblue"),bg = "white",addCoef.col = "black")
par(opar) 
dev.off()
Fst_corr <-　　flattenCorrMatrix(Fst_corr_params$r, Fst_corr_params$P) %>%
　　　　　　　　　select(-p) %>%
　　　　　　　　　mutate(pop_pairs = paste(get("pop_pair1"),"CT",get("pop_pair2"),sep="_")) %>%
　　　　　　　　　select(pop_pairs,cor) %>% 
　　　　　　　　　arrange(pop_pairs) 
Dxy_corr <-　　flattenCorrMatrix(Dxy_corr_params$r, Dxy_corr_params$P) %>% 
　　　　　　　　　select(-p) %>%
　　　　　　　　　mutate(pop_pairs = paste(get("pop_pair1"),"CT",get("pop_pair2"),sep="_")) %>%
　　　　　　　　　select(pop_pairs,cor) %>%
　　　　　　　　　arrange(pop_pairs)
between_CT_peripatric_Fst_Dxy_params <- inner_join(Fst_corr,Dxy_corr,
　　　　　　　　　　　　　　　　　　　　　　　　by="pop_pairs",suffix = c(paste(".","Fst",sep=""),'.Dxy')) %>%
　　　　　　　　　　　　　　　　　　　　　　　　arrange(pop_pairs) 
rm(Fst_corr_params,Dxy_corr_params,Fst_corr,Dxy_corr)

                 #####################Fst~Dxy########### 
ALLpair_Fst_dxy_corr_params <-　Allpop_pair_pixy %>%
                                mutate(pop_pair=paste(get("pop1"),get("pop2"),sep="_"))  %>%
                                select(pop_pair,chromosome,window_pos_1,window_pos_2,summary_statistic,summary_value) %>%
                                inner_join(FastEPRR_result,by=c('chromosome','window_pos_1','window_pos_2'))  %>%  
                                filter(population != "ML") %>%
                                select(-CIL,-CIR,-Rho,-population) %>%
                                distinct(.keep_all = TRUE) %>% 
                                spread(summary_statistic, summary_value) %>% 
                                select(pop_pair,avg_wc_fst,avg_dxy) %>%
                                na.omit() %>% 
                                group_by(pop_pair) %>%
                                dplyr::summarize(COR=cor(avg_wc_fst,avg_dxy,method = "spearman")) %>%
                                dplyr::rename(cor.Dxy_Fst=COR) %>%
                                arrange(pop_pair)  
combined_between_Pop_params <- inner_join(between_populations_params,ALLpair_Fst_dxy_corr_params,by='pop_pair') %>% select(-p.Rho,-p.π) 
rm(ALLpair_Fst_dxy_corr_params,between_populations_params)

               
ALLpair_Fst_corr_params  <-  Allpop_pair_pixy %>%
                            filter(summary_statistic!="avg_dxy")  %>%
                            mutate(pop_pair = paste(get("pop1"),get("pop2"),sep="_")) %>%
                            select(pop_pair,chromosome,window_pos_1,window_pos_2,summary_statistic,summary_value) %>%
                            inner_join(FastEPRR_result,by=c('chromosome','window_pos_1','window_pos_2'))  %>%  
                            filter(population != "ML") %>%
                            select(-CIL,-CIR,-Rho,-population) %>%
                            distinct(.keep_all = TRUE) %>%  
                            spread(pop_pair, summary_value) %>%
                            select(-chromosome,-window_pos_1,-window_pos_2,-summary_statistic) %>%
                            na.omit()  %>% 
                            as.matrix() %>% 
                            rcorr(type = 'spearman')
							
ALLpair_Dxy_corr_params <-  Allpop_pair_pixy %>%
                            filter(summary_statistic!="avg_wc_fst")  %>%
                            mutate(pop_pair = paste(get("pop1"),get("pop2"),sep="_")) %>%
                            select(pop_pair,chromosome,window_pos_1,window_pos_2,summary_statistic,summary_value) %>%
                            inner_join(FastEPRR_result,by=c('chromosome','window_pos_1','window_pos_2'))  %>%  
                            filter(population != "ML") %>%
                            select(-CIL,-CIR,-Rho,-population) %>%
                            distinct(.keep_all = TRUE) %>% 
                            spread(pop_pair, summary_value) %>%
                            select(-chromosome,-window_pos_1,-window_pos_2,-summary_statistic) %>%
                            na.omit()  %>% 
                            as.matrix() %>% 
                            rcorr(type = 'spearman')  

     ####################################FstI~FstJ,DxyI~DxyJ######################################
Fst_no_pseudo_replicated <- flattenCorrMatrix(ALLpair_Fst_corr_params$r,ALLpair_Fst_corr_params$P) %>%
　　　　　　　　　　　　　dplyr::filter(
    　　　　　　　　　　　　　　　　　　　!(
      　　　　　　　　　　　　　　　　　　(grepl("AB", pop_pair1) & grepl("AB", pop_pair2))  |
      　　　　　　　　　　　　　　　　　　(grepl("BM", pop_pair1) & grepl("BM", pop_pair2))  |
      　　　　　　　　　　　　　　　　　　(grepl("QH", pop_pair1) & grepl("QH", pop_pair2))  |
      　　　　　　　　　　　　　　　　　　(grepl("NE", pop_pair1) & grepl("NE", pop_pair2))  |
      　　　　　　　　　　　　　　　　　　(grepl("HN", pop_pair1) & grepl("HN", pop_pair2))  |
      　　　　　　　　　　　　　　　　　　(grepl("TW", pop_pair1) & grepl("TW", pop_pair2))  |
      　　　　　　　　　　　　　　　　　　(grepl("Central", pop_pair1) & grepl("Central", pop_pair2))
   　　　　　　　　　　　　　　　　　　　  )     # no pseudo_replicated_pair  
  　　　　　　　　　　　　　　　　　　　　 )  %>%
  　　　　　　　　　　　 select(-p) %>%
  　　　　　　　　　　　 arrange(pop_pair1, pop_pair2) 
Dxy_no_pseudo_replicated <- flattenCorrMatrix(ALLpair_Dxy_corr_params$r, ALLpair_Dxy_corr_params$P) %>%
　　　　　　　　　　　　　　　　　dplyr::filter(
    　　　　　　　　　　　　　　　　　　　　　　　!(
        　　　　　　　　　　　　　　　　　　　　　(grepl("AB", pop_pair1) & grepl("AB", pop_pair2))  |
        　　　　　　　　　　　　　　　　　　　　　(grepl("BM", pop_pair1) & grepl("BM", pop_pair2))  |
        　　　　　　　　　　　　　　　　　　　　　(grepl("QH", pop_pair1) & grepl("QH", pop_pair2))  |
        　　　　　　　　　　　　　　　　　　　　　(grepl("NE", pop_pair1) & grepl("NE", pop_pair2))  |
        　　　　　　　　　　　　　　　　　　　　　(grepl("HN", pop_pair1) & grepl("HN", pop_pair2))  |
       　　　　　　　　　　　　　　　　　　　　　 (grepl("TW", pop_pair1) & grepl("TW", pop_pair2))  |
        　　　　　　　　　　　　　　　　　　　　　(grepl("Central", pop_pair1) & grepl("Central", pop_pair2))
    　　　　　　　　　　　　　　　　　　　　　　　　)     # no pseudo_replicated_pair  
  　　　　　　　　　　　　　　　　　　　　　　　　)  %>%
  　　　　　　　　　　　　　 　select(-p) %>%
  　　　　　　　　　　　　　 　arrange(pop_pair1, pop_pair2) 
exclude_pseudo_replicated_params <- inner_join(Fst_no_pseudo_replicated,Dxy_no_pseudo_replicated,
                                               by=c('pop_pair1','pop_pair2'),suffix = c(paste(".","Fst",sep=""),'.Dxy')) %>%
  　　　　　　　　　　　　　　　　　　　　　　　　　　　arrange(pop_pair1, pop_pair2) %>% 
　　　　　　　　　　　　　　　　　　　　　　　　　　　　dplyr::rename("FstI~FstJ"=cor.Fst,"DxyI~DxyJ"=cor.Dxy) 
rm(ALLpair_Fst_corr_params,ALLpair_Dxy_corr_params,Fst_no_pseudo_replicated,Dxy_no_pseudo_replicated)        

           ######################################Rhoi~FstJ############################################
Fst_Rho <- Allpop_pair_pixy %>%
           filter(summary_statistic == "avg_wc_fst")  %>%
           mutate(pop_pair=paste(get("pop1"),get("pop2"),sep="_"))  %>%
           select(pop_pair,chromosome,window_pos_1,window_pos_2,summary_value) %>%
           dplyr::rename(Fst=summary_value)%>%
           inner_join(FastEPRR_result,by=c('chromosome','window_pos_1','window_pos_2'))  %>%  
           select(-CIL,-CIR) %>%
           filter(population != "ML") %>%
           filter(
                !(
                   (grepl("AB", pop_pair) & grepl("AB", population))  |
                   (grepl("BM", pop_pair) & grepl("BM", population))  |
                   (grepl("QH", pop_pair) & grepl("QH", population))  |
                   (grepl("NE", pop_pair) & grepl("NE", population))  |
                   (grepl("HN", pop_pair) & grepl("HN", population))  |
                   (grepl("TW", pop_pair) & grepl("TW", population))  |
                   (grepl("Central", pop_pair) & grepl("Central", population))
                 ))     
for (target_pop in c("AB","BM","QH","HN","TW","NE","Central")){
　　　assign(paste(target_pop,"Rho_Fst",sep="_"), Fst_Rho %>%
           filter(population==target_pop) %>%
           na.omit()  %>%
           group_by(pop_pair) %>%
           dplyr::summarize(single_pop_Rho=cor(Fst,Rho,method = "spearman")) %>%
           `names<-`(c("pop_pair",paste(target_pop,"Rho",sep="."))))  
}  
Rho_anotherPair_Fst <- full_join(AB_Rho_Fst,BM_Rho_Fst,by="pop_pair") %>%
                       full_join(QH_Rho_Fst,by="pop_pair") %>%
                       full_join(NE_Rho_Fst,by="pop_pair") %>%
                       full_join(HN_Rho_Fst,by="pop_pair") %>%
                       full_join(TW_Rho_Fst,by="pop_pair") %>%
　　　　　　　　　　　　　　full_join(Central_Rho_Fst,by="pop_pair") %>%
                       arrange(pop_pair)
rm(Fst_Rho,QH_Rho_Fst,AB_Rho_Fst,BM_Rho_Fst,HN_Rho_Fst,TW_Rho_Fst,NE_Rho_Fst,Central_Rho_Fst)
         
               ##########################################Rhoi~DxyJ################################         
Dxy_Rho <- Allpop_pair_pixy %>%
           filter(summary_statistic == "avg_dxy")  %>%
           mutate(pop_pair=paste(get("pop1"),get("pop2"),sep = "_"))  %>%
           select(pop_pair,chromosome,window_pos_1,window_pos_2,summary_value) %>%
           dplyr::rename(Dxy=summary_value)%>%
           inner_join(FastEPRR_result,by=c('chromosome','window_pos_1','window_pos_2'))  %>%  
           select(-CIL,-CIR) %>%
           filter(population != "ML") %>%
           filter(
                !(
                   (grepl("AB", pop_pair) & grepl("AB", population))  |
                   (grepl("BM", pop_pair) & grepl("BM", population))  |
                   (grepl("QH", pop_pair) & grepl("QH", population))  |
                   (grepl("NE", pop_pair) & grepl("NE", population))  |
                   (grepl("HN", pop_pair) & grepl("HN", population))  |
                   (grepl("TW", pop_pair) & grepl("TW", population))  |
                   (grepl("Central", pop_pair) & grepl("Central", population))
                 ))     
for (target_pop in c("AB","BM","QH","HN","TW","NE","Central")){
  assign(paste(target_pop,"Rho_Dxy",sep="_"), Dxy_Rho %>%
           filter(population==target_pop) %>%
           na.omit()  %>%
           group_by(pop_pair) %>%
           dplyr::summarize(single_pop_Rho=cor(Dxy,Rho,method = "spearman")) %>%
           `names<-`(c("pop_pair",paste(target_pop,"Rho",sep="."))))  
}  
Rho_anotherPair_Dxy <- full_join(AB_Rho_Dxy,BM_Rho_Dxy,by="pop_pair") %>%
                       full_join(QH_Rho_Dxy,by="pop_pair") %>%
                       full_join(NE_Rho_Dxy,by="pop_pair") %>%
                       full_join(HN_Rho_Dxy,by="pop_pair") %>%
                       full_join(TW_Rho_Dxy,by="pop_pair") %>%
　　　　　　　　　　　 full_join(Central_Rho_Dxy,by="pop_pair") %>%
                       arrange(pop_pair)
combined_reuslt <- full_join(Rho_anotherPair_Fst,Rho_anotherPair_Dxy,
                              by="pop_pair",suffix = c(paste(".","Fst",sep=""),'.Dxy')) %>%
                   full_join(combined_between_Pop_params,by="pop_pair") %>%
                   select(pop_pair,cor.Rho,cor.π,cor.Dxy_Fst,
                          AB.Rho.Fst,BM.Rho.Fst,QH.Rho.Fst,NE.Rho.Fst,HN.Rho.Fst,TW.Rho.Fst,Central.Rho.Fst,
                          AB.Rho.Dxy,BM.Rho.Dxy,QH.Rho.Dxy,NE.Rho.Dxy,HN.Rho.Dxy,TW.Rho.Dxy,Central.Rho.Dxy) %>%  
                   dplyr::rename("Rhoi~Rhoj"=cor.Rho,"πi~πj"=cor.π,"FstI~DxyI"=cor.Dxy_Fst,
                                 "Rho_AB~Fst"=AB.Rho.Fst,"Rho_BM~Fst"=BM.Rho.Fst,"Rho_QH~Fst"=QH.Rho.Fst,"Rho_NE~Fst"=NE.Rho.Fst,
                                 "Rho_HN~Fst"=HN.Rho.Fst,"Rho_TW~Fst"=TW.Rho.Fst,"Rho_Central~Fst"=Central.Rho.Fst,
                                 "Rho_AB~Dxy"=AB.Rho.Dxy,"Rho_BM~Dxy"=BM.Rho.Dxy,"Rho_QH~Dxy"=QH.Rho.Dxy,"Rho_NE~Dxy"=NE.Rho.Dxy,
                                 "Rho_HN~Dxy"=HN.Rho.Dxy,"Rho_TW~Dxy"=TW.Rho.Dxy,"Rho_Central~Dxy"=Central.Rho.Dxy)  
final_summary_reuslt <- full_join(combined_reuslt,combined_SamPair_Fst_π_Dxy_corr,by="pop_pair")
write.table(final_summary_reuslt ,file="combined_summary.corr.tab",quote=F,row.names=F,sep="\t")
rm(Dxy_Rho,QH_Rho_Dxy,AB_Rho_Dxy,BM_Rho_Dxy,HN_Rho_Dxy,TW_Rho_Dxy,NE_Rho_Dxy,Central_Rho_Dxy,Rho_anotherPair_Fst,Rho_anotherPair_Dxy,combined_between_Pop_params,combined_reuslt)

                 #######################violin plot###############
tab1 <- exclude_pseudo_replicated_params$`FstI~FstJ`  
tab2 <- exclude_pseudo_replicated_params$`DxyI~DxyJ`
tab3 <- final_summary_reuslt$`Rhoi~Rhoj`             
tab4 <- final_summary_reuslt$`πi~πj`
tab5 <- final_summary_reuslt$`FstI~DxyI`
tab6 <- (final_summary_reuslt %>%  select(5:11)  %>% gather(AnotherPair, "Rhoi~FSTJ") %>% select("Rhoi~FSTJ") %>% na.omit())$`Rhoi~FSTJ`
tab7 <- (final_summary_reuslt %>%  select(12:18) %>% gather(AnotherPair, "Rhoi~DxyJ") %>% select("Rhoi~DxyJ") %>% na.omit())$`Rhoi~DxyJ`
tab8 <- (final_summary_reuslt %>%  select(19:25) %>% gather(AnotherPair, "πi~FSTI") %>% select("πi~FSTI") %>% na.omit())$`πi~FSTI`
tab9 <- (final_summary_reuslt %>%  select(26:32) %>% gather(AnotherPair, "πi~DxyI") %>% select("πi~DxyI") %>% na.omit())$`πi~DxyI`
 
max_ln <- max(c(length(tab1),length(tab2), length(tab3),length(tab4),length(tab5),length(tab6),length(tab7),length(tab8),length(tab9)))
 
exclude_pseudo_replicated_Corr  <- data.frame(`FstI~FstJ` = c(tab1,rep(NA, max_ln - length(tab1))),
                                              `DxyI~DxyJ` = c(tab2,rep(NA, max_ln - length(tab2))),
                                              `Rhoi~Rhoj` = c(tab3,rep(NA, max_ln - length(tab3))),
                                              `πi~πj` = c(tab4,rep(NA, max_ln - length(tab4))),
                                              `FstI~DxyI` = c(tab5,rep(NA, max_ln - length(tab5))),
                                              `Rhoi~FSTJ` = c(tab6,rep(NA, max_ln - length(tab6))),
                                              `Rhoi~DxyJ` = c(tab7,rep(NA, max_ln - length(tab7))),
                                              `πi/j~FSTI` = c(tab8,rep(NA, max_ln - length(tab8))),
                                              `πi/j~DxyI` = c(tab9,rep(NA, max_ln - length(tab9))),
                                              check.names=FALSE) 
params_corr <-  exclude_pseudo_replicated_Corr %>%
                gather(summ_corr,cor_val) %>%
                na.omit() %>%
                ggplot(aes(summ_corr,cor_val)) +
                geom_violin(aes(fill = summ_corr)) +
                coord_flip()+
                theme_classic() +
　　　　　　　　　 theme(legend.position = "none") +
                ylim(-0.5,1) +
                scale_fill_brewer(palette = "Pastel1")
ggsave("Firure2A_violin_plot.pdf",width=3,height = 3)


