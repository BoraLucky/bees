library(tidyverse)
library(ggplot2)
library(dplyr)
#library(gg.gap)

args <- commandArgs(TRUE)
window_fst.txt <- args[1] 
window_pi.txt <- args[2]  
fst_site.txt <- args[3] 
Population <- args[4] 

fst_window <-read.table(window_fst.txt,header=T)
pi_window <-read.table(window_pi.txt,header= T)
fst_window<-na.omit(fst_window) 
pi_window<-na.omit(pi_window) 

fst_window_mean <- mean(fst_window$WEIGHTED_FST)
fst_window_sd <- sd(fst_window$WEIGHTED_FST) 
pi_window_mean <- mean(pi_window$PI)
pi_window_sd <- sd(pi_window$PI)
result <- data.frame(Population,fst_window_mean,fst_window_sd,pi_window_mean,pi_window_sd)
names(result) <-c("pop","fst_window_mean","fst_window_sd","pi_window_mean","pi_window_sd")


fst_site <-read.table(fst_site.txt,header= T)
fst_site <-na.omit(fst_site)
binsize <- diff(range(fst_site$WEIR_AND_COCKERHAM_FST))/50  
FST_distribution_His <- fst_site %>% 
                    ggplot(aes(x=WEIR_AND_COCKERHAM_FST)) + 
                    geom_histogram(fill="#FF8C64", colour="grey50", size=.2,binwidth=binsize)+
                    theme_classic()+xlab("Fst")+ylab("Number of SNPs")+
                    scale_x_continuous(limits = c(0,1),breaks = c(0,0.2,0.4,0.6,0.8,1))
ggsave(FST_distribution_His,filename=paste(Population,"fst_His_destribution.pdf",sep=""),width = 4,height = 3)

		    
perFst <- data.frame(fst_site$WEIR_AND_COCKERHAM_FST)
names(perFst) <- Population
write.table(perFst,file=paste(Population,"_PerSite_wcFst.txt",sep=""),quote=F,row.names=F,sep="\t")
