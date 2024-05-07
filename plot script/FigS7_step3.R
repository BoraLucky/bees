suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
library(ggplot2)
library(gg.gap)

args <- commandArgs(TRUE)
input_path <- args[1]
out_path  <- args[2]

setwd(input_path)
myfiles <- list.files(pattern = "*_island_candidate.table") 

for (i in 1:length(myfiles)){
  assign(paste("pop",i,sep=""),strsplit(myfiles[i],split="_")[[1]][1])
}     


for (i in 1:length(myfiles)){
  assign(paste(get(paste("pop",i,sep="")),"_table",sep=""),read.table(myfiles[i],header=T))
}    

################################islands##############################
all_pop_stat <- data.frame()

for (m in paste("pop",1:length(myfiles),sep="")){
     assign(paste(get(m),"_summarise",sep=""),get(paste(get(m),"_table",sep="")) %>% 
           group_by(island_lenth) %>%
           summarise(count_sales = n()) %>%
           mutate(population=get(m)))
     all_pop_stat=rbind(all_pop_stat,get(paste(get(m),"_summarise",sep="")))
}  


p <-  all_pop_stat %>% mutate(island=island_lenth/1000)  %>%
      arrange(island,desc(population)) %>% 
      mutate(island=paste(island,"K",sep="") %>% fct_infreq())  %>%
      ggplot(aes(x=island,y=count_sales, fill=population)) +
      geom_bar(stat = "identity",position="dodge", colour="black") +
      scale_fill_brewer(palette="Pastel1") +
      scale_x_discrete(drop = FALSE)+
      theme_classic()+
      xlab("size distributions") + 
      ylab("number") +
      theme(axis.title.x = element_text(face="bold",size=15))+
      theme(axis.title.y = element_text(face="bold",size=15))+
      theme(legend.position=c(0.5,1), legend.justification=c(0.5,1))
ggsave(paste(out_path,"islands_distributions_base.pdf",sep="/"),p) 


picture <- gg.gap(plot = p,
           segments = c(25, 120),  
           tick_width = 10,       
           rel_heights = c(0.1, 0.2, 0.1),  
           ylim = c(0,160))      
ggsave(paste(out_path,"islands_distributions_subject.pdf",sep="/"),picture) 


