#########Looking overlap outlier window ############


suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(Repitools))
library(stringr)

args <- commandArgs(TRUE)
input_path <- args[1]
out_path  <- args[2]


setwd(input_path)
myfiles <- list.files(pattern = "*_island_candidate.table")  #list all file

for (i in 1:length(myfiles)){
  assign(paste("pop",i,sep=""),strsplit(myfiles[i],split="_")[[1]][1])
}    


for (i in 1:length(myfiles)){
  assign(paste(get(paste("pop",i,sep="")),"_table",sep=""),read.table(myfiles[i],header=T))
}   


for (m in paste("pop",1:length(myfiles),sep="")){
  assign(paste(get(m),"_grange",sep=""), 
         get(paste(get(m),"_table",sep="")) %>% 
           select(-c("island_lenth","window_number")) %>%
           rename(c("chromosome"="seqnames","window_pos_1"="start", "window_pos_2"="end")) %>%
           group_by(seqnames) %>%
           as_granges())
}   


for (j in 1:(length(myfiles)-1)){
  for (k in (j+1):length(myfiles)){
    head_pop <- get(paste("pop",j,sep=""))
    next_pop <- get(paste("pop",k,sep=""))
    head_pop_grange <- get(paste(get(paste("pop",j,sep="")),"_grange",sep=""))
    next_pop_grange <- get(paste(get(paste("pop",k,sep="")),"_grange",sep=""))
    assign(paste(head_pop,next_pop,"ovalap",sep="_"),
           get(paste(get(paste("pop",j,sep="")),"_table",sep="")) %>% 
             mutate(Overlaps=countOverlaps(head_pop_grange,next_pop_grange)) %>%
             filter(Overlaps==1)%>% mutate(overlap_window="Yes") %>%
             select(-Overlaps))
    assign(paste(next_pop,head_pop,"ovalap",sep="_"),
           get(paste(get(paste("pop",k,sep="")),"_table",sep="")) %>% 
             mutate(Overlaps=countOverlaps(next_pop_grange,head_pop_grange)) %>%
             filter(Overlaps==1) %>% mutate(overlap_window="Yes") %>%
             select(-Overlaps))
    write.table(get(paste(head_pop,next_pop,"ovalap",sep="_")),
                file=paste(out_path,paste(head_pop,next_pop,"ovalap_outliner_windows.table",sep="_"),sep = "/"),
                quote=F,row.names=F,sep="\t")
    write.table(get(paste(next_pop,head_pop,"ovalap",sep="_")),
                file=paste(out_path,paste(next_pop,head_pop,"ovalap_outliner_windows.table",sep="_"),sep = "/"),
                quote=F,row.names=F,sep="\t")
  }
}




