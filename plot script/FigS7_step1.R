library(tidyr)
library(dplyr)
args <- commandArgs(TRUE)
outlier_table <- args[1]  
population <- args[2]     #pop name


            ###################################################################################

all_outlier <- read.table(outlier_table,sep='\t',header =T)

##########################adjacent window(overlap=1)####################################
window_pos_1_lead <- all_outlier  %>%  group_by(chromosome)  %>%
                     select(c(1,2,3))  %>%  
                     mutate(lead_win_pos_1 = lead(window_pos_1, 1, order_by = window_pos_1)) %>% 
                     mutate(overlap = lead_win_pos_1 - window_pos_2) 

island_index_up <- c(which(window_pos_1_lead$overlap == 1))
island_index_down  <- c(which(window_pos_1_lead$overlap == 1)+1) 
island_windows_index <- unique(sort(c(island_index_up,island_index_down)))
all_windows_index <- 1:nrow(window_pos_1_lead)
island_candidate_temp1 <- window_pos_1_lead[island_windows_index,] #adjacent window
other_windows_index  <-  setdiff(all_windows_index,island_windows_index)
other_windows  <- window_pos_1_lead[other_windows_index,] %>% 
                      select(-c("lead_win_pos_1","overlap"))  %>%
                      mutate(island_lenth = window_pos_2 - window_pos_1+1) %>%
                       mutate(window_number = (window_pos_2 - window_pos_1+1)/10000)
					   
############################repeat values index ###########################################
runs <- rle(island_candidate_temp1$overlap)
myruns = which(runs$lengths >= 2)
runs.lengths.cumsum = cumsum(runs$lengths)
ends = runs.lengths.cumsum[myruns]
newindex = ifelse(myruns>1, myruns-1, 0)
starts = runs.lengths.cumsum[newindex] + 1
if (0 %in% newindex) starts = c(1,starts)
Consecutive_Repeats <- data.frame(starts,ends)

##################################the first index of repeat values ######################
x <- c() 
for (i in (1:nrow(Consecutive_Repeats))){
   x[[i]] <- c((Consecutive_Repeats[i,]$starts+1):(Consecutive_Repeats[i,]$ends))
}
need_remove_index <- as.vector(unlist(x))
island_candidate_temp2 <- island_candidate_temp1[-need_remove_index,] 
##############################island windows###########################################################
island_candidate_temp3  <- island_candidate_temp2  %>% 
                          group_by(chromosome)   %>% 
                          mutate(lead_window_pos_2 = lead(window_pos_2, 1, order_by = window_pos_1)) %>%
                          filter(overlap==1)  %>%
                          select(-c("window_pos_2","lead_win_pos_1","overlap")) %>%
                          mutate(island_lenth = lead_window_pos_2 - window_pos_1+1) %>%
                          mutate(window_number = (lead_window_pos_2 - window_pos_1+1)/10000)
colnames(island_candidate_temp3)[3] <-"window_pos_2"   

island_candidate <- rbind(other_windows,island_candidate_temp3) %>% arrange(chromosome,window_pos_1)
write.table(island_candidate,file=paste(population,"island_candidate.table",sep="_"),quote=F,row.names=F,sep="\t")


