#!/home/chengennong/anaconda3/envs/R/bin/R
library("tidyverse")

#CRISPRa
count_dir <- "/home/chengennong/project/ylp/Science2022/screening/mageck_count"

aA_count_file <- "CRISPRa.CalabreseSetA.count.txt"
aA_raw_counts_file <- paste0(count_dir,"/",aA_count_file)
aA_raw_counts_table <- read.table(aA_raw_counts_file,sep="\t",header=1)
colnames(aA_raw_counts_table)  <- colnames(aA_raw_counts_table) %>% str_replace("SetA","Set")
#aA_raw_nor_counts_table <- aA_raw_counts_table %>% mutate_at(vars(-sgRNA,-Gene), funs(./sum(.)*1e6))

aB_count_file <- "CRISPRa.CalabreseSetB.count.txt"
aB_raw_counts_file <- paste0(count_dir,"/",aB_count_file)
aB_raw_counts_table <- read.table(aB_raw_counts_file,sep="\t",header=1)
colnames(aB_raw_counts_table)  <- colnames(aB_raw_counts_table) %>% str_replace("SetB","Set")
#aB_raw_nor_counts_table <- aB_raw_counts_table %>% mutate_at(vars(-sgRNA,-Gene), funs(./sum(.)*1e6))

#a_nor_counts_table <- rbind(aA_raw_nor_counts_table,aB_raw_nor_counts_table)
a_nor_counts_table <- rbind(aA_raw_counts_table,aB_raw_counts_table) %>% mutate_at(vars(-sgRNA,-Gene), ~(./sum(.)*1e6))
write.table(a_nor_counts_table,"/home/chengennong/project/ylp/Science2022/screening/mageck_count/CRISPRa.CalabreseSet.normalized_count.txt",row.names=F,quote=F,sep="\t")
write.table((a_nor_counts_table %>% filter(Gene=="NO-TARGET") %>%select(sgRNA)),
    "/home/chengennong/project/ylp/Science2022/screening/mageck_count/CRISPRa.CalabreseSet.no_target.txt",row.names=F,quote=F,sep="\t")


#CRISPRi
count_dir <- "/home/chengennong/project/ylp/Science2022/screening/mageck_count"

iA_count_file <- "CRISPRi.DolcettoSetA.count.txt"
iA_raw_counts_file <- paste0(count_dir,"/",iA_count_file)
iA_raw_counts_table <- read.table(iA_raw_counts_file,sep="\t",header=1)
#iA_raw_nor_counts_table <- iA_raw_counts_table %>% mutate_at(vars(-sgRNA,-Gene), funs(./sum(.)*1e6))
colnames(iA_raw_counts_table)  <- colnames(iA_raw_counts_table) %>% str_replace("SetA","Set")

iB_count_file <- "CRISPRi.DolcettoSetB.count.txt"
iB_raw_counts_file <- paste0(count_dir,"/",iB_count_file)
iB_raw_counts_table <- read.table(iB_raw_counts_file,sep="\t",header=1)
#iB_raw_nor_counts_table <- iB_raw_counts_table %>% mutate_at(vars(-sgRNA,-Gene), funs(./sum(.)*1e6))
colnames(iB_raw_counts_table)  <- colnames(iB_raw_counts_table) %>% str_replace("SetB","Set")

#i_nor_counts_table <- rbind(iA_raw_nor_counts_table,iB_raw_nor_counts_table)
i_nor_counts_table <- rbind(iA_raw_counts_table,iB_raw_counts_table) %>% mutate_at(vars(-sgRNA,-Gene), ~(./sum(.)*1e6))

write.table(i_nor_counts_table,"/home/chengennong/project/ylp/Science2022/screening/mageck_count/CRISPRi.DolcettoSet.normalized_count.txt",row.names=F,quote=F,sep="\t")
write.table((i_nor_counts_table %>% filter(Gene=="NO-TARGET") %>%select(sgRNA)),
    "/home/chengennong/project/ylp/Science2022/screening/mageck_count/CRISPRi.DolcettoSet.no_target.txt",row.names=F,quote=F,sep="\t")
