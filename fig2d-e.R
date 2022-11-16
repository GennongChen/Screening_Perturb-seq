#/home/wucheng/miniconda/envs/R4.1.0/lib/R/bin/R

library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(ggrepel)
library(colorspace)


nfkbs_path_gene <- read.table("/home/chengennong/code-manual/vscode/ylp/Science2022/sample_info/nfkb_regulators.txt",sep="\t") %>% pull(V1)
tcrs_path_gene <- read.table("/home/chengennong/code-manual/vscode/ylp/Science2022/sample_info/hsa04660.txt",sep="\t") %>% separate(V2,c("gene","des"),sep=";") %>% pull(gene)

test_dir <- "/home/chengennong/project/ylp/Science2022/screening/mageck_test"

files <- c("CalabreseSet_IL2.sgrna_summary.txt","CalabreseSet_IFNG.sgrna_summary.txt","DolcettoSet_IL2.sgrna_summary.txt","DolcettoSet_IFNG.sgrna_summary.txt")
gene_files <- c("CalabreseSet_IL2.gene_summary.txt","CalabreseSet_IFNG.gene_summary.txt","DolcettoSet_IL2.gene_summary.txt","DolcettoSet_IFNG.gene_summary.txt")
markers <- c("IL2","IFNG","IL2","IFNG")
strategies <- c("CRISPRa","CRISPRa","CRISPRi","CRISPRi")
ColPaired <- brewer.pal(12, "Paired"); red <- ColPaired[6]; blue <- ColPaired[2]; green <- ColPaired[4]; yellow <- "#E69F00"; pink <- "#D499B9"
###IL2 CRISPRa 2d
marker <- markers[1]
strategy <- strategies[1]
print(marker)
a_marker_raw_gene_tests_file <- paste0(test_dir,"/",gene_files[1])
a_marker_raw_gene_tests_table <- read.table(a_marker_raw_gene_tests_file,sep="\t",header=1) %>% rename(sgrna_gene = id) %>% mutate(Zscore=scale(neg.lfc))
a_marker_raw_tests_file <- paste0(test_dir,"/",files[1])
a_marker_raw_tests_table <- read.table(a_marker_raw_tests_file,sep="\t",header=1) 
a_strategy_sgrna_gene_tests_table <- a_marker_raw_tests_table %>%  # filter(Gene=="CD4")%>%
    mutate(marker=marker) %>% 
    separate("sgrna",c("sgrna_gene","sgrna_seq","sgrna_donor"),sep="_") %>% 
    inner_join(a_marker_raw_gene_tests_table)  %>% 
    group_by(sgrna_gene,marker,pos.fdr,neg.fdr,pos.rank,neg.rank,Zscore)  %>% 
    mutate(med_LFC=median(LFC)) %>% 
    group_by(sgrna_gene,marker,sgrna_donor,pos.fdr,neg.fdr,pos.rank,neg.rank,med_LFC,Zscore)  %>%
    summarise(zscore_donor_LFC=(LFC-mean(LFC))/sd(LFC), med_zscore_donor_LFC=median((LFC-mean(LFC))/sd(LFC))) %>%
    mutate(
    type = case_when(
      #abs(LFC) <= 0.5 | (pos.fdr >= 0.05 & neg.fdr >= 0.05)  ~ "Control",
      med_LFC > 0.5 & pos.fdr < 0.05        ~ "Positive regulator",
      med_LFC < -0.5 & neg.fdr < 0.05        ~ "Negative regulator",
      TRUE                      ~ "Not a hit"
    ))   %>%  select(-zscore_donor_LFC) %>% unique %>% 
    pivot_wider(names_from = sgrna_donor, values_from = med_zscore_donor_LFC) %>%
    mutate(label = ifelse(pos.rank < 16 | neg.rank < 16, sgrna_gene, "")) %>% 
    mutate(type=factor(type,levels=c("Positive regulator","Not a hit","Negative regulator"))) %>% 
    mutate(mean_med_zscore_donor_LFC=(r0+r1)/2) %>% 
    mutate(strategy=strategy)
    #mutate(label = ifelse(sgrna_gene %in% sgrna_gene_highlight, sgrna_gene, ""))
df_list_IL2 <- list()
df_list_IL2[[strategy]] <- a_strategy_sgrna_gene_tests_table 

###IL2 CRISPRi 2d
marker <- markers[3]
strategy <- strategies[3]
print(marker)
a_marker_raw_gene_tests_file <- paste0(test_dir,"/",gene_files[3])
a_marker_raw_gene_tests_table <- read.table(a_marker_raw_gene_tests_file,sep="\t",header=1) %>% rename(sgrna_gene = id) %>% mutate(Zscore=scale(neg.lfc))
a_marker_raw_tests_file <- paste0(test_dir,"/",files[3])
a_marker_raw_tests_table <- read.table(a_marker_raw_tests_file,sep="\t",header=1)
a_strategy_sgrna_gene_tests_table <- a_marker_raw_tests_table %>% 
    mutate(marker=marker) %>% 
    separate("sgrna",c("sgrna_gene","sgrna_seq","sgrna_donor"),sep="_") %>% 
    inner_join(a_marker_raw_gene_tests_table) %>% 
    group_by(sgrna_gene,marker,pos.fdr,neg.fdr,pos.rank,neg.rank,Zscore)  %>% 
    mutate(med_LFC=median(LFC)) %>% 
    group_by(sgrna_gene,marker,sgrna_donor,pos.fdr,neg.fdr,pos.rank,neg.rank,med_LFC,Zscore)  %>% 
    summarise(zscore_donor_LFC=(LFC-mean(LFC))/sd(LFC), med_zscore_donor_LFC=median((LFC-mean(LFC))/sd(LFC))) %>% 
    mutate(
    type = case_when(
      #abs(LFC) <= 0.5 | (pos.fdr >= 0.05 & neg.fdr >= 0.05)  ~ "Control",
      med_LFC > 0.5 & pos.fdr < 0.05        ~ "Negative regulator",
      med_LFC < -0.5 & neg.fdr < 0.05        ~ "Positive regulator",
      TRUE                      ~ "Not a hit"
    ))   %>%  select(-zscore_donor_LFC) %>% unique %>% 
    pivot_wider(names_from = sgrna_donor, values_from = med_zscore_donor_LFC) %>% 
    mutate(label = ifelse(pos.rank < 16 | neg.rank < 16, sgrna_gene, "")) %>% 
    mutate(type=factor(type,levels=c("Positive regulator","Not a hit","Negative regulator"))) %>% 
    mutate(mean_med_zscore_donor_LFC=(r0+r1)/2) %>% 
    mutate(strategy=strategy)
    #mutate(label = ifelse(sgrna_gene %in% sgrna_gene_highlight, sgrna_gene, ""))

df_list_IL2[[strategy]] <- a_strategy_sgrna_gene_tests_table

a_2strategy_sgrna_gene_tests_fea_table0 <- Reduce(full_join,df_list_IL2)  %>% ungroup %>% #filter(sgrna_gene=="IL2") %>%
    select(sgrna_gene,type,marker,strategy) %>% mutate(type=as.character(type)) %>% 
    pivot_wider(names_from = strategy, values_from = type) %>% 
    mutate(type = case_when(
          (str_detect("Not a hit",CRISPRa) & str_detect("Not a hit",CRISPRi))  ~ "Not a hit",
          (str_detect("Positive regulator|Negative regulator",CRISPRa) & str_detect("Not a hit",CRISPRi))  ~ "CRISPRa Only",
          (str_detect("Positive regulator|Negative regulator",CRISPRi) & str_detect("Not a hit",CRISPRa))  ~ "CRISPRi Only",
          (str_detect("Positive regulator|Negative regulator",CRISPRa) & str_detect("Positive regulator|Negative regulator",CRISPRi))  ~ "CRISPRi and CRISPRa",
          TRUE  ~ "other")) %>% 
    mutate(type = ifelse(!sgrna_gene %in% tcrs_path_gene, "Not KEGG TCR pathway", type)) %>% 
    select(sgrna_gene, type)

a_2strategy_sgrna_gene_tests_fea_table <- Reduce(full_join,df_list_IL2)  %>% ungroup %>% #filter(sgrna_gene=="IL2") %>%
    select(sgrna_gene,marker,strategy,Zscore) %>% 
    inner_join(a_2strategy_sgrna_gene_tests_fea_table0) %>% 
    pivot_wider(names_from = strategy, values_from = Zscore) %>% 
    mutate(type=as.character(type)) %>% 
    mutate(label = ifelse(str_detect(type,"CRISPR"), sgrna_gene, ""))%>% 
    mutate(type=factor(type,levels=c("CRISPRa Only","CRISPRi Only","CRISPRi and CRISPRa","Not a hit","Not KEGG TCR pathway"))) # %>% as.data.frame %>% head

p_2d_total <-ggplot(a_2strategy_sgrna_gene_tests_fea_table, aes(-CRISPRi, CRISPRa)) +
    geom_point(aes(color = type, fill = type),size = 1, alpha = 0.5, shape = 21)  +
    #geom_density_2d(color="black") +
    geom_vline(aes(xintercept = 0), colour = "grey",linetype=2) +
    geom_hline(aes(yintercept = 0), colour = "grey",linetype=2) +
    geom_text_repel(aes(label = label, color = type),size = 9/.pt, point.padding = 0.1, box.padding = 0.6,,max.overlaps = 1000,seed = 7654) +
    #scale_x_continuous(limits = c(-5,5),breaks=seq(-5, 5, 2.5)) +
    theme_bw()+
    labs(x = "CRISPRi IL-2 Screen Z-score(IL2hi/IL2lo)",y = "CRISPRa IL-2 Screen Z-score(IFNGhi/IFNGlo)", title="IL-2 Screens:\nKEGG T CELL RECEPTOR SIGNALING PATHWAY") +
    scale_fill_manual(values = c(pink,green,red,"black","grey")) +
    scale_color_manual(values = darken(c(pink,green,red,"black","grey"), 0.3))+
    theme(legend.position='right',
          panel.grid =element_blank()) #, ## 删去网格线 panel.border = element_blank())+ ## 删去外层边框
ggsave(p_2d_total,filename="/home/chengennong/code-manual/vscode/ylp/Science2022/fig/main/2d.pdf",width=7,height=6)

#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################

#a_strategy_sgrna_gene_tests_table %>% filter(sgrna_gene=="CD4")  %>% as.data.frame
#Reduce(full_join,df_list_IL2)  %>% filter(sgrna_gene=="CD4")  %>% as.data.frame
#a_2strategy_sgrna_gene_tests_fea_table %>% filter(label!="")

###IFNG CRISPRa 2f
marker <- markers[2]
strategy <- strategies[2]
print(marker)
a_marker_raw_gene_tests_file <- paste0(test_dir,"/",gene_files[2])
a_marker_raw_gene_tests_table <- read.table(a_marker_raw_gene_tests_file,sep="\t",header=1) %>% rename(sgrna_gene = id) %>% mutate(Zscore=scale(neg.lfc))
a_marker_raw_tests_file <- paste0(test_dir,"/",files[2])
a_marker_raw_tests_table <- read.table(a_marker_raw_tests_file,sep="\t",header=1)
a_strategy_sgrna_gene_tests_table <- a_marker_raw_tests_table %>% 
    mutate(marker=marker) %>% 
    separate("sgrna",c("sgrna_gene","sgrna_seq","sgrna_donor"),sep="_") %>% 
    inner_join(a_marker_raw_gene_tests_table) %>% 
    group_by(sgrna_gene,marker,pos.fdr,neg.fdr,pos.rank,neg.rank,Zscore)  %>% 
    mutate(med_LFC=median(LFC)) %>% 
    group_by(sgrna_gene,marker,sgrna_donor,pos.fdr,neg.fdr,pos.rank,neg.rank,med_LFC,Zscore)  %>% 
    summarise(zscore_donor_LFC=(LFC-mean(LFC))/sd(LFC), med_zscore_donor_LFC=median((LFC-mean(LFC))/sd(LFC))) %>% 
    mutate(
    type = case_when(
      #abs(LFC) <= 0.5 | (pos.fdr >= 0.05 & neg.fdr >= 0.05)  ~ "Control",
      med_LFC > 0.5 & pos.fdr < 0.05        ~ "Positive regulator",
      med_LFC < -0.5 & neg.fdr < 0.05        ~ "Negative regulator",
      TRUE                      ~ "Not a hit"
    ))   %>%  select(-zscore_donor_LFC) %>% unique %>% 
    pivot_wider(names_from = sgrna_donor, values_from = med_zscore_donor_LFC) %>% 
    mutate(label = ifelse(pos.rank < 16 | neg.rank < 16, sgrna_gene, "")) %>% 
    mutate(type=factor(type,levels=c("Positive regulator","Not a hit","Negative regulator"))) %>% 
    mutate(mean_med_zscore_donor_LFC=(r0+r1)/2) %>% 
    mutate(strategy=strategy)
    #mutate(label = ifelse(sgrna_gene %in% sgrna_gene_highlight, sgrna_gene, ""))
df_list_IFNG <- list()
df_list_IFNG[[strategy]] <- a_strategy_sgrna_gene_tests_table 

###IFNG CRISPRi 2e
marker <- markers[4]
strategy <- strategies[4]
print(marker)
a_marker_raw_gene_tests_file <- paste0(test_dir,"/",gene_files[4])
a_marker_raw_gene_tests_table <- read.table(a_marker_raw_gene_tests_file,sep="\t",header=1) %>% rename(sgrna_gene = id) %>% mutate(Zscore=scale(neg.lfc))
a_marker_raw_tests_file <- paste0(test_dir,"/",files[4])
a_marker_raw_tests_table <- read.table(a_marker_raw_tests_file,sep="\t",header=1)
a_strategy_sgrna_gene_tests_table <- a_marker_raw_tests_table %>% 
    mutate(marker=marker) %>% 
    separate("sgrna",c("sgrna_gene","sgrna_seq","sgrna_donor"),sep="_") %>% 
    inner_join(a_marker_raw_gene_tests_table) %>% 
    group_by(sgrna_gene,marker,pos.fdr,neg.fdr,pos.rank,neg.rank,Zscore)  %>% 
    mutate(med_LFC=median(LFC)) %>% 
    group_by(sgrna_gene,marker,sgrna_donor,pos.fdr,neg.fdr,pos.rank,neg.rank,med_LFC,Zscore)  %>% 
    summarise(zscore_donor_LFC=(LFC-mean(LFC))/sd(LFC), med_zscore_donor_LFC=median((LFC-mean(LFC))/sd(LFC))) %>% 
    mutate(
    type = case_when(
      #abs(LFC) <= 0.5 | (pos.fdr >= 0.05 & neg.fdr >= 0.05)  ~ "Control",
      med_LFC > 0.5 & pos.fdr < 0.05        ~ "Negative regulator",
      med_LFC < -0.5 & neg.fdr < 0.05        ~ "Positive regulator",
      TRUE                      ~ "Not a hit"
    ))   %>%  select(-zscore_donor_LFC) %>% unique %>% 
    pivot_wider(names_from = sgrna_donor, values_from = med_zscore_donor_LFC) %>% 
    mutate(label = ifelse(pos.rank < 16 | neg.rank < 16, sgrna_gene, "")) %>% 
    mutate(type=factor(type,levels=c("Positive regulator","Not a hit","Negative regulator"))) %>% 
    mutate(mean_med_zscore_donor_LFC=(r0+r1)/2) %>% 
    mutate(strategy=strategy)
    #mutate(label = ifelse(sgrna_gene %in% sgrna_gene_highlight, sgrna_gene, ""))

df_list_IFNG[[strategy]] <- a_strategy_sgrna_gene_tests_table 

a_2strategy_sgrna_gene_tests_fea_table0 <- Reduce(full_join,df_list_IFNG)  %>% ungroup %>% #filter(sgrna_gene=="CCL4") %>%
    select(sgrna_gene,type,marker,strategy) %>% mutate(type=as.character(type)) %>% 
    pivot_wider(names_from = strategy, values_from = type) %>% 
    mutate(type = case_when(
          (str_detect("Not a hit",CRISPRa) & str_detect("Not a hit",CRISPRi))  ~ "Not a Hit",
          (str_detect("Positive regulator|Negative regulator",CRISPRa) & str_detect("Not a hit",CRISPRi))  ~ "CRISPRa Only",
          (str_detect("Positive regulator|Negative regulator",CRISPRi) & str_detect("Not a hit",CRISPRa))  ~ "CRISPRi Only",
          (str_detect("Positive regulator|Negative regulator",CRISPRa) & str_detect("Positive regulator|Negative regulator",CRISPRi))  ~ "CRISPRi and CRISPRa",
          TRUE  ~ "other")) %>% 
    mutate(type = ifelse(!sgrna_gene %in% nfkbs_path_gene, "Not NF-KB pathway Regulators", type)) %>% 
    select(sgrna_gene, type)

a_2strategy_sgrna_gene_tests_fea_table <- Reduce(full_join,df_list_IFNG)  %>% ungroup %>% #filter(sgrna_gene=="IL2") %>%
    select(sgrna_gene,marker,Zscore,strategy) %>% 
    inner_join(a_2strategy_sgrna_gene_tests_fea_table0) %>% 
    pivot_wider(names_from = strategy, values_from = Zscore) %>% 
    mutate(type=as.character(type))  %>% 
    mutate(label = ifelse(!str_detect(type,"Regulators"), sgrna_gene, ""))%>% 
    mutate(type=factor(type,levels=c("CRISPRa Only","CRISPRi Only","CRISPRi and CRISPRa","Not a Hit","Not NF-KB pathway Regulators"))) 

#a_2strategy_sgrna_gene_tests_fea_table %>% filter(is.na(type))
#a_2strategy_sgrna_gene_tests_fea_table0 %>% filter(sgrna_gene=="CCL4")
#Reduce(full_join,df_list_IFNG) %>% filter(sgrna_gene=="CCL4") %>% as.data.frame

p_2e_total <-ggplot(a_2strategy_sgrna_gene_tests_fea_table%>% filter(!is.na(type)), aes(-CRISPRi, CRISPRa)) +
    geom_point(aes(color = type, fill = type),size = 1, alpha = 0.5, shape = 21)  +
    #geom_density_2d(color="black") +
    geom_vline(aes(xintercept = 0), colour = "grey",linetype=2) +
    geom_hline(aes(yintercept = 0), colour = "grey",linetype=2) +
    geom_text_repel(aes(label = label, color = type),size = 9/.pt, point.padding = 0.1, box.padding = 0.6,max.overlaps = 10000,seed = 7654) +
    #scale_x_continuous(limits = c(-5,5),breaks=seq(-5, 5, 2.5)) +
    theme_bw()+
    labs(x = "CRISPRi IL-2 Screen Z-score(IL2hi/IL2lo)",y = "CRISPRa IL-2 Screen Z-score(IFNGhi/IFNGlo)", title="IFN-G Screens: NF-KB Pathway Regulators") +
    scale_fill_manual(values = c(pink,green,red,"black","grey")) +
    scale_color_manual(values = darken(c(pink,green,red,"black","grey"), 0.3))+
    theme(legend.position='right',
          panel.grid =element_blank()) #, ## 删去网格线 panel.border = element_blank())+ ## 删去外层边框
ggsave(p_2e_total,filename="/home/chengennong/code-manual/vscode/ylp/Science2022/fig/main/2e.pdf",width=7,height=6)
