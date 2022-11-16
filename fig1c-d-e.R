#cda R
#/home/wucheng/miniconda/envs/R4.1.0/lib/R/bin/R
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(ggrepel)
library(colorspace)

test_dir <- "/home/chengennong/project/ylp/Science2022/screening/mageck_test"

files <- c("CalabreseSet_IL2.sgrna_summary.txt","CalabreseSet_IFNG.sgrna_summary.txt")
gene_files <- c("CalabreseSet_IL2.gene_summary.txt","CalabreseSet_IFNG.gene_summary.txt")
markers <- c("IL2","IFNG")
ColPaired <- brewer.pal(12, "Paired"); red <- ColPaired[6]; blue <- ColPaired[2]; green <- ColPaired[4]; yellow <- "#E69F00"
###IL2 1c
marker <- markers[1]
print(marker)
a_marker_raw_gene_tests_file <- paste0(test_dir,"/",gene_files[1])
a_marker_raw_gene_tests_table <- read.table(a_marker_raw_gene_tests_file,sep="\t",header=1) %>% rename(sgrna_gene = id)
a_marker_raw_tests_file <- paste0(test_dir,"/",files[1])
a_marker_raw_tests_table <- read.table(a_marker_raw_tests_file,sep="\t",header=1)
a_marker_sgrna_gene_tests_table <- a_marker_raw_tests_table %>% 
    mutate(marker=marker) %>% 
    separate("sgrna",c("sgrna_gene","sgrna_seq","sgrna_donor"),sep="_") %>% 
    inner_join(a_marker_raw_gene_tests_table) %>% 
    group_by(sgrna_gene,marker,pos.fdr,neg.fdr,pos.rank,neg.rank)  %>% 
    mutate(med_LFC=median(LFC)) %>% 
    group_by(sgrna_gene,marker,sgrna_donor,pos.fdr,neg.fdr,pos.rank,neg.rank,med_LFC)  %>% 
    summarise(mean_donor_LFC=mean(LFC)) %>% 
    mutate(
    type = case_when(
      #abs(LFC) <= 0.5 | (pos.fdr >= 0.05 & neg.fdr >= 0.05)  ~ "Control",
      med_LFC > 0.5 & pos.fdr < 0.05        ~ "Positive Hit",
      med_LFC < -0.5 & neg.fdr < 0.05        ~ "Negative Hit",
      TRUE                      ~ "Not a Hit"
    ))   %>% 
    pivot_wider(names_from = sgrna_donor, values_from = mean_donor_LFC) %>% 
    mutate(label = ifelse(pos.rank < 16 | neg.rank < 16, sgrna_gene, "")) %>% 
    mutate(type=factor(type,levels=c("Not a Hit","Positive Hit","Negative Hit"))) %>% 
    mutate(mean_donor_LFC=(r0+r1)/2)
    #mutate(label = ifelse(sgrna_gene %in% sgrna_gene_highlight, sgrna_gene, ""))

p_1c_total <-ggplot(a_marker_sgrna_gene_tests_table, aes(r0, r1)) +
    geom_point(aes(color = type, fill = type),size = 1, alpha = 0.5, shape = 21)  +
    geom_density_2d(color="black") +
    geom_vline(aes(xintercept = 0), colour = "grey",linetype=2) +
    geom_hline(aes(yintercept = 0), colour = "grey",linetype=2) +
    geom_text_repel(aes(label = label, color = type),size = 9/.pt, point.padding = 0.1, box.padding = 0.6,,max.overlaps = 1000,seed = 7654) +
    #scale_x_continuous(limits = c(-5,5),breaks=seq(-5, 5, 2.5)) +
    theme_bw()+
    labs(x = "Donor 1, log2FoldChange",y = "Donor 2, log2FoldChange", title="IL-2 CRISPRa Screen",fill="Screen Hit") +
    scale_fill_manual(values = c("grey",red,blue)) +
    scale_color_manual(values = darken(c("grey",red,blue), 0.3))+
    theme(legend.position='right',
          panel.grid =element_blank()) #, ## 删去网格线 panel.border = element_blank())+ ## 删去外层边框
ggsave(p_1c_total,filename="/home/chengennong/code-manual/vscode/ylp/Science2022/fig/main/1c.pdf",width=7,height=6)

df_list <- list()
df_list[[marker]] <- a_marker_sgrna_gene_tests_table 
###IFNG 1d
marker <- markers[2]
print(marker)
a_marker_raw_gene_tests_file <- paste0(test_dir,"/",gene_files[2])
a_marker_raw_gene_tests_table <- read.table(a_marker_raw_gene_tests_file,sep="\t",header=1) %>% rename(sgrna_gene = id)
a_marker_raw_tests_file <- paste0(test_dir,"/",files[2])
a_marker_raw_tests_table <- read.table(a_marker_raw_tests_file,sep="\t",header=1)
a_marker_sgrna_gene_tests_table <- a_marker_raw_tests_table %>% 
    mutate(marker=marker) %>% 
    separate("sgrna",c("sgrna_gene","sgrna_seq","sgrna_donor"),sep="_") %>% 
    inner_join(a_marker_raw_gene_tests_table) %>% 
    group_by(sgrna_gene,marker,pos.fdr,neg.fdr,pos.rank,neg.rank)  %>% 
    mutate(med_LFC=median(LFC)) %>% 
    group_by(sgrna_gene,marker,sgrna_donor,pos.fdr,neg.fdr,pos.rank,neg.rank,med_LFC)  %>% 
    summarise(mean_donor_LFC=mean(LFC)) %>% 
    mutate(
    type = case_when(
      #abs(LFC) <= 0.5 | (pos.fdr >= 0.05 & neg.fdr >= 0.05)  ~ "Control",
      med_LFC > 0.5 & pos.fdr < 0.05        ~ "Positive Hit",
      med_LFC < -0.5 & neg.fdr < 0.05        ~ "Negative Hit",
      TRUE                      ~ "Not a Hit"
    ))   %>% 
    pivot_wider(names_from = sgrna_donor, values_from = mean_donor_LFC) %>% 
    mutate(label = ifelse(pos.rank < 16 | neg.rank < 16, sgrna_gene, "")) %>% 
    mutate(type=factor(type,levels=c("Not a Hit","Positive Hit","Negative Hit"))) %>% 
    mutate(mean_donor_LFC=(r0+r1)/2)
    #mutate(label = ifelse(sgrna_gene %in% sgrna_gene_highlight, sgrna_gene, ""))

p_1d_total <-ggplot(a_marker_sgrna_gene_tests_table, aes(r0, r1)) +
    geom_point(aes(color = type, fill = type),size = 1, alpha = 0.5, shape = 21)  +
    geom_density_2d(color="black") +
    geom_vline(aes(xintercept = 0), colour = "grey",linetype=2) +
    geom_hline(aes(yintercept = 0), colour = "grey",linetype=2) +
    geom_text_repel(aes(label = label, color = type),size = 9/.pt, point.padding = 0.1, box.padding = 0.6,,max.overlaps = 1000,seed = 7654) +
    #scale_x_continuous(limits = c(-5,5),breaks=seq(-5, 5, 2.5)) +
    theme_bw()+
    labs(x = "Donor 1, log2FoldChange",y = "Donor 2, log2FoldChange", title="IFN-G CRISPRa Screen",fill="Screen Hit") +
    scale_fill_manual(values = c("grey",red,blue)) +
    scale_color_manual(values = darken(c("grey",red,blue), 0.3))+
    theme(legend.position='right',
          panel.grid =element_blank()) #, ## 删去网格线 panel.border = element_blank())+ ## 删去外层边框
ggsave(p_1d_total,filename="/home/chengennong/code-manual/vscode/ylp/Science2022/fig/main/1d.pdf",width=7,height=6)

df_list[[marker]] <- a_marker_sgrna_gene_tests_table 

###1e
names(df_list)

a_2marker_sgrna_gene_tests_fea_table0 <- Reduce(full_join,df_list)  %>% ungroup %>% 
    select(sgrna_gene,type,marker) %>% mutate(type=as.character(type)) %>% 
    pivot_wider(names_from = marker, values_from = type) %>% 
    mutate(type = case_when(
          (str_detect("Not a Hit",IL2) & str_detect("Not a Hit",IFNG))  ~ "Not a Hit",
          (str_detect("Positive Hit|Negative Hit",IL2) & str_detect("Not a Hit",IFNG))  ~ "IL-2 Only",
          (str_detect("Positive Hit|Negative Hit",IFNG) & str_detect("Not a Hit",IL2))  ~ "IFN-G Only",
          (str_detect("Positive Hit|Negative Hit",IL2) & str_detect("Positive Hit|Negative Hit",IFNG))  ~ "Both Screens",
          TRUE  ~ "other")) %>% 
    select(sgrna_gene, type)

a_2marker_sgrna_gene_tests_fea_table <- Reduce(full_join,df_list)  %>% ungroup %>% 
    select(sgrna_gene,marker,mean_donor_LFC) %>% 
    inner_join(a_2marker_sgrna_gene_tests_fea_table0) %>% 
    pivot_wider(names_from = marker, values_from = mean_donor_LFC) %>% 
    mutate(type=factor(type,levels=c("Both Screens","IFN-G Only","IL-2 Only","Not a Hit"))) %>% 
    mutate(label = ifelse(abs(IL2) > 1.9 | abs(IFNG) > 1.9, sgrna_gene, ""))

a_2marker_sgrna_gene_tests_fea_table %>% filter(sgrna_gene=="MAP4K1")


p_1e_total <-ggplot(a_2marker_sgrna_gene_tests_fea_table, aes(IL2, IFNG)) +
    geom_point(aes(color = type, fill = type),size = 1, alpha = 0.5, shape = 21)  +
    geom_density_2d(color="black") +
    geom_vline(aes(xintercept = 0), colour = "grey",linetype=2) +
    geom_hline(aes(yintercept = 0), colour = "grey",linetype=2) +
    geom_text_repel(aes(label = label, color = type),size = 9/.pt, point.padding = 0.1, box.padding = 0.6,,max.overlaps = 1000,seed = 7654) +
    #scale_x_continuous(limits = c(-5,5),breaks=seq(-5, 5, 2.5)) +
    theme_bw()+
    labs(x = "log2FoldChange(IL2hi/IL2lo)",y = "log2FoldChange(IFNGhi/IFNGlo)", title="CRISPRa Screens",fill="Screen Hit") +
    scale_fill_manual(values = c(blue,red,green,"grey")) +
    scale_color_manual(values = darken(c(blue,red,green,"grey"), 0.3))+
    theme(legend.position='right',
          panel.grid =element_blank()) #, ## 删去网格线 panel.border = element_blank())+ ## 删去外层边框
ggsave(p_1e_total,filename="/home/chengennong/code-manual/vscode/ylp/Science2022/fig/main/1e.pdf",width=7,height=6)
