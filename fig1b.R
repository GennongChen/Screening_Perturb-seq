#cda R
#/home/wucheng/miniconda/envs/R4.1.0/lib/R/bin/R
library(tidyverse)
library(ggpubr)
library(RColorBrewer)

test_dir <- "/home/chengennong/project/ylp/Science2022/screening/mageck_test"

files <- c("CalabreseSet_IL2.sgrna_summary.txt","CalabreseSet_IFNG.sgrna_summary.txt")
markers <- c("IL2","IFNG")
gene_sh <- c("NO-TARGET", "IL2", "IFNG", "CD28", "VAV1", "SLA2", "MAP4K1")
ColPaired <- brewer.pal(12, "Paired"); red <- ColPaired[6]; blue <- ColPaired[2]; green <- ColPaired[4]; yellow <- "#E69F00"#ColPaired[8]
p_1b_list = list()
for (n_file in 1:length(files)){
    marker <- markers[n_file]
    print(marker)
    a_marker_raw_tests_file <- paste0(test_dir,"/",files[n_file])
    a_marker_raw_tests_table <- read.table(a_marker_raw_tests_file,sep="\t",header=1)
    a_marker_sgrna_tests_table <- a_marker_raw_tests_table %>% 
        mutate(marker=marker) %>% 
        filter(str_detect(Gene, "^IL2$|^IFNG$|^CD28$|^VAV1$|^SLA2$|^MAP4K1$|^NO-TARGET$")) %>% 
        separate("sgrna",c("sgrna_gene","sgrna_seq","sgrna_donor"),sep="_") %>% 
        mutate(
        type = case_when(
          str_detect("IL2|IFNG",Gene) ~ "Control",
          str_detect("CD28|VAV1",Gene)        ~ "Positive",
          str_detect("SLA2|MAP4K1",Gene)        ~ "Negative",
          str_detect("NO-TARGET",Gene)        ~ "NO-TARGET",
          TRUE                      ~ "other"
        ) 
      )   %>% mutate(sgrna_gene=factor(sgrna_gene,levels=gene_sh)) %>% 
        mutate(type=factor(type,levels=c("NO-TARGET","Control","Positive","Negative")))
    a_marker_sgrna_lfc_mean_table <- a_marker_sgrna_tests_table %>% 
        group_by(sgrna_gene,sgrna_seq,type,marker) %>% 
        summarise(LFC=mean(LFC))
    flowcy <- paste0("log2FoldChange(",marker,"hi","/",marker,"lo)")
    p_1b_density <-ggplot(a_marker_sgrna_lfc_mean_table, aes(LFC)) +
        geom_density(adjust=1.5, alpha=0.9,fill="grey") +
        scale_x_continuous(limits = c(-5,5),breaks=seq(-5, 5, 2.5)) +
        theme_bw()+
        theme(legend.position='none',
              panel.grid =element_blank())+
        labs(x = flowcy,y = "") 
    flowcy <- paste0("log2FoldChange(",marker,"hi","/",marker,"lo)")
    p_1b_vline <-ggplot(a_marker_sgrna_lfc_mean_table, aes(LFC, colour = type)) +
        geom_vline(aes(xintercept = LFC, colour = type),size=2, alpha=0.3) +
        facet_grid(sgrna_gene ~ .)+
        scale_x_continuous(limits = c(-5,5),breaks=seq(-5, 5, 2.5)) +
        scale_color_manual(values = c("grey",yellow,red,blue)) +
        theme_bw()+
        theme(legend.position='none',
              panel.grid =element_blank(),
              strip.background = element_rect(color="white", fill="white", size=1.5))+
        labs(x = flowcy,y = "") 
    
    p_1b <- ggarrange(p_1b_density,p_1b_vline,nrow=2, heights = c(1,2.3))
    p_1b_list[[marker]] <- p_1b
}


p_1b_total <- ggarrange(plotlist=p_1b_list,ncol=2)
ggsave(p_1b_total,filename="/home/chengennong/code-manual/vscode/ylp/Science2022/fig/main/1b.pdf",width=6,height=8)

