
library(tidyverse)


ptb_gene <- read.table("/home/chengennong/code-manual/vscode/ylp/Science2022/sample_info/perturb_gene.txt",sep="\t",header=1)
ptb_sgRNA <- read.table("/home/chengennong/project/ylp/Science2022/screening/mageck_count/CRISPRa.CalabreseSet.normalized_count.txt",sep="\t",header=1)

ref_csv <- ptb_gene %>% 
    inner_join(ptb_sgRNA %>% select(sgRNA,Gene)) %>% 
    separate(sgRNA,c("sgRNA_gene","sgRNA_seq"),sep="_") %>% 
    group_by(Gene) %>% 
    mutate(sgRNA_id=rank(sgRNA_seq)) %>% 
    mutate(id=paste(sgRNA_gene,sgRNA_id,sep="-")) %>%
    mutate(name=id,read="R2",pattern="(BC)GTTTAAGAGCTATG",sequence=sgRNA_seq,feature_type="CRISPR Guide Capture") %>% 
    select(id,name,read,pattern,sequence,feature_type)

write.table(ref_csv,"/home/chengennong/code-manual/vscode/ylp/Science2022/sample_info/feature_ref.csv",sep=",",row.names=F,col.names=T,quote=F)

