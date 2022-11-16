
#/home/wucheng/miniconda/envs/R4.1.0/lib/R/bin/R
library(Seurat)
library(tidyverse)
library(stringr)
library(ggplot2)
library(patchwork)

setwd("/home/chengennong/code-manual/vscode/ylp/Science2022/fig/main")
res_dir="/home/chengennong/project/ylp/Science2022/perturb/final_10x_outs"
res_dir="/home/chengennong/project/single_cell_fyz/outs_dir/B1/outs/filtered_feature_bc_matrix"
res_dir="/home/chengennong/project/ylp/Science2022/perturb/aggr/stim/outs/count/filtered_feature_bc_matrix"
crispr_file="/home/chengennong/project/ylp/Science2022/perturb/aggr/stim/outs/count/crispr_analysis/protospacer_calls_per_cell.csv"
crispr_table <- read.csv(crispr_file)
filtered_crispr_table  <- crispr_table %>% 
    filter(num_features==1) %>% 
    filter(as.numeric(num_umis) > 5) %>% 
    #separate(feature_call,c("gene","guide_id"),sep="-",remove=F, extra = "merge", fill = "left")
    mutate(gene=str_replace(feature_call,"-[0-9]+",""))

#################
# Load the stim dataset
#################
stim.data <- Read10X(data.dir = res_dir)
# Initialize the Seurat object with the raw (non-normalized data).
#stim <- CreateSeuratObject(counts = stim.data, project = "stim", min.cells = 3, min.features = 200)
stim <- CreateSeuratObject(counts = stim.data$`Gene Expression`, project = "stim", min.cells = 3, min.features = 200)

#################
# QC
#################
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
stim[["percent.mt"]] <- PercentageFeatureSet(stim, pattern = "^MT-")
stim <- subset(stim, subset = nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 25)
stim <- subset(stim, cells = filtered_crispr_table$cell_barcode)
dim(stim)
stim <- stim[,1:100]
pdf("/home/chengennong/code-manual/vscode/ylp/Science2022/fig/perturb_su/qc.pdf", width=8,height=4)
plot1 <- FeatureScatter(stim, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(stim, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

#################
# Normalize
#################

stim <- SCTransform(stim, assay = "RNA", vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
# Scaling the data
###Remove CellCycleScoring
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
stim <- CellCycleScoring(stim s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
stim <- ScaleData(stim, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), features = rownames(stim))

#################
# Feature selection
#################
# Identification of highly variable features (feature selection)
stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(stim), 10)
# plot variable features with and without labels
pdf("/home/chengennong/code-manual/vscode/ylp/Science2022/fig/perturb_su/high_features.pdf", width=8,height=4)
plot1 <- VariableFeaturePlot(stim)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
dev.off()


#################
# Cut dimension
#################
# Perform linear dimensional reduction
stim <- RunPCA(stim, features = VariableFeatures(object = stim))
print(stim[["pca"]], dims = 1:5, nfeatures = 5)

str(stim)
stim@assays$SCT@scale.data %>% dim
stim@reductions$umap %>% dim
stim@commands$FindClusters
#################
# Clustering
#################
##0.4 for stim 0.5 for resting
stim <- FindNeighbors(stim, dims = 1:20)
stim <- FindClusters(stim, resolution = 0.8, algorithm = 3)
stim <- RunUMAP(stim, dims = 1:20)
stim <- BuildClusterTree(object = stim)

#################
# Find marker
#################
stim.markers <- FindAllMarkers(stim, only.pos = TRUE, logfc.threshold = 0.25)
top50 <- stim.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)




#Save RDS
#saveRDS(stim, file = "/home/chengennong/project/ylp/Science2022/perturb/rds/stim_perturb.rds")
#stim <-readRDS(file = "/home/chengennong/project/ylp/Science2022/perturb/rds/stim_perturb.rds")


#Cluster-G
#pdf("/home/chengennong/code-manual/vscode/ylp/Science2022/fig/perturb_su/cluster.pdf",width=6,height=5)
pdf("/home/chengennong/code-manual/vscode/ylp/Science2022/fig/main/4g.pdf",width=6,height=5)
#ElbowPlot(stim)
#VlnPlot(stim, c("nCount_RNA", "nFeature_RNA", "percent.mt"))
#DimPlot(stim, reduction = "pca")
DimPlot(stim, reduction = "umap", raster =T, label = T)
dev.off()


#CD4/CD8-C
#pdf("/home/chengennong/code-manual/vscode/ylp/Science2022/fig/perturb_su/cluster_cd4cd8.pdf",width=5,height=4)
pdf("/home/chengennong/code-manual/vscode/ylp/Science2022/fig/main/4d.pdf",width=6,height=5)
#stim[["CD4_CD8"]]  <- (log2(stim@assays$SCT@scale.data  %>% as.data.frame%>% filter(rownames(.)=="CD4")+1)-log2(colMeans(stim@assays$SCT@scale.data  %>% as.data.frame%>% filter(rownames(.)=="CD8A"|rownames(.)=="CD8B"))+1)) %>% as.matrix %>% t

ColDark2 <- brewer.pal(5, "Dark2"); red <- ColPaired[6]; blue <- ColPaired[2]; green <- ColPaired[4]; yellow <- "#E69F00"
col=colorRampPalette(colors = c(blue, 'white', green))(10)
color.palette=viridis
ggplot(stim[[]], aes(CD4_CD8)) +
        geom_density(adjust=1.5, alpha=0.9,fill="grey") +
        geom_vline(aes(xintercept = -0.9), colour = "grey",size=2, alpha=1,line_type=2) +
        geom_vline(aes(xintercept = 1.4), colour = "grey",size=2, alpha=1,line_type=2) +
        scale_x_continuous(limits = c(-4,4),breaks=seq(-2.5, 0, 2.5)) +
        theme_bw()+
        theme(legend.position='none',panel.grid =element_blank())+
        labs(title="log2(CD4/CD8)")

FeaturePlot(stim,features = "CD4_CD8",reduction = "umap",pt.size = 1, raster=TRUE,label = T, cols = c("navy","firebrick3"), min.cutoff = -4.3, max.cutoff = 4.3,) +
    scale_colour_gradientn(colors = col) +
    labs(title="log2(CD4/CD8)")
dev.off()

#Feature-I
#library(viridis)
#col=viridis(40)
#pdf("/home/chengennong/code-manual/vscode/ylp/Science2022/fig/perturb_su/features.pdf",width=13,height=7)
pdf("/home/chengennong/code-manual/vscode/ylp/Science2022/fig/main/4i.pdf",width=13,height=7)
#FeaturePlot(stim, features = c("IFNG", "IL2", "LTA", "LTB", "IL13", "CCL3"), raster=TRUE, cols = col) 
FeaturePlot(stim, features = c("IFNG", "IL2", "LTA", "LTB", "IL13", "CCL3"), raster=TRUE, label = T,ncol=3, keep.scale="all") &
    scale_color_viridis_c()
dev.off()


#Heatmap-H
#pdf("/home/chengennong/code-manual/vscode/ylp/Science2022/fig/perturb_su/tree.pdf",width=14,height=8)
top50 <- stim.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
heat_top_gene  <- top50 %>% group_by(gene) %>% top_n(n = 1, wt = avg_log2FC)


heat_cg_table <- stim@assays$SCT@scale.data %>% 
    as.data.frame %>%
    filter(rownames(.) %in% heat_top_gene$gene) %>% 
    rownames_to_column("gene") %>% 
    pivot_longer(-gene, names_to = "cell", values_to = "exp") %>% 
    inner_join((stim[[]] %>% select(seurat_clusters) %>% rownames_to_column("cell"))) %>% #filter(gene=="AAK1") %>% filter(cluster=="3")
    group_by(seurat_clusters,gene) %>% 
    summarise(mean_exp=mean(exp)) 

heat_cg_mtx <- heat_cg_table %>% #group_by(seurat_clusters) %>% arrange(seurat_clusters,-mean_exp)%>% 
    pivot_wider(names_from = "seurat_clusters", values_from = "mean_exp")  %>% 
    column_to_rownames("gene") %>% as.matrix  %>% t %>% base::scale(.) %>% t

gene_order <- heat_top_gene %>% group_by(cluster) %>% arrange(cluster,avg_log2FC) %>% pull(gene)
library(ComplexHeatmap)
coul <- colorRampPalette(brewer.pal(8, "PiYG"))(6)
coul <- brewer.pal(8, "PiYG")
pdf("/home/chengennong/code-manual/vscode/ylp/Science2022/fig/main/4h.pdf",width=7,height=13)
Heatmap(heat_cg_mtx, name = "Z-score", column_title = "Cluster", column_title_side = c("bottom"),
    column_order = order(as.numeric(colnames(heat_cg_mtx))),
    row_order = gene_order,
    show_row_names = F)
dev.off()


#Density-J
sgRNA_den_set <- c("NO-TARGET","MAP4K1","VAV1","TNFRSF1A","IL1R1","FOXQ1","GATA3","TBX21")
sg_set_table <- filtered_crispr_table %>% filter(gene %in% sgRNA_den_set) 
sgRNA_stim <- subset(stim, cells = (sg_set_table %>% pull(cell_barcode)))

umap_coor_gene <- sgRNA_stim@reductions$umap@cell.embeddings %>% as.data.frame %>% rownames_to_column("cell_barcode") %>% inner_join(sg_set_table %>% select(cell_barcode,gene)) %>% 
    rbind(sgRNA_stim@reductions$umap@cell.embeddings %>% as.data.frame %>% rownames_to_column("cell_barcode") %>% mutate(gene="Perturbed cells"))

pdf("/home/chengennong/code-manual/vscode/ylp/Science2022/fig/main/4j.pdf",width=12,height=10)
ggplot(umap_coor_gene , aes(x=UMAP_1, y=UMAP_2) ) +
  geom_density_2d(aes(color = ..level..)) +
  scale_color_viridis_c() +
  facet_wrap(vars(gene),scales="fixed") +
  theme_bw()+
  theme(panel.grid =element_blank()) + 
  #theme(panel.border = element_blank()) + 
  theme(axis.line = element_blank()) +
  theme(legend.position = 'top', 
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 11),
        axis.text = element_text(color = "black",size = 11))+
  theme(
    strip.background = element_blank(),
    #  color = "white", fill = "white"),
    panel.grid = element_blank())
dev.off()

str(sgRNA_stim@reductions$umap@cell.embeddings) %>% head
cell_barcode
pdf("/home/chengennong/code-manual/vscode/ylp/Science2022/fig/main/4j.pdf",width=6,height=5)
FeaturePlot(stim, features = "percent.mt", raster=TRUE, label = T,ncol=1, keep.scale="all") 
str(stim[[]])
dev.off()
##############################################################################################################################

#su
pdf("/home/chengennong/code-manual/vscode/ylp/Science2022/fig/perturb_su/heatmap.pdf",width=13,height=7)
DoHeatmap(stim,top50$gene,size=3)
dev.off()