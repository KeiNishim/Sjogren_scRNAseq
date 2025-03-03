library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Seurat)

#Read mono_ind as s.sub
s.sub <- DietSeurat(s.int, assays="integrated", dimreducs = c("pca","umap"), graphs=c("integrated_snn", "integrated_nn"))
pbmc_small_sce <- as.SingleCellExperiment(s.sub)
traj_milo <- Milo(pbmc_small_sce)
plotUMAP(traj_milo, colour_by="sample", point_size=0.1)
traj_milo <- buildGraph(traj_milo, k = 10, d = 30, reduced.dim="PCA")
traj_milo <- makeNhoods(traj_milo, prop = 0.05, k = 10, d=30,reduced_dim="PCA", refined = TRUE)
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), sample="sample")
traj_milo <- calcNhoodDistance(traj_milo, d=30, reduced.dim="PCA")
saveRDS(traj_milo, "sjs_milo.rds")

traj_design <- data.frame(colData(traj_milo))[,c("sample", "disease")]
traj_design <- distinct(traj_design)
age <- data.frame(c(38,62,71,34,22,55,51,50,81,35,53,54))
#SLE,HD class : 
names(age) <- "age"
traj_design <- cbind(traj_design, age)
rownames(traj_design) <- traj_design$sample
traj_design$class <- c("B","C","B","C","C","B","C","C","B","C","C","C")

da_results <- testNhoods(traj_milo,reduced.dim="PCA", design = ~ class, design.df = traj_design)
da_results %>%
  arrange(- SpatialFDR) %>%
  head()
traj_milo <- buildNhoodGraph(traj_milo)
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
ggplot(da_results, aes(logFC, -log10(SpatialFDR)))+ 
  geom_point() +
  geom_hline(yintercept = 1)

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(traj_milo, da_results,size_range=c(0.3,3), layout="UMAP",alpha=1)+
  scale_fill_gradient2(low='blue', mid='white', high="red", midpoint = -2)
nh_graph_pl
ggsave('Fig.2C.png', nh_graph_pl +plot_layout(guides="collect"), width=7.5, height=7.5)

da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "seurat_clusters")
da_results$seurat_clusters <- factor(da_results$seurat_clusters,
                                     levels= rev(c("CD4 TCM","CD4 TRM","TFH","Treg","CD8 CTL","CD8 TRM","B cell","PC","CD69 LLPC",
                                                   "FOS LLPC","ENAM LLPC","MARCO Macrophage","LYVE1 Macrophage","CX3CR1 Macrophage",
                                                   "Monocyte","Neutrophil","Mast","Fibrotic Fibroblast","Fibrotic Fibroblast Str",
                                                   "Inflammatory Fibroblast","Inflammatory Fibroblast Str","FAP Fibroblast",
                                                   "Serous Acinar","Mucous Acinar","Duct","Myoepithelia")))
pDa <- plotDAbeeswarm(da_results, group.by = "seurat_clusters", alpha=1)+geom_boxplot(outlier.shape = NA)+
  scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" )+
  geom_hline(yintercept=0, linetype="dashed")
pDa
ggsave('Fig.2D.png', pDa, width=7.5, height=7.5)

median(da_results[da_results$seurat_clusters=="CD14 Mono_Activated", "logFC"])
median(da_results[da_results$seurat_clusters=="CD14 Mono_ISG", "logFC"])
median(da_results[da_results$seurat_clusters=="CD14 Mono_VCAN", "logFC"])
median(da_results[da_results$seurat_clusters=="CD14 Mono_HLA", "logFC"])
median(da_results[da_results$seurat_clusters=="CD16 Mono", "logFC"])
median(da_results[da_results$seurat_clusters=="cDC", "logFC"])
