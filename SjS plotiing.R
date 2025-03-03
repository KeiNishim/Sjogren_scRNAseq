library(cowplot)
library(dplyr)
library(ggplot2)
library(patchwork)
library(BiocManager)
library(data.table)
library(tidyr)
library(monocle3)
library(SeuratWrappers)
library(data.table)
library(png)
img <- readPNG("Bechedt/Figure/Fig.4d.png")
plot(1:2, type='n', xlab="", ylab="", xlim=c(0, 1), ylim=c(0, 1))
rasterImage(img, 0, 0, 1, 1)

all <- FindAllMarkers(subset(s.int,downsample=20000),only.pos = TRUE,
                          logfc.threshold = 0.25, min.pct=0.25, max.cells.per.ident=5000)

write.csv(all, "Genemarkers.csv")

top10 <- all %>% group_by(cluster) %>% top_n(n = -4, wt = p_val_adj)
top10 <- top10 %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
X <- unique(top10['gene'])
Y <- apply(X,2,rev)

pp <- DotPlot(s.int,features=Y,scale.by="size", scale=T,col.max=5, col.min=-5)+
   theme(axis.text.x = element_text(angle=90,size = 8,vjust=0.5,hjust=1))+ 
  theme(axis.text.y = element_text(size = 7.5))+scale_size(range=c(0.5,3))+
  scale_colour_gradient2(low="blue", mid="white", high="red", midpoint=0)+coord_flip()
pp

DotPlot(s.int,features=c("CXCR5","CXCL13"),scale.by="size",group.by="seurat_clusters", scale=T,col.max=5, col.min=-5)+
  theme(axis.text.x = element_text(angle=90,size = 8,vjust=0.5,hjust=1))+ 
  theme(axis.text.y = element_text(size = 7.5))+scale_size(range=c(0.5,3))+
  scale_colour_gradient2(low="blue", mid="white", high="red", midpoint=0)+coord_flip()

DotPlot(s.int,features=c("TGFB1"),scale.by="size",group.by="seurat_clusters",  scale=T,col.max=5, col.min=-5)+
  theme(axis.text.x = element_text(angle=90,size = 8,vjust=0.5,hjust=1))+ 
  theme(axis.text.y = element_text(size = 7.5))+scale_size(range=c(0.5,3))+
  scale_colour_gradient2(low="blue", mid="white", high="red", midpoint=0)+coord_flip()


DotPlot(s.int,features=c("CXCR5","CXCL13"),scale.by="size", group.by="sample",scale=T,col.max=5, col.min=-5)+
  theme(axis.text.x = element_text(angle=90,size = 8,vjust=0.5,hjust=1))+ 
  theme(axis.text.y = element_text(size = 7.5))+scale_size(range=c(0.5,3))+
  scale_colour_gradient2(low="blue", mid="white", high="red", midpoint=0)+coord_flip()

DotPlot(s.int,features=c("TGFB1"),scale.by="size", group.by="sample",scale=T,col.max=5, col.min=-5)+
  theme(axis.text.x = element_text(angle=90,size = 8,vjust=0.5,hjust=1))+ 
  theme(axis.text.y = element_text(size = 7.5))+scale_size(range=c(0.5,3))+
  scale_colour_gradient2(low="blue", mid="white", high="red", midpoint=0)+coord_flip()

VlnPlot(s.int, "CXCL12",group.by="seurat_clusters", pt.size=0)
Idents(s.int) <- "seurat_clusters"

#For integrated anlysis for fibroblast
s.int$combined_group <- paste(s.int$group...5, s.int$group...7, sep = "_")
Idents(s.int) <- "combined_group"
table(s.int$combined_group)
s.int <- RenameIdents(s.int, "Sjogren syndrome_minor salivary gland"="Sjogren syndrome_minor salivary gland","dry eye syndrome_minor salivary gland"="dry eye syndrome_minor salivary gland",
                      "idiopathic pulmonary fibrosis_lung"="idiopathic pulmonary fibrosis_lung",
                      "interstitial lung disease_lung"="interstitial lung disease_lung","non-specific interstitial pneumonia_lung"="non-specific interstitial pneumonia_lung",
                      "rheumatoid lung disease_lung"="rheumatoid lung disease_lung",
                      "normal_lung"="normal_lung","rheumatoid lung disease_lung"="rheumatoid lung disease_lung","osteoarthritis_synovial membrane of synovial joint"="osteoarthritis_synovial membrane of synovial joint",
                      "rheumatoid arthritis_synovial membrane of synovial joint"="rheumatoid arthritis_synovial membrane of synovial joint",
                      "ulcerative colitis_large intestine"="ulcerative colitis_large intestine","normal_large intestine"="normal_large intestine")
s.int$combined_group <- Idents(s.int)
Idents(s.int) <- "combined_group"
DimPlot(subset(s.int,downsample=5000), group.by="combined_group", shuffle=T)

Idents(s.int) <- "group...7"
DimPlot(subset(s.int,downsample=7500), group.by="group...7", shuffle=T)

VlnPlot(s.int, "CXCL12",group.by="combined_group", pt.size=0)
