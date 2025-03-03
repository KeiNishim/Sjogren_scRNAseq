library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(patchwork)
library(BiocManager)
library(metap)
library(data.table)
library(BiocGenerics)
library("org.Hs.eg.db")
library("hgu133plus2.db")
library(jetset)
library(icellnet)
library(gridExtra)
library(devtools)
library(Connectome)
connectome.genes <- union(Connectome::ncomms8866_human$Ligand.ApprovedSymbol,Connectome::ncomms8866_human$Receptor.ApprovedSymbol)
genes <- connectome.genes[connectome.genes %in% rownames(s.sub)]

Idents(s.int) <- "disease"

s.int1 <- subset(s.int, idents=c("SSA"))
s.int2 <- subset(s.int, idents=c("LOW","CEN"))
Idents(s.int) <- "seurat_clusters3"
Idents(s.int1) <- "seurat_clusters3"
Idents(s.int2) <- "seurat_clusters3"
DimPlot(s.int)

s.con <-  CreateConnectome(subset(s.int,downsample=5000, idents="pDC",invert=T), min.z=0.3, species='human', p.values=T, calculate.DOR=T, assay="RNA")
s.con1 <- CreateConnectome(subset(s.int1,downsample=5000, idents="pDC",invert=T), min.z=0.3, species='human', p.values=T, calculate.DOR=T, assay="RNA")
s.con2 <- CreateConnectome(subset(s.int2,downsample=5000), min.z=0.3, species='human', p.values=T, calculate.DOR=T, assay="RNA")

#####
ggplot(s.con, aes(x=percent.target))+geom_density()
ggplot(s.con, aes(x=percent.source))+geom_density()

s2.con1 <- FilterConnectome(s.con1, min.z=0.3, min.pct=0.05, remove.na=T)
s2.con2 <- FilterConnectome(s.con2, min.z=0.3, min.pct=0.05, remove.na=T)

# ????????????????????????????????????????????????????????????????????????
all_cell_types <- unique(c(levels(s.int1$seurat_clusters3), levels(s.int2$seurat_clusters3), levels(s.int$seurat_clusters3)))

# ??????????????????????????????cell_types?????????????????????
colors <- rainbow(length(all_cell_types))
names(colors) <- all_cell_types

s2.con <- FilterConnectome(s.con, min.z=0.3, max.p=1E-150,  min.pct=0.05, remove.na=T)
CircosPlot(s2.con,min.z=0.3, lab.cex=0.8, cols.use=colors, modes.include=c("CC","CXCL", "Interleukin"), small.gap=1.5)


s2.con <- FilterConnectome(s.con, min.z=0.3,max.p=1E-321,   min.pct=0.15, remove.na=T)
CircosPlot(s2.con,min.z=0.66, lab.cex=0.8, cols.use=colors, modes.include=c("UNCAT"), small.gap=1.5)

s2.con1 <- FilterConnectome(s.con1, min.z=0.3, max.p=1E-150,  min.pct=0.05, remove.na=T)
CircosPlot(s2.con1,min.z=0.3, lab.cex=0.8, cols.use=colors, modes.include=c("CC","CXCL", "Interleukin"), small.gap=1.5)

s2.con2 <- FilterConnectome(s.con2, min.z=0.3, max.p=1E-150,  min.pct=0.05, remove.na=T)
CircosPlot(s2.con2,min.z=0.3, lab.cex=0.8, cols.use=colors, modes.include=c("CC","CXCL", "Interleukin"), small.gap=1.5)

s2.con1 <- FilterConnectome(s.con1, min.z=0.3, max.p=1E-150,  min.pct=0.05, remove.na=T, modes.include=c("CC","CXCL", "Interleukin","UNCAT"))
s2.con2 <- FilterConnectome(s.con1, min.z=0.3, max.p=1E-150,  min.pct=0.05, remove.na=T, modes.include=c("CC","CXCL", "Interleukin","UNCAT"))

s.conX <- DifferentialConnectome(s.con2, s.con1, min.pct=0.05)
s2.conX <- FilterConnectome(s.conX, modes.include=c("CC","CXCL", "Interleukin"))

DifferentialScoringPlot(s2.conX, min.score=1, min.pct=0.1)
#DiffEdgeDotPlot(s2.conX, min.score=1, min.pct=0.1)

s2.conX <- FilterConnectome(s.conX, modes.include=c("CC","CXCL", "Interleukin", "UNCAT"))
diff.up.up <- subset(s2.conX,ligand.norm.lfc > 0 & recept.norm.lfc > 0 )
CircosDiff(diff.up.up,min.score = 2,min.pct = 0.1,lab.cex=0.8, cols.use=colors,small.gap=1.5)

NetworkPlot(s.con, weight.attribute='weight_sc',include.all.nodes = F)
Centrality(s.con, modes.include=NULL, min.z=NULL, weight.attribute='weight_sc', group.by='mode')


#####################################################
s.con <-  CreateConnectome(subset(s.sub,downsample=5000), min.z=0.3, species='human', p.values=T, calculate.DOR=T, assay="RNA")
s2.con <- FilterConnectome(s.con, min.z=0.3, max.p=1E-100,  min.pct=0.05, remove.na=T)

# ????????????????????????????????????????????????????????????????????????
all_cell_types <- unique(c(levels(s.sub$seurat_clusters), levels(s.sub$seurat_clusters)))

# ??????????????????????????????cell_types?????????????????????
colors <- rainbow(length(all_cell_types))
names(colors) <- all_cell_types
CircosPlot(s2.con,min.z=0.5, lab.cex=0.8, cols.use=colors, small.gap=1.5, modes.include=c("CC","CXCL", "Interleukins","UNCAT","TNF"))
    
s.sub <- subset(s.int, idents=c("Macrophage_TREM2", "Macrophage_MRC1"))
VlnPlot(subset(s.sub,downsample=785), c("IL1B","TNF","IL10","TGFB1"), pt.size=0, ncol=2)                     

VlnPlot(subset(s.sub,downsample=785), c("LYVE1","RNASE1","CD163","MRC1","TREM2","PLA2G7"), pt.size=0, ncol=2) 
        