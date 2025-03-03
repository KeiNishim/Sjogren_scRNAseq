library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(patchwork)
library(data.table)
library(tidyr)
library(dplyr)
library(magrittr)
library(spatstat)

devtools::install_github(repo = "alikhuseynov/seurat", ref = "feat/vizgen")

vizgen.obj <- ReadVizgen(data.dir = "/Data_here/202403131602_HuSalivaryGland-VS153-SlideB-S3_VMSC10802/region_1/",
                         type="centroids")
vizgen.obj <- LoadVizgen(data.dir = "/Data_here/202403131602_HuSalivaryGland-VS153-SlideB-S3_VMSC10802/region_1/")
                         
file2read <- "/Data_here/202403131602_HuSalivaryGland-VS153-SlideB-S3_VMSC10802/region_0/cell_boundaries.parquet"
test_parq <- sfarrow::st_read_parquet(file2read)
test_parq

vizgen.obj <- SCTransform(vizgen.obj, assay = "Vizgen", clip.range = c(-10, 10))
vizgen.obj <- RunPCA(vizgen.obj, npcs = 30, features = rownames(vizgen.obj))
vizgen.obj <- RunUMAP(vizgen.obj, dims = 1:30)
vizgen.obj <- FindNeighbors(vizgen.obj, reduction = "pca", dims = 1:30)
vizgen.obj <- FindClusters(vizgen.obj, resolution = 0.3)

DimPlot(viz_seurat, group.by="seurat_clusters3")
