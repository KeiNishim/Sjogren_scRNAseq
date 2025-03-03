library(dbscan)
library(scclusteval)
library(ComplexHeatmap)

mat <- viz_seurat@meta.data[,7:8]
mat <- as.matrix(mat)
all.equal(Cells(viz_seurat),rownames(mat))
mat <- mat[Cells(viz_seurat),]

eps <- 50
nn <- frNN(x=mat, eps=eps)
nn$id %>% head(n=3)
nn_df <- nn$id %>% stack()
cluster_ids <- viz_seurat$seurat_clusters1 %>% unname()
nn_df$cluster_id<- cluster_ids[nn_df$values]
nn_df$cluster_id<- factor(nn_df$cluster_id)
nn_count<- nn_df %>%
  group_by(ind) %>%
  count(cluster_id, .drop = FALSE)
nn_count<- nn_count %>%
  tidyr::pivot_wider(names_from = cluster_id, values_from = n)
nn_mat<- nn_count[,-1] %>% as.matrix()
rownames(nn_mat)<- nn_count$ind
set.seed(123)

# I cheated here to make k=13 so I can compare it with the Louvian clustering later
k_means_res<- kmeans(nn_mat, centers = 10)
k_means_id<- k_means_res$cluster %>%
  tibble::enframe(name = "cell_id", value = "kmeans_cluster")
head(k_means_id)
# add it to the metadata solot for nn_obj below
k_means_df<- as.data.frame(k_means_id)
rownames(k_means_df)<- k_means_id$cell_id
nn_obj<- CreateSeuratObject(counts = t(nn_mat),  min.features = 5)
nn_obj<- SCTransform(nn_obj, vst.flavor = "v2")
nn_obj <- RunPCA(nn_obj, npcs = 15, features = rownames(nn_obj))
ElbowPlot(nn_obj)
nn_obj <- FindNeighbors(nn_obj, reduction = "pca", dims = 1:11)
nn_obj <- RunUMAP(nn_obj, dims = 1:11)
nn_obj <- FindClusters(nn_obj, resolution = 0.1)
DimPlot(viz_seurat, group.by="SCT_snn_res.0.1")


nn_obj$SCT_snn_res.0.1 <- Idents(nn_obj)

nn_obj<- AddMetaData(nn_obj, metadata = k_means_df)
nn_obj@meta.data %>%
  head()
table(nn_obj$kmeans_cluster, nn_obj$seurat_clusters)

PairWiseJaccardSetsHeatmap(Idents(nn_obj), k_means_res$cluster, best_match = TRUE)
nn_meta<- nn_obj@meta.data %>%
  select(cell_id, SCT_snn_res.0.1, kmeans_cluster)
viz_seurat <- AddMetaData(viz_seurat, nn_meta)

# calculate the average abundance of cell types per niche
average_cell_type_abundance<- AverageExpression(
  nn_obj,
  assays = NULL,
  features = rownames(nn_obj),
  return.seurat = FALSE,
  group.by = "SCT_snn_res.0.1",
)

scCustomize::Stacked_VlnPlot(nn_obj, features= rownames(nn_obj),x_lab_rotate = 90)+ coord_flip()
  theme(axis.text.x = element_text(size = 1))
cell_fun = function(j, i, x, y, width, height, fill) {
  grid::grid.rect(x = x, y = y, width = width *0.99, 
                  height = height *0.99,
                  gp = grid::gpar(col = "grey", 
                                  fill = fill, lty = 1, lwd = 0.5))
}

col_fun=circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

Heatmap(t(scale(t(as.matrix(average_cell_type_abundance$SCT)))),
        show_row_dend =FALSE,
        show_column_dend = FALSE,
        rect_gp = grid::gpar(type = "none"),
        cell_fun = cell_fun,
        col = col_fun,
        column_names_rot = 45)

mat2<- table(viz_seurat$SCT_snn_res.0.1, viz_seurat$seurat_clusters) 
mat2
mat2 <- mat2[,-10]
Heatmap(t(scale(as.matrix(mat2))),
        cluster_columns=T,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        rect_gp = grid::gpar(type = "none"),
        cell_fun = cell_fun)

