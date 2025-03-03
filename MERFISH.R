```{r setup, include=FALSE}

knitr::opts_chunk$set(echo = TRUE)
working_directory = './' # change to the directory your data is in

# Setting the working directory for the rest of the file see 2
knitr::opts_knit$set(root.dir = working_directory)
```

# Set up the environment

```{r install_packages, message=FALSE, warning = FALSE, results = 'hide',eval = FALSE}

# Package that has many helpful tools for plotting and data manipulation
install.packages("tidyverse", version = "1.3.0", repos='http://cran.us.r-project.org')
# Has features useful for checking data
install.packages("useful", version = "1.2.6", repos='http://cran.us.r-project.org')
# Main package to work with single cell data
install.packages("Seurat", version = "4.3.0", repos='http://cran.us.r-project.org')
# Package with color palletes for custom visualizations
install.packages("viridis", version = "0.6.3", repos='http://cran.us.r-project.org')
# Package that handle categorical variables
install.packages("forcats", version = "1.0.0", repos='http://cran.us.r-project.org')
# fast reading of large csv files
install.packages('data.table', version = "1.14.8", repos='http://cran.us.r-project.org')
```

```{r load libraries, message=FALSE}

library(tidyverse)
library(Seurat)
library(useful)
library(viridis)
library(forcats)

```


# Read the cell-by-gene.csv and cell-metadata.csv and create a Seurat object

To set up the working directory (where the files are stored) follow this [link](http://www.sthda.com/english/wiki/running-rstudio-and-setting-up-your-working-directory-easy-r-programming)
Section: 'Change your working directory'


```{r set data directories}
cell_by_gene_csv = 'SlideC_1/cell_by_gene.csv'
meta_cell_csv = 'SlideC_1/cell_metadata.csv'
```


```{r read cell-by-gene and cell-metadata tables}

#read cell by gene table
#file_dir = "./"
#files = list.files(file_dir, pattern = ".csv$")

cell_gene_table <- data.table::fread(cell_by_gene_csv,
                                     data.table = F) # returns a data.frame instead

cell_names <- cell_gene_table[,1] #get cell IDs

##set cell IDs as rownames
row.names(cell_gene_table) <- cell_names
cell_gene_table <- cell_gene_table[, -1] # Remove cell IDs

##remove "Blanks"
cell_gene_table <- cell_gene_table[,!grepl('Blank.',
                                           colnames(cell_gene_table))]

# Transposing the table per convention
cell_gene_table <- t(cell_gene_table)

# Looking at the dimensions of the table
dim(cell_gene_table)
print(paste0('We have ',dim(cell_gene_table)[1], ' features (genes) in ',
             dim(cell_gene_table)[2], ' samples (cells).'))

corner(cell_gene_table) #take a look

#read metadata
cell_meta <- read.csv(meta_cell_csv, 
                      header = T, 
                      row.names = as.character(cell_names), #directly reuse cell IDs
                      stringsAsFactors = F)


dim(cell_meta)
corner(cell_meta) #take a look

```

```{r create a Seurat object}

viz_seurat <- CreateSeuratObject(counts = cell_gene_table, 
                                 meta.data = cell_meta, 
                                 assay = "Vizgen")
#Rename Identity!!!!!!!
viz_seurat <- AddMetaData(viz_seurat, "C1", col.name="sample") 

#Add centroid fields
centroids <- data.frame(x = cell_meta$center_x,
                        y = cell_meta$center_y,
                        cell = cell_names,
                        row.names = as.character(cell_names),
                        stringsAsFactors = F)

cents <- CreateCentroids(centroids)

coords <- CreateFOV(coords = cents, type = "centroids", assay = "Vizgen")

viz_seurat[["viz_fov"]] <- coords
viz_seurat[["viz_fov"]]

viz_seurat

```


## Downsize the object to save time

**For a real analysis, skip this code block.** ↓ ↓ ↓ ↓ ↓ ↓ 
This block downsize the object to save time for tuning parameters.
```{r downsize} 

#set.seed(123)
#viz_seurat <- viz_seurat[, sample(colnames(viz_seurat), size =length(colnames(viz_seurat))/10, replace=F)]
#viz_seurat

```
**For the complete experiment analysis, skip downsizing and resume from next code chunk.** ↓ ↓ ↓ ↓ ↓ ↓
\
\

# Single-ell Analysis

## Preprocessing
### Filter cells based on minimum gene expression counts and volume
Set up total gene counts and volume thresholds to filter cells that don't pass the criteria:

**barcode_count_threshold** A minimal barcode_count_threshold value removes low quality cells and poorly segmented cells.

**volume_upper_threshold** and **volume_lower_threshold** The unit for cell volume is expressed in microns. Users can utilize the Vizualizer to estimate the volume by multiplying the area and thickness.\
Please note that the default Cellpose segmentation algorithm has set the minimum cell volume to 100. Thus, any cells smaller than 100 micron^3 was not included in the cell-by-gene table.

```{r, warning=FALSE}

print(paste0("Cell number before filtering = ", length(rownames(viz_seurat@meta.data))))

barcode_count_threshold <- 5
volume_upper_threshold <- 2500
volume_lower_threshold <- 100

# Plotting a histogram of our volume values with the cutoffs visualized
ggplot(viz_seurat@meta.data, aes(x = volume)) +
  geom_histogram(bin = 30) + # plotting a histogram of the volume
  theme_bw() + # a basic black and white them for the plot
  scale_x_log10() + # log scaling the x-axis because volume has a long right tail
  geom_vline(xintercept = c(volume_lower_threshold,volume_upper_threshold),
             colour = 'red',
             linetype = 2) # plotting the cutoffs we're using

# Similar plot but for barcode counts
ggplot(viz_seurat@meta.data, aes(x = nCount_Vizgen)) +
  geom_histogram(bins = 30) +
  theme_bw() +
  scale_x_log10() +
  geom_vline(xintercept = barcode_count_threshold,
             colour = 'red',
             linetype = 2)

```


```{r}
viz_seurat <- subset(viz_seurat, subset = nCount_Vizgen > barcode_count_threshold & 
                       volume < volume_upper_threshold & 
                       volume > volume_lower_threshold)
print(paste0("Cell number after filtering = ", length(rownames(viz_seurat@meta.data))))

```

### Normalize barcode counts by cell volume

To account for variations in cell volume due to different cutting angles and areas, we recommend normalizing transcript counts to the median cell volume. This normalization helps address the tendency of larger volume cells to have more transcripts detected.

```{r normalize to volume}
# Calculate volume factor
volume_factor <- viz_seurat@meta.data$volume / median(viz_seurat@meta.data$volume)

# Divide expression matrix by volume factor
viz_seurat[["Vizgen"]]@layers[["counts"]] <- viz_seurat[["Vizgen"]]@layers[["counts"]] / volume_factor
saveRDS(viz_seurat, file = "../viz_C1.rds")

```


### Normalize barcode counts by total counts over all genes, log transform and scale gene counts
**NormalizeData()** Use LogNormalize method to globally normalizes the feature expression measurements for each cell by the total expression.

**ScaleData()** Adjusts the gene expression values to have a mean of 0 and a variance of 1 so to remove cell-specific effects.

```{r normalization, message=FALSE}

viz_seurat <- NormalizeData(object = viz_seurat,
                            normalization.method = "LogNormalize")
all.genes <- rownames(viz_seurat)
viz_seurat <- ScaleData(viz_seurat,
                        features = all.genes)

```

## Dimension reduction by Principal component analysis (PCA)
PCA reduces the high-dimensional gene expression data to a lower-dimensional space while retaining the most informative variation.

In the context of analyzing scRNA-seq data, where the entire transcriptome is measured, it is common to pre-filter the genes based on their variability using the **FindVariableFeatures()** function. This ensures that only the highly variable genes are included in the PCA analysis.

In Seurat's PBMC [vignette](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html), the **RunPCA()** function is used with the **features = VariableFeatures()** argument to specify that only the highly variable genes should be considered for PCA. This helps to focus on the genes that exhibit the most variation across cells.

However, in the case of the MERFISH panel, all genes in the panel are selected as "highly variable" genes during the experimental design. Therefore, all genes can be used for PCA analysis, and the **features = all.genes** argument is set accordingly. This ensures that all genes are included in the PCA analysis, as they are considered informative for capturing the variation within the MERFISH dataset.

```{r PCA}

viz_seurat<- RunPCA(viz_seurat,
                    npcs = 30,
                    features = all.genes,
                    verbose = F) # not including text output from the function
ElbowPlot(viz_seurat)

```

PCs are arranged based on their contribution to the variance of the data. The number of PCs selected for subsequent analysis should be determined by identifying the point where the 'elbow' occurs, indicating a significant decrease in variance.

It's important to consider a trade-off when choosing the number of PCs. Including more PCs will capture a greater amount of variation in the data. However, late PCs tend to contribute less to the overall variance and can slow down the subsequent analysis process.


## Connect neighbors and embed on UMAP
**FindNeighbors()** Identify the nearest neighbors of each cell in a reduced-dimensional space (PCA here). It helps to construct a shared nearest neighbor graph that captures the local relationships between cells.

**FindClusters()** Perform clustering analysis on the cells using the previously identified neighbors. It assigns cells to clusters based on their similarity in the reduced-dimensional space.

- **resolution:** Control the granularity of the clustering. Higher values result in more and smaller clusters, while lower values lead to fewer and larger clusters.\

**RunUMAP()** Perform a nonlinear dimensionality reduction using UMAP algorithm. We suggest using default min_dist and spread parameter for balanced dispersal of points.

- **dims:** The dimensions (principal components) to be used for the UMAP projection. In this case, 1:30 indicates that the first 30 principal components (obtained from a previous PCA step) will be used for the UMAP calculation
```{r cluster and UMPA, message=FALSE, warning=FALSE}

viz_seurat<- FindNeighbors(viz_seurat,
                           reduction = "pca",
                           dims = 1:19)
viz_seurat<- FindClusters(viz_seurat,
                          resolution =4,graph.name="Vizgen_snn",
                          verbose = T)
viz_seurat<- RunUMAP(viz_seurat,min.dist=0.2,
                     dims = 1:19,
                     verbose = T)

```

# Visualize clusters on UMAP and spatial distribution
```{r UMAP projection}

DimPlot(viz_seurat, reduction = "umap", shuffle=T)

```

```{r spatial projection}

ImageDimPlot(viz_seurat, fov = "viz_fov", cols = "polychrome", axes = TRUE) + 
  theme(panel.grid = element_blank()) # Remove grid lines

```

## Visualize specific cluster:

```{r}

ImageDimPlot(viz_seurat, fov = "viz_fov",
             cols = "Yellow",
             cells = WhichCells(viz_seurat, idents = 0),
             axes = TRUE) + 
  theme(panel.grid = element_blank()) # Remove grid lines

```


# Annotate clusters based on differentially expressed genes

## Find cluster markers
**FindAllMarkers()** Calculate statistical significance and log fold change for each gene and cluster comparison using methods such as Wilcoxon rank sum test (default) or likelihood ratio test.

- **only.pos:** A logical parameter (default is TRUE) that specifies whether to only consider markers that are positive enriched. Setting only.pos = TRUE will exclude markers that are not expressed in other cluster.

- **min.pct:** The minimum percentage of cells within a cluster that must express a gene for it to be considered a marker. By default, min.pct = 0.25 means that a gene must be expressed in at least 25% of the cells in a cluster to be considered a marker.
- **logfc.threshold:** The threshold for log fold change to consider a gene as a marker. By default, logfc.threshold = 0.25 means that a gene must have a log fold change of at least 0.25 (or -0.25 if downregulated) to be considered a marker.

```{r, message=FALSE}

viz.markers <- FindAllMarkers(viz_seurat,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)
viz.markers %>%
  group_by(cluster) %>%
  slice_max(n = 4, order_by = avg_log2FC)

```

## Plot ranking of genes using dotplot
To gain insights into the genes that are most significantly associated with each cell group and visualize the fraction of cells within each group that express those genes.

In **slice_max**, you can set **n = ** to adjust the number of top-ranked genes to include in the dotplot.
```{r, message=FALSE}

markers <- viz.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC) %>% 
  pull(gene) %>% 
  unique()

DotPlot(viz_seurat, features = markers) +
  geom_point(aes(size=pct.exp),
             shape = 21,
             colour="black",
             stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21,
                                             colour="black",
                                             fill="white"))) +
  ylab('Clusters') +
  xlab('Marker Genes') +
  coord_flip() + # flipping axes for easier reading of gene names
  theme(axis.text.x = element_text(size = 8, angle=90),   # Set font size for x-axis
        axis.text.y = element_text(size = 8))   # Set font size for y-axis

```


## Assign cell type identity based on the above gene enrichment analysis
- Before running the following chunk, follow the format **seurat_clusters == number ~ "name"**, to annotate individual clusters.
This syntax will ensure correct assignment.
- Two or more clusters can be assigned the same name.
- Any un-annotated clusters will be assigned as 'unknown cell type' by default
- For example, the hashed seurat_clusters == 7 is not annotated and their cluster name will be assigned as 'unknown cell type' in the following steps.
```{r}

viz_seurat@meta.data <- viz_seurat@meta.data %>%
  mutate(manual_annot = "unknown cell type") %>%
  mutate(manual_annot = case_when(
    seurat_clusters == 0 ~ "Mucous",
    seurat_clusters == 1 ~ "LLPC",
    seurat_clusters == 2 ~ "Serous",
    seurat_clusters == 3 ~ "EC",
    seurat_clusters == 4 ~ "Macrophage_MRC1",
    seurat_clusters == 5 ~ "Mucous",
    seurat_clusters == 6 ~ "CD8T",
    seurat_clusters == 7 ~ "CD4T",
    seurat_clusters == 8 ~ "Serous",
    seurat_clusters == 9 ~ "Serous",
    seurat_clusters == 10 ~ "Mast",
    seurat_clusters == 11 ~ "Fibro",
    seurat_clusters == 12 ~ "Duct",
    seurat_clusters == 13 ~ "Serous",
    seurat_clusters == 14 ~ "CD4T",
    seurat_clusters == 15 ~ "PC",
    seurat_clusters == 16 ~ "Acinar_Prol",
    seurat_clusters == 17 ~ "EC_Act",
    seurat_clusters == 18 ~ "Fibro_Act",
    seurat_clusters == 19 ~ "Myoepithelial",
    seurat_clusters == 20 ~ "Fibro_Act",
    seurat_clusters == 21 ~ "Serous",
    seurat_clusters == 22 ~ "Myoepithelial",
    seurat_clusters == 23 ~ "NK",
    seurat_clusters == 24 ~ "Fibro",
    seurat_clusters == 25 ~ "Serous",
    seurat_clusters == 26 ~ "Fibro_Act",
    seurat_clusters == 27 ~ "Fibro",
    seurat_clusters == 28 ~ "EC",
    seurat_clusters == 29 ~ "EC",
    seurat_clusters == 30 ~ "Fibro",
    seurat_clusters == 31 ~ "Mast_Act",
    seurat_clusters == 32 ~ "LLPC",
    seurat_clusters == 33 ~ "Unknown",
    seurat_clusters == 34 ~ "Macrophage_CX3CR1",
    seurat_clusters == 35 ~ "Serous",
    seurat_clusters == 36 ~ "Fibro_Prol",
    seurat_clusters == 37 ~ "Unknown",
    seurat_clusters == 38 ~ "Unknown",
    seurat_clusters == 39 ~ "Duct_Prol",
    seurat_clusters == 40 ~ "Fibro",
    seurat_clusters == 41 ~ "Unknown",
    seurat_clusters == 42 ~ "Unknown",
    seurat_clusters == 43 ~ "Unknown",
    TRUE ~ manual_annot
  ))

```



## Project annotated clusters on UMAP and spatial plots
```{r project annotated clusters on UMAP}

DimPlot(viz_seurat, group.by = "manual_annot")

```

```{r project annotated clusters on spatial map}

ImageDimPlot(viz_seurat, fov = "viz_fov", group.by = "manual_annot", axes = TRUE) + 
  theme(panel.grid = element_blank()) # Remove grid lines

```

# Export the cluster and UMAP information for further exploring on the MERSCOPE Vizualizer
Under `select(viz_seurat@meta.data, "manual_annot")`, here we set the parameter to **"manual_annot"** to export the annotated cluster names.\ 
You can also set **"seurat_clusters"** to export the original cluster information.

```{r export for the Vizualizer}

cell_cat <- merge(select(viz_seurat@meta.data,  NULL),
                  select(viz_seurat@meta.data, "seurat_clusters1"),
                  by = "row.names",
                  all.x = TRUE)
#This step might take a while

## rename group NA as "unclustered (cells don't pass threshold)"
# Check if column exists and is not NULL
if (!is.null(cell_cat$seurat_clusters1) && "seurat_clusters1" %in% names(cell_cat)) {
  cell_cat$seurat_clusters1 <- fct_na_value_to_level(cell_cat$seurat_clusters1,
                                                     level = "unclustered (cells don't pass threshold)")
} else {
  # Handle the case when column doesn't exist or is NULL
  # You can choose to display a warning or perform other actions
  warning("Column 'manual_annot' doesn't exist or is NULL.")
}
write.csv(cell_cat, file = "seurat_clusterfix.csv", row.names = FALSE, quote = FALSE)
head(cell_cat)

cell_cat <- merge(select(viz_seurat@meta.data,  NULL),
                  select(viz_seurat@meta.data, "seurat_niche1"),
                  by = "row.names",
                  all.x = TRUE)
#This step might take a while

## rename group NA as "unclustered (cells don't pass threshold)"
# Check if column exists and is not NULL
if (!is.null(cell_cat$seurat_niche1) && "seurat_niche1" %in% names(cell_cat)) {
  cell_cat$seurat_niche1 <- fct_na_value_to_level(cell_cat$seurat_niche1,
                                                  level = "unclustered (cells don't pass threshold)")
} else {
  # Handle the case when column doesn't exist or is NULL
  # You can choose to display a warning or perform other actions
  warning("Column 'manual_annot' doesn't exist or is NULL.")
}
write.csv(cell_cat, file = "seurat_nichefix.csv", row.names = FALSE, quote = FALSE)
head(cell_cat)


umap_df <- viz_seurat@reductions[["umap"]]@cell.embeddings
umap_df <- merge(select(viz_seurat@meta.data, NULL), umap_df, by = "row.names", all.x = TRUE)
names(umap_df) <- c("EntityID", "umap_1", "umap_2")
write.csv(umap_df, file = "manual_cluster_umap_coordnates.csv", row.names = FALSE, quote = FALSE)
head(umap_df)

```


# Save the annotated Seurat object

```{r save MERFISH object}

saveRDS(viz_seurat, file = "../viz_all.rds")
```

# Information about the current R session

```{r session information}

sessionInfo()

```