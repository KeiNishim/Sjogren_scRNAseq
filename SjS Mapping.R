library(Seurat)
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

DimPlot(s.int)
FeaturePlot(s.int, "CXCL12")

FeaturePlot(s.int, "CR1", split.by="disease")
VlnPlot(s.sub, "IL32")

s.sub2$sample <- as.character(s.sub2$sample)
s.sub2$clusters <- as.character(s.sub2$clusters)
library(stats)
cell_counts <- table(s.sub2$sample, s.sub2$clusters)
proportion_data <- as.data.frame.matrix(cell_counts)
proportion_data$sample <- rownames(proportion_data)
proportion_data <- proportion_data %>%
  mutate(total = rowSums(select(., -sample))) %>%
  mutate(across(-c(sample, total), ~./total)) %>%
  select(-total) %>%
  pivot_longer(-sample, names_to = "cluster", values_to = "proportion")

sample_disease <- s.sub2@meta.data %>%
  select(sample, disease) %>%
  distinct()

proportion_data <- proportion_data %>%
  left_join(sample_disease, by = "sample")
wilcox_results <- proportion_data %>%
  group_by(cluster) %>%
  summarise(
    p_value = wilcox.test(proportion ~ disease)$p.value,
    .groups = 'drop'
  )
wilcox_results$p_adj <- p.adjust(wilcox_results$p_value, method = "fdr")
print(wilcox_results, n=40)

