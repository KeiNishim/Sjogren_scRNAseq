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
library(circlize)

Sys.setenv("http_proxy" = "http://170547:Kn01202584@www-proxy.chugai-pharm.co.jp:80")
Sys.setenv("https_proxy" = "http://170547:Kn01202584@www-proxy.chugai-pharm.co.jp:80")

library(nichenetr)
library(tidyverse)

Idents(s.int) <- "seurat_clusters3"
s.int <- RenameIdents(s.int,
                      "CD4TRM"="CD4", "Treg"="CD4","Serous"="Serous Acinar","Mucous"="Mucous Acinar",
                      "TFH"="CD4","Fibroblast"="Fibroblast",
                      "Fibroblast_EGR1"="Fibroblast","IGG_PC"="PC",
                      "IGA_PC"="PC","PC_dim"="PC","TFH"="TFH",
                      "CD8CTL"="CD8","CD8TRM"="CD8",
                      "Macrophage_MARCO"="Macrophage","Monocyte"="Macrophage",
                      "B cell"="Bcell",  "CD4Naive"="CD4", "CD8Naive"="CD8","Lymp_Prolif"="CD4"    #CD4 naive not??
                      )
DimPlot(s.int)
s.int$seurat_clusters3 <- Idents(s.int)
DimPlot(s.int, group.by="disease")
DefaultAssay(s.int) <- "SCT"
DefaultAssay(s.sub) <- "SCT"
Ligand =c("CD4","CD8","PC","Bcell", "PC","Duct","Macrophage","Mast", "cDC","Neutrophil","Serous Acinar","Mucous Acinar", "Fibroblast")

Idents(s.int) <- "seurat_clusters3"
s.int <- subset(s.int, idents=Ligand)
s.sub <- subset(s.int, ident="Bcell")
s.sub <- subset(s.int, ident="CD4")

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))

Receptor = c("SSA","CEN")
Idents(s.int) <- "disease"
Idents(s.sub) <- "disease"
#Defaultt Assay is RNA

expressed_genes_sender = rownames(FindMarkers(s.int, ident.1=Receptor,
                                              only.pos =TRUE, min.pct=0.01,logfc.threshold = 0.2, max.cells.per.ident=2000))
expressed_genes_receiver = rownames(FindMarkers(s.sub, ident.1=Receptor,only.pos = FALSE, min.pct=0.01,logfc.threshold = -Inf, max.cells.per.ident=2000)) #0.001 for Bcell, 0.01 for CD4
interested_genes_receiver <- rownames(top_n(
  FindMarkers(s.sub, ident.1=Receptor,logfc.threshold=0.1, min.pct=0.01, only.pos=T),n = -450, wt = p_val)#n=350 for CD4, 450 for Bcell
  )
geneset_oi = interested_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)] 
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
ligands = lr_network %>% pull(from) %>% unique()
ligands <-intersect(colnames(ligand_target_matrix),ligands)#
expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_ligands
receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)
lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities %>% arrange(-pearson) 
best_upstream_ligands = ligand_activities %>% top_n(15, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
best_upstream_ligands

library(dplyr)#not use raster or select
vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(auroc) %>% arrange(auroc) %>% as.matrix(ncol = 1)


(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank()))  

# show histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="red", linetype="dashed", linewidth=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
active_ligand_target_links_df <- na.omit(active_ligand_target_links_df)
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix)
nrow(active_ligand_target_links_df)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
order_targets <- intersect(order_targets, rownames(active_ligand_target_links))
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

#Circus output
Idents(s.int) <- "seurat_clusters3"
avg_expression_ligands = AverageExpression(s.int, features = best_upstream_ligands, assay="RNA")
sender_ligand_assignment = avg_expression_ligands$RNA %>% apply(1, function(ligand_expression){
  ligand_expression > (ligand_expression %>% mean() + ligand_expression %>% sd()*2)
}) %>% t()
sender_ligand_assignment = sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x) > 0})
names(sender_ligand_assignment)
all_assigned_ligands = sender_ligand_assignment %>% lapply(function(x){names(x)}) %>% unlist()
unique_ligands = all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
general_ligands = best_upstream_ligands %>% setdiff(unique_ligands)

CD4_specific_ligands = sender_ligand_assignment$CD4 %>% names() %>% setdiff(general_ligands)
CD8_specific_ligands = sender_ligand_assignment$CD8 %>% names() %>% setdiff(general_ligands)
Bcell_specific_ligands = sender_ligand_assignment$Bcell %>% names() %>% setdiff(general_ligands)
PC_specific_ligands = sender_ligand_assignment$PC %>% names() %>% setdiff(general_ligands)
Neutrophil_specific_ligands = sender_ligand_assignment$Neutrophil %>% names() %>% setdiff(general_ligands)
Duct_specific_ligands = sender_ligand_assignment$Duct %>% names() %>% setdiff(general_ligands)
Mast_specific_ligands = sender_ligand_assignment$Mast %>% names() %>% setdiff(general_ligands)
Serous_specific_ligands = sender_ligand_assignment$'Serous Acinar' %>% names() %>% setdiff(general_ligands)

ligand_type_indication_df = tibble(
  ligand_type = c(
                  rep("CD4-specific", times = CD4_specific_ligands %>% length()),
                  rep("CD8-specific", times = CD8_specific_ligands %>% length()),
                  rep("Bcell-specific", times = Bcell_specific_ligands %>% length()),
                  rep("PC-specific", times = PC_specific_ligands %>% length()),
                  rep("Neutrophil-specific", times = Neutrophil_specific_ligands %>% length()),
                  rep("Duct-specific", times = Duct_specific_ligands %>% length()),
                  rep("Mast-specific", times = Mast_specific_ligands %>% length()),
                  rep("Serous-specific", times = Serous_specific_ligands %>% length()),
                  rep("General", times = general_ligands %>% length())),
  
  ligand = c(CD4_specific_ligands,CD8_specific_ligands,
             Bcell_specific_ligands,PC_specific_ligands, Neutrophil_specific_ligands,
             Duct_specific_ligands, Mast_specific_ligands,Serous_specific_ligands,  general_ligands
             ))


active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "PC") %>% inner_join(ligand_type_indication_df) # if you want ot make circos plots for multiple gene sets, combine the different data frames and differentiate which target belongs to which gene set via the target type
#active_ligand_target_links_df <-active_ligand_target_links_df[-41,]
active_ligand_target_links_df

cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0)
active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)
ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())
circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

grid_col_ligand =c("General" = "lawngreen",
                   "CD4-specific" = "darkgreen",
                   "CD8-specific" = "steelblue4",
                   "PC-specific" = "steelblue1",
                   "Bcell-specific"="orange4",
                   "Mast-specific"="violetred",
                   "cDC-specific"="lightpink4",
                   "Serous-specific"="khaki")
grid_col_target =c(
  "PC" = "tomato")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
library(dplyr)
links_circle = circos_links %>% select(ligand,target, weight)
ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)
grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 0.8-weight) %>% .$transparency 
target_order = circos_links$target %>% unique()
ligand_order = c(CD4_specific_ligands,CD8_specific_ligands,PC_specific_ligands,
                 Bcell_specific_ligands,Mast_specific_ligands,Duct_specific_ligands,Serous_specific_ligands, general_ligands
  
) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)
width_same_cell_same_ligand_type = 2
width_different_cell = 1
width_ligand_target = 1
width_same_cell_same_target_type = 2
table(circos_links[,"ligand_type"])

gaps = c(
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CD4-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CD8-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Bcell-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "PC-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Mast-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Serous-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "PC") %>% distinct(target) %>% nrow() -1)),width_ligand_target
)
circos.clear()
circos.par(gap.degree = gaps)


chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, 
             link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, 
             diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", 
             link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) 
legend(x = 1.2, y = 1, 
       legend = names(grid_col), 
       fill = grid_col, 
       title = "Cells",
       cex = 0.8, xpd=TRUE)
#Save plot image
