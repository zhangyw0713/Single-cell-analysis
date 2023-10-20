



# (after the metacells filtering)
sce_tumor <- readr::read_rds("/cluster/home/projects/mef2d/output/sce_tumor.rds")
Idents(sce_tumor) <- "cluster_final"
sce_tumor2 <-  subset(x = sce_tumor, subset = cluster_final %in% 1:18)
sce_tumor2 <- RunUMAP(sce_tumor2, dims = 1:10)
DimPlot(sce_tumor2, group.by = 'orig.ident', reduction = 'umap')

DoHeatmap(sce_tumor2, features = c(gene03, "MME"))

# first, prepare the UMAP
# extensively tested that harmony and seurat rpca will harm the structure and bury real differences
# so we won't do batch correction, but reduce the dims to avoid high-distinct clusters

ancol4 <- readr::read_rds("/cluster/home/projects/mef2d/data/test/ancol4.rds")
fil_anno <- sce_tumor_fil@meta.data[, c("cell_types_all", "cell_types_broad")]

quick_anno <- function(sce, anno){
  temp <- sce@meta.data
  temp <- left_join(temp, data.frame(metacell = 0:823, anno = as.character(anno)))
  sce@meta.data$anno <- temp$anno
  sce
}
#sce_tumor2$anno <- NULL

sce_tumor2 <- quick_anno(sce_tumor2, ancol4$cluster)
sce_tumor2@meta.data$cell_types_all <- sce_tumor_fil@meta.data$cell_types_all
sce_tumor2@meta.data$cell_types_broad <- sce_tumor_fil@meta.data$cell_types_broad
sce_tumor2 <- RunUMAP(sce_tumor2, dims = 1:10)



DimPlot(sce_tumor2, group.by = "anno", reduction = "umap")
DimPlot(sce_tumor2, group.by = "cell_types_all", reduction = "umap", cols = cluster_colors_all)

meta_others <- unique(sce_tumor2$metacell[sce_tumor2$cell_types_broad %in% c("preB", "proB")])



meta <- sce_tumor_fil@meta.data
meta0 <- meta %>% group_by(metacell) %>% mutate(perc = max(table(orig.ident)) / sum(table(orig.ident)))
meta00 <- meta0 %>% group_by(metacell) %>% summarise(id = names(table(orig.ident)[which.max(table(orig.ident))]))
mcxxx <- unique(meta0$metacell[meta0$perc < 0.5])
meta00$id[meta00$metacell %in% mcxxx] <- NA
meta00 <- meta00[order(meta00$metacell), ]

pht6 <- readr::read_rds("/cluster/home/projects/mef2d/code/plotht.rds")
gene001 <- pht6$tree_row$labels[pht6$tree_row$order]

gene002 <- gene001[-c(1:3, 5, 6 , 8, 11, 14, 26:28, 43:60, 77:82, 84:85)]
gene003x <- setdiff(gene002, c("MS4A4A", "IRF8", "BCL11B", "IL7R", "S100A9", "IGHG1", "MLLT3", "CD3G",
                              "SPI1", "HBD", "CD68"))

rna_mat <- sce_tumor2[["RNA"]]
rnat <- t(rna_mat[])
D <- data.frame(metacell = meta$metacell, rnat)
rna <- D %>% group_by(metacell) %>%
  summarise_at(vars(-group_cols()), sum)
umis_rna <- rna %>% as.data.frame() %>% column_to_rownames(var = "metacell") %>% t() 
rownames(umis_rna) <- sub("\\.", "-", rownames(umis_rna))
log_rna <- log2(1 + umis_rna)
log_rna <- t(log_rna)
median_log_rna <- apply(log_rna, 2, median)
relative_log_rna <- sweep(log_rna, 2, median_log_rna)



# ====== plot =======

# rep
mcn_order <- pht6x$tree_col$labels[pht6x$tree_col$order]
mcn_colclust <- cutree(pht6x$tree_col, k = 16)
mcn_cco <- mcn_colclust[mcn_order]
annotation_col3$cluster <- mcn_colclust
mcn_order <- pht6x$tree_col$labels[pht6x$tree_col$order]
mcn_colclust <- cutree(pht6x$tree_col, k = 30)
mcn_cco <- mcn_colclust[mcn_order]
annotation_col3$cluster_30 <- mcn_colclust

annotation_col3$cluster <- ancol4$cluster

# anno
annotation_col3 <- unique(meta[, c("metacell", "cell_types_all", "pre_anno_pred",
                                   "labels_leu1_pred", "cell_types_pred")])
annotation_col3 <- annotation_col3[order(annotation_col3$metacell), ]
annotation_col3$id <- meta00$id
annotation_col4 <- annotation_col3 %>% remove_rownames() %>% column_to_rownames(var = "metacell")

# colors
cluster_colors2 <- c("#ee5d5a", "#34956c", "#6b6498", "#71a3a2", "#c88978", "#cc9b32", "#53096a", "#03bb26", 
                              RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(8, "Pastel1"), 
                              RColorBrewer::brewer.pal(8, "Set1"))[1:16]
names(cluster_colors2) <- unique(annotation_col3$cluster)
cluster_colors20 <- c("#ee5d5a", "#34956c", "#6b6498", "#71a3a2", "#c88978", "#cc9b32", "#53096a", "#03bb26", 
                              RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(8, "Pastel1"), 
                              RColorBrewer::brewer.pal(8, "Set1"))[1:20]
names(cluster_colors20) <- unique(annotation_col3$cluster_20)
cluster_colors30 <- c("#ee5d5a", "#34956c", "#6b6498", "#71a3a2", "#c88978", "#cc9b32", "#53096a", "#03bb26", 
                               RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(8, "Pastel1"), 
                               RColorBrewer::brewer.pal(8, "Set1"))[1:30]
names(cluster_colors30) <- unique(annotation_col3$cluster_30)


leu1_colors <- c("seagreen", "yellow3", "red", "brown3", "ivory", "lightgreen", "lightblue", "orange")
names(leu1_colors) <- sort(na.omit(unique(annotation_col3$labels_leu1_pred)))
pre_anno_colors <- c("magenta", "red", "blue", "darkgreen", "orange", "chartreuse", "limegreen", "cyan", "lightgreen")
names(pre_anno_colors) <- sort(na.omit(unique(annotation_col3$pre_anno_pred)))
cell_type_colors <- c("brown3", "red", "green", "purple", "orange", "limegreen", RColorBrewer::brewer.pal(4, "Paired"))
names(cell_type_colors) <- sort(na.omit(unique(annotation_col3$cell_types_pred)))
cluster_colors_all <- c("#ee5d5a", "#34956c", "#6b6498", "#71a3a2", "#c88978", "#cc9b32", "#53096a", "#03bb26", 
                                 RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(3, "Pastel1"))
names(cluster_colors_all) <- unique(annotation_col3$cell_types_all)
id_colors <- c("#c88978", "#cc9b32", "gray", RColorBrewer::brewer.pal(3, "Pastel1"))
names(id_colors) <- unique(annotation_col3$id)

# 
annotation_colors4 <- list(cluster = cluster_colors2,
                           cluster_20 = cluster_colors20,
                           cluster_30 = cluster_colors30,
                           labels_leu1_pred = leu1_colors,
                           pre_anno_pred = pre_anno_colors,
                           cell_types_pred = cell_type_colors,
                           cell_types_all = cluster_colors_all,
                           id = id_colors)
# set bar
breaks <- pracma::interp1(0:7, c(-3, -2, -1, 0, 1, 2, 3, 4), 0:140/20)
colors <- colorRampPalette(c('darkblue', 'blue', 'lightblue', 'white', '#ffcccb', 'red', 'darkred', 'darkred'))(141)
options(repr.plot.width = 25, repr.plot.height = 13)


gene003y <- readr::read_rds("/cluster/home/projects/mef2d/code/gene03.rds")

pht6data <- t(relative_log_rna[, c(gene003y)])

pht6x <- pheatmap::pheatmap(
  pht6data,
  # cluster_rows = F,
  # clustering_method = "mcquitty",
  treeheight_col= 90,
  treeheight_row=0,
  cellwidth=1,
  cellheight=10,
  show_rownames=TRUE,
  show_colnames=FALSE,  
  main='Metacell Clusters - BALL anno Genes',
  color=colors,
  breaks=breaks,
  #annotation_row=annotation_row3,
  annotation_col=annotation_col4,
  annotation_colors=annotation_colors4,
  annotation_legend=TRUE,
  legend=TRUE
)

pdf("/cluster/home/projects/mef2d/output/heatmap_tumor_666.pdf", width = 15, height = 28)
pht6x
dev.off()


fin <- data.frame(cluster = unique(mcn_cco), 
                  cell_types_all = c("erythroid_cell", "myeloid", "CLP", "T", "NK", 
                                     "NK", "matureB", "matureB", "matureB", "ery_like_proB", 
                                     "proB_EBF1p", "proB_SOX4p", "proB_SOX4p", "pre_proB", "immatureB", 
                                     "proB_CD83p", "proB_LEF1p", "proB_LEF1p", "proB_LEF1p", "preB_CD74p", 
                                     "preB_CD74p", "preB_CD74p", "immatureB", "proB_BCL11Bp", "preB_II_EBF1p", 
                                     "pre_proB", "preB_II_ATXN1p", "preB_HLADRp", "preB_BLKp", "preB_BLKp"),
                  cell_types_broad = c("erythroid_cell", "myeloid", "CLP", "NK_T", "NK_T", 
                                     "NK_T", "matureB", "matureB", "matureB", "proB", 
                                     "proB", "proB", "proB", "pre_proB", "immatureB", 
                                     "proB", "proB", "proB", "proB", "preB_I", 
                                     "preB_I", "preB_I", "immatureB", "proB", "preB_II", 
                                     "pre_proB", "preB_II", "preB_I", "preB_I", "preB_I"))
fin_meta <- data.frame(metacell = names(mcn_cco), cluster = mcn_cco)
fin_meta <- left_join(fin_meta, fin)
fin_meta$metacell <- as.numeric(fin_meta$metacell)
fin_meta <- fin_meta[order(fin_meta$metacell), ]

quick_anno <- function(sce, anno){
  temp <- sce@meta.data
  temp <- left_join(temp, data.frame(metacell = 0:823, anno = as.character(anno)))
  sce@meta.data$cell_types_final_broad <- temp$anno
  sce
}
sce_tumor2$cell_types_final_broad <- NULL
sce_tumor2 <- quick_anno(sce_tumor2, fin_meta$cell_types_broad)


quick_anno <- function(sce, anno){
  temp <- sce@meta.data
  temp <- left_join(temp, data.frame(metacell = 0:823, anno = as.character(anno)))
  sce@meta.data$cell_types_final <- temp$anno
  sce
}
sce_tumor2$cell_types_final <- NULL
sce_tumor2 <- quick_anno(sce_tumor2, fin_meta$cell_types_all)



sce_tumor2$cell_types_final_broad[sce_tumor2$cell_types_all %in% "proB_IGHDp" &
                                    sce_tumor2$cell_types_final %in% "preB_II_ATXN1p"] <- "proB"




readr::write_rds(sce_tumor2, "/cluster/home/projects/mef2d/data/sce_tumor_final.rds")


facet(DimPlot(sce_tumor2, group.by = 'cell_types_final_broad', reduction = 'umap') , facet.by = "cell_types_final_broad")
facet(DimPlot(sce_tumor2, group.by = 'cell_types_broad', reduction = 'umap') , facet.by = "cell_types_broad")



facet(DimPlot(sce_tumor2, group.by = 'anno', reduction = 'umap') , facet.by = "anno")

Idents(sce_tumor2) <- "cell_types_final_broad"
DotPlot(object = sce_tumor2, features = gene03) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

Idents(sce_tumor2) <- "cell_types_final"
DotPlot(object = sce_tumor2, features = gene01) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))


DimPlot(sce_tumor2, group.by = 'cell_types_final_broad', reduction = 'umap', label = T)

DimPlot(sce_tumor2, group.by = 'cell_types_final', reduction = 'umap', label = T)




DimPlot(sce_tumor2, group.by = 'orig.ident', reduction = 'umap', label = T)


head(sce_tumor_fil@meta.data)






