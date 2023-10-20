pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "igraph",
          "Seurat", "anndata", "stats", "pheatmap", "pracma", "SeuratDisk"
)
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}



sce_tumor <- readr::read_rds("/cluster/home/projects/mef2d/output/sce_tumor.rds")
Idents(sce_tumor) <- "cluster_final"
sce_tumor2 <-  subset(x = sce_tumor, subset = cluster_final %in% 1:18)
sce_tumor_fil <- FindNeighbors(sce_tumor_fil, dims = 1:30, verbose = F) 
sce_tumor_fil <- FindClusters(sce_tumor_fil, verbose = F)
sce_tumor2 <- RunUMAP(sce_tumor2, dims = 1:10)
DimPlot(sce_tumor2, group.by = 'orig.ident', reduction = 'umap')
DoHeatmap(sce_tumor2, features = c(gene03, "MME"))


sce_tumor_fil <- subset(x = sce_tumor, subset = cluster_final %in% 1:18)
Idents(sce_tumor_fil) <- "orig.ident"

FeaturePlot(sce_tumor2, features = c("cnv1", "EBF1", "EIF4EBP1", "CD19", "ATXN1", "DNTT", "MKI67", "MS4A1",
                                     "CD22", "KIT", "FLT3", "CD34", "MEF2D", "HDAC9", "STMN1"))

dat_sce <- sce_tumor_fil@assays$RNA@data["IGKC", ][dat_sce["IGKC", ] > 50] <- 50
sce_tumor_fil@assays$RNA@data["HBA1", ][dat_sce["HBA1", ] > 50] <- 50
sce_tumor_fil@assays$RNA@data["HBA2", ][dat_sce["HBA2", ] > 50] <- 50

sct_sce <- sce_tumor_fil@assays$SCT@data



sce_tumor_fil <- FindNeighbors(sce_tumor_fil, dims = 1:30, verbose = F) 
sce_tumor_fil <- FindClusters(sce_tumor_fil, verbose = F)
sce_tumor_fil <- RunUMAP(sce_tumor_fil, dims = 1:10)


facet(DimPlot(sce_tumor_fil, group.by = 'temp_anno', reduction = 'umap') , facet.by = "temp_anno")''
facet(DimPlot(sce_tumor_fil, group.by = 'seurat_clusters', reduction = 'umap') , facet.by = "seurat_clusters")

DimPlot(sce_tumor_fil, group.by = 'orig.ident', reduction = 'umap')
Idents(sce_tumor_fil) <- "orig.ident"
DoHeatmap(sce_tumor_fil, features = c(gene03, "MME", "CD19"))



sce_TCF3@assays$RNA@data



sce_TCF3 <- subset(x = sce_tumor_fil, subset = orig.ident %in% "E1")
sce_MEF2D <- subset(x = sce_tumor_fil, subset = orig.ident %notin% "E1")



# after naming the clusters
# save plots both from TCF3 and MEF2D to one 
# first, overall cellchat
dat_tcf <- sce_TCF3@assays$RNA@data





sce_tumor_fil2 <- sce_tumor2
DefaultAssay(sce_tumor_fil2) <- "SCT"
sce_tumor_fil2@meta.data$group2 <- ifelse(sce_tumor_fil2@meta.data$orig.ident %in% "B069_Dx_CD19", "B069", "others")

sce_tumor_fil2 <- sce_tumor_fil2 %>% 
  RunHarmony("group2", plot_convergence = TRUE)



sce_tumor_fil2 <- sce_tumor_fil2 %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()


sce_tumor_fil2 <- sce_tumor2 %>% 
  RunUMAP(reduction = "harmony", dims = 1:30)

DimPlot(sce_tumor_fil2, reduction = "umap", label = T)


facet(DimPlot(sce_tumor_fil2, reduction = "umap", group.by = "cell_types_broad", label = T), facet.by = "cell_types_broad")
DimPlot(sce_tumor_fil2, reduction = "umap", group.by = "cluster_anno", label = T)
DimPlot(sce_tumor_fil2, reduction = "umap", group.by = "cluster_final_anno", label = T)
facet(DimPlot(sce_tumor_fil2, reduction = "umap", group.by = "cell_types", label = T), facet.by = "cell_types")
facet(DimPlot(sce_tumor_fil2, reduction = "umap", group.by = "cell_types_final", label = T), facet.by = "cell_types_final")


DimPlot(sce_tumor_fil2, reduction = "umap", group.by = "orig.ident")


sce_tumor_fil2 <- quick_anno(sce_tumor_fil2, ancol4$cluster)
sce_tumor_fil2@meta.data$cell_types_all <- sce_tumor_fil@meta.data$cell_types_all
sce_tumor_fil2@meta.data$cell_types_broad <- sce_tumor_fil@meta.data$cell_types_broad
Idents(sce_tumor_fil2) <- "anno"
Idents(sce_tumor_fil2) <- "orig.ident"

DimPlot(sce_tumor_fil2, reduction = "umap")#, group.by = "orig.ident", pt.size = .1, split.by = 'orig.ident')


ancol4 <- readr::read_rds("/cluster/home/projects/mef2d/data/test/ancol4.rds")

quick_anno <- function(sce, anno){
  temp <- sce@meta.data
  temp <- left_join(temp, data.frame(metacell = 0:823, anno = as.character(anno)))
  sce@meta.data$anno <- temp$anno
  sce
}
sce_tumor_fil2$anno <- NULL
sce_tumor_fil2 <- quick_anno(sce_tumor_fil2, ancol4$cluster)
sce_tumor2 <- quick_anno(sce_tumor2, ancol4$cluster)
sce_tumor2@meta.data$cell_types_all <- sce_tumor_fil@meta.data$cell_types_all
sce_tumor2@meta.data$cell_types_broad <- sce_tumor_fil@meta.data$cell_types_broad

Idents(sce_tumor2) <- "anno"
DotPlot(object = sce_tumor2, features = c(gene03, "MME")) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

Idents(sce_tumor2) <- "cell_types_broad"
DotPlot(object = sce_tumor2, features = c(gene03, "MME")) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))



Idents(sce) <- "orig.ident"
sce_hd <- subset(sce, `orig.ident` %notin% c("B069_Dx_CD19", "M2", "M3", "M4", "E1"))
Idents(sce_hd) <- "cell_types"
sce_hd <- subset(sce, `cell_types` %in% names(table(sce_hd$cell_types)))
Idents(sce_hd) <- "cell_types"
DotPlot(object = sce_hd, features = c(gene03, "MME")) + 
    theme(axis.text.x = element_text(angle = 45, hjust=1))



sce_tumor_fil2 <- RunUMAP(sce_tumor_fil2, dims = 1:30, verbose = F)
sce_tumor_fil <- FindNeighbors(sce_tumor_fil, dims = 1:30, verbose = F) 
sce_tumor_fil <- FindClusters(sce_tumor_fil, verbose = F)




c("CLP", "")




facet(DimPlot(sce_tumor2, group.by = 'anno', reduction = 'umap'), facet.by = "seurat_clusters")
facet(DimPlot(sce_tumor_fil2, group.by = 'cell_types_all', reduction = 'umap'), facet.by = "cell_types_all")
facet(DimPlot(sce_tumor_fil2, group.by = 'cell_types_broad', reduction = 'umap'), facet.by = "cell_types_broad")
facet(DimPlot(sce_tumor_fil2, group.by = 'seurat_clusters', reduction = 'umap'), facet.by = "seurat_clusters")
facet(DimPlot(sce_tumor_fil2, group.by = 'cn_cluster', reduction = 'umap'), facet.by = "cn_cluster")



quick_pht <- function(data, cmethod = "hclust"){
  ht <- pheatmap::pheatmap(
    data,
    # cluster_rows = F,
    clustering_method = cmethod,
    treeheight_col = 30, treeheight_row = 0,
    cellwidth = 1, cellheight = 10,
    show_rownames = TRUE, show_colnames = FALSE,  
    main = 'Metacell Clusters - BALL anno Genes',
    color = colors, breaks = breaks,
    annotation_row = annotation_row3,
    annotation_col = annocol4,
    annotation_colors = annotation_colors4,
    annotation_legend = TRUE,
    legend = TRUE
  )
  ht
}

pht6data <- pht5data[gene02, ]
annocol4 <- annotation_col3[, c("pre_anno_pred", "cluster_anno")]

cluster_anno_4 <- c("#ee5d5a", "#34956c", "#6b6498", "#71a3a2", "#c88978", "#cc9b32", "#53096a", 
                             RColorBrewer::brewer.pal(7, "Set2"), RColorBrewer::brewer.pal(4, "Pastel1"))
names(cluster_anno_4) <- sort(na.omit(unique(annocol4$cluster_anno)))
pre_anno_4 <- c("magenta", "red", "green2", "green", "orange", "chartreuse", "limegreen", "cyan", "lightblue")
names(pre_anno_4) <- sort(na.omit(unique(annocol4$pre_anno_pred)))





c7_order <- pht7$tree_col$labels[pht7$tree_col$order]
c7_colclust <- cutree(pht7$tree_col, k = 18)
c7_colclust_15 <- cutree(pht7$tree_col, k = 15)




annocol4$cluster <- c7_colclust
cluster_4 <- cluster_anno_4
names(cluster_4) <- unique(annocol4$cluster)

annocol4$cluster15 <- c7_colclust_15
cluster_154 <- cluster_anno_4[1:15]
names(cluster_154) <- unique(annocol4$cluster_15)

annotation_colors4 <- list(cluster_anno = cluster_anno_4,
                           pre_anno_pred = pre_anno_4,
                           cluster = cluster_4, 
                           cluster15 = cluster_154)

pht7 <- quick_pht(pht6data, "ward.D2")

pdf("/cluster/home/projects/mef2d/output/test_heatmap7ww_tumor.pdf", width = 18, height = 22)
pht7
dev.off()

pht7i <- quick_pht(pht6data, "complete")

pdf("/cluster/home/projects/mef2d/output/test_heatmap7cp_tumor.pdf", width = 18, height = 22)
pht7i
dev.off()

unique(c7_colclust[c7_order])
#[1]  9  8 15 14  2  1    16 17  6 18 10 12     7 13 11  5  3  4






readr::write_rds(annocol4, "/cluster/home/projects/mef2d/output/annocol.rds")
annocol4 <- readr::read_rds("/cluster/home/projects/mef2d/output/annocol.rds")



  meta <- sce_tumor_fil@meta.data
  dat_join <- data.frame(metacell = 0:823, ncluster = annocol4$cluster)
  meta <- left_join(meta, dat_join)
  sce_tumor_fil@meta.data$ncluster <- meta$ncluster



pl <- DimPlot(sce_tumor_fil, group.by = 'ncluster', reduction = 'umap')
plf <- facet(pl, facet.by = "ncluster")

Idents(sce_tumor_fil) <- "ncluster"
DotPlot(object = sce_tumor_fil, features = gene01) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))



cluster_anno <- c("T", "NK", "proB_IGHDp", "preB_ATXN1p", "CLP_myeloid_HSC", "matureB", 
                  "preB_HLADRp", "preB_CD74p", "erythroid", 
                  "Ery_like_proB", "proB_BCL11Bp", "proB_PAX5p", 
                  "proB_I", "proB_VPREB1p", "proB_CD83p",
                  "proB_MKI67p", "preB_BLKp", "proB_EBF1h")
cell_types_broad <- c("NK_T", "NK_T", "proB", "preB", "CLP_myeloid", "matureB", 
                   "preB", "preB", "erythroid", 
                   "proB", "proB", "proB", 
                   "proB", "proB", "proB",
                   "proB", "preB", "proB")
data_join_fin <- data.frame(ncluster = c(9,  8, 15, 14,  2,  1, 16, 17, 6, 18, 10, 12,
                                         7, 13, 11,  5,  3,  4),
                            temp_anno = cluster_anno, cell_types_broad = cell_types_broad)


meta <- left_join(sce_tumor_fil@meta.data, data_join_fin)
meta$cell_types_all <- meta$temp_anno
meta$cell_types_broad[meta$cell_types_broad %in% "CLP_myeloid"] <- meta$cluster_anno[meta$cell_types_broad %in% "CLP_myeloid"]
meta$cell_types_all[meta$cell_types_all %in% "CLP_myeloid_HSC"] <- meta$cluster_anno[meta$cell_types_all %in% "CLP_myeloid_HSC"]
sce_tumor_fil@meta.data$cell_types_all <- meta$cell_types_all
sce_tumor_fil@meta.data$cell_types_broad <- meta$cell_types_broad

readr::write_rds(sce_tumor_fil, "/cluster/home/projects/mef2d/data/sce_tumor_fil.rds")
                 
                 
meta <- sce_tumor_fil@meta.data
                 
gene03 <- c("DNTT", "IGHM", "CD24", "VPREB1", "LEF1", "FLT3", "EIF4EBP1", 
            "EBF1", "PAX5",  "BLK", "SOX4", "STMN1", "CD19", "CD79A", "HLA-DRA", 
            "CD74", "TNFRSF13C", "IGKC", "MS4A1", "LTB", "MEIS1", "ATP8B4", "MLLT3",
            "CD99", "CST3", "NKG7", "CD3D", "CD3G", "SPI1", "CD68", "S100A9", 
            "LILRB2", "HBA1", "HBA2", "CA1", "ALAS2")
Idents(sce_tumor_fil) <- "cell_types_broad"


DotPlot(object = sce_tumor_fil, features = gene03) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))



sce_tumor_fil@meta.data$cell_types_broad[sce_tumor_fil@meta.data$cell_types_all %in% "proB_I"] <- "proB_I"
Idents(sce_tumor_fil) <- "cell_types_broad"
Idents(sce_tumor_fil) <- relevel(Idents(sce_tumor_fil), "proB_I")
Idents(sce_tumor_fil) <- "cell_types_all"
DotPlot(object = sce_tumor_fil, features = gene03) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

DimPlot(sce_tumor, group.by = 'orig.ident', reduction = 'umap')








