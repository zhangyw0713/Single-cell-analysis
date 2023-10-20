pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse",
          "SCopeLoomR", "AUCell", "SCENIC", "KernSmooth", "BiocParallel", "Seurat",
          "plotly", "ComplexHeatmap", "data.table", "grid")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

# read data
sce_tumor <- readr::read_rds("~/projects/mef2d/analysis/zwn/human/tenx/seurat/merge/merge.rds")

# retrieve information from loom file
scenicLoomPath <- "/cluster/home/projects/mef2d/analysis/zwn/home/SCENIC/sample_SCENIC_all.loom"
loom <- open_loom(scenicLoomPath)
exprMat <- get_dgem(loom)
exprMat_log <- log2(exprMat + 1) # Better if it is logged/normalized
regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
regulonAucThresholds <- get_regulon_thresholds(loom)
close_loom(loom)

# rss
cells <- colnames(exprMat)
cellClusters <- data.frame(cell_type = sce_tumor$cell_types)
rss <- calcRSS(AUC = getAUC(regulonAUC), cellAnnotation = cellClusters[colnames(regulonAUC), "cell_type"])
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)
plotRSS_oneSet(rss, setName = "myeloid")















# regulon separated
cellsPerCluster <- split(rownames(cellClusters), cellClusters[, "cell_type"]) 
regulonAUCx <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]
# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUCx)[, cells]))
# Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale = T))
# plot:
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[70:120, ], name = "Regulon activity",
                                   row_names_gp = grid::gpar(fontsize=6))) # row font size
regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)] # to save the clustered regulons for later


# check top regulators
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators$CellType <- factor(as.character(topRegulators$CellType))
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0), ]
dim(topRegulators)




# ====== separate mef2d and tcf3 samples ======

mef_id <- colnames(sce_tumor)[sce_tumor$orig.ident %in% c("M2", "M3", "M4")]
regulonAUC_mef <- regulonAUC[, match(mef_id, colnames(regulonAUC))]
expr_mef <- exprMat[, mef_id]

cells_mef <- colnames(expr_mef)
cellClusters_mef <- data.frame(cell_type = sce_tumor$cell_types_all[match(mef_id, colnames(sce_tumor))], 
                               cell_type_broad = sce_tumor$cell_types_broad[match(mef_id, colnames(sce_tumor))])
rss <- calcRSS(AUC = getAUC(regulonAUC_mef), cellAnnotation = cellClusters_mef[colnames(regulonAUC_mef), "cell_type_broad"])
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)

# just like this...















# regulonAUC find diff-----------------------------------------------------
regulonAUC_merge <- regulonAUC[, colnames(sce_tumor)]
scenic_data <- CreateSeuratObject(regulonAUC_merge@assays@data$AUC)
scenic_data@meta.data$label <- sce_tumor@meta.data$cell_types

Idents(scenic_data) <- 'label'
marker_tf <- NormalizeData(scenic_data) %>% FindAllMarkers(assay = "RNA")
write_csv(marker_tf,file=glue("{workdir}/marker_tf.csv"))

maker_tf$gene <- str_replace_all(maker_tf$gene,"[(+)]","")
maker_tf_list <- maker_tf %>% filter(avg_log2FC >0.8 ) %>% split(.$cluster)

clone_deg <- read_rds('/cluster/home/ylxie_jh/projects/leukemia/analysis/weinazhang/human/infercnv/subcluster/DGEs.rds')
clone_deg_list <- clone_deg %>% filter(avg_log2FC > 0.8) %>% split(.$cluster)
index <- unique(maker_tf$cluster)
overlap_gene <- map(index,function(x){
  return(intersect(maker_tf_list[[x]]$gene,
                   clone_deg_list[[x]]$gene))
})
names(overlap_gene) <- index
write_rds(overlap_gene,file=glue("{workdir}/overlap.rds"))

# Binarized regulon activity ----------------------------------------------

top10 <- marker_tf %>% group_by(cluster) %>% top_n(10, wt = avg_log2FC)
cellsPerCluster <- split(rownames(sce_tumor@meta.data), sce_tumor@meta.data[, "cell_types"])
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]
# Calculate average expression:
regulonActivity_byCellType <- sapply(cellsPerCluster,
                                     function(cells) rowMeans(getAUC(regulonAUC)[, cells]))
# Scale expression:
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale = T))

pdf(file = glue("{plot_dir}/tf_cluster.pdf"),width = 8,height = 12)
hm <- ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[unique(top10$gene),], name="Regulon activity",
                              row_names_gp=grid::gpar(fontsize=6),cluster_columns = T) # row font size
draw(hm)
dev.off()





