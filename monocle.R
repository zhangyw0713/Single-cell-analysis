pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "igraph",
          "Seurat", "CytoTRACE", "reticulate", "monocle", "numDeriv"#, "monocle3"
)
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

# read in data
mef2d <- readr::read_rds("/cluster/home/projects/mef2d/data/merged_4_tumors_new.rds")
sct_assay <- GetAssayData(mef2d)
rna_assay <- GetAssayData(mef2d, assay = "RNA", slot = "data")
meta_mef2d <- mef2d@meta.data
# subset to remove cells
Idents(mef2d) <- "cell_types_broad"
mef2d_sub <- subset(mef2d, `cell_types_broad` %notin% c("NK_T", "erythroid_cell", "myeloid"))
# set colors
cell_color <- c("HSC_MPP" = "indianred1", "CLP" = "indianred3", "pre_proB" = "darkseagreen1",
                "proB" = "lightgreen", "preB" = "darkolivegreen1",  "preB_I" = "yellowgreen",
                "preB_II" = "chartreuse4", "immatureB" = "#3CB371", "matureB" = "darkgreen",
                "erythroid_cell" = "lightpink2", "NK_T" = "#9575aa", "myeloid" = "deeppink4", "tumor1" = "grey",
                "tumor_cell1" = "yellow4", "tumor_cell2" = "aquamarine2", "tumor_cell3" = "lightseagreen")
samples <- unique(meta_mef2d$orig.ident)
sample_color <- c("#ee5d5a", "#9f2f6a", "#34956c", "#6b6498","#b29bc9", "#71a3a2", 
                           "#c88978",  "#a15217", "#ce662f", "#cc9b32", "#53096a", "#6569c2", "#66676c")[1:length(samples)]
names(sample_color) <- samples
                           
# output umap dimplot
psub1 <- DimPlot(mef2d_sub, group.by = "cell_types_broad", reduction = 'umap', label = T) + NoLegend() + 
  scale_color_manual(values = cell_color)

ggsave(psub1, filename = "/cluster/home/projects/temp/mef2d_umap_allsamples.pdf", width = 5, height = 4)

psubsep1 <- DimPlot(mef2d_sub, split.by = "cell_types_broad", reduction = 'umap', label = F) + NoLegend() + 
  scale_color_manual(values = cell_color)

ggsave(psubsep1, filename = "/cluster/home/projects/temp/mef2d_umap_all_sep.pdf", width = 12, height = 3)




# ====== CytoTRACE ==========

cyt_res_mef2d_sub <- CytoTRACE(as.matrix(mef2d_sub[["RNA"]]@counts), ncores = 30)

# plot the CytoTRACE results
#pdf("/cluster/home/projects/temp/plotCytoTRACE_umap_sub.pdf", width = 12, height = 4)
#plotCytoTRACE(cyt_res_mef2d_sub, colors = NULL, emb = mef2d_sub@reductions$umap@cell.embeddings, 
#              phenotype = cell_types_broad)
#dev.off()

# prepare plot data
cell_type <- as.character(mef2d_sub@meta.data$cell_types_broad)
orig_ident <- as.character(mef2d_sub@meta.data$orig.ident)
names(cell_type) <-rownames(mef2d_sub@meta.data)
names(orig_ident) <- rownames(mef2d_sub@meta.data)

meta_sub <- mef2d_sub@meta.data
emb <- as.matrix(mef2d_sub@reductions$umap@cell.embeddings)
mat <- t(cyt_res_mef2d_sub$exprMatrix)[rownames(emb), ]
cyto <- cyt_res_mef2d_sub$CytoTRACE[rownames(emb)]

# cyto_plot_datfm <- data.frame(emb, CytoTRACE_score = cyto, HDAC9 = mat[, "HDAC9"], 
#                               MEF2D = mat[, "MEF2D"], meta_sub)
# 
# 
# # CytoTRACE plot by ggplot2
# temp_color <- RColorBrewer::brewer.pal(11, "Spectral")
# temp_color[6] <- "gold"
# rbPal <- colorRampPalette(temp_color)
#   
# p_cyt <- ggplot(cyto_plot_datfm, aes(x = UMAP_1, y = UMAP_2, color = CytoTRACE_score)) + geom_point(size = .5) + 
#     ggpubr::theme_pubr() + scale_colour_gradientn(name = "CytoTRACE score", colours = rev(rbPal(50)),
#                                                   guide = ggplot2::guide_colourbar(ticks.colour = "black", ticks.linewidth = 1, frame.colour = "black"),
#                                                   breaks = seq(0, 1, 0.2), labels = c("0.0 (More diff.)", 0.2, 0.4, 0.6, 0.8, "1.0 (Less diff.)")) +
#     theme(legend.text = ggplot2::element_text(size = 6), legend.title = ggplot2::element_text(size = 7), 
#           plot.title = ggplot2::element_text(size = 10, hjust = 0.5), axis.title.x = ggplot2::element_text(size = 8), 
#           axis.title.y = ggplot2::element_text(size = 7), axis.text = ggplot2::element_text(size = 6), 
#           legend.position = "right", plot.margin = ggplot2::unit(c(0.5, 1, 0.5, 1), "cm"))
# 


# ======= monocle2 ========

"reducedDimW<-" <- function (cds, value) {
  stopifnot(is(cds, "CellDataSet"))
  cds@reducedDimW <- value
  validObject(cds)
  cds
}
"reducedDimK<-" <- function (cds, value) {
  stopifnot(is(cds, "CellDataSet"))
  cds@reducedDimK <- value
  validObject(cds)
  cds
}


# construct cds in monocle2 style
mef2d_sub_data <- GetAssayData(mef2d_sub, assay = "RNA", slot = 'counts')
fd <- data.frame(gene_short_name = row.names(mef2d_sub_data), row.names = row.names(mef2d_sub_data))
pd <- new('AnnotatedDataFrame', data = meta_sub)
fd1 <- new('AnnotatedDataFrame', data = fd)

cds0 <- monocle::newCellDataSet(mef2d_sub_data,
                                phenoData = pd,
                                featureData = fd1,
                                lowerDetectionLimit = 0.5,
                                expressionFamily = VGAM::negbinomial.size())

cds0 <- cds0 %>% estimateSizeFactors() %>% estimateDispersions()

# select DE genes among celltypes
cds1 <- monocle::detectGenes(cds0, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds1), num_cells_expressed >= 0.1 * dim(mef2d_sub_data)[2])) 
# if monocle3 has been loaded, detach it: detach("package:monocle3",unload = TRUE)
diff_test_res <- monocle::differentialGeneTest(cds1[expressed_genes, ],
                                               fullModelFormulaStr = "~cell_types_broad", cores = 1)

deg <- subset(diff_test_res, qval < 0.001)
deg <- deg[order(deg$qval, decreasing = F), ]
#ordering_genes <- row.names(deg[1:2000,])
ordering_genes <- row.names(deg[1:1000, ])


# orderingfilter cds object for plot
cds2 <- monocle::setOrderingFilter(cds1[expressed_genes, ], ordering_genes)
readr::write_rds(cds2, "/cluster/home/projects/mef2d/analysis/cds_ordered_m2x.rds")

cds_test2 <- monocle::reduceDimension(cds2, max_components = 2, method = 'DDRTree')
readr::write_rds(cds_test2, "/cluster/home/projects/mef2d/analysis/cds_ordered_m2x.rds")

#cds_test2 <- readr::read_rds("/cluster/home/projects/mef2d/analysis/cds_ordered_m2x.rds")
#source("/cluster/home/projects/mef2d/code/order_cells.R")

cds_test2 <- orderCells(cds_test2)
readr::write_rds(cds_test2, "/cluster/home/projects/mef2d/analysis/cds_ordered_m2xx.rds")

# plot
pm1 <- monocle::plot_cell_trajectory(cds_test2, color_by = "Pseudotime", size = 1, show_backbone = TRUE)
pm2 <- monocle::plot_cell_trajectory(cds_test2, color_by = "cell_types_broad", size = 1, show_backbone = TRUE) +
  scale_color_discrete(cell_color)
pm_all <- pm1 + pm2 + plot_layout(nrow = 1)

pdf("/cluster/home/projects/temp/trajectory_ddrtree_monocle2_0321.pdf", width = 10, height = 5)
pm_all
dev.off()


table(pData(cds_test2)[, c("cell_types_broad", "State")])
# manually change state in pdata?
cds_test3 <- orderCells(cds_test2, root_state = 8)

new_color <- cell_color[match(levels(as.factor(pData(cds_test3)$cell_types_broad)), names(cell_color))]

readr::write_rds(cds_test3, "/cluster/home/projects/mef2d/output/ddrtree_object_tumor.rds")

# reordered plot
pm1 <- monocle::plot_cell_trajectory(cds_test3, color_by = "Pseudotime", size = 1, show_backbone = TRUE)
pm2 <- monocle::plot_cell_trajectory(cds_test3, color_by = "cell_types_broad", size = 1, show_backbone = TRUE) +
  scale_color_discrete(new_color) + theme(legend.title = element_blank())
pm_all <- pm1 + pm2 + plot_layout(nrow = 1)

pdf("/cluster/home/projects/temp/trajectory_ddrtree_monocle2_0321_ordered.pdf", width = 10, height = 5)
pm_all
dev.off()
