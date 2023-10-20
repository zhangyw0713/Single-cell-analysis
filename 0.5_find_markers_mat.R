pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "SummarizedExperiment", "Matrix", "Seurat")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}



# ===== find markers for adult bone marrow dataset


# merits all the same as below, fill up later











# ==== find markers for leukemia dataset from PMID:30827681 =======

ref_mat <- readr::read_rds("/cluster/home/projects/mef2d/data/ref_cells_lsc_mat.rds")
ref_anno <- readr::read_rds("/cluster/home/projects/mef2d/data/ref_cells_lsc_anno.rds")
ref_mat <- as(as.matrix(ref_mat), "dgCMatrix")
rownames(ref_anno) <- ref_anno$Cell
leu_ref <- CreateSeuratObject(counts = ref_mat, project = "leu1_anno", meta.data = ref_anno,
                              min.cells = 3, min.features = 200)
Idents(leu_ref) <- leu_ref@meta.data$CellType
leu_ref <- NormalizeData(leu_ref)
all_markers_leu <- FindAllMarkers(leu_ref)
readr::write_rds(all_markers_leu, "/cluster/home/projects/mef2d/data/public/rds/all_markers_leu_1.rds")


