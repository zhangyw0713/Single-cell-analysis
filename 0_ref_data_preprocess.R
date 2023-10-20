pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "loomR", "SummarizedExperiment", "Matrix", "anndata")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

# 0 This script is for generating singlecell data matrix and annotation files
rds_fn <- "/cluster/home/projects/mef2d/data/public/rds"

# ====== adult bone marrow ======

bm_1 <- readr::read_csv("/cluster/home/projects/mef2d/data/public/CensusImmune-BoneMarrow-10x_cell_type_2020-03-12.csv")
l_bm_1 <- connect(filename = "/cluster/home/projects/mef2d/data/public/1M-immune-human-immune-10XV2.loom",
                  mode = "r+", skip.validate = T)

mat_bm_1 <- l_bm_1[["matrix"]][, ]
coln_bm1 <- l_bm_1[["col_attrs/cell_names"]][]
rown_bm1 <- l_bm_1[["row_attrs/gene_names"]][]
colnames(mat_bm_1) <- rown_bm1
rownames(mat_bm_1) <- coln_bm1

bm_1 <- as.data.frame(bm_1)
colnames(bm_1)[1:5] <- c("document_id", "biomaterial_id", "cell_identity", "cell_ontology", "cell_label")
toj1 <- data.frame(biomaterial_id = strsplit(l_bm_1[["attrs/input_name"]][], ", ")[[1]],
                   order_index = 0:(length(unique(bm_1$document_id)) - 1))
bm_1 <- left_join(bm_1, toj1)
bm_1$cname <- paste0(bm_1$barcode, "-", bm_1$order_index)

mat_bm_1_sub <- mat_bm_1[bm_1$cname, ]
mat_bm_1_sub <- Matrix::Matrix(mat_bm_1_sub, sparse=TRUE)   # 64G rds object to 3G

readr::write_rds(mat_bm_1_sub, glue::glue("{rds_fn}/adultbm_cells_matrix_sub.rds"))
readr::write_rds(bm_1, glue::glue("{rds_fn}/adultbm_cells_metadata.rds"))


# ======= cordblood ========

bld_1 <- readr::read_csv("/cluster/home/projects/mef2d/data/public/CensusImmune-CordBlood-10x_cell_type_2020-03-12.csv")
l_bld_1 <- connect(filename = "/cluster/home/projects/mef2d/data/public/1M-immune-human-blood-10XV2.loom",
                   mode = "r+", skip.validate = T)

mat_bld_1 <- l_bld_1[["matrix"]][, ]
coln_bld1 <- l_bld_1[["col_attrs/cell_names"]][]
rown_bld1 <- l_bld_1[["row_attrs/gene_names"]][]
colnames(mat_bld_1) <- rown_bld1
rownames(mat_bld_1) <- coln_bld1

bld_1 <- as.data.frame(bld_1)
colnames(bld_1)[1:5] <- c("document_id", "biomaterial_id", "cell_identity", "cell_ontology", "cell_label")
toj2 <- data.frame(biomaterial_id = strsplit(l_bld_1[["attrs/input_name"]][], ", ")[[1]],
                   order_index = 0:(length(unique(bld_1$document_id)) - 1))
bld_1 <- left_join(bld_1, toj2)
bld_1$cname <- paste0(bld_1$barcode, "-", bld_1$order_index)

mat_bld_1_sub <- mat_bld_1[bld_1$cname, ]
mat_bld_1_sub <- Matrix::Matrix(mat_bld_1_sub, sparse=TRUE)

readr::write_rds(mat_bld_1_sub, glue::glue("{rds_fn}/cordblood_cells_matrix_sub.rds"))
readr::write_rds(bld_1, glue::glue("{rds_fn}/cordblood_cells_metadata.rds"))


# ====== fetal bone marrow ======

# the original data was stored in anndata format, 
# switch to a conda environment with anndata module installed (e.g. scanpy-1.8.2)
h5ad1 <- read_h5ad("/cluster/home/projects/mef2d/data/public/fig1b_fbm_scaled_gex_updated_dr_20210104.h5ad")
fbm_mat <- as(as.matrix(h5ad1$X), "dgCMatrix")
readr::write_rds(fbm_mat, "/cluster/home/projects/mef2d/data/h5ad1_mat_dgc.rds")
readr::write_rds(h5ad1$var, "/cluster/home/projects/mef2d/data/h5ad1_var.rds")
readr::write_rds(h5ad1$obs, "/cluster/home/projects/mef2d/data/h5ad1_obs.rds")


# ====== leukemia dataset from PMID:30827681 ========
# more data annotations
data_fn <- "/cluster/home/projects/mef2d/data"
ref_loc <- "/cluster/home/projects/mef2d/data/public/GSE116256_RAW"
ref_files <- list.files(ref_loc)
ref_annos <- ref_files[grepl("anno", ref_files)]
ref_mats <- ref_files[grepl("dem", ref_files)]
#ref_nanopores <- ref_files[grepl("nanopore", ref_files)]
ref_samples <- sub(".dem.txt", "", ref_mats)

test <- lapply(glue::glue("{ref_loc}/{ref_mats}"), read_delim)
for (i in 1:length(test)){
  names(test)[i] <- ref_samples[i]
  test[[i]] <- test[[i]] %>% as.data.frame() %>% column_to_rownames("Gene")
}
ref_mat <- bind_cols(test)
readr::write_rds(ref_mat, glue::glue("{/data_fn}/ref_cells_lsc_mat.rds"))

test <- lapply(glue::glue("{ref_loc}/{ref_annos}"), read_delim)
ref_anno <- as.data.frame(bind_rows(test))
readr::write_rds(ref_anno, glue::glue("{data_fn}/ref_cells_lsc_anno.rds"))







