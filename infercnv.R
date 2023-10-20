pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "cowplot",
          "Seurat", "vroom", "harmony", "SingleR", "scuttle", "Matrix", "infercnv")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}



obj <- readr::read_rds("/cluster/home/projects/mef2d/data/erdata_no_harmony.rds")
sce_hd <- readr::read_rds("/cluster/home/jhuang/projects/mef2d/analysis/zwn/human/tenx/seurat/merge/sce_hd.rds")
setwd("/cluster/home/projects/mef2d/output/2")

sce_hd_ref <- sce_hd[,sce_hd$cell_types %in% c("pre_proB", "proB", "preB_I", "preB_II", "immatureB", "matureB")]
idx <- sample(1:18989, 6000) %>% sort
sce_hd_ref <- sce_hd_ref[, idx]
sce_hd_ref[["cell_types"]] <- str_c("HD_", sce_hd_ref$cell_types)
sce_ball_obs <- obj[, obj$cell_types %in% c("matureB", "immatureB","preB_I", "preB_II", "proB", "pre_proB")]
idx2 <- sample(1:27300, 21000) %>% sort
sce_ball_obs_e <- sce_ball_obs[, idx2]

#E1
sce_infercnv_e <- merge(sce_hd_ref, sce_ball_obs_e)

# sample of E1 ------------------------------------------------------------
output_dir <- "./E" %>% checkdir()
gtf_fn <- "/cluster/home/ylxie_jh/share/annotation/infercnv/gene_info.txt"
object <- sce_infercnv_e

#创建gene_info
gene_info <- readr::read_table(gtf_fn, col_names = FALSE) %>%
  dplyr::transmute(gene_name=X4, seqnames=X1, start=X2, end=X3) %>%
  dplyr::filter(gene_name %in% rownames(object)) %>%
  dplyr::filter(!duplicated(gene_name))
readr::write_tsv(gene_info, file = glue("{output_dir}/gene_info.txt"), col_names = FALSE)

#创建exp_file
exp <-  GetAssayData(object[["RNA"]], slot = 'count') %>% as.data.frame()
exp <- exp[gene_info$gene_name, ]
readr::write_tsv(exp %>% rownames_to_column(), file = glue("{output_dir}/exp_file.txt"), col_names = TRUE)

#创建group_info
ref_group_name <- c("HD_matureB", "HD_immatureB", "HD_preB_II", "HD_preB_I", "HD_proB", "HD_pre_proB", "matureB")

groupinfo <- data.frame(cellID = names(object$cell_types), cellType = object$cell_types)
readr::write_tsv(groupinfo, glue("{output_dir}/groupinfo.txt"), col_names = FALSE)

#run infercnv
result_dir <- glue::glue("{output_dir}/result") %>% checkdir()
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = glue("{output_dir}/exp_file.txt"),
                                     annotations_file = glue("{output_dir}/groupinfo.txt"),
                                     gene_order_file = glue("{output_dir}/gene_info.txt"),
                                     delim = '\t',
                                     ref_group_names = ref_group_name,
                                     chr_exclude = c("chrM"))

readr::write_rds(infercnv_obj, glue::glue("{output_dir}/infercnv_obj.rds"))

#random_trees
infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff = 0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                              out_dir = result_dir,
                              output_format = "pdf",
                              cluster_by_groups = T,
                              analysis_mode='subclusters',
                              tumor_subcluster_partition_method = "random_trees",
                              denoise=T,
                              tumor_subcluster_pval=0.05,
                              num_threads = 40,
                              HMM= FALSE)
write_rds(infercnv_obj, file = "infercnv_obj1.rds")


#infercnv_groups <- read.table("/cluster/home/projects/mef2d/output/2/E/result/infercnv.observation_groupings.txt")



# re-order the output plot

order1 <- rev(c("pre_proB", "proB", "preB_I", "preB_II", "immatureB"))
order2 <- c("HD_pre_proB", "HD_proB", "HD_preB_I", "HD_preB_II", "HD_immatureB", "HD_matureB", "matureB")

infercnv_obj_e <- readr::read_rds("/cluster/home/projects/mef2d/analysis/zwn/home/infercnv/infercnv_obj_E.rds")

infercnv_obj_e@reference_grouped_cell_indices <- infercnv_obj_e@reference_grouped_cell_indices[order2]
infercnv_obj_e@observation_grouped_cell_indices <- infercnv_obj_e@observation_grouped_cell_indices[order1]

plot_cnv(infercnv_obj_e,
         cluster_by_groups=TRUE,
         cluster_references=TRUE,
         out_dir=output_dir,
         output_format = "pdf",
         output_filename="rev_chr_infercnv_E"
)









