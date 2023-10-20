pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "igraph",
          "Seurat", "Asgard", "celldex", "cmapR", "pracma", "SeuratDisk",
          "cowplot", "edgeR"
)
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

# set paths 
lincs_dir <- "/cluster/home/jhuang/reference/scell/ASGARD/LINCS"
output_dir <- "/cluster/home/projects/mef2d/analysis/zwn/home/Asgard"
setwd(output_dir)
set.seed(2023)

PrepareReference(cell.info = glue("{lincs_dir}/GSE70138_Broad_LINCS_cell_info.txt"),
                 gene.info = glue("{lincs_dir}/GSE70138_Broad_LINCS_gene_info.txt"),
                 GSE70138.sig.info = glue("{lincs_dir}/GSE70138_Broad_LINCS_sig_info.txt"),
                 GSE92742.sig.info = glue("{lincs_dir}/GSE92742_Broad_LINCS_sig_info.txt"),
                 GSE70138.gctx = glue("{lincs_dir}/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx"),
                 GSE92742.gctx = glue("{lincs_dir}/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"),
                 Output.Dir = "output_dir"
)


# subset from all to avoid merging
sce_tumor <- readr::read_rds("~/projects/mef2d/analysis/zwn/human/tenx/seurat/merge/merge.rds")
sce_hd <- readr::read_rds("~/projects/mef2d/analysis/zwn/human/tenx/seurat/merge/sce_hd.rds")
sce_all <- readr::read_rds("~/projects/mef2d/analysis/zwn/human/tenx/seurat/merge/all_merged.rds")
new_cell_id <- str_sub(colnames(sce_all), 1 + str_locate(colnames(sce_all), ":")[, 1] / 2)
sce_all <- RenameCells(sce_all, new.names = new_cell_id)
sce_all <- sce_all[, c(colnames(sce_tumor), colnames(sce_hd))]

# set annotation
sce_all$orig.ident[sce_all$orig.ident %in% "Healthy"] <- sce_hd$orig.ident
sce_all@meta.data$cell_types <- c(sce_tumor$cell_types, sce_hd$cell_types)

# set basic 
case_type = "MEF2D"
case = c("M2", "M3", "M4")
target_cell = setdiff(unique(sce_hd$cell_types), c("erythroid_cell", "remove"))
control = setdiff(unique(sce_all$orig.ident), c("E1", "M2", "M3", "M4"))

# separate
mef_id <- colnames(sce_tumor)[sce_tumor$orig.ident %in% case]
tcf_id <- colnames(sce_tumor)[sce_tumor$orig.ident %in% "E1"]
hd_id <- colnames(sce_hd)
sce_all_mef <- sce_all[, c(mef_id, hd_id)]
sce_all_tcf <- sce_all[, c(tcf_id, hd_id)]
sce_all_mef@meta.data$type <- ifelse(sce_all_mef@meta.data$orig.ident %in% case, "MEF2D", "normal")
sce_all_tcf@meta.data$type <- ifelse(sce_all_tcf@meta.data$orig.ident %in% "E1", "TCF3", "normal")



# data <- as.data.frame(sce@meta.data)
# celltypes <- unique(data$celltype)
# final_table <- data.frame()
# for(i in c("normal", case_type)){
#   sub_data <- subset(data, type == i)
#   severity <- unique(sub_data$type)
#   celltype_freq <- round(100 * table(sub_data$celltype) / nrow(sub_data), 2)
#   celltype_freq <- celltype_freq[celltypes]
#   final_table <- rbind(final_table, celltype_freq)
# }
# colnames(final_table) <- celltypes
# rownames(final_table) <- c("normal", case_type)
# 


# functions calculating differential expressed genes
# limma
limma_diff <- function(sce, case, control){
  DefaultAssay(sce) <- "RNA"
  Idents(sce) <- "cell_types"
  min.cells <- 3
  gene_list <- list()
  c_names <- NULL
  for (i in unique(sce@meta.data$cell_types)) {
    c_cells <- subset(sce, cell_types == i)
    Idents(c_cells) <- "type"
    samples <- c_cells@meta.data
    controlsample <- row.names(subset(samples, orig.ident %in% control))
    casesample <- row.names(subset(samples, orig.ident %in% case))
    if(length(controlsample) > min.cells & length(casesample) > min.cells){
      expr <- c_cells@assays$RNA@data
      expr <- expr[, c(casesample, controlsample)]
      new_sample <- data.frame(sample_id = c(casesample, controlsample),
                               type = c(rep("case", length(casesample)), rep("control", length(controlsample))))
      row.names(new_sample) <- paste(new_sample$sample_id, row.names(new_sample), sep = "_")
      bad <- which(rowSums(expr > 0) < 3)
      if(length(bad) > 0){
        expr <- expr[-bad, ]
      }
      mm <- model.matrix(~0 + type, data = new_sample)
      fit <- lmFit(expr, mm)
      contr <- makeContrasts(typecase - typecontrol, levels = colnames(coef(fit)))
      tmp <- eBayes(contrasts.fit(fit, contrasts = contr))
      c_data <- topTable(tmp, sort.by = "P", n = nrow(tmp))
      c_data_for_drug <- data.frame(row.names = row.names(c_data), score = c_data$t, 
                                    adj.P.Val = c_data$adj.P.Val, P.Value = c_data$P.Value)
      gene_list[[i]] <- c_data_for_drug
      c_names <- c(c_names, i)
    }
  }
  names(gene_list) <- c_names
  gene_list
}

# Seurat
seurat_diff <- function(sce, case_type){
  gene_list <- list()
  c_names <- NULL
  Idents(sce) <- "cell_types"
  for (i in unique(sce@meta.data$cell_types)){
    c_cells <- subset(sce, cell_types == i)
    Idents(c_cells) <- "type"
    if (length(table(c_cells$type)) > 1) {
      c_data <- FindMarkers(c_cells, ident.1 = case_type, ident.2 = "normal")
      c_data_for_drug <- data.frame(row.names = row.names(c_data), score = c_data$avg_log2FC,
                                    adj.P.Val = c_data$p_val_adj, P.Value = c_data$p_val)
      gene_list[[i]] <- c_data_for_drug
      c_names <- c(c_names, i)
    }
  }
  names(gene_list) <- c_names
  gene_list
}





# run the functions
# ====== MEF2D =========

gene_list_mef_limma <- limma_diff(sce_all_mef, case, control)
gene_list_mef_seurat <- seurat_diff(sce_all_mef, case_type)

readr::write_rds(gene_list_mef_limma, file = "MEF2D_genelist_limma.rds")
readr::write_rds(gene_list_mef_seurat, file = "MEF2D_genelist_seurat.rds")


# Drug repurposing
for (i in c("limma", "seurat")){
  gene_list <- readr::read_rds(paste0("MEF2D_genelist_", i, ".rds"))
  my_gene_info <- read.table(file = "haematopoietic-and-lymphoid-tissue_gene_info.txt", sep = "\t", header = T, quote = "")
  my_drug_info <- read.table(file = "haematopoietic-and-lymphoid-tissue_drug_info.txt", sep = "\t", header = T, quote = "")
  cmap.ref.profiles <- GetDrugRef(drug.response.path = 'haematopoietic-and-lymphoid-tissue_rankMatrix.txt',
                                  probe.to.genes = my_gene_info, drug.info = my_drug_info)
  Drug.ident.res <- GetDrug(gene.data = gene_list, drug.ref.profiles = cmap.ref.profiles, 
                            repurposing.unit = "drug", connectivity = "negative", drug.type = "FDA")
  readr::write_rds(Drug.ident.res, paste0("MEF2D_drugs_FDA_", i, ".rds"))
}


# drug score
sce_all_mef@meta.data$sample <- sce_all_mef$orig.ident
sce_all_mef@meta.data$celltype <- sce_all_mef$cell_types
GSE92742.gctx.path <- glue::glue("{lincs_dir}/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx")
GSE70138.gctx.path <- glue::glue("{lincs_dir}/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx")
Tissue <- "haematopoietic and lymphoid tissue"
for(i in c("limma", "seurat")){
  gene_list <- read_rds(paste0("MEF2D_genelist_", i, ".rds"))
  Drug.ident.res <- read_rds(paste0("MEF2D_drugs_FDA_", i, ".rds"))  
  Drug.score <- DrugScore(SC.integrated = sce_all_mef,
                          Gene.data = gene_list,
                          Cell.type = NULL,
                          Drug.data = Drug.ident.res,
                          FDA.drug.only = T,
                          Case = case,
                          Tissue = Tissue,
                          GSE92742.gctx = GSE92742.gctx.path,
                          GSE70138.gctx = GSE70138.gctx.path)
  readr::write_rds(Drug.score, file = paste0("MEF2D_drugscore_FDA_", i, ".rds"))
}


# ============== TCF3 =============

case_type = "TCF3"
case <- "E1"

gene_list_tcf_limma <- limma_diff(sce_all_tcf, case, control)
gene_list_tcf_seurat <- seurat_diff(sce_all_tcf, case_type)

readr::write_rds(gene_list_tcf_limma, file = "TCF3_genelist_limma.rds")
readr::write_rds(gene_list_tcf_seurat, file = "TCF3_genelist_seurat.rds")


# Drug repurposing
for (i in c("limma", "seurat")){
  gene_list <- readr::read_rds(paste0("TCF3_genelist_", i, ".rds"))
  my_gene_info <- read.table(file = "haematopoietic-and-lymphoid-tissue_gene_info.txt", sep = "\t", header = T, quote = "")
  my_drug_info <- read.table(file = "haematopoietic-and-lymphoid-tissue_drug_info.txt", sep = "\t", header = T, quote = "")
  cmap.ref.profiles <- GetDrugRef(drug.response.path = 'haematopoietic-and-lymphoid-tissue_rankMatrix.txt',
                                  probe.to.genes = my_gene_info, drug.info = my_drug_info)
  Drug.ident.res <- GetDrug(gene.data = gene_list, drug.ref.profiles = cmap.ref.profiles, 
                            repurposing.unit = "drug", connectivity = "negative", drug.type = "FDA")
  readr::write_rds(Drug.ident.res, paste0("TCF3_drugs_FDA_", i, ".rds"))
}


# drug score
sce_all_tcf@meta.data$sample <- sce_all_tcf$orig.ident
sce_all_tcf@meta.data$celltype <- sce_all_tcf$cell_types
for(i in c("limma", "seurat")){
  gene_list <- read_rds(paste0("TCF3_genelist_", i, ".rds"))
  Drug.ident.res <- read_rds(paste0("TCF3_drugs_FDA_", i, ".rds"))  
  Drug.score <- DrugScore(SC.integrated = sce_all_tcf,
                          Gene.data = gene_list,
                          Cell.type = NULL,
                          Drug.data = Drug.ident.res,
                          FDA.drug.only = T,
                          Case = case,
                          Tissue = Tissue,
                          GSE92742.gctx = GSE92742.gctx.path,
                          GSE70138.gctx = GSE70138.gctx.path)
  readr::write_rds(Drug.score, file = paste0("TCF3_drugscore_FDA_", i, ".rds"))
}



# ======= output figures ===========

output_plot <- function(type, method){}







Drug.score <- readr::read_rds(paste0("TCF3_drugscore_FDA_", "seurat", ".rds"))

#Drug score plot
Score.list<-data.frame(Patient = "E1", Drug = row.names(Drug.score),
                       DrugScore = Drug.score$Drug.therapeutic.score,
                       Pvalue = Drug.score$P.value, FDR = Drug.score$FDR)

Score.list <- Score.list[Score.list$Pvalue < 0.03, ]
Score.list$Drug <- Hmisc::capitalize(Score.list$Drug)
Drug.order <- Score.list$Drug[order(Score.list$DrugScore, decreasing = F)]
Score.list$Drug <- factor(Score.list$Drug, levels = Drug.order)

pdf(file = "plots/Drug_individual_TCF3_seurat.pdf", width = 3, height = 3)
ggplot(Score.list, aes(x = Patient, y = Drug, size = DrugScore, color = Patient)) +
  geom_point(alpha = 1) +
  scale_size(name = "DrugScore") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+ guides(color = FALSE)
dev.off()








# 

Drug.score <- readr::read_rds(paste0("MEF2D_drugscore_FDA_", "limma", ".rds"))
Combined.Drug.score <- readr::read_rds(paste0("MEF2D_drugscore_FDA_", "limma", ".rds"))
Gene.list <- readRDS("MEF2D_genelist_limma.rds")
Drug.ident.res <- readRDS("MEF2D_drugs_FDA_limma.rds")


Drug.score <- readr::read_rds(paste0("MEF2D_drugscore_FDA_", "seurat", ".rds"))
Combined.Drug.score <- readr::read_rds(paste0("MEF2D_drugscore_FDA_", "seurat", ".rds"))
Gene.list <- readRDS("MEF2D_genelist_seurat.rds")
Drug.ident.res <- readRDS("MEF2D_drugs_FDA_seurat.rds")


#Drug score plot
Score.list<-data.frame(Patient = "Overall", Drug = row.names(Drug.score),
                       DrugScore = Drug.score$Drug.therapeutic.score,
                       Pvalue = Drug.score$P.value, FDR = Drug.score$FDR)

j <- 0
for(i in case){
  j <- j + 1
  Drug.score <- DrugScore(SC.integrated = sce_all_mef,
                          Gene.data = Gene.list,
                          Cell.type = NULL,
                          Drug.data = Drug.ident.res,
                          FDA.drug.only = T,
                          Case = i,
                          Tissue = Tissue,
                          GSE92742.gctx = GSE92742.gctx.path,
                          GSE70138.gctx = GSE70138.gctx.path)
  Temp = data.frame(Patient = paste0("M", j + 1), Drug = row.names(Drug.score), 
                    DrugScore = Drug.score$Drug.therapeutic.score,Pvalue=Drug.score$P.value,FDR=Drug.score$FDR)
  Score.list = rbind(Score.list,Temp)
}

library(Hmisc)
Score.list$Drug<-capitalize(Score.list$Drug)
sig.list<-subset(Score.list,DrugScore>quantile(Score.list$DrugScore, 0.9,na.rm=T))
Drug.order<-row.names(Combined.Drug.score)[order(Combined.Drug.score$Drug.therapeutic.score,decreasing = F)]
Drug.order<-capitalize(Drug.order)
sig.list$Drug<-factor(sig.list$Drug,levels = Drug.order)

pdf(file = "Drug_individual_MEF2D_seurat.pdf",width = 3,height = 3)
ggplot(sig.list,aes(x=Patient, y=Drug, size=DrugScore, color=Patient)) +
  geom_point(alpha=1) +
  scale_size(name="DrugScore")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+ guides(color=FALSE)
dev.off()
