pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse",
          "CellChat", "ggplot2", "patchwork", "ggalluvial", "svglite", "Seurat")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

work_dir <- c('/cluster/home/projects/mef2d/analysis/zwn/home/cellchat/')
project_dir <- file.path(work_dir,"tumor_cell_types")
setwd(project_dir)

sce_tumor2 <- readr::read_rds("/cluster/home/projects/mef2d/data/sce_tumor_final.rds")

sce <- readr::read_rds("~/projects/mef2d/analysis/zwn/human/tenx/seurat/merged/mef2d.rds")

sce_hd <- subset(sce, `orig.ident` %notin% c("B069_Dx_CD19", "M2", "M3", "M4", "E1"))
sce_hd <- subset(sce_hd, `cell_types` %in% names(table(sce_hd$cell_types)))

dat_tumor <- sce_tumor2@assays$RNA@data
meta_tumor <- data.frame(all = sce_tumor2$cell_types_all, broad = sce_tumor2$cell_types_final_broad)


sce_hd$cell_types[sce_hd$cell_types == "imatureB"] <- "immatureB"


dat_hd <- sce_hd@assays$RNA@data
meta_hd <- data.frame(broad = sce_hd$cell_types)

cc_tumor <- createCellChat(object = dat_tumor, meta = meta_tumor, group.by = "broad")
cc_tumor@DB <- CellChatDB.human ## get DB data
cc_tumor <- subsetData(cc_tumor) ## filter for signaling genes
cc_tumor_all <- createCellChat(object = dat_tumor, meta = meta_tumor, group.by = "all")
cc_tumor_all@DB <- CellChatDB.human ## get DB data
cc_tumor_all <- subsetData(cc_tumor_all) ## filter for signaling genes


cc_hd <- createCellChat(object = dat_hd, meta = meta_hd, group.by = "broad")
cc_hd@DB <- CellChatDB.human ## get DB data
cc_hd <- subsetData(cc_hd) ## filter for signaling genes





# calculate cellchat
future::plan("multicore", workers = 20) ## do parallel
cc_tumor <- identifyOverExpressedGenes(cc_tumor)
cc_tumor <- identifyOverExpressedInteractions(cc_tumor)
cc_tumor <- projectData(cc_tumor, PPI.human)  ## PPI.human 4875x4875
cc_tumor <- computeCommunProb(cc_tumor)
cc_tumor <- filterCommunication(cc_tumor, min.cells = 10) ## filter groups under 10 cells

cc_hd <- identifyOverExpressedGenes(cc_hd)
cc_hd <- identifyOverExpressedInteractions(cc_hd)
cc_hd <- projectData(cc_hd, PPI.human)  ## PPI.human 4875x4875
cc_hd <- computeCommunProb(cc_hd)
cc_hd <- filterCommunication(cc_hd, min.cells = 10) ## filter groups under 10 cells

cc_tumor <- computeCommunProbPathway(cc_tumor)
cc_tumor <- aggregateNet(cc_tumor)
cc_hd <- computeCommunProbPathway(cc_hd)
cc_hd <- aggregateNet(cc_hd)

readr::write_rds(cc_tumor, "/cluster/home/projects/mef2d/output/cellchat_tumor_broad.rds")
readr::write_rds(cc_hd, "/cluster/home/projects/mef2d/output/cellchat_healthy.rds")


cc_tumor_all <- identifyOverExpressedGenes(cc_tumor_all)
cc_tumor_all <- identifyOverExpressedInteractions(cc_tumor_all)
cc_tumor_all <- projectData(cc_tumor_all, PPI.human)  ## PPI.human 4875x4875
cc_tumor_all <- computeCommunProb(cc_tumor_all)
cc_tumor_all <- filterCommunication(cc_tumor_all, min.cells = 10) ## filter groups under 10 cells
cc_tumor_all <- computeCommunProbPathway(cc_tumor_all)
cc_tumor_all <- aggregateNet(cc_tumor_all)

readr::write_rds(cc_tumor_all, "/cluster/home/projects/mef2d/output/cellchat_tumor_detailed.rds")









