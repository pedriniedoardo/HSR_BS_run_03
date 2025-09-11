# # LIBRARIES ---------------------------------------------------------------
# library(scater)
# library(Seurat)
# library(tidyverse)
# library(robustbase)
# # library(SeuratData)
# library(patchwork)
# 
# # read in the data --------------------------------------------------------
# # read in the metadata for the samples
# LUT <- read_csv("../../data/LUT_samples.csv")
# 
# # in this case wa are going to use the fix threshold filtered data
# df_file <- dir("../../out/object/") %>%
#   str_subset(pattern = c("datasc_fix_filter_norm_doublet_Singlet.*_W8_18h_myelin_plus_untreated_multiplexed.rds")) %>%
#   data.frame(sample_file = .)
#   # mutate(sample_id = str_remove_all(sample_file,pattern = "datasc_fix_filter_norm_doublet_|_cr_61_5K|_cr_61_SoupX|.rds")) %>%
#   # left_join(.,LUT,by = c("sample_id"="ID"))
# 
# # set them up in a list
# data.list <- lapply(df_file$sample_file, function(x){
#   readRDS(paste0("../../out/object/",x))
# }) %>%
#   setNames(df_file$sample_id)
# 
# # save the list of features for the integration
# features <- SelectIntegrationFeatures(object.list = data.list)
# combined.anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features)
# 
# # this command creates an 'integrated' data assay
# data.combined <- IntegrateData(anchorset = combined.anchors)
# 
# # fix the meta
# # data.combined$treat <- factor(data.combined$treat,levels = c("CTRL","CSF"))
# # data.combined$infection <- factor(data.combined$infection,levels = c("WT","SOX10"))
# 
# # specify that we will perform downstream analysis on the corrected data note that the
# # original unmodified data still resides in the 'RNA' assay
# # # Run the standard workflow for visualization and clustering
# # data.combined <- ScaleData(data.combined, verbose = FALSE)
# # add the cell cycle analysis
# DefaultAssay(data.combined) <- "RNA"
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# data.combined <- CellCycleScoring(data.combined, s.features = s.genes, g2m.features = g2m.genes)
# # head(data.combined@meta.data)
# DefaultAssay(data.combined) <- "integrated"
# 
# # Run the standard workflow for visualization and clustering
# # data.combined <- ScaleData(data.combined,vars.to.regress = c("percent.mt","nCount_RNA"))
# data.combined <- ScaleData(data.combined,vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"))
# # data.combined <- ScaleData(data.combined,vars.to.regress = c("percent.mt","nCount_RNA"))
# # data.combined <- ScaleData(data.combined)
# data.combined <- RunPCA(data.combined, npcs = 30, verbose = FALSE)
# data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:30)
# data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:30)
# data.combined <- FindClusters(data.combined, resolution = 0.2)
# 
# # wrangling ---------------------------------------------------------------
# meta_full <- data.combined@meta.data %>% 
#   rownames_to_column("barcodes") %>% 
#   separate(col = "barcodes",into = c("barcode_id","sample_id"),sep = "_",remove = F) %>% 
#   mutate(sample_id_fix = case_when(sample_id=="1"~"02_60",
#                                    sample_id=="2"~"06_60",
#                                    T~"10_60"))
# 
# # add the sample id
# data.combined$sample_id <- meta_full$sample_id_fix
# DimPlot(data.combined,split.by = "sample_id")
# 
# saveRDS(data.combined,"../../out/object/test_myelin_integration.rds")
