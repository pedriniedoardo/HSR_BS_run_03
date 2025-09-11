# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)

# read in the final object ------------------------------------------------
scobj <- readRDS("../../data/data.combined_NOT_annotated_norm_fix_DoubletSinglet_harmonyMartina.rds")

meta <- read_tsv("../../data/GSE233295_metadata_allsample.tsv")
meta2 <- read_tsv("../../data/data.combined_NOT_annotated_fix_normalized_meta_AzimuthHumanCortexRef_resolution_DoubletSinglet_harmonyMartina.tsv")

list1 <- lapply(meta, function(x){
 x 
})

list2 <- lapply(meta2, function(x){
  x 
})

pmap(list(list1,list2),function(l1,l2){
  identical(l1,l2)
})

# export the normalized reads
GetAssayData(scobj,assay = "RNA",slot = "data") %>%
  saveRDS("../../out/object/normCount_BSrun1run2.rds")

# export the raw reads
GetAssayData(scobj,assay = "RNA",slot = "count") %>%
  saveRDS("../../out/object/rawCount_BSrun1run2.rds")

# pull only the data from the sample requested
scobj@meta.data$sample %>% table()

scobj_subset <- subset(scobj, sample == "CTRL8_SOX10_CTRL")
DimPlot(scobj_subset)

GetAssayData(scobj_subset,assay = "RNA",slot = "count") %>%
  saveRDS("../../out/object/rawCount_BSrun1run2_donor2_CTRL.rds")


scobj_subset2 <- subset(scobj, sample == "RR16_SOX10_CTRL")
DimPlot(scobj_subset2)

GetAssayData(scobj_subset2,assay = "RNA",slot = "count") %>%
  saveRDS("../../out/object/rawCount_BSrun1run2_donor3_CTRL.rds")

scobj_subset3 <- subset(scobj, sample == "CTRL4_SOX10_CTRL")
DimPlot(scobj_subset3)

GetAssayData(scobj_subset3,assay = "RNA",slot = "count") %>%
  saveRDS("../../out/object/rawCount_BSrun1run2_donor4_CTRL.rds")

DimPlot(scobj,split.by = "sample")
