# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)

# read in the final object ------------------------------------------------
scobj <- readRDS("../../out/object/sobj_processed_donor.rds")

DimPlot(scobj,label = T,raster = T,group.by = "seurat_clusters")
DimPlot(scobj,label = T,raster = T,group.by = "expertAnno.l1")
DimPlot(scobj,label = T,raster = T,group.by = "treat_full")
DimPlot(scobj,label = T,raster = T,group.by = "expertAnno.l1",split.by = "treat_full")

# wrangling ---------------------------------------------------------------
# explore the full meta
scobj@meta.data

# subset the dataset only for the new samples that Martina indicated
sobj_new <- subset(scobj,subset = orig.ident %in% c("hBS_CTR4_MG",
                                                    "hBS_RR16_MG",
                                                    "hBS_RR24_MG",
                                                    "hBS_RR25_MG",
                                                    "W8_24h_CSF-controls_plus_untreated_multiplexed",
                                                    "W8_24h_CSF-MS_plus_untreated_multiplexed",
                                                    "W8_24h_cytokines",
                                                    "W8_48h_CSF-MS_multiplexed"))

DimPlot(sobj_new,label = T,raster = T,group.by = "expertAnno.l1")
DimPlot(sobj_new,label = T,raster = T,group.by = "expertAnno.l1",split.by = "treat_full")

dim(sobj_new)
# export the metadata for the new samples
df_meta <- sobj_new@meta.data %>%
  rownames_to_column("barcodes")

df_meta %>%
  write_tsv("../../out/table/metadata_GEOsubmission.tsv")

df_meta %>%
  group_by(orig.ident) %>%
  summarise(n = n()) %>%
  pull(n) %>%
  sum()

scobj@meta.data %>%
  group_by(orig.ident) %>%
  summarise(n = n()) %>% 
  # pull(orig.ident) %>%
  dplyr::filter(orig.ident %in% c("hBS_CTR4_MG",
                                  "hBS_RR16_MG",
                                  "hBS_RR24_MG",
                                  "hBS_RR25_MG",
                                  "W8_24h_CSF-controls_plus_untreated_multiplexed",
                                  "W8_24h_CSF-MS_plus_untreated_multiplexed",
                                  "W8_24h_cytokines",
                                  "W8_48h_CSF-MS_multiplexed")) %>%
  pull(n) %>%
  sum()

# export the normalized reads
GetAssayData(sobj_new,assay = "RNA",slot = "data") %>%
  saveRDS("../../out/object/normCount_GEOsubmission.rds")

# export the raw reads
GetAssayData(sobj_new,assay = "RNA",slot = "count") %>%
  saveRDS("../../out/object/rawCount_GEOsubmission.rds")
