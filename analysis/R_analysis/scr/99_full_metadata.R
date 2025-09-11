# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(scales)
library(ggrepel)
library(cowplot)

# read the data -----------------------------------------------------------
# read in the dataset
data.combined <- readRDS("../../out/object/sobj_processed_donor.rds")
DimPlot(data.combined,label = T,raster = T,group.by = "seurat_clusters")
DimPlot(data.combined,label = T,raster = T,group.by = "expertAnno.l1")
DimPlot(data.combined,label = T,raster = T,group.by = "treat_full")



# calculate the summary statistics ----------------------------------------
# calculate the total numbe of cells and the averange number of features anc counts per condition
test <- data.combined@meta.data %>%
  group_by(treat_full) %>%
  summarise(n_cell = n(),
            avg_features = mean(nFeature_RNA),
            tot_counts = sum(nCount_RNA))
write_tsv(test,"../../out/table/summary_meta.tsv")

test %>% mutate(tot_counts/n_cell)
test %>% mutate(tot_counts/n_cell*2)

test2 <- data.combined@meta.data %>%
  group_by(treat_full,expertAnno.l1) %>%
  summarise(n_cell = n())
write_tsv(test2,"../../out/table/summary_meta_celltype.tsv")

# calculate the number of unique genes ------------------------------------
dim(data.combined)

