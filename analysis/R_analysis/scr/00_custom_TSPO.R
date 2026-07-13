# AIM ---------------------------------------------------------------------
# plot GOI for the collaboration in HUMANITAS

# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(scales)
library(ComplexHeatmap)
library(lemon)
library(finalfit)
library(cowplot)
library(patchwork)
library(Nebulosa)
library(pals)
library(speckle)

# read in the dataset -----------------------------------------------------
data.combined <- readRDS("../../out/object/sobj_processed_donor.rds")
DimPlot(data.combined,label = T,raster = T,group.by = "expertAnno.l1")

# plotting ----------------------------------------------------------------
# Determine the number of clusters you need colors for
# RColorBrewer::display.brewer.all()
num_clusters <- length(unique(data.combined$expertAnno.l1))

# Generate a palette (e.g., "Set1" or "Paired")
# my_colors <- brewer.pal(n = num_clusters, name = "Paired")
color_id <- tableau20(num_clusters)
names(color_id) <- unique(data.combined$expertAnno.l1 %>% unlist())

show_col(color_id)

# Apply to DimPlot
p01 <- DimPlot(data.combined, 
               label = TRUE, 
               raster = TRUE, 
               group.by = "expertAnno.l1", 
               cols = color_id)
ggsave(plot = p01,filename = "../../out/image/00_UMAP_general_paletteUpdate.pdf",width = 6,height = 5)

# wrangling ---------------------------------------------------------------
# Sofia asked to focus only on the Fe cytokine and baseline
data.combined@meta.data$treat_full %>% table()

data.combined_subset <- subset(data.combined,subset = treat_full %in% c("BASELINE","cytokine","Fe"))
data.combined_subset

# make the barplot for the proporiton of cell types
df_meta_summary <- data.combined_subset@meta.data %>%
  # filter(!harmonized_donor2 %in% c("unassigned","doublet")) %>%
  filter(!harmonized_donor2 %in% c("unassigned")) %>%
  group_by(treat,expertAnno.l1) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot)

# save the summary table
df_meta_summary %>%
  write_tsv("../../out/table/00_prop_summary.tsv")

p02 <- df_meta_summary %>%
  ggplot(aes(x=treat,y=prop,fill=expertAnno.l1)) + geom_col() + scale_fill_manual(values = color_id) + theme_cowplot() + theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave(plot = p02,filename = "../../out/image/00_barplot_prop_paletteUpdate.pdf",width = 4,height = 5)

# run the test ------------------------------------------------------------
# run propeller
meta_test <- data.combined_subset@meta.data %>%
  # filter(!harmonized_donor2 %in% c("unassigned","doublet")) %>%
  filter(!harmonized_donor2 %in% c("unassigned"))

table(meta_test$orig.ident,meta_test$harmonized_donor2)
table(meta_test$treat_full,meta_test$harmonized_donor2)

table(paste0(meta_test$treat_full,meta_test$harmonized_donor2))

#
out_diagnosis <- propeller(clusters = meta_test$expertAnno.l1,
                           sample = paste0(meta_test$treat_full,meta_test$harmonized_donor2),
                           group = meta_test$treat_full)

out_diagnosis %>%
  rownames_to_column("expertAnno.l1") %>%
  write_tsv("../../out/table/00_propeller_expertAnno.l1_subsetSofi.tsv")

