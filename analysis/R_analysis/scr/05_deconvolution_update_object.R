# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(scales)
library(ComplexHeatmap)
library(limma)
library(multtest)
library(metap)
library(ggbreak)
library(patchwork)
# library(lemon)
library(future)
library(ggnewscale)

#  read in the data -------------------------------------------------------
test <- readRDS("../../out/object/data.combined_harmonySkipIntegration_AllSoupX_01000_06000_15.rds")
DefaultAssay(test)
DimPlot(test,raster = T,label = T)

# # defirne the pattern in the barcodes
# pattern_remove <- table(test@meta.data$orig.ident) |> 
#   names() |>
#   paste0("_") |> 
#   paste0(collapse = "|")

# extract the full meta from the scobj
meta <- test@meta.data |> 
  rownames_to_column("full_barcode")

# define the table of the sample identity after donor devonvolution
# martina also asked to add another treat varaible that divide the CSF.MS 24 from the CSF.MS 48 treatmentnt
LUT_sample <- data.frame(id_deconvolution = c("hBS_CTR4_MG","W8_24h_CSF-controls_plus_untreated_multiplexed","W8_24h_CSF-MS_plus_untreated_multiplexed","W8_48h_CSF-MS_multiplexed","W8_24h_cytokines","W8_6h_Fe_RSL3","W8_18h_myelin_plus_untreated_multiplexed","W8_48h_TBHP_and_TBHP-dasatinib_multiplexed","hBS_RR16_MG","hBS_RR24_MG","hBS_RR25_MG"),
  id_sample_scRNAseq = c("CTR4_BASELINE_0h","mix_CSF.ctrl_24h","mix_CSF.MS_24h","mix_CSF.MS_48h","mix_cytokine_24h","mix_Fe_6h","mix_myelin_18h","mix_TBHP_48h","RR16_BASELINE_0h","RR24_BASELINE_0h","RR25_BASELINE_0h"),
  id_sample_short = c("CTR4","CSFctrl","CSFms24h","CSFms48h","cytokine","Fe","myelin","TBHP","RR16","RR24","RR25"),
  treat_full = c("BASELINE","CSF.ctrl.24h","CSF.MS.24h","CSF.MS.48h","cytokine","Fe","myelin","TBHP","BASELINE","BASELINE","BASELINE"))

# From the classification approach built the lut for the donor identity. Use a common donor ID
LUT_donor <- data.frame(id_donor_sample = c("donor1.TBHP","donor4.RR24","donor2.CSFms48h","donor0.cytokine","donor2.Fe","donor2.CSFms24h","donor0.CSFctrl","donor1.CSFctrl","donor0.Fe","donor0.RR25","donor1.cytokine","donor1.CSFms48h","donor0.CSFms24h","donor0.TBHP","donor1.myelin","donor1.Fe","donor2.TBHP","donor0.CSFms48h","donor2.cytokine","donor1.CSFms24h","donor5.RR16","donor2.CSFctrl","donor2.CTR4","donor0.myelin"),
                        harmonized_donor = c(rep("donRR24",7),
                                             rep("donRR25",8),
                                             rep("donRR16",9))) |> 
  separate(id_donor_sample,into = c("donor","id_sample_short"),sep = "\\.",remove = F)

# pull the deconvolution results from vireo
# read in one run of genotyping
df_donorIds <- lapply(LUT_sample$id_deconvolution,function(x){
  file <- paste0("../genetic_demux/results/02_vireo/",x,"/donor_ids.tsv")
  read_tsv(file)
}) %>% 
  setNames(LUT_sample$id_sample_short) |> 
  bind_rows(.id = "id_sample_short") |> 
  dplyr::select(id_sample_short,cell,donor_id)

# merge the info
meta_full <- left_join(df_donorIds,LUT_sample,by="id_sample_short") |> 
  # build the join variable
  mutate(id_donor_sample = paste0(donor_id,".",id_sample_short)) |> 
  # join the summary from the LUT_donor
  left_join(LUT_donor,by = c("id_donor_sample","id_sample_short")) |> 
  # build the barcode column
  mutate(full_barcode = paste0(id_deconvolution,"_",cell)) |> 
  # if the hamonyzed donor is missing pull the category from the donor_id imputation
  mutate(harmonized_donor2 = case_when(is.na(harmonized_donor)~donor_id,
                   T~harmonized_donor)) |> 
  # from the pure samples do not allow for more than one donor for the baseline
  mutate(harmonized_donor2 = case_when(id_deconvolution == "hBS_CTR4_MG"~"donRR16",
                                       id_deconvolution == "hBS_RR16_MG"~"donRR16",
                                       id_deconvolution == "hBS_RR24_MG"~"donRR24",
                                       id_deconvolution == "hBS_RR25_MG"~"donRR25",
                                       T~harmonized_donor2))


# add the expert annotation to the meta -----------------------------------
LUT_annotation <- data.frame(seurat_clusters = as.factor(0:17)) |> 
  mutate(expertAnno.l1 = case_when(seurat_clusters %in% c(1,2,9,10,17)~"NEU",
                                   seurat_clusters %in% c(6,13,16)~"MG",
                                   seurat_clusters %in% c(7)~"OLIGO",
                                   seurat_clusters %in% c(5)~"OPC",
                                   seurat_clusters %in% c(8,15)~"PROG",
                                   seurat_clusters %in% c(3,11,14)~"GLIA_IMM",
                                   seurat_clusters %in% c(0,4,12)~"ASTRO"))


# add the full meta to the sc meta table
meta_full2 <- meta |> 
  left_join(meta_full,by = c("full_barcode")) |> 
  left_join(LUT_annotation,by = c("seurat_clusters"))

# confrim the correcto identity in the full dataset
table(meta_full2$expertAnno.l1,meta_full2$seurat_clusters)

# confirm the dimension of the metadata
dim(meta_full2)
dim(meta)

# swap the metadata in the sc object
test@meta.data <- meta_full2 |> 
  mutate(rowname = full_barcode) |> 
  column_to_rownames("rowname")

# count the number of cells per harmonized donors
meta_full2 |> 
  group_by(harmonized_donor2,orig.ident) |> 
  summarise(n = n()) |> 
  pivot_wider(names_from = harmonized_donor2,values_from = n)

# martina wanted to have an idea fo the relaitve proportion of doublet and unassigned per sample
meta_full2 |> 
  group_by(harmonized_donor2,orig.ident) |> 
  summarise(n = n()) |> 
  group_by(orig.ident) |> 
  mutate(tot_sample = sum(n)) |> 
  ungroup() |> 
  mutate(prop = n/tot_sample) |>
  ggplot(aes(x=orig.ident,y=prop,fill=harmonized_donor2))+
  geom_col()+
  theme_void()+
  theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 45), plot.margin = margin(0.5, 0.5, 2, 2, "cm"),axis.text.y = element_text())
ggsave("../../out/image/sobj_processed_donor_proportions.pdf",width = 8,height = 6)

# plot it as a scatter plot
meta_full2 |> 
  group_by(harmonized_donor2,orig.ident) |> 
  summarise(n = n()) |> 
  group_by(orig.ident) |> 
  mutate(tot_sample = sum(n)) |> 
  ungroup() |> 
  mutate(prop = n/tot_sample) |>
  ggplot(aes(x = orig.ident,
             y = prop,col = harmonized_donor2,
             group = harmonized_donor2))+
  geom_point() +
  geom_line() +
  theme_bw() +
  facet_grid(harmonized_donor2~"prop over sample") +
  theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 45),
        plot.margin = margin(0.5, 0.5, 2, 2, "cm"),
        strip.background = element_blank())+
  # scale_y_sqrt(breaks = seq(from=0,to=1,by=0.1))
  scale_y_sqrt(breaks = c(0,0.1,0.2,0.4,0.6,0.8,1))
ggsave("../../out/image/sobj_processed_donor_proportions_scatter.pdf",width = 8,height = 10)

# confirm the swap of the metadata
DimPlot(test,split.by = "id_sample_short",raster = T,group.by = "harmonized_donor2")+facet_wrap(~id_sample_short,nrow=2)
DimPlot(test,split.by = "harmonized_donor2",raster = T,group.by = "id_sample_short")+facet_wrap(~harmonized_donor2,nrow=2)
DimPlot(test,split.by = "harmonized_donor2",raster = T,group.by = "treat_full")+facet_wrap(~harmonized_donor2,nrow=2)

# meta_full2 |> 
#   group_by(id_sample_short,harmonized_donor2) |> 
#   summarise(n = n()) |> 
#   mutate(tot = sum(n)) |> 
#   mutate(prop = n/tot) |> 
#   ggplot(aes(x=id_sample_short,y=prop,fill=harmonized_donor2))+geom_col(position = "dodge")+theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# ggsave("../../out/image/donor_prop.pdf",width = 6,height = 3)
# -------------------------------------------------------------------------
# save the object with the updated matadata
saveRDS(test,"../../out/object/sobj_processed_donor.rds")
