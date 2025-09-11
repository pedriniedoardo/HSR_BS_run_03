# try the same with propeller on the same data
# renv::install("phipsonlab/speckle")
# renv::install("statmod")
library(speckle)
library(limma)
library(statmod)
library(cowplot)
library(ggrepel)
library(finalfit)
library(Seurat)
library(tidyverse)

# read in the object ------------------------------------------------------
scobj <- readRDS("../../out/object/sobj_processed_donor.rds")

DimPlot(scobj,label = T,raster = T,group.by = "seurat_clusters")
DimPlot(scobj,label = T,raster = T,group.by = "expertAnno.l1")
DimPlot(scobj,label = T,raster = T,group.by = "treat_full")
DimPlot(scobj,label = T,raster = T,group.by = "expertAnno.l1",split.by = "treat_full")

# wrangling ---------------------------------------------------------------
# martina wanted to focus on the 
scobj_subset <- subset(scobj, treat_full %in% c("BASELINE","CSF.ctrl.24h","CSF.MS.24h","CSF.MS.48h","cytokine"))
DimPlot(scobj_subset,label = T,raster = T,group.by = "expertAnno.l1")
DimPlot(scobj_subset,label = T,raster = F,group.by = "expertAnno.l1",split.by = "treat_full")

# run the test ------------------------------------------------------------
# run propeller
meta_test <- scobj_subset@meta.data
# filter(treat != "CSF.MS_RAPA")

table(meta_test$orig.ident,meta_test$harmonized_donor2)
table(meta_test$treat_full,meta_test$harmonized_donor2)

#
out_diagnosis <- propeller(clusters = meta_test$expertAnno.l1,
                           sample = paste0(meta_test$treat_full,meta_test$harmonized_donor2),
                           group = meta_test$treat_full)

out_diagnosis %>%
  rownames_to_column("expertAnno.l1") %>%
  write_tsv("../../out/table/07_propeller_expertAnno.l1.tsv")

# plotting diagnosis ------------------------------------------------------
# default plot
speckle::plotCellTypeProps(x = scobj_subset,
                           clusters = scobj_subset$expertAnno.l1,
                           sample = scobj_subset$treat_full)+theme_minimal()+theme(panel.grid = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/image/07_plot_propeller_expertAnno.l1.pdf",height = 5,width = 5)

# custom plot
df_summary_diagnosis <- meta_test %>% 
  mutate(group_id = paste0(treat_full,harmonized_donor2)) %>%
  group_by(expertAnno.l1,
           group_id,
           treat_full) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  group_by(group_id) %>% 
  mutate(tot = sum(n),
         prop = n/tot) %>%
  mutate(clone = str_remove_all(group_id,"BASELINE|CSF.ctrl.24h|CSF.MS.24h|CSF.MS.48h|cytokine"))

# plot 01
df_summary_diagnosis %>%
  ggplot(aes(x=treat_full,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),shape=1,alpha =0.7,aes(col=clone))+
  facet_wrap(~expertAnno.l1,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/image/07_plot_propeller_expertAnno.l1_replicates.pdf",width = 10,height = 10)

df_summary_diagnosis %>% filter(expertAnno.l1=="NEU") %>% filter(treat_full == "cytokine")

df_summary_diagnosis %>% filter(treat_full == "cytokine") %>% filter(clone=="unassigned")

df_summary_diagnosis %>% filter(clone=="unassigned") %>% print(n = 30)

# plot 02
df_summary_diagnosis %>%
  ggplot() +
  geom_boxplot(aes(x=RNA_snn_res.0.2,y=prop,color=treat),outlier.shape = NA) +
  geom_point(aes(x=RNA_snn_res.0.2,y=prop,color=treat),position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha=0.7) +
  theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_y_sqrt()
# ggsave("../../out/plot/manualClean/propeller_plot02_diagnosis_cellid.pdf",width = 8,height = 5)
