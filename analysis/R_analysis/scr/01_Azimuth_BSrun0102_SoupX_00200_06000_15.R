# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)

# read in the dataset -----------------------------------------------------
# # read in the list of singlets
# list_query <- readRDS("../../out/object/list_datasc_fix_filter_norm_doubletSinglet_SoupX_00200_06000_15.rds")
# 
# list_diet <- lapply(list_query,function(x){
#   # track the process
#   print(x)
#   DefaultAssay(x) <- "RNA"
#   object_diet <- DietSeurat(object = x, assays = "RNA")
#   return(object_diet)
# })
# # save the diet versino of the objects
# saveRDS(list_diet,"../../out/object/diet.list_datasc_fix_filter_norm_doubletSinglet_SoupX_00200_06000_15.rds")

# available datasets ------------------------------------------------------
# load the reference that was built previously
reference <- LoadReference(path = "../../data/ref_BS_run_01_02/")

# # save the reference umap
# df_point <- reference$map@reductions$refUMAP@cell.embeddings %>%
#   data.frame() %>%
#   rownames_to_column("barcode")
# # write_tsv(df_point,"../../out/table/azimuth_ref_BSrun0102_CoordUMAP.tsv")
# 
# # save the meta of the ref
# df_meta <- reference$map@meta.data %>%
#   data.frame() %>%
#   rownames_to_column("barcode")
# # write_tsv(df_meta,"../../out/table/azimuth_ref_BSrun0102_Metadata.tsv")

# run azimuth -------------------------------------------------------------
list_diet <- readRDS("../../out/object/diet.list_datasc_fix_filter_norm_doubletSinglet_SoupX_00200_06000_15.rds")

query <- list_diet[[1]]
head(query@meta.data)

list_azimuth <- pmap(list(list_diet,names(list_diet)), function(x,name){
  # track the progress
  print(name)
  # The RunAzimuth function can take a Seurat object as input
  query <- RunAzimuth(query,reference = "../../data/ref_BS_run_01_02/")
  return(query)
})

saveRDS(list_azimuth,"../../out/object/list_out_Azimuth_BSrun0102_SoupX_00200_06000_15.rds")

# plotting the results ----------------------------------------------------
# read in the result of Azimuth annotation
list_azimuth <- readRDS("../../out/object/list_out_Azimuth_BSrun0102_SoupX_00200_06000_15.rds")


query <- list_azimuth[[1]]
# the metadata have all the info of interest
head(query@meta.data)

DimPlot(query, group.by = "predicted.annotation.l1", label = TRUE, label.size = 3) + NoLegend()
DimPlot(query, group.by = "predicted.annotation.l1", label = TRUE, label.size = 3,reduction = "ref.umap") + NoLegend()
DimPlot(query, group.by = "predicted.annotation.l1", label = TRUE, label.size = 3,split.by = "orig.ident") + NoLegend()

# # save the UMAP coordinates of the new reference
# query@reductions$ref.umap@cell.embeddings %>% 
#   data.frame() %>% 
#   rownames_to_column("barcode")
# # write_tsv("out/table/data.combined_NOT_annotated_fix_normalized_CoordUMAP_AzimuthHumanCortexRef_resolution_DoubletSinglet_harmonyMartina.tsv")
# 
# # save the new annotation from azimuth
# query@meta.data %>% 
#   data.frame() %>% 
#   rownames_to_column("barcode")
# # write_tsv("out/table/data.combined_NOT_annotated_fix_normalized_meta_AzimuthHumanCortexRef_resolution_DoubletSinglet_harmonyMartina.tsv")

# costume plots threshold 0.75 --------------------------------------------
# costume visualization
UMAP_query <- query@reductions$ref.umap@cell.embeddings %>% 
  data.frame() %>% 
  rownames_to_column("barcode")
meta_query <- query@meta.data %>% 
  data.frame() %>% 
  rownames_to_column("barcode") %>%
  mutate(robust_score_subclass = case_when(predicted.annotation.l1.score>0.75 & mapping.score>0.75 ~ predicted.annotation.l1,
                                           T~"uncertain"))

# define the levels of the cluster variable
level_annotation <- meta_query %>%
  group_by(predicted.annotation.l1) %>%
  summarise(med = median(predicted.annotation.l1.score)) %>%
  mutate(predicted.annotation.l1 = fct_reorder(predicted.annotation.l1, med,.desc = T)) %>%
  pull(predicted.annotation.l1) %>%
  levels()

# for each assignment what is the distribution of the scores
meta_query %>%
  mutate(predicted.annotation.l1 = factor(predicted.annotation.l1, levels = level_annotation)) %>%
  ggplot(aes(x=predicted.annotation.l1,y=predicted.annotation.l1.score))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitter(width = 0.2),alpha=0.01,size=0.1)+
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 1,angle = 45)) + 
  geom_hline(yintercept = 0.75,col="red")
# ggsave("out/image/data.combined_annotated_norm_fix_resolution_DoubletSinglet_human_motorcortexv1.0.0_subclass_score_075_harmonyMartina.pdf",height = 4,width = 4)

meta_query %>%
  mutate(predicted.annotation.l1 = factor(predicted.annotation.l1, levels = level_annotation)) %>%
  ggplot(aes(y=predicted.annotation.l1,x=predicted.annotation.l1.score))+
  ggridges::geom_density_ridges()+
  theme_bw() + 
  geom_vline(xintercept = 0.75,col="red",linetype="dashed")
# ggsave("out/image/data.combined_annotated_norm_fix_resolution_DoubletSinglet_human_motorcortexv1.0.0_subclass_score_ridges_075_harmonyMartina.pdf",height = 5,width = 4)

# identifyt he most likely assignment for each seurat cluster
# first using all the subcluster annotation, not filtering for threshold of scores
prop_table_subclass <- meta_query %>%
  group_by(seurat_clusters,predicted.annotation.l1) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters,predicted.annotation.l1,prop) %>%
  pivot_wider(names_from = seurat_clusters,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.annotation.l1")

# pdf("out/image/data.combined_annotated_norm_fix_resolution_DoubletSinglet_human_motorcortexv1.0.0_subclass_heatmapAll_harmonyMartina.pdf",height = 3,width = 5)
Heatmap(prop_table_subclass,
        name = "prop", 
        column_title = "subclass score all",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
# dev.off()

meta_query %>%
  filter(robust_score_subclass != "uncertain") %>%
  group_by(seurat_clusters,predicted.annotation.l1) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot)

# first using all the subcluster annotation, filtering for threshold of scores
prop_table_subclass_filter <- meta_query %>%
  # filter(robust_score_subclass != "uncertain") %>%
  filter(robust_score_subclass != "uncertain") %>%
  group_by(seurat_clusters,predicted.annotation.l1) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters,predicted.annotation.l1,prop) %>%
  pivot_wider(names_from = seurat_clusters,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.annotation.l1") %>%
  as.matrix()

# pdf("out/image/data.combined_annotated_norm_fix_resolution_DoubletSinglet_human_motorcortexv1.0.0_subclass_heatmapFilter_075_harmonyMartina.pdf",height = 3,width = 5)
Heatmap(prop_table_subclass_filter,
        name = "prop", 
        column_title = "subclass score high confidence",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
# dev.off()

meta_query %>%
  filter(robust_score_subclass != "uncertain") %>%
  group_by(seurat_clusters,predicted.annotation.l1) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters,predicted.annotation.l1,prop) %>%
  pivot_wider(names_from = seurat_clusters,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.annotation.l1")

# add the meta to the coordinates
data_query <- left_join(UMAP_query,meta_query,"barcode")
# add the groupign suggested by Martina
# data_query$group_martina <- case_when(data_query$orig.ident %in% c("pool_MOCK_CTRL")~"MOCK",
#                                       data_query$orig.ident %in% c("CTRL8_WT_CSF","CTRL8_WT_CTRL ")~"CTRL",
#                                       T~"SOX10")

# divide the dataset into uncertain and not
data_query_unc <- data_query %>%
  filter(robust_score_subclass == "uncertain")
#
data_query_defined <- data_query %>%
  filter(robust_score_subclass != "uncertain")

# average the position of the clusters
data_query_avg <- data_query_defined %>% group_by(robust_score_subclass) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)
data_query_avg2 <- data_query_defined %>% group_by(seurat_clusters) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# build the plot using both info
ggplot(label= TRUE)+
  geom_point(data = data_query_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_query_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score_subclass),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_query_avg,aes(x = UMAP_1,y = UMAP_2,label = robust_score_subclass),col="black")+theme_bw()
# ggsave("out/image/data.combined_annotated_norm_fix_resolution_DoubletSinglet_human_motorcortexv1.0.0_subclass_UMAP_075_harmonyMartina.pdf",width = 5,height = 3)

# build the plot using both info
ggplot(label= TRUE)+
  geom_point(data = data_query_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_query_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score_subclass),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_query_avg,aes(x = UMAP_1,y = UMAP_2,label = robust_score_subclass),col="black")+theme_bw()+facet_grid(~orig.ident)+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
# ggsave("out/image/data.combined_annotated_norm_fix_resolution_DoubletSinglet_human_motorcortexv1.0.0_subclass_split_UMAP_075_harmonyMartina.pdf",width = 20,height = 3)

# costume visualization
UMAP_ref <- reference$map@reductions$refUMAP@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column("barcode")

UMAP_ref2 <- query@reductions$ref.umap@cell.embeddings %>% 
  data.frame() %>% 
  rownames_to_column("barcode")

meta_ref <- reference$map@meta.data %>%
  data.frame() %>%
  rownames_to_column("barcode")

meta_query <- query@meta.data %>% 
  data.frame() %>% 
  rownames_to_column("barcode") %>%
  mutate(robust_score_subclass = case_when(predicted.annotation.l1.score>0.75 & mapping.score>0.75 ~ predicted.annotation.l1,
                                           T~"uncertain"))

# add the meta to the coordinates
data_ref <- left_join(UMAP_ref,meta_ref,"barcode")
data_ref2 <- left_join(UMAP_ref2,meta_query,"barcode") %>%
  mutate(seurat_clusters=factor(seurat_clusters))

# average the position of the clusters
data_ref_avg <- data_ref %>% group_by(annotation.l1) %>% select(refumap_1,refumap_2) %>% summarize_all(mean)

# average the position of the clusters
data_ref_avg2 <- data_ref2 %>%
  filter(robust_score_subclass != "uncertain") %>%
  group_by(seurat_clusters) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# build the plot using both info
ggplot(label= TRUE)+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref2 %>%
               filter(robust_score_subclass != "uncertain"),aes(x = UMAP_1,y = UMAP_2, col = seurat_clusters),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=2)) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_ref_avg2,aes(x = UMAP_1,y = UMAP_2,label = seurat_clusters),col="black")+theme_bw()
# ggsave("out/image/human_motorcortexv1.0.0_subclass_data.combined_annotated_norm_fix_resolution_UMAP_075_harmonyMartina.pdf",width = 9,height = 6)

ggplot(label= TRUE)+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref2 %>%
               filter(robust_score_subclass != "uncertain"),aes(x = UMAP_1,y = UMAP_2, col = seurat_clusters),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=2)) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_ref_avg2,aes(x = UMAP_1,y = UMAP_2,label = seurat_clusters),col="black")+theme_bw()+facet_wrap(~orig.ident)+theme(strip.background = element_blank())
# ggsave("out/image/human_motorcortexv1.0.0_subclass_data.combined_annotated_norm_fix_resolution_UMAP_075_split_harmonyMartina.pdf",width = 20,height = 22)

# average the position of the clusters
data_ref_avg2subclass <- data_ref2 %>%
  filter(robust_score_subclass != "uncertain") %>%
  group_by(predicted.annotation.l1) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# build the plot using both info
ggplot(label= TRUE)+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref2 %>%
               filter(robust_score_subclass != "uncertain"),aes(x = UMAP_1,y = UMAP_2, col = predicted.annotation.l1),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=1)) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_ref_avg2subclass,aes(x = UMAP_1,y = UMAP_2,label = predicted.annotation.l1),col="black")+theme_bw()
# ggsave("out/image/human_motorcortexv1.0.0_subclass_data.combined_annotated_norm_fix_resolution_UMAP_075_subclass_harmonyMartina.pdf",width = 8,height = 6)

ggplot(label= TRUE)+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref2 %>%
               filter(robust_score_subclass != "uncertain"),aes(x = UMAP_1,y = UMAP_2, col = predicted.annotation.l1),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=1)) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_ref_avg2subclass,aes(x = UMAP_1,y = UMAP_2,label = predicted.annotation.l1),col="black")+theme_bw()+facet_wrap(~orig.ident)+theme(strip.background = element_blank())
# ggsave("out/image/human_motorcortexv1.0.0_subclass_data.combined_annotated_norm_fix_resolution_UMAP_075_subclass_split_harmonyMartina.pdf",width = 20,height = 21)
