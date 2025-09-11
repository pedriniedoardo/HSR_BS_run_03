# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(SeuratData)
library(ggridges)
library(ComplexHeatmap)

# read in the data --------------------------------------------------------
# # in this case I want to use the Leng 2021 dataset as reference to map the cells of martina's dataset
# # read in the reference dataset for the leng dataset
# ref_BS <- readRDS("/beegfs/scratch/ric.cosr/pedrini.edoardo/scRNAseq_Brainsphere_Absinta/BS_drop_RR25CTRL/data/data.combined_NOT_annotated_norm_fix_DoubletSinglet_harmonyMartina.rds")
# # add the model to the obejct if it is missing
# ref_BS <- RunUMAP(object = ref_BS,reduction = "harmony",dims = 1:10, return.model = TRUE)
# 
# # save ccordinates fro reference of the UMAP and metadata
# ref_BS_meta <- left_join(ref_BS@reductions$umap@cell.embeddings %>%
#                            data.frame() %>%
#                            rownames_to_column("barcodes"),
#                          ref_BS@meta.data %>%
#                            data.frame() %>%
#                            rownames_to_column("barcodes"),by="barcodes")
# 
# write_tsv(ref_BS_meta,file = "../../out/table/ref_BS_meta.tsv")
# DimPlot(ref_BS,label = T)
# # save the reference object as defined above
# saveRDS(ref_BS,"../../out/object/ref_BS_labelTransfer.rds")

ref_BS <- readRDS("../../out/object/ref_BS_labelTransfer.rds")

# loop the label transfer procedure ---------------------------------------
# read in the list of singlets
list_query <- readRDS("../../out/object/list_datasc_fix_filter_norm_doubletSinglet_SoupX_01000_06000_15.rds")

# obj <- list_query[[1]]
# obj_name <- names(list_query)[1]
list_out <- pmap(list(list_query,names(list_query)),function(obj,obj_name){
  print(obj_name)
  
  # read in the query dataset
  query <- obj
  # fix the meta to add the general cell annotation
  query$seurat_clusters_fix1 <- paste0("cluster_",query$seurat_clusters)
  
  # -------------------------------------------------------------------------
  # to work the reference should be the integrated dataset
  DefaultAssay(ref_BS)<-"integrated"
  DefaultAssay(query) <- "RNA"
  test.anchors_0 <- FindTransferAnchors(reference = ref_BS,
                                        query = query,
                                        dims = 1:30,
                                        reference.reduction = "pca")
  
  predictions_0 <- TransferData(anchorset = test.anchors_0,
                                refdata = ref_BS$seurat_clusters,
                                dims = 1:30)
  
  # add the predictions ot the query dataset
  query_transfer <- AddMetaData(query, metadata = predictions_0)
  
  # -------------------------------------------------------------------------
  # add the meta to the coordinates
  data_transfer <- left_join(
    # get the coordinates of the query dataset
    query_transfer@reductions$umap@cell.embeddings %>%
      data.frame() %>%
      rownames_to_column(var = "barcodes"),
    
    # get the metadata from the query dataset, after adding the prediciton information. use a threshold for the score
    query_transfer@meta.data %>%
      rownames_to_column(var = "barcodes") %>%
      # this is based on martina's table, I need to ask her how she provided the accuracy of the imputation
      mutate(robust_score = case_when(prediction.score.max>0.70~predicted.id,
                                      T ~ "uncertain"))
    ,"barcodes")
  
  # Unimodal UMAP Projection ------------------------------------------------
  # ref_BS <- RunUMAP(object = ref_BS,reduction = "harmony",dims = 1:10, return.model = TRUE)
  query_transfer <- MapQuery(anchorset = test.anchors_0,
                             reference = ref_BS,
                             query = query_transfer,
                             refdata = list(celltype = "seurat_clusters"),
                             reference.reduction = "pca", reduction.model = "umap")
  
  # add the meta to the coordinates
  data2_transfer <- left_join(
    # get the coordinates of the query dataset
    query_transfer@reductions$ref.umap@cell.embeddings %>%
      data.frame() %>%
      rownames_to_column(var = "barcodes"),
    
    # get the metadata from the query dataset, after adding the prediciton information. use a threshold for the score
    query_transfer@meta.data %>%
      rownames_to_column(var = "barcodes") %>%
      # this is based on martina's table, I need to ask her how she provided the accuracy of the imputation
      mutate(robust_score = case_when(prediction.score.max>0.70~predicted.id,
                                      T ~ "uncertain"))
    ,"barcodes")
  
  # save the metadata with the new clsuter informations
  return(data2_transfer)
}) %>% 
  setNames(names(list_query))

saveRDS(list_out,"../../out/object/list_out_label_transfer_BSrun0102_SoupX_01000_06000_15.rds")

# # read in the reference metadata for the run 01 and run 02
# df_ref <- read_tsv(file = "../../out/table/ref_BS_meta.tsv")
# 
# data2_transfer %>%
#   ggplot(label= TRUE) +
#   # reference layer
#   geom_point(data = df_ref %>% dplyr::select(UMAP_1,UMAP_2),aes(x=UMAP_1,y=UMAP_2),col="gray",alpha=0.3,size=0.1)+
#   # transfer layer
#   geom_point(data = data2_transfer,aes(x = refUMAP_1,y = refUMAP_2, col = predicted.id),size=0.3,alpha=0.8) +
#   # labs(color= "Clusters") +
#   theme_bw() +
#   guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))+theme(strip.background = element_blank())
