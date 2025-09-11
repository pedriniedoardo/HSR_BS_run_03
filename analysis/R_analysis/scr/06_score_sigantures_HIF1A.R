# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(scales)
library(ggrepel)
library(cowplot)
library(msigdbr)
library(UpSetR)

# read the data -----------------------------------------------------------
# read in the dataset
data.combined <- readRDS("../../out/object/sobj_processed_donor.rds")

# 
msigdbr_collections() %>%
  print(n=30)

#
gene_sets <- msigdbr(species = "Homo sapiens", category = "C2")
gene_sets %>%
  pull(gs_name) %>%
  unique() %>%
  str_subset(pattern = "HYPOXIA")
# REACTOME_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:REACTOME")
# get all the C2 terms
# cgp_gene_sets <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CGP")
# C2_gene_sets = msigdbr(species = "Mus musculus", category = "C2")
# h_gene_sets <- msigdbr(species = "Mus musculus", category = "H")
head(gene_sets)

# format in order to be accepted by GSEA
# subset only the signatures of interest
gene_sets_subset <- gene_sets %>%
  filter(gs_name %in% c("REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA",
                        "REACTOME_REGULATION_OF_GENE_EXPRESSION_BY_HYPOXIA_INDUCIBLE_FACTOR",
                        "BIOCARTA_P53HYPOXIA_PATHWAY",
                        "MANALO_HYPOXIA_UP"))

# generate the specific object
pathways <- split(x = gene_sets_subset$gene_symbol, f = gene_sets_subset$gs_name)
list_sig <- pmap(list(names(pathways),pathways),function(name,gene){
  df <- data.frame(Pathwway = name,Genes = gene)  
  return(df)
  }) %>%
  setNames(names(pathways))

# wrangling ---------------------------------------------------------------
# add the new classification to the metadata
meta <- data.combined@meta.data %>%
  rownames_to_column("barcodes")

# meta_full <- left_join(meta,LUT,by=c("official_id"))
meta_full <- meta

# add to the original dataset
# data.combined$treat_full <- meta_full$treat_full

# score the siganture in the UMAP -----------------------------------------
# run the enrichment for the signature. do it on the UMAP using the score siganatures
DefaultAssay(data.combined) <- "RNA"
# x <- "REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA"

# run the snippet over the whole signatures
lapply(names(list_sig),function(x){
  # extract the dataframe of genes in the signature
  signature.genes.df <- list_sig[[x]]
  
  # pull the genes
  signature.genes <- signature.genes.df %>%
    pull(Genes) %>%
    unique()
  
  # score the module
  data.combined <- AddModuleScore(data.combined,
                                  features = list(signature.genes),
                                  name="signature_score")
  
  # confirm the addition of the score for the module
  # data.combined@meta.data
  df_meta <- data.combined@meta.data %>%
    rownames_to_column("barcode") %>% 
    mutate(signature = x) 
  # mutate(treat_full = factor(treat_full,levels = c("control cortex","myelinated cortex","demyelinated cortex")))
  
  # save the table with the scores
  df_meta %>% 
    write_tsv(paste0("../../out/table/modules_HIF1A/06_Module_score_",x,".tsv"))
  
  # save the UMAP coordinates
  df_UMAP <- data.combined@reductions$umap@cell.embeddings %>%
    data.frame() %>%
    rownames_to_column("barcode")
  
  # data2 <- left_join(df_UMAP,df_meta,"barcode")
  # data2_avg <- data2 %>% group_by(expertAnno.l1) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)
  data2 <- left_join(df_UMAP,df_meta,"barcode")
  data2_avg <- data2 %>% group_by(expertAnno.l1) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)
  
  # data2 %>%
  #   arrange(signature_score1) %>%
  #   # mutate(gene = "Ptx3") %>%
  #   ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2,col = signature_score1),alpha = 0.5,size = 0.2) +
  #   geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = expertAnno.l1)) +
  #   facet_grid(~pathology) + theme_bw() + 
  #   theme(strip.background = element_blank()) +
  #   # scale_color_gradientn(colours = c("blue","gray","red"))
  #   scale_color_gradientn(colours = viridis::turbo(10))
  # # scale_color_gradientn(colours = c("blue","gray","red"), 
  # #                       values = rescale(c(-0.1,0,0.58)),
  # #                       guide = "colorbar", limits=c(-0.1,0.58))
  # ggsave(paste0("out/image/modules/UMAP_score_",x,".pdf"),width = 12,height = 3)
  
  data2 %>%
    arrange(signature_score1) %>%
    # mutate(gene = "Ptx3") %>%
    ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2,col = signature_score1),alpha = 0.5,size = 0.2) +
    geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = expertAnno.l1)) +
    facet_grid(~treat_full) + theme_bw() + 
    theme(strip.background = element_blank()) +
    # scale_color_gradientn(colours = c("blue","gray","red"))
    scale_color_gradientn(colours = viridis::turbo(10))
  ggsave(paste0("../../out/image/modules_HIF1A/06_UMAP_score_",x,".pdf"),width = 18,height = 3)
  ggsave(paste0("../../out/image/modules_HIF1A/06_UMAP_score_",x,".png"),width = 18,height = 3)
  
  data2 %>%
    arrange(signature_score1) %>%
    # mutate(gene = "Ptx3") %>%
    ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2,col = signature_score1),alpha = 0.5,size = 0.2) +
    geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = expertAnno.l1)) +
    facet_wrap(treat_full~harmonized_donor2) + theme_bw() + 
    theme(strip.background = element_blank()) +
    # scale_color_gradientn(colours = c("blue","gray","red"))
    scale_color_gradientn(colours = viridis::turbo(10))
  ggsave(paste0("../../out/image/modules_HIF1A/06_UMAP_score_split_",x,".pdf"),width = 32,height = 30)
  ggsave(paste0("../../out/image/modules_HIF1A/06_UMAP_score_split_",x,".png"),width = 32,height = 30)
  
  # # plot the score as distribution
  # data2 %>%
  #   # mutate(BraakStage=as.factor(BraakStage)) %>% 
  #   ggplot(aes(x=signature_score1,fill=disease))+geom_density(alpha=0.5)+
  #   facet_grid(origin~expertAnno.l1,scales = "free")+
  #   theme_bw()+
  #   theme(strip.background = element_blank(), 
  #         panel.border = element_rect(colour = "black", fill = NA))+
  #   scale_fill_manual(values = c("blue","red"))
  # ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/HVG4000/dist_score_",x,"_2.pdf"),width = 12,height =9)
  
  # plot the score as distribution but as ridges
  data2 %>%
    # mutate(BraakStage=as.factor(BraakStage)) %>% 
    ggplot(aes(x=signature_score1,y=treat_full,fill=treat_full))+
    ggridges::geom_density_ridges(alpha=0.5)+
    facet_wrap(~expertAnno.l1,scales = "free")+
    theme_bw()+
    theme(strip.background = element_blank(), 
          panel.border = element_rect(colour = "black", fill = NA))
  # scale_fill_manual(values = c("green","yellow","red"))
  ggsave(paste0("../../out/image/modules_HIF1A/06_dist_score_ridges_",x,".pdf"),width = 12,height =9)
  
  # # plot the distribution of the score per cluster
  # data2 %>%
  #   # mutate(BraakStage=as.factor(BraakStage)) %>% 
  #   arrange(signature_score1) %>%
  #   # mutate(gene = "Ptx3") %>%
  #   mutate(name_sample = paste0(treat_full,"_",orig.ident)) %>% 
  #   mutate(rank = rank(treat_full,ties.method = "min")) %>% 
  #   mutate(name_sample = fct_reorder(name_sample, rank,.desc = F)) %>% 
  #   
  #   ggplot(aes(x=name_sample,y = signature_score1,fill=treat_full)) + 
  #   geom_violin() +
  #   geom_boxplot(width=0.1,outlier.shape = NA) +
  #   facet_wrap(~expertAnno.l1) + theme_bw() + 
  #   # scale_fill_manual(values = c("green","yellow","red"))+
  #   theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))
  # ggsave(paste0("../../out/image/ManualClean/modules_SENESCENCE/HVG4000/violin_score_split_",x,".pdf"),width = 25,height = 9)
  
  # plot the distribution of the score per cluster
  data2 %>%
    # mutate(BraakStage=as.factor(BraakStage)) %>% 
    arrange(signature_score1) %>%
    # mutate(gene = "Ptx3") %>%
    ggplot(aes(x=treat_full,y = signature_score1,fill=treat_full)) + 
    geom_violin() +
    geom_boxplot(width=0.1,outlier.shape = NA) +
    facet_wrap(~expertAnno.l1) + theme_bw() + 
    # scale_fill_manual(values = c("green","yellow","red"))+
    theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))
  ggsave(paste0("../../out/image/modules_SENESCENCE/06_violin_score_",x,".pdf"),width = 9,height = 9)
  
  data2 %>%
    arrange(signature_score1) %>%
    # mutate(treat_full = factor(treat_full,levels = c("CX_Ctrl","CX_Mye","CX_Demye","WM_Ctrl","WM_NAWM","WM_CA","WM_CI","WM_Core"))) %>% 
    # mutate(gene = "Ptx3") %>%
    ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2,col = signature_score1),alpha = 0.5,size = 0.2) +
    # geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = expertAnno.l1)) +
    facet_wrap(~treat_full,nrow = 3) + theme_void() + 
    theme(strip.background = element_blank()) +
    scale_color_gradientn("sig score",
                          colours = c("gray","orange","red"),
                          oob = scales::squish,limits = c(0.1,0.6))
    # scale_color_gradientn("sig score",
    #                       colours = c("gray","orange","red"))
  # scale_color_gradientn(colours = c("gray","yellow","red"))
  # scale_color_gradientn(colours = viridis::turbo(10))
  ggsave(paste0("../../out/image/modules_HIF1A/06_UMAP_score_",x,"_tailored1.pdf"),width = 11,height = 10)
  ggsave(paste0("../../out/image/modules_HIF1A/06_UMAP_score_",x,"_tailored1.png"),width = 11,height = 10,bg="white")
  
  data2 %>%
    arrange(signature_score1) %>%
    # mutate(treat_full = factor(treat_full,levels = c("CX_Ctrl","CX_Mye","CX_Demye","WM_Ctrl","WM_NAWM","WM_CA","WM_CI","WM_Core"))) %>% 
    # mutate(gene = "Ptx3") %>%
    ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2,col = signature_score1),alpha = 0.5,size = 0.2) +
    # geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = expertAnno.l1)) +
    facet_wrap(~treat_full,nrow = 3) + theme_void() + 
    theme(strip.background = element_blank()) +
    scale_color_gradientn("sig score",
                          colours = viridis::turbo(10),limits = c(-0.1,0.6),oob = scales::squish)
  ggsave(paste0("../../out/image/modules_HIF1A/06_UMAP_score_",x,"_tailored2.pdf"),width = 11,height = 10)
  ggsave(paste0("../../out/image/modules_HIF1A/06_UMAP_score_",x,"_tailored2.png"),width = 11,height = 10,bg="white")
})

# tailored plotting -------------------------------------------------------
# x <- "senmayo"
# 
# signature.genes.df <- list_sig[[x]]
# 
# # pull the genes
# signature.genes <- signature.genes.df %>%
#   pull(Genes) %>%
#   unique()
# 
# # score the module
# data.combined <- AddModuleScore(data.combined,
#                                 features = list(signature.genes),
#                                 name="signature_score")
# 
# # confirm the addition of the score for the module
# # data.combined@meta.data
# df_meta <- data.combined@meta.data %>%
#   rownames_to_column("barcode") %>% 
#   mutate(signature = x) 
# # mutate(treat_full = factor(treat_full,levels = c("control cortex","myelinated cortex","demyelinated cortex")))
# 
# # save the UMAP coordinates
# df_UMAP <- data.combined@reductions$umap@cell.embeddings %>%
#   data.frame() %>%
#   rownames_to_column("barcode")
# 
# # data2 <- left_join(df_UMAP,df_meta,"barcode")
# # data2_avg <- data2 %>% group_by(expertAnno.l1) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)
# data2 <- left_join(df_UMAP,df_meta,"barcode")
# data2_avg <- data2 %>% group_by(expertAnno.l1) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)
# 
# # data2 %>%
# #   arrange(signature_score1) %>%
# #   # mutate(gene = "Ptx3") %>%
# #   ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2,col = signature_score1),alpha = 0.5,size = 0.2) +
# #   geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = expertAnno.l1)) +
# #   facet_grid(~pathology) + theme_bw() + 
# #   theme(strip.background = element_blank()) +
# #   # scale_color_gradientn(colours = c("blue","gray","red"))
# #   scale_color_gradientn(colours = viridis::turbo(10))
# # # scale_color_gradientn(colours = c("blue","gray","red"), 
# # #                       values = rescale(c(-0.1,0,0.58)),
# # #                       guide = "colorbar", limits=c(-0.1,0.58))
# # ggsave(paste0("out/image/modules/UMAP_score_",x,".pdf"),width = 12,height = 3)
# 
# # data2 %>%
# #   arrange(signature_score1) %>%
# #   mutate(treat_full = factor(treat_full,levels = c("CX_Ctrl","CX_Mye","CX_Demye","WM_Ctrl","WM_NAWM","WM_CA","WM_CI","WM_Core"))) %>% 
# #   # mutate(gene = "Ptx3") %>%
# #   ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2,col = signature_score1),alpha = 0.5,size = 0.2) +
# #   geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = expertAnno.l1)) +
# #   facet_wrap(~treat_full,nrow = 3) + theme_void() + 
# #   theme(strip.background = element_blank()) +
# #   scale_color_gradientn("sig score",
# #                         colours = viridis::turbo(10),
# #                         values = rescale(c(-0.1,0,0.3)),
# #                         oob = scales::squish,limits = c(-0.1,0.3))
# 
# data2 %>%
#   arrange(signature_score1) %>%
#   mutate(treat_full = factor(treat_full,levels = c("CX_Ctrl","CX_Mye","CX_Demye","WM_Ctrl","WM_NAWM","WM_CA","WM_CI","WM_Core"))) %>% 
#   # mutate(gene = "Ptx3") %>%
#   ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2,col = signature_score1),alpha = 0.5,size = 0.2) +
#   # geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = expertAnno.l1)) +
#   facet_wrap(~treat_full,nrow = 3) + theme_void() + 
#   theme(strip.background = element_blank()) +
#   scale_color_gradientn("sig score",
#                         colours = viridis::turbo(10),limits = c(-0.1,0.3),oob = scales::squish)
# ggsave(paste0("../../out/image/revision/modules_SENESCENCE/119_UMAP_score_",x,"_tailored2.pdf"),width = 11,height = 10)
# ggsave(paste0("../../out/image/revision/modules_SENESCENCE/119_UMAP_score_",x,"_tailored2.png"),width = 11,height = 10,bg="white")

# run Upset plot for the signatures ---------------------------------------
# load the siganture file
list_sig

# shortlist the signatures and rename them
list_sig_shortlist <- list_sig %>%
  bind_rows(.id = "signature") %>%
  # filter(signature %in% c("Induces","Inhibits","SAEPHIA_CURATED_SASP","senmayo")) %>% 
  # mutate(signature2 = case_when(signature == "Induces"~"CellAge Induces",
  #                               signature == "Inhibits"~"CellAge Inhibits",
  #                               signature == "SAEPHIA_CURATED_SASP"~"Curated SASP",
  #                               signature == "senmayo"~"Senmayo")) %>%
  # split(f = .$signature2) %>%
  split(f = .$signature) %>%
  lapply(function(x){
    x %>% pull(Genes) %>%
      unique()
  })


library(UpSetR)
pdf("../../out/image/modules_HIF1A/06_upset_signatures_shortlist_HIF1A.pdf",width = 10,height = 6,onefile = F)
upset(fromList(list_sig_shortlist), order.by = "freq",nsets = 100,nintersects = 100)
dev.off()

# pull the intersections
str(list_sig_shortlist)

df1 <- lapply(list_sig_shortlist,function(x){
  data.frame(gene = x)
}) %>% 
  bind_rows(.id = "path")

head(df1)
df2 <- data.frame(gene=unique(unlist(list_sig_shortlist)))

head(df2)

# now loop through each individual gene and pick the list of all the intersections they belong to
df_int <- lapply(df2$gene,function(x){
  # pull the name of the intersections
  intersection <- df1 %>% 
    dplyr::filter(gene==x) %>% 
    arrange(path) %>% 
    pull("path") %>% 
    paste0(collapse = "|")
  
  # build the dataframe
  data.frame(gene = x,int = intersection)
}) %>% 
  bind_rows()

head(df_int,n=20)

# save the intersection table
df_int %>%
  write_tsv("../../out/table/modules_HIF1A/06_upset_signatures_shortlist_intersection_HIF1A.tsv")

# confirm the data and the list are congruent
df_int %>% 
  group_by(int) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n))
