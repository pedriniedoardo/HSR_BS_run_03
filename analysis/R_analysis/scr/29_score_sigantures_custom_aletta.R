# AIM ---------------------------------------------------------------------
# this is the updated script for the generation of the modules scores for aletta's custom signatures

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

# load the siganture file
list_sig <- readRDS("../../out/object/28_custom_signatures.rds")

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
# x <- "ferroptosis_all"

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
    write_tsv(paste0("../../out/table/modules_custom/29_Module_score_",x,".tsv"))
  
  # save the UMAP coordinates
  df_UMAP <- data.combined@reductions$umap@cell.embeddings %>%
    data.frame() %>%
    rownames_to_column("barcode")
  
  # data2 <- left_join(df_UMAP,df_meta,"barcode")
  # data2_avg <- data2 %>% group_by(expertAnno.l1) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)
  data2 <- left_join(df_UMAP,df_meta,"barcode")
  data2_avg <- data2 %>% group_by(expertAnno.l1) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)
  
  data2 %>%
    arrange(signature_score1) %>%
    # mutate(gene = "Ptx3") %>%
    ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2,col = signature_score1),alpha = 0.5,size = 0.2) +
    geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = expertAnno.l1)) +
    facet_grid(~treat_full) + theme_bw() + 
    theme(strip.background = element_blank()) +
    # scale_color_gradientn(colours = c("blue","gray","red"))
    scale_color_gradientn(colours = viridis::turbo(10))
  ggsave(paste0("../../out/image/modules_custom/29_UMAP_score_",x,".pdf"),width = 18,height = 3)
  ggsave(paste0("../../out/image/modules_custom/29_UMAP_score_",x,".png"),width = 18,height = 3)
  
  data2 %>%
    arrange(signature_score1) %>%
    # mutate(gene = "Ptx3") %>%
    ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2,col = signature_score1),alpha = 0.5,size = 0.2) +
    geom_text_repel(data = data2_avg,aes(x = UMAP_1, y = UMAP_2,label = expertAnno.l1)) +
    facet_wrap(treat_full~harmonized_donor2) + theme_bw() + 
    theme(strip.background = element_blank()) +
    # scale_color_gradientn(colours = c("blue","gray","red"))
    scale_color_gradientn(colours = viridis::turbo(10))
  ggsave(paste0("../../out/image/modules_custom/29_UMAP_score_split_",x,".pdf"),width = 32,height = 30)
  ggsave(paste0("../../out/image/modules_custom/29_UMAP_score_split_",x,".png"),width = 32,height = 30)
  
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
  ggsave(paste0("../../out/image/modules_custom/29_dist_score_ridges_",x,".pdf"),width = 12,height =9)
  
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
  ggsave(paste0("../../out/image/modules_custom/29_violin_score_",x,".pdf"),width = 9,height = 9)
  
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
  ggsave(paste0("../../out/image/modules_custom/29_UMAP_score_",x,"_tailored1.pdf"),width = 11,height = 10)
  ggsave(paste0("../../out/image/modules_custom/29_UMAP_score_",x,"_tailored1.png"),width = 11,height = 10,bg="white")
  
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
  ggsave(paste0("../../out/image/modules_custom/29_UMAP_score_",x,"_tailored2.pdf"),width = 11,height = 10)
  ggsave(paste0("../../out/image/modules_custom/29_UMAP_score_",x,"_tailored2.png"),width = 11,height = 10,bg="white")
})

# custom plotting ---------------------------------------------------------

# load the module score for the ferroptosis_all custom signature
df_modules <- read_tsv("../../out/table/modules_custom/29_Module_score_ferroptosis_all.tsv") %>%
  mutate(treat_full = factor(treat_full,levels = c("BASELINE","CSF.ctrl.24h","CSF.MS.24h","CSF.MS.48h","cytokine","Fe","myelin","TBHP")))

# define the threshold per cell type for the signature
df_90 <- df_modules %>% 
  filter(signature %in% "ferroptosis_all",
         treat_full %in% "BASELINE") %>% 
  group_by(signature,expertAnno.l1) %>% 
  summarise(thr = quantile(signature_score1,prob=0.90))

# show the distributino of the score per cell type and condition
df_modules %>%
  # focus only on MG and include only baseline Fe and myelin 
  filter(expertAnno.l1 %in% c("MG"),
         treat_full %in% c("BASELINE","Fe","myelin","cytokine")) %>%
  # mutate(BraakStage=as.factor(BraakStage)) %>% 
  ggplot(aes(x=signature_score1,y=treat_full))+
  ggridges::geom_density_ridges(alpha=0.5)+
  facet_grid(signature~expertAnno.l1,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA)) +
  # add the threshold for the specific cell type
  geom_vline(data = df_90 %>% filter(expertAnno.l1 %in% c("MG")), aes(xintercept=thr),linetype = "dashed",col="red")
# scale_fill_manual(values = c("green","yellow","red"))
ggsave(paste0("../../out/image/modules_custom/29_dist_score_ridges_ferroptosis_all_tailored.pdf"),width = 5,height =5)
