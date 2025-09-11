# AIM ---------------------------------------------------------------------
# run GSEA on the table of DE from pbulk analysis
# this is a task for Aletta's project, she is interested in exploring the effect of Fe and Myelin treatments in the MG subset. I am using the GO:BP annotaiton
# use as ranking the signed pvalue as recommended here:
# https://github.com/crazyhottommy/compbio_tutorials/blob/main/scripts/07_gene_set_enrichment_RNAseq.Rmd

# library -----------------------------------------------------------------
library(tidyverse)
library(fgsea)
# install.packages('msigdbdf', repos = 'https://igordot.r-universe.dev')
# library(msigdbdf)
library(msigdbr)
library(GSEABase)
library(patchwork)

# prepare the dataset with all the annoration needed ---------------------- 
results <- read_tsv("../../out/table/DE_treatvsBASELINE_pseudobulk_MG_shr.tsv") %>%
  split(f = .$conditionVsBASELINE)

# read in the list of terms provided by Aletta
TOI <- read_csv("../../data/GOterms_iron_myelin.csv")

# GSEA -------------------------------------------------------------------- 
# use the FC dataset to create the ranked list of genes 
# Symbol or Entrez? 
# x <- results$`30_DE_pseudobulk_filterExp_shr_IMMUNE`
list_ranks <- lapply(results, function(x){
  
  x <- dplyr::filter(x,!is.na(symbol)) %>%
    dplyr::filter(!is.na(pvalue)) %>%
    group_by(symbol) %>%
    # average pvalue in case of duplicated genenames
    summarise(pvalue = mean(pvalue),
              log2FC = mean(log2FoldChange)) %>%
    ungroup() %>%
    # arrange(pvalue) %>%
    # change the inf to big numbers
    # mutate(signed_rank_stats = sign(log2FC) * -log10(pvalue)) %>%
    mutate(negative_log10pvalue = -log10(pvalue)) %>%
    # multiply the 1000 by the log2FC to break the ranking ties
    # mutate(negative_log10pvalue = ifelse(is.infinite(negative_log10pvalue), 1000 * log2FC, negative_log10pvalue)) %>%
    mutate(negative_log10pvalue = ifelse(is.infinite(negative_log10pvalue), 1000, negative_log10pvalue)) %>%
    mutate(signed_rank_stats = sign(log2FC) * negative_log10pvalue)
  
  
  ranks <- setNames(x$signed_rank_stats, x$symbol)
  ranks
}) 
glimpse(list_ranks)

# score all the signatures in MsigDB from C2 category ---------------------
# pull all the annotations
msigdbr_collections() %>% print(n=30)

# extrac the annotion of interest
gene_sets <- msigdbr(species = "Homo sapiens", category = "C5",subcategory = "GO:BP")
head(gene_sets)

# check that the terms provided by Aletta are present in the current dataset
geneset_summary <- gene_sets %>%
  group_by(gs_exact_source,gs_name) %>%
  summarise(n = n())

geneset_summary

# check that the terms provided are present in the current dataset
TOI %>%
  left_join(geneset_summary,by = c("GO_id" = "gs_exact_source")) %>%
  mutate(test = !is.na(n)) %>%
  print(n = 30)

# not all the terms are present

# format in order to be accepted by GSEA
pathways <- split(x = gene_sets$gene_symbol, f = gene_sets$gs_name)
# head(pathways)

# RUN GSEA ----------------------------------------------------------------
# add a seed to fix the GSEA result
set.seed(123)

df_tables_GSEA_all <- lapply(list_ranks, function(x){
  fgsea(pathways, x, minSize=10, maxSize=500)  
}) %>%
  bind_rows(.id = "dataset") %>% 
  # the ladingEdge columns has to be re-arranged in order to save the file as a table (originally is a list) 
  mutate(leadingEdge = unlist(lapply(.$leadingEdge, function(x){
    paste0(x,collapse = "|")
  }))) %>%
  arrange(padj,-abs(NES)) 

dim(df_tables_GSEA_all)

head(df_tables_GSEA_all,n=20) 

# save the whole table
df_tables_GSEA_all %>%
  write_tsv(paste0("../../out/table/23_df_table_GSEA_GOBP_rankSigLogPval_","MG",".tsv"))

# check senescent terms
df_tables_GSEA_all %>%
  filter(str_detect(pathway,"SEN"))

# COLLAPSE REDUNDANT ------------------------------------------------------
# collapsing the similar pathways 

# split the dataset per type
list_tables_GSEA_all <- split(df_tables_GSEA_all,f = df_tables_GSEA_all$dataset)

names(list_tables_GSEA_all)

list_collapsedPathways <- lapply(names(list_tables_GSEA_all),function(x){
  collapsePathways(list_tables_GSEA_all[[x]], pathways, list_ranks[[x]])
}) %>%
  setNames(names(list_tables_GSEA_all))

str(list_collapsedPathways)

list_mainPathways <- pmap(list(list_tables_GSEA_all,list_collapsedPathways),function(x,y){
  x %>%
    dplyr::filter(pathway %in% y$mainPathways) %>%
    arrange(padj,-abs(NES)) %>%
    pull(pathway) 
})

str(list_mainPathways)

# save list of non redundant terms
# chackt the order of the names is the same
sum(!names(list_tables_GSEA_all) == names(list_mainPathways))

# filter only the non redundant fro each comparison
df_tables_GSEA_all_non_redundant <- 
  pmap(list(list_tables_GSEA_all,list_mainPathways),function(x,y){
    x %>%
      dplyr::filter(pathway %in% y)
  }) %>%
  bind_rows()

# save the table
df_tables_GSEA_all_non_redundant %>%
  write_tsv(file = paste0("../../out/table/23_df_table_GSEA_GOBP_nonredundant_rankSigLogPval_","MG",".tsv"))

test <- df_tables_GSEA_all_non_redundant %>%
  group_by(dataset) %>%
  top_n(wt = padj*(-1),n = 5)

# test plot to show the main terms in each dataset
library(ggrepel)
df_tables_GSEA_all_non_redundant %>%
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = "KEGG_") %>%
           str_sub(start = 1,end = 35)) %>%
  # mutate(min_log10_padj = -log10(padj)) %>%
  ggplot(aes(y = -log(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() +
  theme(strip.background = element_blank())+
  geom_text_repel(size = 2,box.padding = 0.5,segment.alpha = 0.6,max.overlaps = 10)+
  geom_hline(yintercept = -log(0.05),col="gray",linetype="dashed")
ggsave(paste0("../../out/image/23_GSEA_unbiased_GOBP_nonredundant_rankSigLogPval_","MG",".pdf"),width = 15,height = 15)

# library(ggrepel)
df_tables_GSEA_all %>%
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = "KEGG_") %>%
           str_sub(start = 1,end = 35)) %>%
  # mutate(min_log10_padj = -log10(padj)) %>%
  ggplot(aes(y = -log(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() +
  theme(strip.background = element_blank())+
  geom_text_repel(size = 2,box.padding = 0.5,segment.alpha = 0.6,max.overlaps = 10)+
  geom_hline(yintercept = -log(0.05),col="gray",linetype="dashed")
ggsave(paste0("../../out/image/23_GSEA_unbiased_GOBP_rankSigLogPval_","MG",".pdf"),width = 15,height = 15)

# # PLOT PROFILE ------------------------------------------------------------ 
# # library("patchwork")
# plot_pathway <- lapply(list_mainPathways$res_GMPvsRESEARCH_shr[1:9],function(x){ 
#   name <- str_sub(x,start = 1,end = -4) 
#   # pdf(file = paste0("HUVEC_BMP9_24h_",name,".pdf"),width = 3,height = 3) 
#   plotEnrichment(pathways[[x]], list_ranks$res_GMPvsRESEARCH_shr) + labs(title = name) 
# })
# 
# wrap_plots(plot_pathway)
# ggsave(filename = "../../out/image/GSEA_plot_profile_GMPvsRESEARCH_shr_nonredundant_KEGG.pdf",width = 15,height = 10)
# # wrap_plots(plot_pathway) + ggsave(filename = "image/GSEA_plot_profile_CGP.pdf",width = 15,height = 10)
# # wrap_plots(plot_pathway) + ggsave(filename = "image/GSEA_plot_profile_H.pdf",width = 15,height = 10)


# compare the two results from the two ranking systems --------------------
#
df_tables_GSEA_rankPvalue <- read_tsv(paste0("../../out/table/23_df_table_GSEA_GOBP_rankSigLogPval_MG.tsv")) %>%
  dplyr::select(dataset,pathway,NES,padj,leadingEdge)

df_tables_GSEA_rankLogFC <- read_tsv(paste0("../../out/table/23_df_table_GSEA_GOBP_ranklogfc_MG.tsv")) %>%
  dplyr::select(dataset,pathway,NES,padj,leadingEdge)

# join the two tables
test <- df_tables_GSEA_rankPvalue %>%
  left_join(df_tables_GSEA_rankLogFC,by = c("dataset","pathway"),suffix = c(".pvalue",".logfc")) %>%
  # add more annotations to the dataset
  left_join(geneset_summary,by = c("pathway" = "gs_name"))

# the two ranking methods do not produce exactly the same results, therefore there might be interest in checking both outputs.
test %>%  
  ggplot(aes(x=NES.pvalue, y = NES.logfc)) + geom_point(shape = 1,alpha = 0.2) + theme_bw()

# check the stats from the subset fo terms provided by Aletta
test %>%
  filter(gs_exact_source %in% TOI$GO_id) %>%
  filter(dataset %in% c("myelin_shr","Fe_shr")) %>%
  dplyr::select(dataset,pathway,NES.pvalue,padj.pvalue,leadingEdge.pvalue) %>%
  arrange(padj.pvalue)

test %>%
  filter(gs_exact_source %in% TOI$GO_id) %>%
  filter(dataset %in% c("myelin_shr","Fe_shr")) %>%
  dplyr::select(dataset,pathway,NES.logfc,padj.logfc,leadingEdge.logfc) %>%
  arrange(padj.logfc)

# plot the top 20 for Fe an Myelin
list_plot <- lapply(c("myelin_shr","Fe_shr"),function(dat){
  global_min <- min(-log(c(test %>%
                             filter(dataset == dat) %>%
                             slice_max(abs(NES.logfc),n = 20) %>%
                             pull(padj.logfc),
                           test %>%
                             filter(dataset == dat) %>%
                             slice_max(abs(NES.pvalue),n = 20) %>%
                             pull(padj.pvalue))))
  global_max <- max(-log(c(test %>%
                             filter(dataset == dat) %>%
                             slice_max(abs(NES.logfc),n = 20) %>%
                             pull(padj.logfc),
                           test %>%
                             filter(dataset == dat) %>%
                             slice_max(abs(NES.pvalue),n = 20) %>%
                             pull(padj.pvalue))))
  
  
  p1 <- test %>%
    filter(dataset == dat) %>%
    slice_max(abs(NES.pvalue),n = 20) %>%
    mutate(Term = str_remove_all(pathway,pattern = "GOBP_") %>% str_sub(start = 1,end = 50)) %>%
    mutate(Term = fct_reorder(Term, NES.pvalue,.desc = F)) %>%
    mutate(direction = factor(sign(NES.pvalue),labels = c(1,-1),levels = c(1,-1))) %>%
    ggplot(aes(x = NES.pvalue,y=Term)) + 
    geom_point(aes(size = -log(padj.pvalue),col=direction)) +
    theme_bw() +
    geom_vline(xintercept = 0,col="gray",linetype = "dashed") +
    scale_color_manual(values = c("red","blue")) +
    scale_size_continuous(limits = c(global_min, global_max)) +
    coord_cartesian(xlim = c(-3,3)) + 
    ggtitle(paste0("GSEA ", dat," signed NegLogPval ranking"))
  
  p2 <- test %>%
    filter(dataset == dat) %>%
    slice_max(abs(NES.logfc),n = 20) %>%
    mutate(Term = str_remove_all(pathway,pattern = "GOBP_") %>% str_sub(start = 1,end = 50)) %>%
    mutate(Term = fct_reorder(Term, NES.logfc,.desc = F)) %>%
    mutate(direction = factor(sign(NES.logfc),labels = c(1,-1),levels = c(1,-1))) %>%
    ggplot(aes(x = NES.logfc,y=Term)) + 
    geom_point(aes(size = -log(padj.logfc),col=direction)) +
    theme_bw() +
    geom_vline(xintercept = 0,col="gray",linetype = "dashed") +
    scale_color_manual(values = c("red","blue")) +
    scale_size_continuous(limits = c(global_min, global_max)) +
    coord_cartesian(xlim = c(-3,3)) + 
    ggtitle(paste0("GSEA ", dat," logFC ranking"))
  
  p1 / p2
  
})

wrap_plots(list_plot)
ggsave("../../out/image/23_GSEA_top20_GOBP_Fe_Myelin.pdf",width = 20,height = 10)

# non rudundant -----------------------------------------------------------
#
df_tables_GSEA_rankPvalue2 <- read_tsv(paste0("../../out/table/23_df_table_GSEA_GOBP_nonredundant_rankSigLogPval_MG.tsv")) %>%
  dplyr::select(dataset,pathway,NES,padj,leadingEdge)

df_tables_GSEA_rankLogFC2 <- read_tsv(paste0("../../out/table/23_df_table_GSEA_GOBP_nonredundant_ranklogfc_MG.tsv")) %>%
  dplyr::select(dataset,pathway,NES,padj,leadingEdge)

# join the two tables
test2 <- df_tables_GSEA_rankPvalue2 %>%
  left_join(df_tables_GSEA_rankLogFC2,by = c("dataset","pathway"),suffix = c(".pvalue",".logfc")) %>%
  # add more annotations to the dataset
  left_join(geneset_summary,by = c("pathway" = "gs_name"))

# the two ranking methods do not produce exactly the same results, therefore there might be interest in checking both outputs.
test2 %>%  
  ggplot(aes(x=NES.pvalue, y = NES.logfc)) + geom_point(shape = 1,alpha = 0.2) + theme_bw()

# check the stats from the subset fo terms provided by Aletta
test2 %>%
  filter(gs_exact_source %in% TOI$GO_id) %>%
  filter(dataset %in% c("myelin_shr","Fe_shr")) %>%
  dplyr::select(dataset,pathway,NES.pvalue,padj.pvalue,leadingEdge.pvalue) %>%
  arrange(padj.pvalue)

test2 %>%
  filter(gs_exact_source %in% TOI$GO_id) %>%
  filter(dataset %in% c("myelin_shr","Fe_shr")) %>%
  dplyr::select(dataset,pathway,NES.logfc,padj.logfc,leadingEdge.logfc) %>%
  arrange(padj.logfc)

read_tsv(paste0("../../out/table/23_df_table_GSEA_GOBP_nonredundant_ranklogfc_MG.tsv")) %>%
  # filter(pathway %in% c("GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS"))
  filter(pathway %in% c("GOBP_RIBOSOME_BIOGENESIS"))

# plot the top 20 for Fe an Myelin
list_plot2 <- lapply(c("myelin_shr","Fe_shr"),function(dat){
  global_min <- min(-log(c(test2 %>%
                             filter(dataset == dat) %>%
                             slice_max(abs(NES.logfc),n = 20) %>%
                             pull(padj.logfc),
                           test2 %>%
                             filter(dataset == dat) %>%
                             slice_max(abs(NES.pvalue),n = 20) %>%
                             pull(padj.pvalue))))
  global_max <- max(-log(c(test2 %>%
                             filter(dataset == dat) %>%
                             slice_max(abs(NES.logfc),n = 20) %>%
                             pull(padj.logfc),
                           test2 %>%
                             filter(dataset == dat) %>%
                             slice_max(abs(NES.pvalue),n = 20) %>%
                             pull(padj.pvalue))))
  
  
  p1 <- test2 %>%
    filter(dataset == dat) %>%
    slice_max(abs(NES.pvalue),n = 20) %>%
    mutate(Term = str_remove_all(pathway,pattern = "GOBP_") %>% str_sub(start = 1,end = 50)) %>%
    mutate(Term = fct_reorder(Term, NES.pvalue,.desc = F)) %>%
    mutate(direction = factor(sign(NES.pvalue),labels = c(1,-1),levels = c(1,-1))) %>%
    ggplot(aes(x = NES.pvalue,y=Term)) + 
    geom_point(aes(size = -log(padj.pvalue),col=direction)) +
    theme_bw() +
    geom_vline(xintercept = 0,col="gray",linetype = "dashed") +
    scale_color_manual(values = c("red","blue")) +
    scale_size_continuous(limits = c(global_min, global_max)) +
    coord_cartesian(xlim = c(-3,3)) + 
    ggtitle(paste0("GSEA ", dat," signed NegLogPval ranking"))
  
  p2 <- test2 %>%
    filter(dataset == dat) %>%
    slice_max(abs(NES.logfc),n = 20) %>%
    mutate(Term = str_remove_all(pathway,pattern = "GOBP_") %>% str_sub(start = 1,end = 50)) %>%
    mutate(Term = fct_reorder(Term, NES.logfc,.desc = F)) %>%
    mutate(direction = factor(sign(NES.logfc),labels = c(1,-1),levels = c(1,-1))) %>%
    ggplot(aes(x = NES.logfc,y=Term)) + 
    geom_point(aes(size = -log(padj.logfc),col=direction)) +
    theme_bw() +
    geom_vline(xintercept = 0,col="gray",linetype = "dashed") +
    scale_color_manual(values = c("red","blue")) +
    scale_size_continuous(limits = c(global_min, global_max)) +
    coord_cartesian(xlim = c(-3,3)) + 
    ggtitle(paste0("GSEA ", dat," logFC ranking"))
  
  p1 / p2
  
})

wrap_plots(list_plot2)
# ggsave("../../out/image/23_GSEA_top20_GOBP_Fe_Myelin.pdf",width = 20,height = 10)
