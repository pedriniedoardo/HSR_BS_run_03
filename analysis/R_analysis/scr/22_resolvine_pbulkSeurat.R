# # AIM ---------------------------------------------------------------------
# # the aim of this test is to compare seurat's implementaiton of pseudobulk analysis, to the default process of DGE
# # the reference of the test is presented here: https://satijalab.org/seurat/articles/de_vignette.html
# # this second script run the comparison on the DGE
# 
# # libraries ---------------------------------------------------------------
# library(Seurat)
# library(SeuratData)
# library(tidyverse)
# library(harmony)
# library(ggExtra)
# library(ComplexUpset)
# library(cowplot)
# library(UpSetR)
# 
# # read in the dataset -----------------------------------------------------
# sobj <- readRDS("../../out/object/sobj_processed_donor.rds")
# 
# # filter out the doublet and the unassigned donors
# sobj_filter <- subset(sobj,subset = harmonized_donor2 %in% c("donRR16","donRR24","donRR25"))
# 
# # check the object version
# class(sobj@assays$RNA)
# 
# # load in the genes sets to check
# set_01 <- c("FPR2", "CMKLR1","GPR32","GPR18","GPR37","LGR6","LTB4R","RORA")
# set_02 <- c("ABCC1", "AKR1C3", "ALOX12", "ALOX15", "ALOX5", "CBR1", "CYP1A2", "CYP2C8", "CYP2C9", "CYP2D6", "CYP2E1", "CYP3A4", "CYP4F2", "CYP4F3", "DPEP1", "DPEP2", "DPEP3", "EPHX1", "EPHX2", "EPHX3", "GGT1", "GGT2", "GGT5", "GPX2", "GPX4", "GSTM4", "HPGDS", "LTA4H", "LTC4S", "PTGDS", "PTGES", "PTGES2", "PTGIS", "PTGS1", "PTGS2", "TBXAS1", "TXN")
# 
# # test DGE at pseudobulk --------------------------------------------------
# # pseudobulk the counts based on donor-condition-celltype
# pseudo_sobj <- AggregateExpression(sobj, assays = "RNA", return.seurat = T, group.by = c("treat_full", "harmonized_donor2", "expertAnno.l1"),slot = "counts")
# 
# # confirm the count slot is correctly populated, notice that the slot contains integers
# pseudo_sobj@assays$RNA@counts
# pseudo_sobj@assays$RNA@data
# 
# # add the covariate for the stimulation per cell type
# pseudo_sobj$sample_id <- rownames(pseudo_sobj@meta.data)
# pseudo_sobj$treat <- pseudo_sobj$sample_id %>% str_split(pattern = "_") %>% map(function(x){x[1]}) %>% unlist()
# pseudo_sobj$donor <- pseudo_sobj$sample_id %>% str_split(pattern = "_") %>% map(function(x){x[2]}) %>% unlist()
# pseudo_sobj$cellid <- pseudo_sobj$sample_id %>% str_split(pattern = "_") %>% map(function(x){x[3]}) %>% unlist()
# pseudo_sobj$celltype.stim <- paste(pseudo_sobj$cellid, pseudo_sobj$treat, sep = "_")
# 
# # explore the dimensionality of the new dataset
# pseudo_sobj@meta.data %>%
#   filter(cellid %in% c("MG")) %>%
#   group_by(celltype.stim) %>%
#   summarise(n = n())
# 
# # run the DGE over the same cell type for stim vs ctrl. This time the unist are the pseudobulks per donor/celltype/stimulus
# Idents(pseudo_sobj) <- "celltype.stim"
# table(pseudo_sobj$treat)
# 
# # treat_id <- "Fe"
# 
# list_bulk_de <- lapply(c("CSF.ctrl.24h","CSF.MS.24h","CSF.MS.48h","cytokine","Fe","myelin","TBHP"),function(treat_id){
#   print(treat_id)
#   
#   bulk.mono.de <- FindMarkers(object = pseudo_sobj, 
#                               ident.1 = paste0("MG_",treat_id), 
#                               ident.2 = "MG_BASELINE",
#                               test.use = "DESeq2")
#   bulk.mono.de %>%
#     rownames_to_column("gene") %>%
#     mutate(test = paste0(treat_id,"_vs_BASELINE")) %>%
#     mutate(cell_id = "MG")
# })
# 
# # make it a df
# df_bulk_de <- list_bulk_de %>%
#   bind_rows()
# 
# head(df_bulk_de)
# 
# # save the table
# df_bulk_de %>%
#   write_tsv("../../out/table/22_DE_pbulkSeurat_MG.tsv")
# 
# # do the same for AST
# list_bulk_de2 <- lapply(c("CSF.ctrl.24h","CSF.MS.24h","CSF.MS.48h","cytokine","Fe","myelin","TBHP"),function(treat_id){
#   print(treat_id)
#   
#   bulk.mono.de <- FindMarkers(object = pseudo_sobj, 
#                               ident.1 = paste0("ASTRO_",treat_id), 
#                               ident.2 = "ASTRO_BASELINE",
#                               test.use = "DESeq2")
#   bulk.mono.de %>%
#     rownames_to_column("gene") %>%
#     mutate(test = paste0(treat_id,"_vs_BASELINE")) %>%
#     mutate(cell_id = "ASTRO")
# })
# 
# # make it a df
# df_bulk_de2 <- list_bulk_de2 %>%
#   bind_rows()
# 
# head(df_bulk_de2)
# 
# # save the table
# df_bulk_de2 %>%
#   write_tsv("../../out/table/22_DE_pbulkSeurat_ASTRO.tsv")
# 
# 
# # volcano -----------------------------------------------------------------
# plot_volcano <- df_bulk_de %>%
#   # add a clor variable in case significant
#   mutate(col=ifelse(((p_val_adj<0.05)&abs(avg_log2FC)>1&!is.na(gene)),yes = 1,no = 0)) %>%
#   filter(!is.na(col))
# 
# plot_volcano %>%
#   ggplot(aes(x=avg_log2FC,y=-log(p_val_adj)))+
#   geom_point(alpha=0.5,shape = 1) +
#   # geom_point(data = plot_volcano2[plot_volcano2$col==0,],aes(x=avg_log2FC,y=-log(p_val_adj),col=factor(col)),alpha=0.05)+
#   # geom_point(data = plot_volcano2[plot_volcano2$col==1,],aes(x=avg_log2FC,y=-log(p_val_adj),col=factor(col)),alpha=0.5)+
#   geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
#   geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
#   # scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
#   # ggrepel::geom_text_repel(
#   #   data = plot_volcano[plot_volcano$col==1,][1:1000,],
#   #   aes(label = symbol),max.overlaps = 1,segment.alpha=0.4,
#   #   size = 2,
#   #   box.padding = unit(0.35, "lines"),
#   #   point.padding = unit(0.3, "lines")) +
#   # ggrepel::geom_text_repel(
#   #   data = plot_volcano %>% 
#   #     group_by(conditionVsCTRL) %>% 
#   #     arrange(padj) %>% 
#   #     dplyr::slice(1:30) %>% 
#   #     dplyr::filter(abs(log2FoldChange)>1) %>% 
#   #     dplyr::filter(padj<0.05),
#   #   aes(label = symbol),segment.alpha=0.4) +
#   facet_wrap(~test)+
#   theme_bw()+
#   theme(strip.background = element_blank())+
#   theme(legend.position = "none")
# # ggsave("../out/plot/102_vulcano_plot_pseudobulk_filterExp.pdf",width = 10,height = 10)
# 
# 
# 
# 
# plot_volcano2 <- df_bulk_de2 %>%
#   # add a clor variable in case significant
#   mutate(col=ifelse(((p_val_adj<0.05)&abs(avg_log2FC)>1&!is.na(gene)),yes = 1,no = 0)) %>%
#   filter(!is.na(col))
# 
# plot_volcano2 %>%
#   ggplot(aes(x=avg_log2FC,y=-log(p_val_adj)))+
#   geom_point(alpha=0.5,shape = 1) +
#   # geom_point(data = plot_volcano2[plot_volcano2$col==0,],aes(x=avg_log2FC,y=-log(p_val_adj),col=factor(col)),alpha=0.05)+
#   # geom_point(data = plot_volcano2[plot_volcano2$col==1,],aes(x=avg_log2FC,y=-log(p_val_adj),col=factor(col)),alpha=0.5)+
#   geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
#   geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
#   # scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
#   # ggrepel::geom_text_repel(
#   #   data = plot_volcano[plot_volcano$col==1,][1:1000,],
#   #   aes(label = symbol),max.overlaps = 1,segment.alpha=0.4,
#   #   size = 2,
#   #   box.padding = unit(0.35, "lines"),
#   #   point.padding = unit(0.3, "lines")) +
#   # ggrepel::geom_text_repel(
#   #   data = plot_volcano %>% 
#   #     group_by(conditionVsCTRL) %>% 
#   #     arrange(padj) %>% 
#   #     dplyr::slice(1:30) %>% 
#   #     dplyr::filter(abs(log2FoldChange)>1) %>% 
#   #     dplyr::filter(padj<0.05),
#   #   aes(label = symbol),segment.alpha=0.4) +
#   facet_wrap(~test)+
#   theme_bw()+
#   theme(strip.background = element_blank())+
#   theme(legend.position = "none")
# # ggsave("../out/plot/102_vulcano_plot_pseudobulk_filterExp.pdf",width = 10,height = 10)
# 
