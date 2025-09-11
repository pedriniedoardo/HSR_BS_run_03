# AIM ---------------------------------------------------------------------
# sample routine for the pseudobulk analysis.
# the following is a readaptation of different resources
# https://github.com/sib-swiss/single-cell-training/
# https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(GGally)
library(cowplot)
library(ComplexHeatmap)
library(scales)
library(circlize)
library(DESeq2)
library(RNAseqQC)
library(limma)
library(ashr)
library(magick)
library(UpSetR)

# read in the data --------------------------------------------------------
# read in the sample dataset
data.combined <- readRDS("../../out/object/sobj_processed_donor.rds")

# check the object version
class(data.combined@assays$RNA)

# load in the genes sets to check
set_01 <- c("FPR2", "CMKLR1","GPR32","GPR18","GPR37","LGR6","LTB4R","RORA")
set_02 <- c("ABCC1", "AKR1C3", "ALOX12", "ALOX15", "ALOX5", "CBR1", "CYP1A2", "CYP2C8", "CYP2C9", "CYP2D6", "CYP2E1", "CYP3A4", "CYP4F2", "CYP4F3", "DPEP1", "DPEP2", "DPEP3", "EPHX1", "EPHX2", "EPHX3", "GGT1", "GGT2", "GGT5", "GPX2", "GPX4", "GSTM4", "HPGDS", "LTA4H", "LTC4S", "PTGDS", "PTGES", "PTGES2", "PTGIS", "PTGS1", "PTGS2", "TBXAS1", "TXN")

# add in one covariate the cell anntation and the stimulation status
data.combined$celltype.stim <- paste(data.combined$expertAnno.l1, data.combined$treat_full, sep = "_")

# set the ident
Idents(data.combined) <- "celltype.stim"

# For this test I would focus on the subset only the presumed cells of interest cells for the test
scobj_subset <- subset(data.combined,subset = expertAnno.l1 %in% c("ASTRO") & harmonized_donor2 %in% c("donRR16","donRR24","donRR25"))

# DimPlot(scobj_subset_subset,label = T,raster = T,group.by = "stim")

# wrangling ---------------------------------------------------------------
# explore the full meta
scobj_subset@meta.data %>%
  group_by(harmonized_donor2,treat_full) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = treat_full,values_from = n)

# aggregate the expression per sample per donor.id and stimulation
cts_sample_all <- AggregateExpression(object = scobj_subset,
                                      group.by = c("treat_full","harmonized_donor2"),
                                      assays = 'RNA',
                                      slot = "counts",
                                      return.seurat = FALSE)

# processing --------------------------------------------------------------
# extract
# 1. Get counts matrix
counts <- cts_sample_all$RNA %>%
  as.data.frame()

# 2. generate sample level metadata
LUT_sample <- scobj_subset@meta.data %>%
  group_by(harmonized_donor2,treat_full,expertAnno.l1,celltype.stim) %>%
  summarise() %>%
  # harmonyze the naming of the donor_id
  mutate(stim_donor = paste0(treat_full,"_",harmonized_donor2))

# match the LUT with the expression matrix
colData <- data.frame(sample.id = colnames(counts)) %>%
  # separate(samples,into = c(c("orig.ident","origin_01","pathology_class","origin_02")),remove = F,sep = "_") %>%
  left_join(LUT_sample,by = c("sample.id" = "stim_donor"))

# save matrix and metadata
saveRDS(counts,file = "../../out/object/22_counts_ASTRO_pBulk.rds")
saveRDS(colData,file = "../../out/object/22_colData_ASTRO_pBulk.rds")

# perform DESeq2 ----------------------------------------------------------
# build the model
treat <- colData$treat_full
block_donor <- colData$harmonized_donor2

# design <- model.matrix(~ block_donor + treat)
design <- model.matrix(~ treat)
colnames(design)[1] <- c("intercept")

# save the disign
saveRDS(design,"../../out/object/22_design_ASTRO_pBulk.rds")

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = design)

# filter low aboundant features
feat_keep <- edgeR::filterByExpr(counts(dds), group = colData$treat_full)
dds_filter <- dds[feat_keep,]

# plot the raw number of reads per sample
colSums(counts(dds_filter)) %>%
  data.frame(tot_reads = .) %>%
  rownames_to_column("sample") %>%
  ggplot(aes(x=sample,y=tot_reads)) + geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))
ggsave("../../out/image/22_counts_ASTRO.pdf",width = 6,height = 4)

# scale the data
vds_filter <- vst(dds_filter, blind = F)

# clustering samples ------------------------------------------------------
# set seed to control random annotation colors
pdf("../../out/image/22_heatmap_SampleCluster_ASTRO_filterExp.pdf",width = 10,height = 6)
set.seed(1)
hm <- plot_sample_clustering(vds_filter,
                             anno_vars = c("harmonized_donor2","treat_full"),
                             distance = "euclidean")
draw(hm,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2, 2, 2, 30), "mm"))
dev.off()

# PCA ---------------------------------------------------------------------
plot_vsd <- plotPCA(vds_filter,
                    intgroup = c("harmonized_donor2","treat_full")) +
  theme_bw()

plot_vsd$data %>%
  ggplot(aes(x=PC1,y=PC2,col=treat_full)) +
  geom_point(size =3,alpha=0.6) +
  # ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vsd$labels[1]) + xlab(plot_vsd$labels[2])
ggsave("../../out/image/22_PCA_pseudobulk_ASTRO_filterExp.pdf",width = 6,height = 4)

# run DE ------------------------------------------------------------------
# run DESeq2
ddsHTSeq_filter <- DESeq(dds_filter)

# Check the coefficients for the comparison
resultsNames(ddsHTSeq_filter)

# save the filtered object
saveRDS(ddsHTSeq_filter,"../../out/object/22_ddsHTSeq_ASTRO_pseudobulk_filterExp.rds")

# print the contrast
resultsNames(ddsHTSeq_filter)

# save the resut of the contrast of interest. notice that is also possible to set the alpha value for the significance (adj-p value)
contrast <- makeContrasts(ASTRO_Cytokine_vs_BASELINE = treatcytokine,
                          ASTRO_CSFCtrl24_vs_BASELINE = treatCSF.ctrl.24h,
                          ASTRO_CSFMS24_vs_BASELINE = treatCSF.MS.24h,
                          ASTRO_CSFMS48_vs_BASELINE = treatCSF.MS.48h,
                          ASTRO_Fe_vs_BASELINE = treatFe,
                          ASTRO_Myelin_vs_BASELINE = treatmyelin,
                          ASTRO_TBHP_vs_BASELINE = treatTBHP,
                          levels = design)

# pull the results table
res_cytokine <- results(ddsHTSeq_filter, contrast=contrast[,"ASTRO_Cytokine_vs_BASELINE"],alpha = 0.05)
res_CSFCtrl24 <- results(ddsHTSeq_filter, contrast=contrast[,"ASTRO_CSFCtrl24_vs_BASELINE"],alpha = 0.05)
res_CSFMS24 <- results(ddsHTSeq_filter, contrast=contrast[,"ASTRO_CSFMS24_vs_BASELINE"],alpha = 0.05)
res_CSFMS48 <- results(ddsHTSeq_filter, contrast=contrast[,"ASTRO_CSFMS48_vs_BASELINE"],alpha = 0.05)
res_Fe <- results(ddsHTSeq_filter, contrast=contrast[,"ASTRO_Fe_vs_BASELINE"],alpha = 0.05)
res_Myelin <- results(ddsHTSeq_filter, contrast=contrast[,"ASTRO_Myelin_vs_BASELINE"],alpha = 0.05)
res_TBHP <- results(ddsHTSeq_filter, contrast=contrast[,"ASTRO_TBHP_vs_BASELINE"],alpha = 0.05)

summary(res_cytokine)
summary(res_CSFCtrl24)
summary(res_CSFMS24)
summary(res_CSFMS48)
summary(res_Fe)
summary(res_Myelin)
summary(res_TBHP)

# add the gene symbols
df_res <-
  list(res_cytokine = res_cytokine,
       res_CSFCtrl24 = res_CSFCtrl24,
       res_CSFMS24 = res_CSFMS24,
       res_CSFMS48 = res_CSFMS48,
       res_Fe = res_Fe,
       res_Myelin = res_Myelin,
       res_TBHP = res_TBHP) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue)
    # left_join(LUT_gene,by = c("symbol"="gene_name"))
    df
  }) %>%
  bind_rows(.id = "conditionVsBASELINE")

# save the table
df_res %>%
  write_tsv("../../out/table/22_DE_pseudobulk_ASTRO_filterExp.tsv")

df_res %>%
  filter(symbol %in% set_01)

df_res %>%
  filter(symbol %in% set_02)

# shrink ------------------------------------------------------------------
res_cytokine_shr <- lfcShrink(ddsHTSeq_filter, res = res_cytokine, type = "ashr")
res_CSFCtrl24_shr <- lfcShrink(ddsHTSeq_filter, res = res_CSFCtrl24, type = "ashr")
res_CSFMS24_shr <- lfcShrink(ddsHTSeq_filter, res = res_CSFMS24, type = "ashr")
res_CSFMS48_shr <- lfcShrink(ddsHTSeq_filter, res = res_CSFMS48, type = "ashr")
res_Fe_shr <- lfcShrink(ddsHTSeq_filter, res = res_Fe, type = "ashr")
res_Myelin_shr <- lfcShrink(ddsHTSeq_filter, res = res_Myelin, type = "ashr")
res_TBHP_shr <- lfcShrink(ddsHTSeq_filter, res = res_TBHP, type = "ashr")

summary(res_cytokine_shr)
summary(res_CSFCtrl24_shr)
summary(res_CSFMS24_shr)
summary(res_CSFMS48_shr)
summary(res_Fe_shr)
summary(res_Myelin_shr)
summary(res_TBHP_shr)

# save the table of DEGs
df_res_shr <-
  list(res_cytokine_shr = res_cytokine_shr,
       res_CSFCtrl24_shr = res_CSFCtrl24_shr,
       res_CSFMS24_shr = res_CSFMS24_shr,
       res_CSFMS48_shr = res_CSFMS48_shr,
       res_Fe_shr = res_Fe_shr,
       res_Myelin_shr = res_Myelin_shr,
       res_TBHP_shr = res_TBHP_shr) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue)
    # left_join(LUT_gene,by = c("symbol"="gene_name"))
    df
  }) %>%
  bind_rows(.id = "conditionVsBASELINE")

# save the table
df_res_shr %>%
  write_tsv("../../out/table/22_DE_pseudobulk_ASTRO_filterExp_shr.tsv")

df_res_shr %>%
  filter(symbol %in% set_01)

df_res_shr %>%
  filter(symbol %in% set_02)

# Another useful diagnostic plot is the histogram of the p values (figure below). This plot is best formed by excluding genes with very small counts, which otherwise generate spikes in the histogram.
df_res %>%
  data.frame()%>%
  dplyr::filter(baseMean>1) %>%
  ggplot(aes(x=pvalue))+geom_histogram(breaks = 0:20/20) +
  facet_wrap(~conditionVsBASELINE)+
  theme_bw()+
  theme(strip.background = element_blank())
ggsave("../../out/image/22_histogram_pvalue_pseudobulk_ASTRO_filterExp.pdf",width = 6,height = 4)

# volcano -----------------------------------------------------------------
# add the info of the genename
# focus only on the genes in the set
plot_volcano <- df_res %>%
  filter(symbol %in% c(set_01,set_02)) %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1&!is.na(symbol)),yes = 1,no = 0)) %>%
  filter(!is.na(col))

plot_volcano %>%
  ggplot(aes(x=log2FoldChange,y=-log(padj)))+
  # geom_point()
  geom_point(data = df_res,aes(x=log2FoldChange,y=-log(padj)),col="gray",alpha=0.3,shape = 1) +
  geom_point(data = plot_volcano[plot_volcano$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = plot_volcano[plot_volcano$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  # ggrepel::geom_text_repel(
  #   data = plot_volcano[plot_volcano$col==1,][1:1000,],
  #   aes(label = symbol),max.overlaps = 1,segment.alpha=0.4,
  #   size = 2,
  #   box.padding = unit(0.35, "lines"),
  #   point.padding = unit(0.3, "lines")) +
  ggrepel::geom_text_repel(
    data = plot_volcano %>% 
      group_by(conditionVsBASELINE) %>% 
      arrange(padj) %>% 
      # dplyr::slice(1:30) %>% 
      dplyr::filter(abs(log2FoldChange)>1) %>% 
      dplyr::filter(padj<0.05),
    aes(label = symbol),segment.alpha=0.4) +
  facet_wrap(~conditionVsBASELINE)+
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(legend.position = "none")
ggsave("../../out/image/22_vulcano_plot_pseudobulk_ASTRO_filterExp.pdf",width = 10,height = 10)

#
plot_volcano_shr <- df_res_shr %>%
  filter(symbol %in% c(set_01,set_02)) %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1&!is.na(symbol)),yes = 1,no = 0)) %>%
  filter(!is.na(col))

plot_volcano_shr %>%
  ggplot(aes(x=log2FoldChange,y=-log(padj)))+
  geom_point(data = df_res_shr,aes(x=log2FoldChange,y=-log(padj)),col="gray",alpha=0.3,shape = 1) +
  geom_point(data = plot_volcano_shr[plot_volcano_shr$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = plot_volcano_shr[plot_volcano_shr$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  ggrepel::geom_text_repel(
    data = plot_volcano_shr %>% 
      group_by(conditionVsBASELINE) %>%
      arrange(padj) %>%
      # dplyr::slice(1:10)%>%
      filter(abs(log2FoldChange)>1),
    aes(label = symbol),segment.alpha=0.4) +
  facet_wrap(~conditionVsBASELINE)+
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(legend.position = "none")
ggsave("../../out/plot/22_vulcano_plot_pseudobulk_ASTRO_filterExp_shr.pdf",width = 10,height = 10)

# MA plot -----------------------------------------------------------------
df_res %>%
  filter(!is.na(padj)) %>%
  data.frame()%>%
  mutate(color = ifelse(padj<0.05,yes = 1,no = 0))%>%
  mutate(color = factor(ifelse(is.na(color),yes = 0,no = color)))%>%
  arrange(color) %>%
  ggplot(aes(baseMean,y = log2FoldChange,col=color)) + geom_point(alpha=0.2) +
  scale_x_log10() + scale_color_manual(values = c("gray","red")) +
  facet_wrap(~conditionVsBASELINE)+
  theme_bw()+
  theme(strip.background = element_blank()) +
  theme(legend.position = "none")
# ggsave("../out/plot/22_MA_plot_pseudobulk_filterExp.pdf",width = 10,height = 10)

# shr
df_res_shr %>%
  filter(!is.na(padj)) %>%
  data.frame()%>%
  mutate(color = ifelse(padj<0.05,yes = 1,no = 0))%>%
  mutate(color = factor(ifelse(is.na(color),yes = 0,no = color)))%>%
  arrange(color) %>%
  ggplot(aes(baseMean,y = log2FoldChange,col=color)) + geom_point(alpha=0.2) +
  scale_x_log10() + scale_color_manual(values = c("gray","red")) +
  facet_wrap(~conditionVsBASELINE)+
  theme_bw()+
  theme(strip.background = element_blank()) +
  theme(legend.position = "none")
# ggsave("../out/plot/101_MA_plot_pseudobulk_filterExp_shr.pdf",width = 10,height = 10)

# heatmaps ----------------------------------------------------------------
# pull the scaled values
mat_filter <- assay(vds_filter) %>%
  as.data.frame() %>%
  # dplyr::select(contains(c("_0_","_6_"))) %>%
  as.matrix()

# the DEGs plot
DEG_2 <- df_res_shr %>%
  as.data.frame()%>%
  filter(symbol %in% c(set_01,set_02)) %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1),yes = 1,no = 0)) %>%
  dplyr::filter(col==1) %>%
  pull(symbol) %>%
  unique()

# df_res_shr %>%
#   as.data.frame()%>%
#   filter(symbol %in% c(set_01,set_02)) %>%
#   # add a clor variable in case significant
#   mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1),yes = 1,no = 0)) %>%
#   filter(symbol %in% c("ABCC1","EPHX1"))

# mat_filter <- assay(vds_filter) %>%
#   data.frame() %>%
#   # dplyr::select(contains(c("_0_","_6_"))) %>%
#   as.matrix()

mat_shr <- mat_filter[rownames(vds_filter) %in% DEG_2, ]
mat2_shr <- (mat_shr - rowMeans(mat_shr))/rowSds(mat_shr,useNames = TRUE)
#
meta_sample <- data.frame(colname = colnames(mat2_shr)) %>%
  left_join(colData,by=c("colname"="sample.id"))

# make the column of the matrix more readable
colnames(mat2_shr) <- meta_sample$colname

column_ha_shr <- HeatmapAnnotation(treat = meta_sample$treat_full,
                                   # gender = meta_sample_ENDO$sex,
                                   # facility = meta_sample_ENDO$facility,
                                   col = list(treat = c("BASELINE" = "blue",
                                                        "CSF.ctrl.24h" = "gray",
                                                        "CSF.MS.24h" = "gray30",
                                                        "CSF.MS.48h" = "black",
                                                        "cytokine" = "purple",
                                                        "Fe" = "orange",
                                                        "myelin" = "cyan",
                                                        "TBHP" = "yellow")))

ht2_shr <- Heatmap(mat2_shr, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                   name = "exp",
                   column_title = "test_shr",
                   # row_names_gp = gpar(fontsize = 3),
                   top_annotation = column_ha_shr,show_row_names = T,
                   # cluster_rows = F,
                   # right_annotation = row_ha,
                   # row_split = rep(c(1,2,3,4),c(2,3,4,7))
)
pdf("../../out/image/22_heatmap_DEG_plot_pseudobulk_ASTRO_filterExp_shr.pdf",width = 7,height = 5)
draw(ht2_shr,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
dev.off()

# -------------------------------------------------------------------------
# Martina suggested to produce an heatmap only fro BASELINE vs csf ms or BASELINE vs cytokines

# BASELINE vs cytokines
mat_shr_cytokines <- mat_filter[rownames(vds_filter) %in% c(set_01,set_02), str_subset(colnames(mat_filter),pattern = c("BASELINE|cytokine"))]
mat2_shr_cytokines <- (mat_shr_cytokines - rowMeans(mat_shr_cytokines))/rowSds(mat_shr_cytokines,useNames = TRUE)

#
meta_sample_cytokines <- data.frame(colname = colnames(mat2_shr_cytokines)) %>%
  left_join(colData,by=c("colname"="sample.id"))

# make the column of the matrix more readable
colnames(mat2_shr_cytokines) <- meta_sample_cytokines$colname

column_ha_shr_cytokines <- HeatmapAnnotation(treat = meta_sample_cytokines$treat_full,
                                   # gender = meta_sample_ENDO$sex,
                                   # facility = meta_sample_ENDO$facility,
                                   col = list(treat = c("BASELINE" = "blue",
                                                        "cytokine" = "purple")))

ht2_shr_cytokines <- Heatmap(mat2_shr_cytokines, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                   name = "exp",
                   column_title = "test_shr",
                   # row_names_gp = gpar(fontsize = 3),
                   top_annotation = column_ha_shr_cytokines,show_row_names = T,
                   # cluster_rows = F,
                   # right_annotation = row_ha,
                   # row_split = rep(c(1,2,3,4),c(2,3,4,7))
)
pdf("../../out/image/22_heatmap_DEG_plot_pseudobulk_ASTRO_cytokinesALL_filterExp_shr.pdf",width = 5,height = 7)
draw(ht2_shr_cytokines,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
dev.off()

# plot only the significant ones for the comparison
# the DEGs plot
DEG_2_cytokines <- df_res_shr %>%
  as.data.frame()%>%
  filter(symbol %in% c(set_01,set_02)) %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1),yes = 1,no = 0)) %>%
  dplyr::filter(conditionVsBASELINE %in% "res_cytokine_shr") %>%
  dplyr::filter(col==1) %>%
  pull(symbol) %>%
  unique()

ht2_shr_cytokines_sig <- Heatmap(mat2_shr_cytokines[DEG_2_cytokines,], show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                             name = "exp",
                             column_title = "test_shr",
                             # row_names_gp = gpar(fontsize = 3),
                             top_annotation = column_ha_shr_cytokines,show_row_names = T,
                             # cluster_rows = F,
                             # right_annotation = row_ha,
                             # row_split = rep(c(1,2,3,4),c(2,3,4,7))
)
pdf("../../out/image/22_heatmap_DEG_plot_pseudobulk_ASTRO_cytokinesSIG_filterExp_shr.pdf",width = 5,height = 5)
draw(ht2_shr_cytokines_sig,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
dev.off()

# -------------------------------------------------------------------------
# BASELINE vs CSFMS
mat_shr_csfms <- mat_filter[rownames(vds_filter) %in% c(set_01,set_02), str_subset(colnames(mat_filter),pattern = c("BASELINE|CSF.MS"))]
mat2_shr_csfms <- (mat_shr_csfms - rowMeans(mat_shr_csfms))/rowSds(mat_shr_csfms,useNames = TRUE)

#
meta_sample_csfms <- data.frame(colname = colnames(mat2_shr_csfms)) %>%
  left_join(colData,by=c("colname"="sample.id"))

# make the column of the matrix more readable
colnames(mat2_shr_csfms) <- meta_sample_csfms$colname

column_ha_shr_csfms <- HeatmapAnnotation(treat = meta_sample_csfms$treat_full,
                                             # gender = meta_sample_ENDO$sex,
                                             # facility = meta_sample_ENDO$facility,
                                             col = list(treat = c("BASELINE" = "blue",
                                                                  "CSF.MS.24h" = "gray30",
                                                                  "CSF.MS.48h" = "black")))

ht2_shr_csfms <- Heatmap(mat2_shr_csfms, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                             name = "exp",
                             column_title = "test_shr",
                             # row_names_gp = gpar(fontsize = 3),
                             top_annotation = column_ha_shr_csfms,show_row_names = T,
                             # cluster_rows = F,
                             # right_annotation = row_ha,
                             # row_split = rep(c(1,2,3,4),c(2,3,4,7))
)
pdf("../../out/image/22_heatmap_DEG_plot_pseudobulk_ASTRO_csfmsALL_filterExp_shr.pdf",width = 5,height = 7)
draw(ht2_shr_csfms,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
dev.off()

# plot only the significant ones for the comparison
# the DEGs plot
DEG_2_csfms <- df_res_shr %>%
  as.data.frame()%>%
  filter(symbol %in% c(set_01,set_02)) %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>0.8),yes = 1,no = 0)) %>%
  dplyr::filter(conditionVsBASELINE %in% c("res_CSFMS24_shr","res_CSFMS48_shr")) %>%
  dplyr::filter(col==1) %>%
  pull(symbol) %>%
  unique()

ht2_shr_csfms_sig <- Heatmap(mat2_shr_csfms[DEG_2_csfms,], show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                                 name = "exp",
                                 column_title = "test_shr",
                                 # row_names_gp = gpar(fontsize = 3),
                                 top_annotation = column_ha_shr_csfms,show_row_names = T,
                                 # cluster_rows = F,
                                 # right_annotation = row_ha,
                                 # row_split = rep(c(1,2,3,4),c(2,3,4,7))
)
pdf("../../out/image/22_heatmap_DEG_plot_pseudobulk_ASTRO_csfmsSIG_filterExp_shr.pdf",width = 5,height = 5)
draw(ht2_shr_csfms_sig,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
dev.off()

# -------------------------------------------------------------------------
# data <- readRDS("../../out/object/ddsHTSeq_filter_DESeq.rds")

lut <- colData(ddsHTSeq_filter) %>%
  data.frame()

# GOI <- c(sig_B$symbol,sig_D$symbol)
GOI <- c(set_01,set_02)
# GOI <- subset_genes

MR <- counts(ddsHTSeq_filter,normalized=T) %>%
  data.frame()%>%
  rownames_to_column("symbol") %>%
  pivot_longer(names_to = "sample",values_to = "exp",-symbol)%>%
  group_by(sample)%>%
  summarise(MR = sum(exp)/10^6)

# plot the data following the methods implemented in the plotCounts funciton from DESeq2
# Normalized counts plus a pseudocount of 0.5 are shown by default.
counts(ddsHTSeq_filter,normalized=T) %>%
  data.frame() %>%
  rownames_to_column("symbol") %>%
  filter(symbol %in% GOI) %>%
  pivot_longer(names_to = "sample",values_to = "count",-symbol) %>%
  # add the milion reads per sample
  left_join(MR,by = "sample") %>%
  left_join(lut,by = c("sample" = "sample.id")) %>%
  mutate(count_norm_adj = count + 0.5) %>%
  mutate(treat = fct_relevel(treat_full,"BASELINE")) %>%
  # mutate(symbol = factor(symbol,levels = c("PROCR","GDF9", "GDF11","TGFB1", "TGFB2", "TGFB3","INHBA", "INHBB","MSTN"))) %>%
  ggplot(aes(x=treat,y = count_norm_adj))+
  geom_boxplot(outlier.shape = NA,linewidth = 1)+
  geom_point(position = position_jitter(width = 0.1),alpha=0.6,size = 3)+
  facet_wrap(~symbol,scales = "free_y",ncol=5)+
  scale_y_continuous(trans = "log1p") +
  theme_bw(base_size = 22,base_rect_size = 2)+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_text(angle = 90,hjust = 1))
ggsave("../../out/image/22_boxplot_ASTRO_GOI.pdf",width = 20,height = 20)

# # upset plot --------------------------------------------------------------
# # read in the table of DEGs
# df_res_shr <- read_tsv("../../out/table/DE_pseudobulk_ASTRO_ctrl_refCX_shr.tsv")
# 
# # build a list of common elements belonging to each set of fegs
# list_DE_up <- df_res_shr %>%
#   split(f = .$conditionVsCX) %>%
#   lapply(function(x){
#     x %>%
#       filter(padj < 0.05,log2FoldChange>1) %>%
#       pull(symbol) %>%
#       unique()
#   })
# 
# glimpse(list_DE_up)
# 
# list_DE_down <- df_res_shr %>%
#   split(f = .$conditionVsCX) %>%
#   lapply(function(x){
#     x %>%
#       filter(padj < 0.05,log2FoldChange<(-1)) %>%
#       pull(symbol) %>%
#       unique()
#   })
# glimpse(list_DE_down)
# 
# # try the upset plot version
# # library(UpSetR)
# pdf("../../out/image/upset_DEG_UP_plot_pseudobulk_ASTRO_refBASELINE_shr.pdf",width = 14,height = 7)
# upset(fromList(list_DE_up), order.by = "freq",nsets = 7)
# dev.off()
# 
# pdf("../../out/image/upset_DEG_DOWN_plot_pseudobulk_ASTRO_refBASELINE_shr.pdf",width = 14,height = 7)
# upset(fromList(list_DE_down), order.by = "freq",nsets = 7)
# dev.off()
# 
# # pull the intersections
# df1_UP <- lapply(list_DE_up,function(x){
#   data.frame(gene = x)
# }) %>%
#   bind_rows(.id = "path")
# 
# df1_DOWN <- lapply(list_DE_down,function(x){
#   data.frame(gene = x)
# }) %>%
#   bind_rows(.id = "path")
# 
# head(df1_UP)
# head(df1_DOWN)
# 
# df2_UP <- data.frame(gene=unique(unlist(list_DE_up)))
# df2_DOWN <- data.frame(gene=unique(unlist(list_DE_down)))
# 
# head(df2_UP)
# head(df2_DOWN)
# 
# df_int_UP <- lapply(df2_UP$gene,function(x){
#   # pull the name of the intersections
#   intersection <- df1_UP %>%
#     dplyr::filter(gene==x) %>%
#     arrange(path) %>%
#     pull("path") %>%
#     paste0(collapse = "|")
# 
#   # build the dataframe
#   data.frame(gene = x,int = intersection)
# }) %>%
#   bind_rows()
# 
# df_int_DOWN <- lapply(df2_DOWN$gene,function(x){
#   # pull the name of the intersections
#   intersection <- df1_DOWN %>%
#     dplyr::filter(gene==x) %>%
#     arrange(path) %>%
#     pull("path") %>%
#     paste0(collapse = "|")
# 
#   # build the dataframe
#   data.frame(gene = x,int = intersection)
# }) %>%
#   bind_rows()
# 
# df_int_UP %>%
#   write_tsv("../../out/table/upset_DEG_UP_plot_pseudobulk_ASTRO_refBASELINE_shr.tsv")
# 
# df_int_DOWN %>%
#   write_tsv("../../out/table/upset_DEG_DOWN_plot_pseudobulk_ASTRO_refBASELINE_shr.tsv")
# 
# head(df_int_UP,n=20)
# head(df_int_DOWN,n=20)
# 
# df_int_UP %>%
#   group_by(int) %>%
#   summarise(n=n()) %>%
#   arrange(desc(n))
# 
# df_int_DOWN %>%
#   group_by(int) %>%
#   summarise(n=n()) %>%
#   arrange(desc(n))

# PLOT DISPERSION ---------------------------------------------------------
pdf("../../out/image/22_dispersion_plot_pseudobulk_ASTRO_filterExp.pdf",width = 5,height = 5)
plotDispEsts(ddsHTSeq_filter)
dev.off()