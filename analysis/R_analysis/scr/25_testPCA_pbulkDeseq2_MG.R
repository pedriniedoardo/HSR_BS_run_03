# AIM ---------------------------------------------------------------------
# aletta suggested to try to make a PCA with only the samples from BASELINE, Fe and Myelin

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

# add in one covariate the cell anntation and the stimulation status
data.combined$celltype.stim <- paste(data.combined$expertAnno.l1, data.combined$treat_full, sep = "_")

# set the ident
Idents(data.combined) <- "celltype.stim"

# For this test I would focus on the subset only the presumed cells of interest cells for the test
scobj_subset <- subset(data.combined,subset = expertAnno.l1 %in% c("MG") & harmonized_donor2 %in% c("donRR16","donRR24","donRR25") & treat_full %in% c("BASELINE","myelin","Fe"))

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
# saveRDS(counts,file = "../../out/object/22_counts_MG_pBulk.rds")
# saveRDS(colData,file = "../../out/object/22_colData_MG_pBulk.rds")

# perform DESeq2 ----------------------------------------------------------
# build the model
treat <- colData$treat_full
block_donor <- colData$harmonized_donor2

# design <- model.matrix(~ block_donor + treat)
design <- model.matrix(~ treat)
colnames(design)[1] <- c("intercept")

# save the disign
# saveRDS(design,"../../out/object/22_design_MG_pBulk.rds")

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
# ggsave("../../out/image/22_counts_MG.pdf",width = 6,height = 4)

# scale the data
vds_filter <- vst(dds_filter, blind = F)

# clustering samples ------------------------------------------------------
# set seed to control random annotation colors
# pdf("../../out/image/22_heatmap_SampleCluster_MG_filterExp.pdf",width = 10,height = 6)
set.seed(1)
hm <- plot_sample_clustering(vds_filter,
                             anno_vars = c("harmonized_donor2","treat_full"),
                             distance = "euclidean")
draw(hm,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2, 2, 2, 30), "mm"))
# dev.off()

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
ggsave("../../out/image/22_PCA_pseudobulk_MG_filterExp.pdf",width = 6,height = 4)
