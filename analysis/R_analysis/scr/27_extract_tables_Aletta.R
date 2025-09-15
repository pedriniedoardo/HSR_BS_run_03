# AIM ---------------------------------------------------------------------
# Aletta suggested to provide the raw and normalized table of expression fo the MG subset across BASELINE, Fe and Myelin treatment

# libraries ---------------------------------------------------------------
library(tidyverse)
library(GSVA)
library(limma)
library(ComplexHeatmap)
# library(org.Hs.eg.db)
library(AnnotationHub)
library(AnnotationDbi)
library(msigdbr)
library(ggrepel)

# read in the data --------------------------------------------------------
# read in the filtered pbulk object for MG subset
ddsHTSeq_filter <- readRDS("../../out/object/22_ddsHTSeq_MG_pseudobulk_filterExp.rds")

# wrangling ---------------------------------------------------------------
# extract the raw expression
exp_01 <- counts(ddsHTSeq_filter,normalized = F) %>%
  data.frame() %>%
  dplyr::select(contains("BASELINE")|contains("Fe")|contains("myelin")) %>%
  as.matrix()

exp_01 %>%
  saveRDS("../../out/object/27_MG_pseudobulk_filterExp_rawCounts.rds")

# extract the normalized expression
exp_02 <- counts(ddsHTSeq_filter,normalized = T) %>%
  data.frame() %>%
  dplyr::select(contains("BASELINE")|contains("Fe")|contains("myelin")) %>%
  as.matrix()

exp_02 %>%
  saveRDS("../../out/object/27_MG_pseudobulk_filterExp_deseq2NormCounts.rds")
