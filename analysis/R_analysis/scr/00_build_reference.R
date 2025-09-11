# libraries ---------------------------------------------------------------
library(Matrix)
library(Seurat)
library(Azimuth)
library(presto)
library(dplyr)
# args <- commandArgs(trailingOnly = TRUE)

# read in the object for which I want to build a reference
ref <- readRDS("/beegfs/scratch/ric.cosr/pedrini.edoardo/scRNAseq_Brainsphere_Absinta/BS_drop_RR25CTRL/data/data.combined_NOT_annotated_norm_fix_DoubletSinglet_harmonyMartina.rds")

DimPlot(ref)

# ref.dir <- "reference/"
# ob.dir <- "seurat_objects/"
# ref <- readRDS(file = args[1])
# annotations <- readRDS(file = args[2])
# Idents(object = ref) <- annotations

Idents(object = ref) <- "seurat_clusters"

# if ("remove" %in% levels(x = ref)) {
#   ref <- subset(x = ref, idents = "remove", invert = TRUE)
#   ref <- RunPCA(object = ref, verbose = FALSE)
# }
# ref$annotation.l1 <- Idents(object = ref)
# needs to have the model
ref_SCT <- SCTransform(ref, method = "glmGamPoi", vars.to.regress = c("percent.mt","nCount_RNA"), verbose = T)
ref2 <- RunUMAP(object = ref_SCT,reduction = "harmony",dims = 1:10, return.model = TRUE)

DimPlot(ref2)
# ref2@reductions$umap@misc$model
full.ref <- ref2

full.ref$annotation.l1 <- Idents(object = full.ref)
colormap <- list(annotation.l1 = CreateColorMap(object = ref2, seed = 2))
colormap[["annotation.l1"]] <- colormap[["annotation.l1"]][sort(x = names(x = colormap[["annotation.l1"]]))]

ref_final <- AzimuthReference(
  object = full.ref,
  refUMAP = "umap",
  refDR = "pca",
  refAssay = "SCT",
  metadata = c("annotation.l1"),
  dims = 1:50,
  k.param = 31,
  colormap = colormap,
  reference.version = "1.0.0"
)

DimPlot(ref_final)

SaveAnnoyIndex(object = ref_final[["refdr.annoy.neighbors"]], file = file.path("../../data/ref_BS_run_01_02/","idx.annoy"))
saveRDS(object = ref_final, file = file.path("../../data/ref_BS_run_01_02/", "ref.Rds"))
saveRDS(object = full.ref, file = file.path("../../data/ref_BS_run_01_02/", "fullref.Rds"))

# test reading azimuth reference ------------------------------------------
# try to load the new reference object created
reference <- LoadReference(path = "../../data/ref_BS_run_01_02/")
DimPlot(reference$plot,group.by = "annotation.l1",label = T)
