# AIM ---------------------------------------------------------------------
# single cell DE on sobj_processed_donor.rds, per cell type (expertAnno.l1), each treatment vs BASELINE (Seurat FindMarkers, accelerated by presto).
# same settings of 07_DE_analysis_refBASELINE.R, extended to all the cell types.
# for each treatment the script also plots the pvalue distribution and the volcano, and summarises the number of DEGs (padj < padj_thr and |avg_log2FC| > lfc_thr, split in up and down):
# 1. barplots of the number of DEGs per cell type and treatment
# 2. "perturbation" score = number of DEGs / number of cells in the comparison. barplot of the score per cell type and UMAP colored by the score.
#    NOTICE: the number of DEGs depends on the statistical power (number of cells), dividing by the number of cells is a rough normalization, interpret the score as a relative ranking between cell types.
# the pseudobulk DE is in 06_DE_pbulk_fullObject.R
# plotting adapted from project_edoardo/240531_scRNAseq_MSSpinal_Absinta/analysis/R_analysis_R45_pixi/scr/05_DE_sc_fullObject.R
# NOTICE: the DE is run on the RNA assay (log normalized data layer).

# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(presto)
library(ggrepel)
library(patchwork)
library(cowplot)

# PARAMETERS ------------------------------------------------------------------
options(stringsAsFactors = FALSE)

var_treat <- "treat_full"
var_celltype <- "expertAnno.l1"
ref <- "BASELINE"
list_treat <- c("CSF.ctrl.24h", "CSF.MS.24h", "CSF.MS.48h", "cytokine", "Fe", "myelin")

# minimum number of cells per side of the comparison, min.pct as in the previous DE
min_cells_sc <- 10
min_pct <- 0.05
# thresholds used to call a gene as DEG in the summary and in the volcano (the full tables are saved without any threshold)
padj_thr <- 0.05
lfc_thr <- 1

# dedicated output subfolders: out/table/05_DE_sc, out/image/05_DE_sc, out/object/05_DE_sc
folder_table <- "../../out/table/05_DE_sc/"
folder_plot <- "../../out/image/05_DE_sc/"
folder_object <- "../../out/object/05_DE_sc/"

# read in the data --------------------------------------------------------
sobj <- readRDS("../../out/object/sobj_processed_donor.rds")
DimPlot(scobj,group.by = var_celltype,label = T,repel = T)

DefaultAssay(sobj) <- "RNA"

meta <- sobj@meta.data %>%
  mutate(barcode = rownames(.))
df_umap <- Embeddings(sobj, "umap") %>%
  as.data.frame() %>%
  set_names(c("UMAP_1", "UMAP_2")) %>%
  rownames_to_column("barcode")
list_celltype <- sort(unique(as.character(meta[[var_celltype]])))

# run DE ------------------------------------------------------------------
# avg_log2FC: log fold-change of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the treatment (ident.1)
# ct <- "ASTRO"
# tr <- "cytokine"
df_sc <- map(list_celltype, function(ct){
  # subset the object to the cell type and use the treatment as identity (FindMarkers needs ident.1/ident.2)
  sobj_ct <- subset(sobj, cells = meta$barcode[meta[[var_celltype]] == ct])
  Idents(sobj_ct) <- var_treat
  n_cells <- table(Idents(sobj_ct))

  map(list_treat, function(tr){
    n_1 <- if(tr %in% names(n_cells)) n_cells[[tr]] else 0
    n_2 <- if(ref %in% names(n_cells)) n_cells[[ref]] else 0
    if(n_1 < min_cells_sc || n_2 < min_cells_sc){
      message("skip ", ct, " ", tr, " vs ", ref, ": ", n_1, " vs ", n_2, " cells")
      return(NULL)
    }
    message("run ", ct, " ", tr, " vs ", ref)
    tryCatch({
      FindMarkers(sobj_ct, assay = "RNA", ident.1 = tr, ident.2 = ref,
                  min.pct = min_pct, logfc.threshold = 0, only.pos = FALSE, verbose = FALSE) %>%
        rownames_to_column("gene") %>%
        mutate(cell_id = ct, ident.1 = tr, ident.2 = ref, n_1 = n_1, n_2 = n_2)
    }, error = function(e){
      message("failed ", ct, " ", tr, ": ", conditionMessage(e))
      NULL
    })
  }) %>%
    bind_rows()
}) %>%
  bind_rows()

# the full table, and one table per treatment
df_sc %>%
  write_tsv(paste0(folder_table, "05_DE_sc_treatvsBASELINE_", var_celltype, ".tsv"))

# save the tables as individula plots
df_sc %>%
  split(f = .$ident.1) %>%
  iwalk(function(df, tr){
    write_tsv(df, paste0(folder_table, "05_DE_sc_", var_celltype, "_", tr, "_vs_", ref, ".tsv"))
  })

# plot the results --------------------------------------------------------
# one pvalue distribution and one volcano per treatment, faceted by cell type
df_sc %>%
  split(f = .$ident.1) %>%
  iwalk(function(df, tr){
    n_ct <- n_distinct(df$cell_id)
    n_col <- ceiling(sqrt(n_ct))
    n_row <- ceiling(n_ct / n_col)
    file_id <- paste0(var_celltype, "_", tr, "_vs_", ref)

    # show the distribution of the pvalues
    p01 <- df %>%
      ggplot(aes(x = p_val)) +
      geom_histogram(breaks = 0:20/20) +
      theme_bw() +
      facet_wrap(~cell_id, ncol = n_col) +
      theme(strip.background = element_blank()) +
      ggtitle(file_id)
    ggsave(paste0(folder_plot, "05_dist_p_value_sc_", file_id, ".pdf"), plot = p01, width = 3 * n_col + 1, height = 2.5 * n_row + 1, limitsize = FALSE)

    # volcano plot
    test_significant <- df %>%
      filter(abs(avg_log2FC) > lfc_thr & p_val_adj < padj_thr)

    # label only the top genes per cell type to keep the plot readable
    test_label <- test_significant %>%
      group_by(cell_id) %>%
      arrange(p_val_adj) %>%
      dplyr::slice(1:10) %>%
      ungroup()

    p02 <- df %>%
      ggplot(aes(x = avg_log2FC, y = -log(p_val_adj))) +
      geom_point(alpha = 0.01) +
      geom_point(data = test_significant, col = "red", alpha = 0.5) +
      geom_vline(xintercept = c(-lfc_thr, lfc_thr), col = "gray", linetype = "dashed") +
      geom_hline(yintercept = (-log(padj_thr)), col = "gray", linetype = "dashed") +
      geom_text_repel(data = test_label, aes(label = gene), size = 2, segment.alpha = 0.4) +
      facet_wrap(~cell_id, ncol = n_col) +
      theme_bw() +
      theme(strip.background = element_blank()) +
      ggtitle(file_id)
    ggsave(paste0(folder_plot, "05_volcano_sc_", file_id, ".pdf"), plot = p02, width = 4 * n_col + 1, height = 3.5 * n_row + 1, limitsize = FALSE)
  })

# summary of the number of DEGs and perturbation score --------------------
# only the tested comparisons are in df_sc. n_1 + n_2 = number of cells in the comparison
df_summary <- df_sc %>%
  group_by(cell_id, ident.1, ident.2, n_1, n_2) %>%
  summarise(n_up = sum(p_val_adj < padj_thr & avg_log2FC > lfc_thr, na.rm = TRUE),
            n_down = sum(p_val_adj < padj_thr & avg_log2FC < -lfc_thr, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(n_cells = n_1 + n_2,
         n_deg = n_up + n_down,
         score = n_deg / n_cells,
         score_up = n_up / n_cells,
         score_down = n_down / n_cells)

df_summary %>%
  write_tsv(paste0(folder_table, "05_DE_sc_summary_", var_celltype, ".tsv"))

df_summary %>%
  print(n = Inf)

# number of DEGs: all the treatments together, up and down in the two directions
df_summary %>%
  select(cell_id, ident.1, n_up, n_down) %>%
  pivot_longer(cols = c(n_up, n_down), names_to = "direction", values_to = "n") %>%
  # mutate(direction = str_remove(direction, "n_"),
  #        n = ifelse(direction == "down", -n, n)) %>%
  mutate(direction = str_remove(direction, "n_"),cell_id = fct_reorder(cell_id, n)) %>%
  ggplot(aes(y = cell_id, x = n, fill = direction)) +
  geom_col() +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = c(up = "firebrick", down = "steelblue")) +
  facet_wrap(~ident.1, scales = "free_y") +
  theme_bw() +
  theme(strip.background = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = paste0("number of DEGs vs ", ref, " (single cell)"),
       subtitle = paste0("padj < ", padj_thr, ", |log2FC| > ", lfc_thr),
       x = NULL, y = "n DEGs")
ggsave(paste0(folder_plot, "05_barplot_nDEG_sc_", var_celltype, "_vs_", ref, ".pdf"), width = 12, height = 8)

# same, by treatment on the x axis, one panel per cell type
df_summary %>%
  mutate(n_deg_label = n_deg) %>%
  ggplot(aes(x = ident.1, y = n_deg)) +
  geom_col(fill = "grey30") +
  geom_text(aes(label = n_deg), vjust = -0.3, size = 2.5) +
  facet_wrap(~cell_id, scales = "free_y") +
  theme_bw() +
  theme(strip.background = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = paste0("number of DEGs vs ", ref, " (single cell)"), x = NULL, y = "n DEGs") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
ggsave(paste0(folder_plot, "05_barplot_nDEG_sc_perCellType_", var_celltype, "_vs_", ref, ".pdf"), width = 12, height = 8)

# one barplot of the score (stacked by direction) and one UMAP colored by the score per treatment
df_summary %>%
  split(f = .$ident.1) %>%
  iwalk(function(df, tr){
    file_id <- paste0(var_celltype, "_", tr, "_vs_", ref)

    p_bar <- df %>%
      mutate(cell_id = fct_reorder(cell_id, score)) %>%
      pivot_longer(cols = c(score_up, score_down), names_to = "direction", values_to = "score_dir") %>%
      mutate(direction = str_remove(direction, "score_")) %>%
      ggplot(aes(x = cell_id, y = score_dir, fill = direction)) +
      geom_col() +
      scale_fill_manual(values = c(up = "firebrick", down = "steelblue")) +
      coord_flip() +
      theme_bw() +
      theme(strip.background = element_blank()) +
      labs(title = paste0("perturbation score ", file_id),
           subtitle = paste0("n DEGs (padj < ", padj_thr, ", |log2FC| > ", lfc_thr, ") / n cells"),
           x = NULL, y = "perturbation score")

    # the cell types not tested are in grey
    p_umap <- df_umap %>%
      left_join(meta %>% select(barcode, cell_id = all_of(var_celltype)), by = "barcode") %>%
      mutate(cell_id = as.character(cell_id)) %>%
      left_join(df %>% select(cell_id, score), by = "cell_id") %>%
      # plot the cells with a score on top of the not tested ones
      arrange(!is.na(score)) %>%
      ggplot(aes(x = UMAP_1, y = UMAP_2, color = score)) +
      geom_point(size = 0.2, shape = 16) +
      scale_color_viridis_c(na.value = "grey85", option = "inferno", direction = -1) +
      theme_cowplot() +
      labs(title = paste0("perturbation score ", file_id), color = "score")

    ggsave(paste0(folder_plot, "05_perturbation_score_", file_id, ".pdf"),
           plot = p_bar | p_umap, width = 16, height = 6, limitsize = FALSE)
  })

message("\nAll done.")
