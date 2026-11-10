# AIM ---------------------------------------------------------------------
# pseudobulk DE (DESeq2) on sobj_processed_donor.rds, per cell type (expertAnno.l1), each treatment vs BASELINE.
# follows analysis/R_analysis/scr/07_DE_analysis_pseudobulk_refBASELINE.R (written for MG), run on every cell type, with the plotting and the low expression filter of
# project_edoardo/240531_scRNAseq_MSSpinal_Absinta/analysis/R_analysis_R45_pixi/scr/06_DE_pbulk_deseq2_fullObject.R
# raw counts summed per id_sample_short x harmonized_donor2 x treat_full x cell type (AggregateExpression), low expressed genes dropped with edgeR::filterByExpr, DESeq2 ~ harmonized_donor2 + treat_full (donor as fixed effect, paired design; BASELINE as reference), ashr shrunk lfc.
# the donor is a fixed effect because every treatment is measured in each donor (myelin only in RR16 and RR25), so the donor differences are removed from the treatment contrasts. It costs n_donor - 1 df.
# for each cell type:
# 1. QC: total reads, sample clustering, PCA (500 most variable genes), PCA panel, PC scores vs sample metadata, sample clustering and PCA before/after removing the donor (limma::removeBatchEffect, visualisation only)
# 2. DE: pvalue histogram, volcano, heatmap of the DEGs, MA plot (unshrunk and shrunk), dispersion, UpSet of the DEGs shared between treatments
# at the end of the script:
# 1. barplots of the number of DEGs per cell type, stacked by direction
# 2. upset plot of the DEGs shared between the single cell (05_DE_fullObject.R) and the pseudobulk analysis, per treatment and cell type.
#    NOTICE: a gene is counted as shared only if it is called in the same direction by the two approaches (the element of the sets is gene|direction). The same thresholds are applied to the single cell table
# 3. correlation of the log2FC and of the adjusted p-value between single cell and pseudobulk, on the genes tested by both approaches
# run 05_DE_fullObject.R first to have the comparison with the single cell, otherwise that part is skipped.

# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(GGally)
library(cowplot)
library(DESeq2)
library(edgeR)
library(ashr)
library(RNAseqQC)
library(ComplexHeatmap)
library(circlize)
library(ComplexUpset)
library(ggrepel)
library(magick)

# PARAMETERS ------------------------------------------------------------------
options(stringsAsFactors = FALSE)

var_celltype <- "expertAnno.l1"
ref <- "BASELINE"
# list_treat <- c("CSF.ctrl.24h", "CSF.MS.24h", "CSF.MS.48h", "cytokine", "Fe", "myelin")
list_treat <- c("cytokine", "Fe", "myelin")

# donors to drop (not assigned to a single donor)
donor_drop <- c("doublet", "unassigned")

# minimum pseudobulk samples per treatment to include it in the model
min_samples <- 2

# thresholds used to call a gene as DEG in the volcano, heatmap, UpSet and summary
padj_thr <- 0.05
lfc_thr <- 1

color_treat <- c("BASELINE" = "blue",
                 "CSF.ctrl.24h" = "gray90",
                 "CSF.MS.24h" = "gray40",
                 "CSF.MS.48h" = "gray10",
                 "cytokine" = "purple",
                 "Fe" = "orange",
                 "myelin" = "cyan")

# dedicated output subfolders: out/table/06_DE_pbulk, out/image/06_DE_pbulk, out/object/06_DE_pbulk
folder_table <- "../../out/table/06_DE_pbulk/"
folder_plot <- "../../out/image/06_DE_pbulk/"
folder_object <- "../../out/object/06_DE_pbulk/"

# single cell table of 05_DE_fullObject.R
file_sc <- paste0("../../out/table/05_DE_sc/", "05_DE_sc_treatvsBASELINE_", var_celltype, ".tsv")

# read in the final object ------------------------------------------------
scobj <- readRDS("../../out/object/sobj_processed_donor.rds")
DimPlot(scobj,group.by = var_celltype,label = T,repel = T)

DefaultAssay(scobj) <- "RNA"

# AggregateExpression replaces "_" with "-" in the group names (GLIA_IMM -> GLIA-IMM): lookup to restore the original cell type names
celltype_obj <- sort(unique(as.character(scobj@meta.data[[var_celltype]])))
restore_celltype <- function(x){
  celltype_obj[match(x, str_replace_all(celltype_obj, "_", "-"))]
}

# number of cells of each pseudobulk sample
df_ncells <- scobj@meta.data %>%
  dplyr::count(id_sample_short, harmonized_donor2, treat_full, cell = .data[[var_celltype]], name = "n_cells") %>%
  mutate(samples = paste(id_sample_short, harmonized_donor2, treat_full, str_replace_all(cell, "_", "-"), sep = "_")) %>%
  select(samples, n_cells)

# aggregate the counts per sample per treatment, donor and cell type
cts_sample_all <- AggregateExpression(object = scobj,
                                      group.by = c("id_sample_short", "harmonized_donor2", "treat_full", var_celltype),
                                      assays = "RNA",
                                      return.seurat = FALSE)
rm(scobj)
invisible(gc())

counts_all <- cts_sample_all$RNA

# sample level metadata, the donors not assigned are removed
colData_all <- data.frame(samples = colnames(counts_all)) %>%
  separate(samples, into = c("id_sample_short", "harmonized_donor2", "treat_full", "cell_id"), remove = FALSE, sep = "_") %>%
  mutate(cell_id = restore_celltype(cell_id)) %>%
  left_join(df_ncells, by = "samples") %>%
  filter(!harmonized_donor2 %in% donor_drop,
         treat_full %in% c(ref, list_treat))

colData_all %>%
  dplyr::count(cell_id, treat_full) %>%
  pivot_wider(names_from = treat_full, values_from = n, values_fill = 0) %>%
  print(n = Inf)

list_celltype <- sort(unique(colData_all$cell_id))

# functions ---------------------------------------------------------------
## sample level QC plots for one cell type: total reads, sample clustering, PCA, panel of PC vs metadata
plot_pseudobulk_qc <- function(dds, vsd, col_data, ct){
  # total reads per sample
  p_reads <- colSums(counts(dds)) %>%
    data.frame(tot_reads = .) %>%
    rownames_to_column("sample") %>%
    ggplot(aes(x = sample, y = tot_reads)) + geom_col() +
    theme_bw() +
    theme(axis.text.x = element_text(hjust = 1, angle = 90)) +
    ggtitle(ct)
  ggsave(paste0(folder_plot, "06_pBulk_totReads_", ct, ".pdf"), plot = p_reads, width = 10, height = 6)

  # clustering of the samples. the seed controls the random annotation colors
  tryCatch({
    pdf(paste0(folder_plot, "06_heatmap_SampleCluster_filterExp_", ct, ".pdf"), width = 12, height = 6)
    set.seed(21)
    hm <- plot_sample_clustering(vsd, anno_vars = c("treat_full", "harmonized_donor2", "id_sample_short", "n_cells", "M_reads"), distance = "euclidean")
    draw(hm, heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(2, 2, 2, 30), "mm"))
    dev.off()
  }, error = function(e){
    message("  sample clustering failed ", ct, ": ", conditionMessage(e))
    while(!is.null(dev.list())) dev.off()
  })

  # PCA on the 500 genes with the highest variance
  rv <- rowVars(assay(vsd), useNames = TRUE)
  select_var <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
  pca <- prcomp(t(assay(vsd)[select_var, ]))
  percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  n_pc <- min(10, ncol(pca$x))

  df_pca <- pca$x[, 1:n_pc, drop = FALSE] %>%
    data.frame() %>%
    rownames_to_column("samples") %>%
    left_join(col_data, by = "samples")

  p_pca1 <- df_pca %>%
    ggplot(aes(x = PC1, y = PC2, col = treat_full, shape = harmonized_donor2)) +
    geom_point(size = 3, alpha = 0.6) +
    scale_x_continuous(expand = expansion(mult = 0.1)) +
    scale_color_manual(values = color_treat) +
    theme_bw() + xlab(paste0("PC1: ", percent_var[1], "% variance")) + ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    ggtitle(ct)
  ggsave(paste0(folder_plot, "06_PCA_pseudobulk_filterExp_", ct, ".pdf"), plot = p_pca1, width = 7, height = 6)

  p_pca2 <- df_pca %>%
    ggplot(aes(x = PC1, y = PC2, col = treat_full, label = samples)) +
    geom_point(size = 3, alpha = 0.6) +
    geom_text_repel(show.legend = FALSE, size = 2.5) +
    scale_x_continuous(expand = expansion(mult = 0.1)) +
    scale_color_manual(values = color_treat) +
    theme_bw() + xlab(paste0("PC1: ", percent_var[1], "% variance")) + ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    ggtitle(ct)
  ggsave(paste0(folder_plot, "06_PCA_pseudobulk_filterExp2_", ct, ".pdf"), plot = p_pca2, width = 8, height = 6)

  # more PC by condition
  p_panel <- df_pca %>%
    ggpairs(columns = paste0("PC", 1:n_pc), ggplot2::aes(colour = treat_full), upper = "blank") +
    theme_bw() + theme(strip.background = element_blank(), axis.text.x = element_text(hjust = 1, angle = 45))
  ggsave(paste0(folder_plot, "06_PCA_pseudobulk_panel_", ct, ".pdf"), plot = p_panel, width = 2 * n_pc, height = 2 * n_pc)

  # PC scores by sample metadata
  df_pc_long <- df_pca %>%
    dplyr::select(samples, all_of(paste0("PC", 1:n_pc))) %>%
    pivot_longer(names_to = "var_2", values_to = "value_2", -samples) %>%
    mutate(var_2 = str_remove(var_2, "PC") %>% str_pad(side = "left", width = 2, pad = "0") %>% paste0("PC", .))

  p_cat <- df_pca %>%
    dplyr::select(samples, treat_full, harmonized_donor2, id_sample_short) %>%
    mutate(across(-samples, as.character)) %>%
    pivot_longer(names_to = "var_1", values_to = "value_1", -samples) %>%
    left_join(df_pc_long, by = "samples", relationship = "many-to-many") %>%
    mutate(comparison = paste0(var_1, "_vs_", var_2)) %>%
    ggplot(aes(x = value_1, y = value_2)) +
    facet_wrap(~comparison, scales = "free", ncol = n_pc) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.2)) +
    theme_bw() + theme(strip.background = element_blank(), axis.text.x = element_text(hjust = 1, angle = 45))
  ggsave(paste0(folder_plot, "06_panel_metadata_PC_pbulk_cat_", ct, ".pdf"), plot = p_cat, width = 3 * n_pc, height = 10, limitsize = FALSE)

  p_cont <- df_pca %>%
    dplyr::select(samples, n_cells, M_reads) %>%
    mutate(across(-samples, ~ suppressWarnings(as.numeric(.x)))) %>%
    pivot_longer(names_to = "var_1", values_to = "value_1", -samples) %>%
    left_join(df_pc_long, by = "samples", relationship = "many-to-many") %>%
    mutate(comparison = paste0(var_1, "_vs_", var_2)) %>%
    ggplot(aes(x = value_1, y = value_2)) +
    facet_wrap(~comparison, scales = "free", ncol = n_pc) +
    geom_smooth(method = "lm") +
    geom_point() +
    theme_bw() + theme(strip.background = element_blank(), axis.text.x = element_text(hjust = 1, angle = 45))
  ggsave(paste0(folder_plot, "06_panel_metadata_PC_pbulk_cont_", ct, ".pdf"), plot = p_cont, width = 3 * n_pc, height = 7, limitsize = FALSE)

  # donor removal, for visualisation only ---------------------------------
  # mimics on the plots what the donor term does in the DE design (DE uses the raw counts, never these values).
  # the design ~ treat_full is passed to protect the effect of the treatment, otherwise the treatment signal is partly removed when the donors are unbalanced between the groups (e.g. myelin has 2 donors).
  # a covariate with a single level can not be removed and is skipped
  cd_vsd <- as.data.frame(colData(vsd))
  if(nlevels(cd_vsd$harmonized_donor2) >= 2){
    tryCatch({
      vsd_batch <- vsd
      assay(vsd_batch) <- limma::removeBatchEffect(assay(vsd),
                                                   batch = cd_vsd$harmonized_donor2,
                                                   design = model.matrix(~ treat_full, data = cd_vsd))

      # clustering of the samples after the removal
      pdf(paste0(folder_plot, "06_heatmap_SampleCluster_donorRemoved_", ct, ".pdf"), width = 12, height = 6)
      set.seed(21)
      hm_batch <- plot_sample_clustering(vsd_batch, anno_vars = c("treat_full", "harmonized_donor2", "id_sample_short", "n_cells", "M_reads"), distance = "euclidean")
      draw(hm_batch, heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(2, 2, 2, 30), "mm"))
      dev.off()

      # PCA before (top) and after (bottom) the removal, coloured by treatment (left) and donor (right)
      plot_pca <- function(v, col_var, title){
        df <- plotPCA(v, intgroup = c("treat_full", "harmonized_donor2"), returnData = TRUE)
        pv <- round(100 * attr(df, "percentVar"))
        p <- ggplot(df, aes(x = PC1, y = PC2, col = .data[[col_var]], shape = harmonized_donor2)) +
          geom_point(size = 3, alpha = 0.7) +
          xlab(paste0("PC1: ", pv[1], "% variance")) + ylab(paste0("PC2: ", pv[2], "% variance")) +
          theme_bw() + ggtitle(title)
        if(col_var == "treat_full") p <- p + scale_color_manual(values = color_treat)
        p
      }

      p_batch <- plot_grid(plot_pca(vsd, "treat_full", paste0(ct, " - vst")),
                           plot_pca(vsd, "harmonized_donor2", paste0(ct, " - vst")),
                           plot_pca(vsd_batch, "treat_full", paste0(ct, " - donor removed")),
                           plot_pca(vsd_batch, "harmonized_donor2", paste0(ct, " - donor removed")),
                           ncol = 2)
      ggsave(paste0(folder_plot, "06_PCA_pseudobulk_donorRemoved_", ct, ".pdf"), plot = p_batch, width = 14, height = 11)
    }, error = function(e){
      message("  donor removal plots failed ", ct, ": ", conditionMessage(e))
      while(!is.null(dev.list())) dev.off()
    })
  }
}

## diagnostic plots of the DE results for one cell type: p-value histogram, volcano, DEG heatmap, MA plot (unshrunk and shrunk), dispersion
plot_de <- function(dds, vsd, df_res, df_res_shr, col_data, n_test, ct){
  # histogram of the p values, excluding genes with very small counts
  p_hist <- df_res %>%
    dplyr::filter(baseMean > 1) %>%
    ggplot(aes(x = pvalue)) + geom_histogram(breaks = 0:20/20) +
    facet_wrap(~conditionVsBASELINE) +
    theme_bw() +
    theme(strip.background = element_blank()) +
    ggtitle(ct)
  ggsave(paste0(folder_plot, "06_histogram_pvalue_pseudobulk_filterExp_", ct, ".pdf"), plot = p_hist, width = 3 * n_test + 1, height = 4)

  # volcano, on the unshrunk (res) and on the shrunk (shr) estimates
  walk2(list(df_res, df_res_shr), c("", "shr_"), function(df, tag){
    title_tag <- ifelse(tag == "", "unshrunk", "shrunk")
    df_plot <- df %>%
      mutate(col = ifelse(((padj < padj_thr) & abs(log2FoldChange) > lfc_thr & !is.na(symbol)), yes = 1, no = 0)) %>%
      filter(!is.na(col))

    p_volcano <- df_plot %>%
      ggplot(aes(x = log2FoldChange, y = -log(padj))) +
      geom_point(data = df_plot[df_plot$col == 0, ], aes(col = factor(col)), alpha = 0.05) +
      geom_point(data = df_plot[df_plot$col == 1, ], aes(col = factor(col)), alpha = 0.5) +
      geom_vline(xintercept = c(-lfc_thr, lfc_thr), col = "red", linetype = "dashed") +
      geom_hline(yintercept = (-log(padj_thr)), col = "red", linetype = "dashed") +
      scale_color_manual(values = c("0" = "black", "1" = "red")) +
      geom_text_repel(
        data = df_plot %>%
          group_by(conditionVsBASELINE) %>%
          arrange(padj) %>%
          dplyr::slice(1:30) %>%
          dplyr::filter(abs(log2FoldChange) > lfc_thr) %>%
          dplyr::filter(padj < padj_thr),
        aes(label = symbol), segment.alpha = 0.4) +
      facet_wrap(~conditionVsBASELINE) +
      theme_bw() +
      theme(strip.background = element_blank(), legend.position = "none") +
      ggtitle(paste0(ct, " - ", title_tag))
    ggsave(paste0(folder_plot, "06_volcano_pseudobulk_filterExp_", tag, ct, ".pdf"), plot = p_volcano, width = 6 * min(n_test, 3) + 1, height = 6 * ceiling(n_test / 3) + 1)
  })

  # heatmap of the DEGs and MA plot, on the unshrunk (res) and on the shrunk (shr) estimates
  walk2(list(df_res, df_res_shr), c("", "shr_"), function(df, tag){
    title_tag <- ifelse(tag == "", "unshrunk", "shrunk")
    DEG_2 <- df %>%
      dplyr::filter(padj < padj_thr, abs(log2FoldChange) > lfc_thr) %>%
      pull(symbol) %>%
      unique()

    if(length(DEG_2) > 1){
      mat <- assay(vsd)[rownames(vsd) %in% DEG_2, , drop = FALSE]
      mat2 <- (mat - rowMeans(mat)) / rowSds(mat, useNames = TRUE)
      mat2 <- mat2[!is.na(rowSums(mat2)), , drop = FALSE]

      column_ha <- HeatmapAnnotation(treat = as.character(col_data[colnames(mat2), "treat_full"]),
                                     donor = as.character(col_data[colnames(mat2), "harmonized_donor2"]),
                                     col = list(treat = color_treat))

      ht <- Heatmap(mat2, show_column_names = TRUE, show_row_dend = FALSE, use_raster = TRUE,
                    name = "exp",
                    column_title = paste0(ct, " - ", title_tag),
                    row_names_gp = gpar(fontsize = ifelse(length(DEG_2) > 100, 3, 8)),
                    top_annotation = column_ha, show_row_names = TRUE)
      pdf(paste0(folder_plot, "06_heatmap_DEG_pseudobulk_filterExp_", tag, ct, ".pdf"), width = 10, height = 8)
      draw(ht, heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(30, 2, 2, 2), "mm"))
      dev.off()
    }

    p_ma <- df %>%
      filter(!is.na(padj)) %>%
      mutate(color = factor(ifelse(padj < padj_thr, yes = 1, no = 0), levels = c(0, 1))) %>%
      arrange(color) %>%
      ggplot(aes(x = baseMean, y = log2FoldChange, col = color)) +
      geom_point(alpha = 0.2) +
      scale_x_log10() +
      scale_color_manual(values = c("0" = "gray", "1" = "red"), drop = FALSE) +
      facet_wrap(~conditionVsBASELINE) +
      theme_bw() +
      theme(strip.background = element_blank(), legend.position = "none") +
      ggtitle(paste0(ct, " - ", title_tag))
    ggsave(paste0(folder_plot, "06_MA_pseudobulk_filterExp_", tag, ct, ".pdf"), plot = p_ma, width = 6 * min(n_test, 3) + 1, height = 6 * ceiling(n_test / 3) + 1)
  })

  # dispersion
  pdf(paste0(folder_plot, "06_dispersion_pseudobulk_filterExp_", ct, ".pdf"), width = 5, height = 5)
  plotDispEsts(dds)
  dev.off()
}

# run the analysis per cell type -----------------------------------------
# ct <- "MG"
list_out <- map(list_celltype, function(ct){
  colData_ct <- colData_all %>%
    filter(cell_id == ct)
  n_rep <- table(colData_ct$treat_full)

  # keep only the treatments with enough samples
  treat_ok <- intersect(list_treat, names(n_rep)[n_rep >= min_samples])
  if(!ref %in% names(n_rep) || n_rep[[ref]] < min_samples || length(treat_ok) == 0){
    message("skip ", ct, ": not enough samples")
    return(NULL)
  }
  message("run ", ct, ": ", paste(treat_ok, collapse = ", "))

  tryCatch({
    counts_ct <- as.matrix(counts_all[, colData_ct %>% filter(treat_full %in% c(ref, treat_ok)) %>% pull(samples)])
    colData_ct <- colData_ct %>%
      filter(treat_full %in% c(ref, treat_ok)) %>%
      mutate(treat_full = factor(treat_full, levels = c(ref, treat_ok)),
             harmonized_donor2 = factor(harmonized_donor2),
             M_reads = colSums(counts_ct)[samples] / 10^6)
    rownames(colData_ct) <- colData_ct$samples

    # perform DESeq2 ----------------------------------------------------------
    dds <- DESeqDataSetFromMatrix(countData = counts_ct,
                                  colData = colData_ct,
                                  design = ~ harmonized_donor2 + treat_full)

    # filter the low expressed genes (using the design of the model)
    dds_filter <- dds[edgeR::filterByExpr(counts(dds), design = model.matrix(design(dds), colData(dds))), ]
    message("  genes kept by filterByExpr: ", nrow(dds_filter), " / ", nrow(dds))

    # scale the data (vst needs > 1000 genes, otherwise use the slower function)
    vsd <- if(nrow(dds_filter) > 1000) vst(dds_filter, blind = TRUE) else varianceStabilizingTransformation(dds_filter, blind = TRUE)

    plot_pseudobulk_qc(dds_filter, vsd, colData_ct, ct)

    # run DE ------------------------------------------------------------------
    dds_de <- DESeq(dds_filter)
    saveRDS(dds_de, paste0(folder_object, "06_ddsHTSeq_pbulk_filterExp_", ct, "_ref", ref, ".rds"))

    # results table of each treatment vs the reference
    list_res <- map(treat_ok, function(tr){
      results(dds_de, contrast = c("treat_full", tr, ref), alpha = padj_thr)
    }) %>%
      set_names(treat_ok)
    walk2(list_res, names(list_res), function(x, tr){ message(ct, " ", tr); summary(x) })

    to_table <- function(l){
      l %>%
        map(function(x){
          x %>%
            data.frame() %>%
            rownames_to_column(var = "symbol") %>%
            arrange(pvalue)
        }) %>%
        bind_rows(.id = "conditionVsBASELINE") %>%
        mutate(cell_id = ct)
    }
    df_res <- to_table(list_res)
    df_res %>%
      write_tsv(paste0(folder_table, "06_DE_pbulk_filterExp_treatvsBASELINE_", ct, ".tsv"))

    # shrink ------------------------------------------------------------------
    list_res_shr <- map(list_res, function(x){
      lfcShrink(dds_de, res = x, type = "ashr")
    })
    df_res_shr <- to_table(list_res_shr)
    df_res_shr %>%
      write_tsv(paste0(folder_table, "06_DE_pbulk_filterExp_treatvsBASELINE_", ct, "_shr.tsv"))

    plot_de(dds_de, vsd, df_res, df_res_shr, colData_ct, length(treat_ok), ct)

    # upset plot --------------------------------------------------------------
    # genes shared between the treatments. the treatments without DEGs are dropped, the plot needs at least two sets
    walk(c("up", "down"), function(direction){
      list_deg <- df_res_shr %>%
        filter(padj < padj_thr, if(direction == "up") log2FoldChange > lfc_thr else log2FoldChange < -lfc_thr) %>%
        split(f = .$conditionVsBASELINE) %>%
        map(~ unique(.x$symbol))
      list_deg <- list_deg[lengths(list_deg) > 0]

      if(length(list_deg) >= 2){
        # ComplexUpset needs a table with one row per gene and one logical column per treatment
        genes_all <- unique(unlist(list_deg))
        df_upset <- map(list_deg, ~ genes_all %in% .x) %>%
          as.data.frame(check.names = FALSE) %>%
          mutate(gene = genes_all, .before = 1)

        p_upset <- ComplexUpset::upset(df_upset,
                                       intersect = names(list_deg),
                                       name = paste0(direction, " DEGs"),
                                       min_size = 1,
                                       width_ratio = 0.2,
                                       sort_intersections_by = "cardinality",
                                       sort_intersections = "descending") +
          ggtitle(paste0(ct, " - ", direction, " DEGs vs ", ref))
        ggsave(paste0(folder_plot, "06_upset_DEG_", toupper(direction), "_pbulk_", ct, "_ref", ref, "_shr.pdf"),
               plot = p_upset, width = 14, height = 7)
      }
    })

    df_res_shr
  }, error = function(e){
    message("failed ", ct, ": ", conditionMessage(e))
    while(!is.null(dev.list())) dev.off()
    NULL
  })
})

# summary of the number of DEGs ------------------------------------------
df_shr_all <- bind_rows(list_out)

df_shr_all %>%
  write_tsv(paste0(folder_table, "06_DE_pbulk_filterExp_treatvsBASELINE_", var_celltype, "_shr.tsv"))

# samples per comparison
df_n_samples <- colData_all %>%
  filter(treat_full %in% list_treat) %>%
  dplyr::count(cell_id, conditionVsBASELINE = treat_full, name = "n_samples_1") %>%
  left_join(colData_all %>%
              filter(treat_full == ref) %>%
              dplyr::count(cell_id, name = "n_samples_2"), by = "cell_id")

df_summary <- df_shr_all %>%
  group_by(cell_id, conditionVsBASELINE) %>%
  summarise(n_up = sum(padj < padj_thr & log2FoldChange > lfc_thr, na.rm = TRUE),
            n_down = sum(padj < padj_thr & log2FoldChange < -lfc_thr, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(n_deg = n_up + n_down) %>%
  left_join(df_n_samples, by = c("cell_id", "conditionVsBASELINE"))

df_summary %>%
  write_tsv(paste0(folder_table, "06_DE_pbulk_summary_", var_celltype, ".tsv"))

df_summary %>%
  print(n = Inf)

# barplot of the number of DEGs: up and down in the two directions, one panel per treatment
df_summary %>%
  pivot_longer(cols = c(n_up, n_down), names_to = "direction", values_to = "n") %>%
  # mutate(direction = str_remove(direction, "n_"),
  #        n = ifelse(direction == "down", -n, n)) %>%
  mutate(direction = str_remove(direction, "n_"),cell_id = fct_reorder(cell_id, n)) %>%
  ggplot(aes(y = cell_id, x = n, fill = direction)) +
  geom_col() +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = c(up = "firebrick", down = "steelblue")) +
  facet_wrap(~conditionVsBASELINE, scales = "free_y") +
  theme_bw() +
  theme(strip.background = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = paste0("number of DEGs vs ", ref, " (pseudobulk, shrunken lfc)"),
       subtitle = paste0("padj < ", padj_thr, ", |log2FC| > ", lfc_thr),
       x = NULL, y = "n DEGs")
ggsave(paste0(folder_plot, "06_barplot_nDEG_pbulk_", var_celltype, "_vs_", ref, ".pdf"), width = 12, height = 8)

# same, by treatment on the x axis, one panel per cell type
df_summary %>%
  ggplot(aes(x = conditionVsBASELINE, y = n_deg)) +
  geom_col(fill = "grey30") +
  geom_text(aes(label = n_deg), vjust = -0.3, size = 2.5) +
  facet_wrap(~cell_id, scales = "free_y") +
  theme_bw() +
  theme(strip.background = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = paste0("number of DEGs vs ", ref, " (pseudobulk, shrunken lfc)"), x = NULL, y = "n DEGs")
ggsave(paste0(folder_plot, "06_barplot_nDEG_pbulk_perCellType_", var_celltype, "_vs_", ref, ".pdf"), width = 12, height = 8)

# single cell vs pseudobulk -----------------------------------------------
if(!file.exists(file_sc)){
  message("no single cell table found (", file_sc, "): skip the comparison with the pseudobulk")
} else {
  key <- c("cell_id", "treat")

  df_sc <- read_tsv(file_sc, show_col_types = FALSE) %>%
    dplyr::rename(treat = ident.1)
  df_pb <- df_shr_all %>%
    dplyr::rename(treat = conditionVsBASELINE, gene = symbol)

  # comparisons tested by both approaches
  groups_both <- inner_join(distinct(df_sc, across(all_of(key))),
                            distinct(df_pb, across(all_of(key))),
                            by = key)

  # DEGs with the same thresholds, element = gene|direction
  df_deg_both <- bind_rows(
    df_sc %>%
      filter(p_val_adj < padj_thr, abs(avg_log2FC) > lfc_thr) %>%
      mutate(direction = ifelse(avg_log2FC > 0, "up", "down"), method = "single_cell") %>%
      select(all_of(key), gene, direction, method),
    df_pb %>%
      filter(padj < padj_thr, abs(log2FoldChange) > lfc_thr) %>%
      mutate(direction = ifelse(log2FoldChange > 0, "up", "down"), method = "pseudobulk") %>%
      select(all_of(key), gene, direction, method)
  ) %>%
    semi_join(groups_both, by = key) %>%
    mutate(element = paste(gene, direction, sep = "|"))

  df_deg_both %>%
    write_tsv(paste0(folder_table, "06_DEG_singleCell_and_pseudobulk.tsv"))

  # overlap table, one row per comparison
  df_overlap <- df_deg_both %>%
    group_by(across(all_of(key))) %>%
    summarise(n_single_cell = sum(method == "single_cell"),
              n_pseudobulk = sum(method == "pseudobulk"),
              n_shared = length(intersect(element[method == "single_cell"], element[method == "pseudobulk"])),
              .groups = "drop") %>%
    right_join(groups_both, by = key) %>%
    replace_na(list(n_single_cell = 0, n_pseudobulk = 0, n_shared = 0)) %>%
    mutate(jaccard = n_shared / pmax(n_single_cell + n_pseudobulk - n_shared, 1))

  df_overlap %>%
    write_tsv(paste0(folder_table, "06_overlap_DEG_singleCell_vs_pseudobulk.tsv"))

  # one pdf per treatment, one upsetplot per cell type (one page each)
  df_deg_both %>%
    split(f = .$treat) %>%
    iwalk(function(df, tr){
      file_id <- paste0(var_celltype, "_", tr, "_vs_", ref)
      pdf(paste0(folder_plot, "06_upset_DEG_singleCell_vs_pseudobulk_", file_id, ".pdf"), width = 6, height = 4)
      walk(sort(unique(df$cell_id)), function(cl){
        list_set <- list(single_cell = df %>% filter(cell_id == cl, method == "single_cell") %>% pull(element),
                         pseudobulk = df %>% filter(cell_id == cl, method == "pseudobulk") %>% pull(element))
        tryCatch({
          comb_mat <- make_comb_mat(list_set)
          # the number of genes is added on top of the intersection bars and at the end of the set bars
          draw(UpSet(comb_mat, comb_order = order(comb_size(comb_mat), decreasing = TRUE),
                     top_annotation = upset_top_annotation(comb_mat, add_numbers = TRUE),
                     right_annotation = upset_right_annotation(comb_mat, add_numbers = TRUE),
                     column_title = paste0(file_id, " - ", cl)))
        }, error = function(e) message("  upset failed ", file_id, " ", cl, ": ", conditionMessage(e)))
      })
      dev.off()
    })

  # correlation single cell vs pseudobulk ------------------------------------
  # genes present in both full tables (no threshold). NOTICE: the single cell table only contains the genes returned by FindMarkers (min.pct filter), so the correlation is computed on the genes tested by both approaches
  # the pseudobulk LFC is the ashr shrunk one. The significance is -log10 of the adjusted p-value, floored to avoid Inf for p_val_adj == 0
  padj_floor <- 1e-300

  df_cor <- df_sc %>%
    semi_join(groups_both, by = key) %>%
    select(all_of(key), gene, avg_log2FC, p_val_adj) %>%
    inner_join(df_pb %>% select(all_of(key), gene, log2FoldChange, padj), by = c(key, "gene")) %>%
    filter(!is.na(avg_log2FC), !is.na(log2FoldChange), !is.na(p_val_adj), !is.na(padj)) %>%
    mutate(neglog10_padj_sc = -log10(pmax(p_val_adj, padj_floor)),
           neglog10_padj_pbulk = -log10(pmax(padj, padj_floor)))

  df_cor_summary <- df_cor %>%
    group_by(across(all_of(key))) %>%
    filter(n() >= 3) %>%
    summarise(n_genes = n(),
              pearson_lfc = cor(avg_log2FC, log2FoldChange, method = "pearson"),
              spearman_lfc = cor(avg_log2FC, log2FoldChange, method = "spearman"),
              pearson_padj = cor(neglog10_padj_sc, neglog10_padj_pbulk, method = "pearson"),
              spearman_padj = cor(neglog10_padj_sc, neglog10_padj_pbulk, method = "spearman"),
              .groups = "drop")

  df_cor_summary %>%
    write_tsv(paste0(folder_table, "06_correlation_singleCell_vs_pseudobulk.tsv"))

  df_cor_summary %>%
    print(n = Inf)

  df_lab_cor <- df_cor_summary %>%
    mutate(lab_lfc = paste0("r = ", round(pearson_lfc, 2), "\nrho = ", round(spearman_lfc, 2), "\nn = ", n_genes),
           lab_padj = paste0("r = ", round(pearson_padj, 2), "\nrho = ", round(spearman_padj, 2), "\nn = ", n_genes))

  # one pair of scatter plots per treatment, faceted by cell type
  df_cor %>%
    split(f = .$treat) %>%
    iwalk(function(df, tr){
      file_id <- paste0(var_celltype, "_", tr, "_vs_", ref)
      df_lab <- df_lab_cor %>%
        semi_join(df %>% distinct(across(all_of(key))), by = key)
      n_ct <- n_distinct(df$cell_id)
      n_col <- ceiling(sqrt(n_ct))
      n_row <- ceiling(n_ct / n_col)

      p_lfc <- df %>%
        ggplot(aes(x = avg_log2FC, y = log2FoldChange)) +
        geom_hline(yintercept = 0, color = "grey70") +
        geom_vline(xintercept = 0, color = "grey70") +
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
        geom_point(alpha = 0.2, size = 0.7) +
        geom_text(data = df_lab, aes(label = lab_lfc), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, size = 3, inherit.aes = FALSE) +
        facet_wrap(~cell_id, scales = "free", ncol = n_col) +
        theme_bw() +
        theme(strip.background = element_blank()) +
        labs(title = paste0("log2FC single cell vs pseudobulk - ", file_id),
             subtitle = "r = Pearson, rho = Spearman",
             x = "single cell avg_log2FC", y = "pseudobulk log2FoldChange (shrunk)")
      ggsave(paste0(folder_plot, "06_scatter_log2FC_singleCell_vs_pseudobulk_", file_id, ".pdf"),
             plot = p_lfc, width = 4 * n_col + 1, height = 3.5 * n_row + 1, limitsize = FALSE)

      p_padj <- df %>%
        ggplot(aes(x = neglog10_padj_sc, y = neglog10_padj_pbulk)) +
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
        geom_point(alpha = 0.2, size = 0.7) +
        geom_text(data = df_lab, aes(label = lab_padj), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, size = 3, inherit.aes = FALSE) +
        facet_wrap(~cell_id, scales = "free", ncol = n_col) +
        theme_bw() +
        theme(strip.background = element_blank()) +
        labs(title = paste0("adjusted p-value single cell vs pseudobulk - ", file_id),
             subtitle = paste0("r = Pearson, rho = Spearman, p-values floored at ", padj_floor),
             x = "single cell -log10(p_val_adj)", y = "pseudobulk -log10(padj)")
      ggsave(paste0(folder_plot, "06_scatter_padj_singleCell_vs_pseudobulk_", file_id, ".pdf"),
             plot = p_padj, width = 4 * n_col + 1, height = 3.5 * n_row + 1, limitsize = FALSE)
    })
}

message("\nAll done.")
