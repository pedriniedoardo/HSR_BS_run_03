# AIM ---------------------------------------------------------------------
# run GSVA on the MG subset of the BS run 03 dataset

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

exp <- counts(ddsHTSeq_filter,normalized = T) %>%
  data.frame() %>%
  dplyr::select(contains("BASELINE")|contains("Fe")|contains("myelin")) %>%
  as.matrix()

# read in the list of terms provided by Aletta
TOI <- read_csv("../../data/GOterms_iron_myelin.csv")

# generate the signature file

gene_sets <- msigdbr(species = "Homo sapiens", category = "C5",subcategory = "GO:BP")
head(gene_sets)


test01_01_summary <- gene_sets %>%
  group_by(gs_exact_source,gs_name) %>%
  summarise(n = n())

#
test01_01_summary

# check that the terms provided are present in the current dataset
TOI %>%
  left_join(test01_01_summary,by = c("GO_id" = "gs_exact_source")) %>%
  mutate(test = !is.na(n)) %>%
  print(n = 30)

# format in order to be accepted by GSEA
pathways <- split(x = gene_sets$gene_symbol, f = gene_sets$gs_exact_source)

# -------------------------------------------------------------------------
# set.seed for reproducibility
set.seed(123)

# perform the GSEA based on the normalized counts
es <- gsva(exp,
           pathways,
           min.sz=5,
           max.sz=500,
           kcdf="Poisson",
           mx.diff=TRUE,
           verbose=FALSE,
           parallel.sz=1)

es_log <- gsva(log(exp+1),
               pathways,
               min.sz=5,
               max.sz=500,
               kcdf="Gaussian",
               mx.diff=TRUE,
               verbose=FALSE,
               parallel.sz=1)

# show correlation between the two estimates ------------------------------
left_join(
  es %>%
    data.frame() %>%
    rownames_to_column("pathway") %>%
    pivot_longer(names_to = "sample",values_to = "estimate",-pathway),
  es_log %>%
    data.frame() %>%
    rownames_to_column("pathway") %>%
    pivot_longer(names_to = "sample",values_to = "estimate_log",-pathway),by = c("pathway","sample")) %>%
  ggplot(aes(x=estimate_log,y=estimate))+geom_point()

# STATISTICAL TESTING -----------------------------------------------------
# run the analysis only on the logged normalized values

# use limma for testing significance
# library(limma)
# adjPvalueCutoff <- 0.001
# logFCcutoff <- log2(2)
# logFCcutoff

# define the factor for the treatmentnt basesd on the colnames
lut <- data.frame(colname = colnames(es)) %>%
  mutate(treat = unlist(str_extract_all(colname,pattern = "BASELINE|Fe|myelin"))) %>%
  mutate(treat = factor(treat,levels = c("BASELINE","Fe","myelin")))

design <- model.matrix(~ lut$treat)
colnames(design)[1] <- c("intercept")
design

# fit <- lmFit(es, design)
fit_log <- lmFit(es_log, design)
# fit <- eBayes(fit)
fit_log <- eBayes(fit_log)

# allGenesets_Fe <- topTable(fit, coef="lut$treatFe", number=Inf)
allGenesets_log_Fe <- topTable(fit_log, coef="lut$treatFe", number=Inf)

# allGenesets_myelin <- topTable(fit, coef="lut$treatmyelin", number=Inf)
allGenesets_log_myelin <- topTable(fit_log, coef="lut$treatmyelin", number=Inf)

# pull the stats from Aletta's list
allGenesets_log_Fe %>%
  rownames_to_column("GO_id") %>%
  left_join(TOI,by = "GO_id") %>%
  dplyr::filter(!is.na(GO_term))

allGenesets_log_myelin %>%
  rownames_to_column("GO_id") %>%
  left_join(TOI,by = "GO_id") %>%
  dplyr::filter(!is.na(GO_term))

# save the tables
allGenesets_log_Fe %>%
  rownames_to_column("GO_id") %>%
  left_join(test01_01_summary,by = c("GO_id"="gs_exact_source")) %>%
  write_tsv("../../out/table/24_df_table_GSVA_GOBP_Fe_MG.tsv")

allGenesets_log_myelin %>%
  rownames_to_column("GO_id") %>%
  left_join(test01_01_summary,by = c("GO_id"="gs_exact_source")) %>%
  write_tsv("../../out/table/24_df_table_GSVA_GOBP_myelin_MG.tsv")

# volcano GSVA ------------------------------------------------------------
df_plot_Fe <- allGenesets_log_Fe %>%
  rownames_to_column("pathway")

df_plot_myelin <- allGenesets_log_myelin %>%
  rownames_to_column("pathway")

df_plot_Fe2 <- df_plot_Fe %>%   
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = "HALLMARK_|KEGG_|GOBP_") %>%
           str_sub(start = 1,end = 35))
  # mutate(color = case_when(pathway %in% c("KEGG_PRIMARY_IMMUNODEFICIENCY",
  #                                       "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
  #                                       "KEGG_VEGF_SIGNALING_PATHWAY",
  #                                       "KEGG_CELL_CYCLE",
  #                                       "KEGG_CELL_ADHESION_MOLECULES_CAMS",
  #                                       "KEGG_ECM_RECEPTOR_INTERACTION")~"red",
  #                        T~"black")) %>%
  # mutate(color = factor(color))

df_plot_myelin2 <- df_plot_myelin %>%   
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = "HALLMARK_|KEGG_|GOBP_") %>%
           str_sub(start = 1,end = 35))
  # mutate(color = case_when(pathway %in% c("KEGG_PRIMARY_IMMUNODEFICIENCY",
  #                                     "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
  #                                     "KEGG_VEGF_SIGNALING_PATHWAY",
  #                                     "KEGG_CELL_CYCLE",
  #                                     "KEGG_CELL_ADHESION_MOLECULES_CAMS",
  #                                     "KEGG_ECM_RECEPTOR_INTERACTION")~"red",
  #                      T~"black")) %>%
  # mutate(color = factor(color))

df_plot_Fe2 %>%
  # ggplot(aes(y = -log10(adj.P.Val),x = logFC,label = pathway2,col=color)) +
  ggplot(aes(y = -log10(adj.P.Val),x = logFC)) +
  geom_point(alpha = 0.2) +
  # geom_point(aes(size = size),alpha = 0.2) +
  # facet_wrap(~dataset) +
  # theme_bw(base_rect_size = 2)+
  theme_bw()+
  # theme(strip.background = element_blank(), 
  #       panel.border = element_rect(colour = "black", fill = NA),
  #       axis.ticks = element_line(colour = "black"),
  #       #axis.ticks.length = unit(.25, "cm")
  #       legend.position = "none"
  # )+
  theme(strip.background = element_blank(),legend.position = "none") +
  geom_hline(yintercept = -log10(0.05),linetype="dashed",col="gray",alpha=0.8)
  # geom_text_repel(data = df_plot_Fe2 %>% filter(color == "black"),aes(-log10(adj.P.Val),x = logFC,label = pathway2,col=color),
  #                 size = 2,
  #                 box.padding = 0.5,
  #                 segment.alpha = 0.5,
  #                 max.overlaps = 10)+
  # geom_text_repel(data = df_plot_Fe2 %>% filter(color == "red"),aes(-log10(adj.P.Val),x = logFC,label = pathway2,col=color),
  #                 size = 2,
  #                 box.padding = 0.5,
  #                 segment.alpha = 0.5,
  #                 max.overlaps = 10,min.segment.length = 0,nudge_x = 0.2,nudge_y = 0.2)+
  # scale_color_manual(values = levels(df_plot$color))
# ggsave("../../out/image/volcano_GSVA_KEGG_GMPVsResearch.pdf",width = 12,height = 6)

df_plot_myelin2 %>%
  # ggplot(aes(y = -log10(adj.P.Val),x = logFC,label = pathway2,col=color)) +
  ggplot(aes(y = -log10(adj.P.Val),x = logFC)) +
  geom_point(alpha = 0.2) +
  # geom_point(aes(size = size),alpha = 0.2) +
  # facet_wrap(~dataset) +
  # theme_bw(base_rect_size = 2)+
  theme_bw()+
  # theme(strip.background = element_blank(), 
  #       panel.border = element_rect(colour = "black", fill = NA),
  #       axis.ticks = element_line(colour = "black"),
  #       #axis.ticks.length = unit(.25, "cm")
  #       legend.position = "none"
  # )+
  theme(strip.background = element_blank(),legend.position = "none") +
  geom_hline(yintercept = -log10(0.05),linetype="dashed",col="gray",alpha=0.8)
# geom_text_repel(data = df_plot_myelin2 %>% filter(color == "black"),aes(-log10(adj.P.Val),x = logFC,label = pathway2,col=color),
#                 size = 2,
#                 box.padding = 0.5,
#                 segment.alpha = 0.5,
#                 max.overlaps = 10)+
# geom_text_repel(data = df_plot_myelin2 %>% filter(color == "red"),aes(-log10(adj.P.Val),x = logFC,label = pathway2,col=color),
#                 size = 2,
#                 box.padding = 0.5,
#                 segment.alpha = 0.5,
#                 max.overlaps = 10,min.segment.length = 0,nudge_x = 0.2,nudge_y = 0.2)+
# scale_color_manual(values = levels(df_plot$color))
# ggsave("../../out/image/volcano_GSVA_KEGG_GMPVsResearch.pdf",width = 12,height = 6)

# plot heatmap ------------------------------------------------------------
# library(ComplexHeatmap)
# define only the significnat terms
DEgeneSets_Fe <- allGenesets_log_Fe %>%
  rownames_to_column("pathway") %>%
  filter(adj.P.Val < 0.05) %>%
  pull(pathway) 
  # # force in the one of interest
  # c("KEGG_PRIMARY_IMMUNODEFICIENCY",
  #   "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
  #   "KEGG_VEGF_SIGNALING_PATHWAY",
  #   "KEGG_CELL_CYCLE",
  #   "KEGG_CELL_ADHESION_MOLECULES_CAMS",
  #   "KEGG_ECM_RECEPTOR_INTERACTION")

DEgeneSets_myelin <- allGenesets_log_myelin %>%
  rownames_to_column("pathway") %>%
  filter(adj.P.Val < 0.05) %>%
  pull(pathway) 
# # force in the one of interest
# c("KEGG_PRIMARY_IMMUNODEFICIENCY",
#   "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
#   "KEGG_VEGF_SIGNALING_PATHWAY",
#   "KEGG_CELL_CYCLE",
#   "KEGG_CELL_ADHESION_MOLECULES_CAMS",
#   "KEGG_ECM_RECEPTOR_INTERACTION")

# mat <- es

# define the matrix for the heatmap
mat_norm_Fe <- es %>%
  data.frame() %>%
  rownames_to_column() %>%
  # plot only the significant terms
  filter(rowname %in% DEgeneSets_Fe) %>%
  # scale the values rowwise
  gather(key = sample,value = exp,-rowname) %>%
  group_by(rowname) %>%
  mutate(norm = (exp - mean(exp))/sd(exp)) %>%
  dplyr::select(-exp) %>%
  spread(key = sample,value = norm) %>%
  column_to_rownames()

mat_norm_myelin <- es %>%
  data.frame() %>%
  rownames_to_column() %>%
  # plot only the significant terms
  filter(rowname %in% DEgeneSets_myelin) %>%
  # scale the values rowwise
  gather(key = sample,value = exp,-rowname) %>%
  group_by(rowname) %>%
  mutate(norm = (exp - mean(exp))/sd(exp)) %>%
  dplyr::select(-exp) %>%
  spread(key = sample,value = norm) %>%
  column_to_rownames()

sample_ordered_Fe <- str_extract(colnames(mat_norm_Fe),pattern = "BASELINE|Fe|myelin")
sample_ordered_myelin <- str_extract(colnames(mat_norm_myelin),pattern = "BASELINE|Fe|myelin")
sample_ordered_Fe
sample_ordered_myelin

# build the annotation object
column_ha_Fe <- HeatmapAnnotation(treat = sample_ordered_Fe,
                                  col = list(treat = c("BASELINE" = "blue", "Fe"="orange","myelin"= "cyan"))) 
column_ha_myelin <- HeatmapAnnotation(treat = sample_ordered_myelin,
                                  col = list(treat = c("BASELINE" = "blue", "Fe"="orange","myelin"= "cyan"))) 

# highlight some terms
# rowname_color <- case_when(rownames(mat_norm)%in%c("KEGG_PRIMARY_IMMUNODEFICIENCY",
#                                                    "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
#                                                    "KEGG_VEGF_SIGNALING_PATHWAY",
#                                                    "KEGG_CELL_CYCLE",
#                                                    "KEGG_CELL_ADHESION_MOLECULES_CAMS",
#                                                    "KEGG_ECM_RECEPTOR_INTERACTION")~"red",
#                            T~"black")

hm_Fe <- Heatmap(mat_norm_Fe,show_row_names = F,
              # add annotation for the columns
              # hide columns labels
              # show_column_names = F,
              # fix width of the lables
              top_annotation = column_ha_Fe,
              # row_names_gp = gpar(col = rowname_color),
              row_names_max_width = max_text_width(
                rownames(mat_norm_Fe),
                gp = gpar(fontsize = 12)
              ))

pdf(file = "../../out/image/24_heatmap_GSVA_GOBP_Fe_es_Log_ZScore.pdf", width = 8, height = 5)
draw(hm_Fe,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2, 2, 2, 100), "mm"))
dev.off()

hm_myelin <- Heatmap(mat_norm_myelin,show_row_names = F,
                 # add annotation for the columns
                 # hide columns labels
                 # show_column_names = F,
                 # fix width of the lables
                 top_annotation = column_ha_myelin,
                 # row_names_gp = gpar(col = rowname_color),
                 row_names_max_width = max_text_width(
                   rownames(mat_norm_myelin),
                   gp = gpar(fontsize = 12)
                 ))

pdf(file = "../../out/image/24_heatmap_GSVA_GOBP_myelin_es_Log_ZScore.pdf", width = 8, height = 5)
draw(hm_myelin,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2, 2, 2, 100), "mm"))
dev.off()
