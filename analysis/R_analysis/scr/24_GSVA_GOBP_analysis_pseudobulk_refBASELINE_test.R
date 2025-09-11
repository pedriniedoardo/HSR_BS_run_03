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
ah <- AnnotationHub()
human_orgdb_record <- query(ah, c("OrgDb", "Homo sapiens"))
# Retrieve the database object using its ID
orgdb <- human_orgdb_record[["AH107059"]]

# explore the keytypes for the database
keytypes(orgdb)

# pull all the genes from the specific terms
test03_02 <- AnnotationDbi::select(orgdb,
                                   keys = unique(TOI$GO_id),
                                   columns = c('SYMBOL'),
                                   keytype = "GOALL") %>%
  # filter the NA genes
  filter(!is.na(SYMBOL))

# summarise the number of genes per term
test03_02_summary <- test03_02 %>%
  group_by(GOALL) %>%
  summarise(n = n())

TOI %>%
  left_join(test03_02_summary,by = c("GO_id" = "GOALL")) %>%
  mutate(test = !is.na(n)) %>%
  print(n = 30)

TOI %>%
  left_join(test03_02_summary,by = c("GO_id" = "GOALL")) %>%
  mutate(test = !is.na(n)) %>%
  filter(n > 5)
  print(n = 30)


pathways <- split(x = test03_02$SYMBOL, f = test03_02$GOALL)

# -------------------------------------------------------------------------
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

fit <- lmFit(es, design)
fit_log <- lmFit(es_log, design)
fit <- eBayes(fit)
fit_log <- eBayes(fit_log)

allGenesets_Fe <- topTable(fit, coef="lut$treatFe", number=Inf)
allGenesets_log_Fe <- topTable(fit_log, coef="lut$treatFe", number=Inf)

allGenesets_myelin <- topTable(fit, coef="lut$treatmyelin", number=Inf)
allGenesets_log_myelin <- topTable(fit_log, coef="lut$treatmyelin", number=Inf)

# volcano GSVA ------------------------------------------------------------
df_plot <- allGenesets_Fe %>%
  rownames_to_column("pathway")

df_plot2 <- df_plot %>%   
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = "HALLMARK_|KEGG_|GOBP_") %>%
           str_sub(start = 1,end = 35))

df_plot2 %>%
  ggplot(aes(y = -log(adj.P.Val),x = logFC,label = pathway2)) +
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
  theme(strip.background = element_blank(),legend.position = "none")+
  geom_text_repel(aes(-log(adj.P.Val),x = logFC,label = pathway2),
                  size = 2,
                  box.padding = 0.5,
                  segment.alpha = 0.5,
                  max.overlaps = 10)+
  geom_hline(yintercept = -log(0.05),linetype="dashed",col="gray",alpha=0.8)
ggsave("../../out/image/volcano_GSVA_KEGG_GMPVsResearch.pdf",width = 12,height = 6)

# plot heatmap ------------------------------------------------------------
# library(ComplexHeatmap)
# define only the significnat terms
DEgeneSets <- allGenesets_Fe %>%
  rownames_to_column("pathway") %>%
  filter(adj.P.Val < 0.05) %>%
  pull(pathway) %>%
  # force in the one of interest
  c("KEGG_PRIMARY_IMMUNODEFICIENCY",
    "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
    "KEGG_VEGF_SIGNALING_PATHWAY",
    "KEGG_CELL_CYCLE",
    "KEGG_CELL_ADHESION_MOLECULES_CAMS",
    "KEGG_ECM_RECEPTOR_INTERACTION")

mat <- es
mat_norm <- es %>%
  data.frame() %>%
  rownames_to_column() %>%
  # plot only the significant terms
  # filter(rowname %in% DEgeneSets) %>%
  # scale the values rowwise
  gather(key = sample,value = exp,-rowname) %>%
  group_by(rowname) %>%
  mutate(norm = (exp - mean(exp))/sd(exp)) %>%
  dplyr::select(-exp) %>%
  spread(key = sample,value = norm) %>%
  column_to_rownames()

sample_ordered <- str_extract(colnames(mat_norm),pattern = "BASELINE|Fe|myelin")
sample_ordered

# build the annotation object
column_ha <- HeatmapAnnotation(treat = sample_ordered,  
                               col = list(treat = c("BASELINE" = "blue", "Fe"="orange","myelin"= "cyan"))) 

# highlight some terms
# rowname_color <- case_when(rownames(mat_norm)%in%c("KEGG_PRIMARY_IMMUNODEFICIENCY",
#                                                    "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
#                                                    "KEGG_VEGF_SIGNALING_PATHWAY",
#                                                    "KEGG_CELL_CYCLE",
#                                                    "KEGG_CELL_ADHESION_MOLECULES_CAMS",
#                                                    "KEGG_ECM_RECEPTOR_INTERACTION")~"red",
#                            T~"black")

hm <- Heatmap(mat_norm,
              # add annotation for the columns
              # hide columns labels
              # show_column_names = F,
              # fix width of the lables
              top_annotation = column_ha,
              # row_names_gp = gpar(col = rowname_color),
              row_names_max_width = max_text_width(
                rownames(mat),
                gp = gpar(fontsize = 12)
              ))

pdf(file = "../../out/image/heatmap_GSVA_KEGG_GMPVsResearch_tailored_es_nonLog_ZScore.pdf", width = 16, height = 5)
draw(hm,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2, 2, 2, 100), "mm"))
dev.off()
