# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(scales)
library(ggrepel)
library(cowplot)

# read the data -----------------------------------------------------------
# read in the dataset
data.combined <- readRDS("../../out/object/sobj_processed_donor.rds")
DimPlot(data.combined,label = T,raster = T,group.by = "seurat_clusters")
DimPlot(data.combined,label = T,raster = T,group.by = "expertAnno.l1")
DimPlot(data.combined,label = T,raster = T,group.by = "treat_full")

# pull the annotation derived from the SIT algorithm (generated from 00_test_SIT.R)
df_SIT <- read_tsv("../../out/table/meta_harmonySkipIntegration_AllSoupX_01000_06000_15_SIT.tsv")

# merge the data
# add also the coordinates for the UMAP
meta_full <- data.combined@meta.data %>%
  rownames_to_column("rowname") %>%
  left_join(df_SIT %>%
              dplyr::select("rowname",
                            "KEGG_CELL_CYCLE","Cell_cycle_arrest_signatureK","REACTOME_CELL_CYCLE","Cell_cycle_arrest_signatureR","WP_CELL_CYCLE","Cell_cycle_arrest_signatureWP",
                            "Senescence_scoreK",
                            "Senescence_scoreWP",
                            "Senescence_scoreR",
                            "Senescence_scoreKQUANTILE",
                            "Senescence_scoreRQUANTILE",
                            "Senescence_scoreWPQUANTILE",
                            "SENEQUANTILE"),by = "rowname") %>%
  left_join(data.combined@reductions$umap@cell.embeddings %>%
              data.frame() %>%
              rownames_to_column("rowname"),by = "rowname")

# update the meta in the seurat object
data.combined@meta.data <- meta_full %>%
  column_to_rownames("rowname")

# filter the cells of interest --------------------------------------------
# filter only the MG cells, also focus on MG cells only
sobj_subset <- subset(data.combined,subset = expertAnno.l1 %in% c("MG") & treat_full != "TBHP")

# run the differential expression -----------------------------------------
Idents(sobj_subset) <- "SENEQUANTILE"
# subset the object to include only the one of interest
# 

# check the oveall count of cells per condition
sobj_subset@meta.data %>%
  group_by(SENEQUANTILE) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot)

# per treatment
sobj_subset@meta.data %>%
  group_by(treat_full,SENEQUANTILE) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot)

# pre treatment/donor
sobj_subset@meta.data %>%
  group_by(treat_full,SENEQUANTILE,harmonized_donor2) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(treat_full,harmonized_donor2) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  print(n=50)

# run DE ------------------------------------------------------------------
# split the data per cell id wise and run DE comparing each conditions vs BASELINE.
test <- RunPresto(sobj_subset, ident.1 = "YES",ident.2 = "NO",min.pct = 0.05, logfc.threshold = 0,return.thresh=1,only.pos = F) %>%
  rownames_to_column("gene") %>%
  mutate(ident.1 = "YES",
         ident.2 = "NO")

# confirm the trend of some top genes
VlnPlot(sobj_subset,features = test$gene[1:3],group.by = "SENEQUANTILE")

test |>
  write_tsv("../../out/table/DE_MG_SITpos_vs_SITneg_minpct5.tsv")

# plot volcano ------------------------------------------------------------
# plot the volcano for the cluster 12 DE
volcano_tot <- read_tsv("../../out/table/DE_MG_SITpos_vs_SITneg_minpct5.tsv") %>%
  mutate(DE_cat = case_when(avg_log2FC > 1 & p_val_adj < 0.01~"up",
                            avg_log2FC < (-1) & p_val_adj < 0.01~"down",
                            T~"no")) %>%
  mutate(test = paste0(ident.1,"_",ident.2)) %>%
  mutate(cell_id = "MG")
# filter(cluster %in% c("IMM"))

ggplot() +
  geom_point(data = volcano_tot%>%filter(DE_cat=="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=1,alpha=0.2)+theme_bw()+
  geom_point(data = volcano_tot%>%filter(DE_cat!="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=1,alpha=0.2,col="red")+
  theme_bw()+
  ggrepel::geom_text_repel(data = volcano_tot%>%filter(DE_cat!="no"),aes(x=avg_log2FC,y=-log10(p_val_adj),label=gene))+
  theme_bw()+
  facet_grid(test~cell_id)+
  theme(strip.background = element_blank())
ggsave("../../out/image/volcano_DE_MG_SITpos_vs_SITneg_minpct5.pdf",width = 10,height = 10)

ggplot() +
  geom_point(data = volcano_tot%>%filter(DE_cat=="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=1,alpha=0.2)+theme_bw()+
  geom_point(data = volcano_tot%>%filter(DE_cat!="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=1,alpha=0.2,col="red")+
  theme_bw()+
  # ggrepel::geom_text_repel(data = volcano_tot%>%filter(DE_cat!="no"),aes(x=avg_log2FC,y=-log10(p_val_adj),label=gene))+
  theme_bw()+
  facet_grid(test~cell_id)+
  theme(strip.background = element_blank())
ggsave("../../out/image/volcano_DE_MG_SITpos_vs_SITneg_minpct5_noLab.pdf",width = 10,height = 10)

# check if TSPO is in any comparision
# volcano_tot |>
#   filter(gene %in% c("TSPO")) |>
#   arrange(desc(avg_log2FC))

# run enrichR on each comparison ------------------------------------------
dbs <- listEnrichrDbs()

# filter fo the db of interest
# dbs %>%
#   filter(str_detect(libraryName,pattern = "Atlas"))
#
# dbs %>%
#   filter(str_detect(libraryName,pattern = "Reactome"))

dbs_db <- c("KEGG_2021_Human","MSigDB_Hallmark_2020","Reactome_2022","HDSigDB_Human_2021","GO_Biological_Process_2023")

# query -------------------------------------------------------------------
# seelct only the clusters with more than 10 genes as degs

# pull the gene names dividing the up regulated from the downregulated by cell_id and group
list_genes_UP <- volcano_tot %>%
  mutate(test2 = paste0(.$cell_id,"|",.$test,"|UP")) %>%
  filter(DE_cat=="up") %>%
  split(f = .$test2) %>%
  map(function(x){
    x %>%
      pull(gene)
  })
# pull the downregulated
list_genes_DOWN <- volcano_tot %>%
  mutate(test2 = paste0(.$cell_id,"|",.$test,"|DOWN")) %>%
  filter(DE_cat=="down") %>%
  split(f = .$test2) %>%
  map(function(x){
    x %>%
      pull(gene)
  })

# define the background
# background <- df_modules$feature

# x <- list_res_tot_UP_filter$`DeMye_vs_Norm|clust_5`
list_enrichr_UP <- lapply(list_genes_UP,function(x){
  genes <- x
  # out_enrich <- enrichr(genes, dbs_db,background = background,include_overlap=T)
  out_enrich <- enrichr(genes, dbs_db)
  
  # filter out the annotations without an output
  filter_out <- lapply(out_enrich,function(x){dim(x)[1]}) %>%
    unlist()
  
  out_enrich[filter_out>0] %>%
    bind_rows(.id = "annotation")
}) %>%
  bind_rows(.id = "comparison")

# run the same for the downregulated
list_enrichr_DOWN <- lapply(list_genes_DOWN,function(x){
  genes <- x
  # out_enrich <- enrichr(genes, dbs_db,background = background,include_overlap=T)
  out_enrich <- enrichr(genes, dbs_db)
  
  # filter out the annotations without an output
  filter_out <- lapply(out_enrich,function(x){dim(x)[1]}) %>%
    unlist()
  
  out_enrich[filter_out>0] %>%
    bind_rows(.id = "annotation")
}) %>%
  bind_rows(.id = "comparison")

# save the tables
list_enrichr_UP %>%
  write_tsv("../../out/table/enrichR_DE_MG_SITpos_vs_SITneg_minpct5_UP.tsv")

list_enrichr_DOWN %>%
  write_tsv("../../out/table/enrichR_DE_MG_SITpos_vs_SITneg_minpct5_DOWN.tsv")

# plot the results --------------------------------------------------------
plot_list_UP <- list_enrichr_UP %>%
  split(f = .$comparison)

# library(scales)
list_plot_UP <- pmap(list(plot_list_UP,names(plot_list_UP)), function(x,y){
  x %>%
    group_by(annotation) %>%
    arrange(P.value) %>%
    dplyr::slice(1:10) %>%
    mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
    mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
    # ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    scale_color_gradientn(colors = c("red","blue"),
                          values = rescale(c(0,1)),
                          limits = c(0,0.2))+
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))+
    ggtitle(y)
  # scale_color_gradient(low = "red",high = "blue")
  
  #ggsave(paste0("image/enrichR_out_",y,".pdf"),width = 7,height = 15)
})

wrap_plots(list_plot_UP)
ggsave("../../out/image/enrichR_DE_MG_SITpos_vs_SITneg_minpct5_UP.pdf",width = 6,height = 12,limitsize = FALSE)

# down genes
plot_list_DOWN <- list_enrichr_DOWN %>%
  split(f = .$comparison)

# library(scales)
list_plot_DOWN <- pmap(list(plot_list_DOWN,names(plot_list_DOWN)), function(x,y){
  x %>%
    group_by(annotation) %>%
    arrange(P.value) %>%
    dplyr::slice(1:10) %>%
    mutate(Term = str_sub(Term,start = 1,end = 30)) %>%
    mutate(Term = fct_reorder(Term, Combined.Score,.desc = F)) %>%
    # ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    ggplot(aes(y=Term,x=Combined.Score,size = Odds.Ratio,col = Adjusted.P.value)) + geom_point() + facet_wrap(~annotation,scales = "free",ncol = 1)+theme_bw() +
    scale_color_gradientn(colors = c("red","blue"),
                          values = rescale(c(0,1)),
                          limits = c(0,0.2))+
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))+
    ggtitle(y)
  # scale_color_gradient(low = "red",high = "blue")
  
  #ggsave(paste0("image/enrichR_out_",y,".pdf"),width = 7,height = 15)
})

wrap_plots(list_plot_DOWN)
ggsave("../../out/image/enrichR_DE_MG_SITpos_vs_SITneg_minpct5_DOWN.pdf",width = 6,height = 12,limitsize = FALSE)

# # confirm some of the trends
# list_enrichr_UP |>
#   filter(comparison == "MG|YES_NO|UP") |>
#   head()
# 
# # JUN seems to be an upregulated genes in MS.48h
# FeaturePlot(object = data.combined,features = "CDKN1A",split.by = "treat_full",order = T)
