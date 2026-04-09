# AIM ---------------------------------------------------------------------
# explore the table of genes queried by Giuliana for her proposal

# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(scales)
library(ComplexHeatmap)
library(lemon)
library(finalfit)
library(cowplot)
library(patchwork)
library(Nebulosa)
library(pals)
library(colorRamp2)
library(RColorBrewer)

# read in the dataset -----------------------------------------------------
data.combined <- readRDS("../../out/object/sobj_processed_donor.rds")
DimPlot(data.combined,label = T,raster = T,group.by = "seurat_clusters")
DimPlot(data.combined,label = T,raster = T,group.by = "expertAnno.l1")
DimPlot(data.combined,label = F,raster = T,split.by = "treat_full",group.by = "expertAnno.l1")
DimPlot(data.combined,label = T,raster = T,group.by = "harmonized_donor2")

# str_subset(rownames(data.combined),pattern = "HIF")
# define the gene of interest GOI
# GOI <- c("Irf7","Ddx58")
# GOI <- c("FPR2", "CMKLR1","GPR32","GPR18","GPR37","LGR6","LTB4R","RORA")
df_GOI <- read_csv("../../data/GOI_Giuliana.csv")
GOI <- df_GOI$gene_symbol %>% unique()

# divide the gene into modules
list_GOI <- df_GOI %>%
  split(f = .$class) %>%
  lapply(function(x){
    x %>%
      pull(gene_symbol)
  })

table(data.combined@meta.data$expertAnno.l1)

# generate the table for the plots ----------------------------------------
# treat the set as a module and score it
data.combined <- AddModuleScore(data.combined,
                                features = list_GOI,
                                name="sig_score")

# rename the modules for simplicity
df_rename <- data.frame(name_old = paste0("sig_score",1:length(list_GOI)),
           name_new = paste0("sig_",names(list_GOI) %>% str_replace_all(pattern = "\\s",replacement = "\\.")))

lookup_long <- df_rename$name_old
names(lookup_long) <- df_rename$name_new
# lookup_long
data.combined@meta.data <- dplyr::rename(data.combined@meta.data,all_of(lookup_long))

# get the metadata from the other object
meta <- data.combined@meta.data %>%
  rownames_to_column(var = "barcodes")

# extrac the expression value
df_exp <- FetchData(data.combined, vars = GOI,slot = "data") |> 
  rownames_to_column("barcodes") |> 
  pivot_longer(names_to = "gene",values_to = "exp",-barcodes) |> 
  # try to min/max normalize the count varaible per gene in order to rescale the difference in terms of expression
  group_by(gene) %>%
  # threshold of positiveness is based on the distriubtion of the expression of the signal in tihs case
  mutate(exp_min_max = ((exp - min(exp))/(max(exp)-min(exp))),
         exp_cat = case_when(exp > 0~"pos",
                             T~"neg")) %>%
  ungroup() %>%
  mutate(exp_fix = exp + rnorm(nrow(.))/100000)

# get the coordinates
UMAP1_df <- data.combined@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column(var = "barcodes")

# generate the dataset for mapping the data in the umamp
dim(UMAP1_df)
dim(df_exp)
dim(meta)

# generate a full table with the gene expression values
df_tot <- purrr::reduce(list(meta,UMAP1_df,df_exp),left_join, by="barcodes")
df_tot_avg <- df_tot %>% group_by(expertAnno.l1) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# generate a smaller table for the plotting of the signature score
df_tot2 <- left_join(UMAP1_df,meta,by = "barcodes")
df_tot2_avg <- df_tot2 %>% group_by(expertAnno.l1) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)

dim(df_tot)
dim(df_tot2)

head(data.combined@meta.data)
table(data.combined$ID)

# plot the average expression per sample. Use the variable cell type per donor as grouping
# data.combined$group <- paste0(data.combined$orig.ident,".",data.combined$cell_type2)
data.combined$group <- paste0(data.combined$treat_full,"-",data.combined$expertAnno.l1,"-",data.combined$harmonized_donor2)
# data.combined$group2 <- paste0(data.combined$orig.ident,".",data.combined$treat,".",data.combined$cell_type2)
Idents(data.combined) <- "group"
DefaultAssay(data.combined) <- "RNA"

average_GOI <- AverageExpression(data.combined,features = GOI,group.by = c("group"))

# plot the average expression per condition Use the variable cell type per donor as grouping
# data.combined$group <- paste0(data.combined$orig.ident,".",data.combined$cell_type2)
data.combined$group2 <- paste0(data.combined$treat_full,"-",data.combined$expertAnno.l1)
# data.combined$group2 <- paste0(data.combined$orig.ident,".",data.combined$treat,".",data.combined$cell_type2)
Idents(data.combined) <- "group2"
DefaultAssay(data.combined) <- "RNA"

average_GOI2 <- AverageExpression(data.combined,features = GOI,group.by = c("group2"))

# wrangling ---------------------------------------------------------------

# try to depict the average expression
df_avg <- average_GOI$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  # mutate(gene = GOI) |> 
  pivot_longer(names_to = "group",values_to = "avg_exp",-gene) %>%
  filter(!str_detect(group,pattern="doublet|unassigned")) |> 
  mutate(treat_full = str_extract(group,pattern = c("BASELINE|CSF.ctrl.24h|CSF.MS.24h|CSF.MS.48h|cytokine|Fe|myelin|TBHP"))) |> 
  mutate(donor = str_extract(group,pattern = c("don\\w+"))) |> 
  mutate(expertAnno.l1 = str_remove_all(group,pattern = c("BASELINE.|CSF.ctrl.24h.|CSF.MS.24h.|CSF.MS.48h.|cytokine.|Fe.|myelin.|TBHP.|.don\\w+"))) |> 
  # add a more general treat annotation
  mutate(treat = case_when(treat_full%in%c("CSF.MS.24h","CSF.MS.48h")~"CSF.MS",
                           treat_full%in%c("CSF.ctrl.24h")~"CSF.ctrl",
                           T~treat_full))

# save the table
# df_avg %>%
#   write_tsv("../../out/table/20_df_avg_cellTypeTreatDonor_Giuliana_BS03.tsv")

df_avg2 <- average_GOI2$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  # mutate(gene = GOI) |> 
  pivot_longer(names_to = "group",values_to = "avg_exp",-gene) %>%
  filter(!str_detect(group,pattern="doublet|unassigned")) |> 
  mutate(treat_full = str_extract(group,pattern = c("BASELINE|CSF.ctrl.24h|CSF.MS.24h|CSF.MS.48h|cytokine|Fe|myelin|TBHP"))) |> 
  mutate(expertAnno.l1 = str_remove_all(group,pattern = c("BASELINE.|CSF.ctrl.24h.|CSF.MS.24h.|CSF.MS.48h.|cytokine.|Fe.|myelin.|TBHP.|.don\\w+"))) |> 
  # add a more general treat annotation
  mutate(treat = case_when(treat_full%in%c("CSF.MS.24h","CSF.MS.48h")~"CSF.MS",
                           treat_full%in%c("CSF.ctrl.24h")~"CSF.ctrl",
                           T~treat_full))

# save the table
# df_avg2 %>%
#   write_tsv("../../out/table/20_df_avg_cellTypeTreat_Giuliana_BS03.tsv")


# DGE ---------------------------------------------------------------------
# check is it is a marker in the cell type markes

Idents(data.combined) <- "expertAnno.l1"
# subset the object to include only the one of interest
# cell_id <- "MG"

# list_test <- lapply(names(table(data.combined$expertAnno.l1)),function(cell_id){
#   # check the progress
#   print(cell_id)
#   test_obj <- subset(data.combined,subset = expertAnno.l1 == cell_id)
#   Idents(test_obj) <- "treat_full"
#   # return all the genes to show the realtive position of TSPO treatment
#   test <- RunPrestoAll(test_obj, min.pct = 0.05, logfc.threshold = 0,return.thresh=1)
#   
#   test |> 
#     mutate(cell_id = cell_id)
# })
# 
# # save a table
# df_test_subset <- list_test |> 
#   bind_rows()

# df_test_subset |> 
#   write_tsv("../../out/table/06_FindAllMarkers_cellID_treat_minpct5_logfcthr01.tsv")

df_test_subset <- read_tsv("../../out/table/06_FindAllMarkers_cellID_treat_minpct5_logfcthr01.tsv")

df_test_subset %>%
  filter(gene %in% GOI) %>%
  mutate(col= case_when(p_val_adj<0.05 & abs(avg_log2FC)>0.5~"sig",T~"non-sig")) %>%
  filter(col == "sig")


df_test02_subset <- read_tsv("../../out/table/20_FindMarkers_cellID_treat_minpct5_logfcthr01_refBASELINE.tsv")

# check the stats for the GOI
df_sig <- df_test02_subset %>%
  filter(gene %in% GOI) %>%
  mutate(col= case_when(p_val_adj<0.05 & abs(avg_log2FC)>0.5~"sig",T~"non-sig"))

GOI_test <- df_sig %>%
  filter(col == "sig") %>%
  # filter(cell_id == "ASTRO") %>%
  filter(str_detect(test,pattern = "CSF.ctrl.24h|CSF.MS.24h|CSF.MS.48h")) %>%
  pull(gene) %>%
  unique()

# df_test02_subset %>%
df_sig %>%
  filter(str_detect(test,pattern = "CSF.ctrl.24h|CSF.MS.24h|CSF.MS.48h")) %>%
  # filter(cell_id == c("ASTRO","MG","NEU","OLIGO","OPC")) %>%
  filter(cell_id == c("ASTRO")) %>%
  filter(col=="sig") %>%
  pull(gene) %>%
  unique()


df_plot_volcano <- df_sig %>%
  filter(str_detect(test,pattern = "CSF.ctrl.24h|CSF.MS.24h|CSF.MS.48h")) %>%
  filter(cell_id %in% c("ASTRO","MG","NEU","OLIGO","OPC"))
  # filter(cell_id == c("ASTRO"))

df_plot_volcano$test %>% table()
df_plot_volcano$cell_id %>% table()

df_plot_volcano %>%
  ggplot(aes(x = avg_log2FC,y = -log(p_val_adj))) +
  # geom_point(col="gray50",alpha=0.2,shape=1) +
  geom_point(aes(col = col)) +
  ggrepel::geom_text_repel(data = df_plot_volcano %>% filter(col=="sig"),aes(label = gene),col="black") +
  facet_grid(cell_id~test) +
  theme_bw() +
  theme(strip.background = element_blank())+
  scale_color_manual(values = c("gray","red")) +
  geom_vline(xintercept = c(0),col="gray",linetype = "dashed")+
  geom_hline(yintercept = -log(0.05),col="gray",linetype = "dashed")
ggsave("../../out/image/32_FindMarkers_cellID_treat_minpct5_logfcthr02_refBASELINE.pdf",width = 17,height = 15)

# heatmap -----------------------------------------------------------------
# define the expression value to use
df_use <- df_avg2

list_hm <- lapply(c("ASTRO","MG","NEU","OLIGO","OPC"), function(cell_type){
  
  # use the cell specific GOI
  GOI_test <- df_sig %>%
    filter(col == "sig") %>%
    filter(cell_id == cell_type) %>%
    filter(str_detect(test,pattern = "CSF.ctrl.24h|CSF.MS.24h|CSF.MS.48h")) %>%
    pull(gene) %>%
    unique()
  
  # plot a sample of genes
  df_exp <- df_use %>%
    filter(expertAnno.l1 == cell_type) %>%
    filter(treat_full %in% c("BASELINE","CSF.ctrl.24h","CSF.MS.24h","CSF.MS.48h")) %>%
    # filter(gene %in% c("GAPDH","VDAC1")) %>%
    filter(gene %in% GOI_test) %>%
    # make the gene scaling
    group_by(gene) %>%
    mutate(avg_exp_scaled = scale(avg_exp)[,1]) %>%
    ungroup()
  
  # confirm the scaling
  df_exp %>%
    group_by(gene) %>%
    summarise(avg = mean(avg_exp_scaled),
              sd = sd(avg_exp_scaled))
  
  mat_exp <- df_exp %>%
    select(group,gene,avg_exp_scaled) %>%
    pivot_wider(names_from = group,values_from = avg_exp_scaled) %>%
    column_to_rownames("gene")
  
  
  # annotations -------------------------------------------------------------
  # build the column annotation object
  LUT_sample <- df_use %>%
    filter(expertAnno.l1 == cell_type) %>%
    filter(treat_full %in% c("BASELINE","CSF.ctrl.24h","CSF.MS.24h","CSF.MS.48h")) %>%
    group_by(group,treat_full,expertAnno.l1) %>%
    summarise(.groups = "drop")
  
  df_sample_ordered <- LUT_sample %>% 
    dplyr::slice(match(colnames(mat_exp),LUT_sample$group))
  
  # make color for the cell type
  color_id2 <- alphabet(length(unique(df_use$expertAnno.l1)))
  # check the colors
  # show_col(color_id2)
  # build the named vector
  names(color_id2) <- unique(df_use$expertAnno.l1)
  
  column_ha <- HeatmapAnnotation(treat = df_sample_ordered$treat_full,
                                 cell_type = df_sample_ordered$expertAnno.l1,
                                 col = list(treat = c("BASELINE" = "blue",
                                                      "CSF.ctrl.24h" = "gray",
                                                      "CSF.MS.24h" = "gray30",
                                                      "CSF.MS.48h" = "black"),
                                            cell_type = color_id2)) 
  
  
  # build the row annotation object
  LUT_gene <- df_GOI %>%
    group_by(class,gene_symbol) %>%
    summarise(.groups = "drop")
  
  df_gene_ordered <- LUT_gene %>% 
    dplyr::slice(match(rownames(mat_exp),LUT_gene$gene_symbol))
  
  # make color for the genes
  color_id3 <- tableau20(length(unique(df_GOI$class)))
  # check the colors
  # show_col(color_id3)
  # build the named vector
  names(color_id3) <- unique(df_GOI$class)
  
  row_ha <- rowAnnotation(class = df_gene_ordered$class, 
                          col = list(class = color_id3)) 
  
  # show_col(viridis::viridis(option = "turbo",n = 10))
  # col_fun <- colorRamp2(breaks = seq(from = -2,to=2,length.out=20),
  #                       colors = viridis::viridis(option = "turbo",n = 20))
  
  # show_col(rev(brewer.pal(n = 11, name = "Spectral")))
  # col_fun <- colorRamp2(breaks = seq(from = -2,to=2,length.out=11),
  #                       colors = rev(brewer.pal(n = 11, name = "Spectral")))
  
  col_fun <- colorRamp2(breaks = c(-1.5,0,1.5),
                        colors = c("blue","white","red"))
  
  
  hm <- Heatmap(mat_exp, 
                name = "exp", 
                column_title = cell_type, 
                top_annotation = column_ha, 
                cluster_rows = T,
                cluster_columns = F,
                right_annotation = row_ha,
                col = col_fun
                # col = viridis::viridis(option = "magma",n = 20)
                
                # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                )
  return(hm)
})

# list_hm
# 
# list_hm_full <- Reduce("+", list_hm)
# 
# draw(list_hm_full, 
#      ht_gap = unit(5, "mm"),             # Space between heatmaps
#      main_heatmap = 1,                   # Use the first HM as the reference for legends
#      merge_legend = TRUE,                # Combine identical legends (like 'exp' and 'treat')
#      heatmap_legend_side = "right", 
#      annotation_legend_side = "right")

list_hm2 <- lapply(list_hm, function(x){
  # 2. Capture them as "grob" objects (graphical objects)
  hm <- grid.grabExpr(draw(x,heatmap_legend_side = "left",annotation_legend_side = "left",merge_legend = TRUE) )
  return(hm)
})

pdf("../../out/image/32_heatmap_test_reinhold.pdf",width = 15,height = 15) 
wrap_plots(list_hm2)
dev.off()
