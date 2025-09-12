# AIM ---------------------------------------------------------------------
# plot LGALS3, LAMP1 and TSPO expression for Aletta's project

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

# read in the dataset -----------------------------------------------------
data.combined <- readRDS("../../out/object/sobj_processed_donor.rds")
DimPlot(data.combined,label = T,raster = T,group.by = "seurat_clusters")
DimPlot(data.combined,label = T,raster = T,group.by = "expertAnno.l1")
DimPlot(data.combined,label = F,raster = T,split.by = "treat_full",group.by = "expertAnno.l1")
DimPlot(data.combined,label = T,raster = T,group.by = "harmonized_donor2")

# str_subset(rownames(data.combined),pattern = "HIF")
# define the gene of interest GOI
# GOI <- c("Irf7","Ddx58")
GOI <- c("LGALS3","LAMP1","TSPO")

table(data.combined@meta.data$expertAnno.l1)

# generate the table for the plots ----------------------------------------
# treat the set as a module and score it
data.combined <- AddModuleScore(data.combined,
                                features = list(GOI),
                                name="signature_score")

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
df_tot <- reduce(list(meta,UMAP1_df,df_exp),left_join, by="barcodes")
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

# expression distribution -------------------------------------------------
# crop the 0 expressing cells
df_exp %>%
  ggplot(aes(x=exp))+geom_histogram()+facet_grid(~gene)+theme_bw()+scale_x_log10()+geom_vline(xintercept = 1,col="red",linetype="dotted")

# keep the 0 expressing cells
df_exp %>%
  ggplot(aes(x=exp))+geom_histogram()+facet_wrap(~gene)+theme_bw()+
  # scale_x_log10()+
  geom_vline(xintercept = 2.5,col="red",linetype="dotted")

# library(scales)
# show_col(c("#4662D7FF","#FABA39FF","#7A0403FF"))

# plotting expression -----------------------------------------------------
# plot the signature score on the UMAP
df_tot2 %>%
  arrange(signature_score1) %>%
  # mutate(gene = "Ptx3") %>%
  ggplot() + geom_point(aes(x = UMAP_1, y = UMAP_2,col = signature_score1),alpha = 0.5,size = 0.2) +
  ggrepel::geom_text_repel(data = df_tot2_avg,aes(x = UMAP_1, y = UMAP_2,label = expertAnno.l1)) +
  # facet_grid(~treat) + 
  theme_cowplot() + 
  theme(strip.background = element_blank()) +
  # scale_color_gradientn(colours = c("blue","gray","red"))
  scale_color_gradientn(colours = viridis::turbo(10))

# plot the density using Nebulosa
plot_density(data.combined,reduction = "umap",features = "signature_score1")+ggplot2::scale_color_viridis_c(option = "turbo")

# by counts
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp)) + geom_point(alpha = 0.5,size = 0.2) +
  facet_wrap(gene~treat_full) +
  theme_cowplot() +
  scale_color_gradient(low = "gray",high = "blue") +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
# ggsave("../../out/image/06_UMAPggplot_annotationConfident_GLP1R_count.pdf",width = 13,height = 12)

# by counts category
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp_cat)) + geom_point(alpha = 0.5,size = 0.2) +
  facet_wrap(~gene) +
  theme_cowplot() +
  scale_color_manual(values = c("gray","blue")) +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())

# do the same using Seurat
FeaturePlot(data.combined,features = GOI,split.by = "treat_full",raster = T,order = T)
# ggsave("../../out/image/06_UMAPSeurat_annotationConfident_GLP1R_count.pdf",width = 25,height = 3)

# do the same using Seurat
FeaturePlot(data.combined,features = GOI,raster = T,order = T)
ggsave("../../out/image/26_UMAPSeurat_annotationConfident_GOI01_Aletta_count.pdf",width = 12,height = 10)

# try plotting gene expression with Nebulosa
GOI_subset <- str_subset(string = GOI,pattern = "GPR32",negate = T)
p01 <- plot_density(data.combined,features = GOI_subset,reduction = "umap")
ggsave(plot = p01,filename = "../../out/image/26_UMAPSeuratNebulosa_annotationConfident_GOI01_Aletta_count.pdf",width = 15,height = 5)

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
df_avg %>%
  write_tsv("../../out/table/26_df_avg_cellTypeTreat_GOI01_Aletta_BS03.tsv")

# plot the average expresison by cell annotation
df_avg |>
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=expertAnno.l1,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(~gene,scales = "free")
ggsave("../../out/image/26_avgExp_annotationConfident_GOI01_Aletta.pdf",width = 15,height = 5)

# plot splitting by treat full
df_avg |>
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=treat_full,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_grid(gene~expertAnno.l1,scales = "free")
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
# ggsave("../../out/image/06_dotplot_annotationConfident_GLP1R_expressionAvg_treatFull.pdf",width = 9,height = 9)
ggsave("../../out/image/26_avgExpSplit_annotationConfident_GOI01_Aletta.pdf",width = 10,height = 10)

# try to keep the same scale
# calculate the median per annotation
df_avg_summary <- df_avg %>%
  group_by(expertAnno.l1) %>%
  summarise(med = median(avg_exp)) %>%
  ungroup() %>%
  mutate(expertAnno.l1 = fct_reorder(expertAnno.l1,desc(med)))

# df_avg %>%
#   mutate(expertAnno.l1 = factor(expertAnno.l1,levels = levels(df_avg_summary$expertAnno.l1))) %>%
#   # ggplot(aes(x=NMDA_time,y=count)) + 
#   ggplot(aes(x=treat,y=avg_exp))+
#   geom_boxplot(outlier.shape = NA)+
#   geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
#   geom_hline(data = df_avg_summary,aes(yintercept = med),col="red",linetype="dashed") +
#   # geom_col()+
#   # facet_wrap(~cell_type2,scales = "free")+
#   theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
#   theme(strip.background = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA))+
#   facet_wrap(~expertAnno.l1,nrow=1)
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
# ggsave("../../out/image/06_dotplot_annotationConfident_GLP1R_expressionAvg_treat_scale.pdf",width = 14,height = 3)

# group the expression by donor id
df_avg %>%
  mutate(expertAnno.l1 = factor(expertAnno.l1,levels = levels(df_avg_summary$expertAnno.l1))) %>%
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=expertAnno.l1,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_hline(data = df_avg_summary,aes(yintercept = med),col="red",linetype="dashed") +
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_grid(gene~donor,scales = "free")
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
# ggsave("../../out/image/06_dotplot_annotationConfident_GLP1R_expressionAvg_donor_scale.pdf",width = 9,height = 3)

# add also the doubtles as donor
# df_avg2 <- average_GOI$RNA %>%
#   data.frame() %>%
#   rownames_to_column("gene") %>%
#   # mutate(gene = GOI) |> 
#   pivot_longer(names_to = "group",values_to = "avg_exp",-gene) %>%
#   filter(!str_detect(group,pattern="unassigned")) |> 
#   mutate(treat_full = str_extract(group,pattern = c("BASELINE|CSF.ctrl.24h|CSF.MS.24h|CSF.MS.48h|cytokine|Fe|myelin|TBHP"))) |> 
#   mutate(donor = str_extract(group,pattern = c("don\\w+|doublet"))) |> 
#   mutate(expertAnno.l1 = str_remove_all(group,pattern = c("BASELINE.|CSF.ctrl.24h.|CSF.MS.24h.|CSF.MS.48h.|cytokine.|Fe.|myelin.|TBHP.|.don\\w+|.doublet"))) |> 
#   # add a more general treat annotation
#   mutate(treat = case_when(treat_full%in%c("CSF.MS.24h","CSF.MS.48h")~"CSF.MS",
#                            treat_full%in%c("CSF.ctrl.24h")~"CSF.ctrl",
#                            T~treat_full))
# 
# # save the table
# df_avg2 %>%
#   write_tsv("../../out/table/06_df_avg_cellTypeTreat_C1Q_C3_BS03.tsv")
#
# # plot splitting by treat full
# df_avg2 |>
#   mutate(donor_id = case_when(donor =="doublet"~"doublet",
#                               T~"donor")) |> 
#   # ggplot(aes(x=NMDA_time,y=count)) + 
#   ggplot(aes(x=treat_full,y=avg_exp))+
#   geom_boxplot(outlier.shape = NA)+
#   geom_point(aes(color=donor_id),position = position_jitter(width = 0.1),alpha = 0.6)+
#   # geom_col()+
#   # facet_wrap(~cell_type2,scales = "free")+
#   theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 45))+
#   theme(strip.background = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA))+scale_color_manual(values = c("black","red"))+
#   # theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
#   # theme(strip.background = element_blank(),
#   #       panel.border = element_rect(colour = "black", fill = NA))+
#   facet_grid(gene~expertAnno.l1,scales = "free")
# # scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
# # ggsave("../../out/image/06_dotplot_annotationConfident_GLP1R_expressionAvg_treatFull2.pdf",width = 9,height = 9)

# -------------------------------------------------------------------------
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

# check if GLP1R is in there
# the gene is probably not passign eh threshold of expression at the single cell level
df_test_subset |> 
  filter(gene %in% GOI) |> 
  ggplot(aes(x = avg_log2FC,y = -log(p_val_adj),col=cluster,label=cluster))+geom_point()+
  facet_grid(gene~cell_id,scales = "free")+theme_bw()+theme(strip.background = element_blank())+
  geom_vline(xintercept = c(0),linetype="dashed",col="gray")+
  geom_vline(xintercept = c(-0.5,0.5),linetype="dashed",col="red")+
  geom_hline(yintercept = c(-log(0.05)),linetype="dashed",col="red")+
  ggrepel::geom_text_repel(force = 100,segment.alpha=0.5,max.overlaps = 5)
# ggsave("../../out/image/06_DE_annotationConfident_GLP1R_sc.pdf",width = 12,height = 9)

# run DE using BASELINE as reference
# cell_id <- "MG"
# list_test02 <- lapply(names(table(data.combined$expertAnno.l1)),function(cell_id){
#   # check the progress
#   print(cell_id)
#   test_obj <- subset(data.combined,subset = expertAnno.l1 == cell_id)
#   Idents(test_obj) <- "treat_full"
#   
#   # return all the genes to show the realtive position of TSPO treatment
#   test_cond <- lapply(c("myelin","CSF.MS.24h","CSF.ctrl.24h","cytokine","CSF.MS.48h","TBHP","Fe"), function(cond){
#     print(cond)
#     df_cond <- RunPresto(test_obj, ident.1 = cond,ident.2 = "BASELINE", min.pct = 0.05, logfc.threshold = 0,return.thresh=1)
#     
#     out <- df_cond %>%
#       mutate(test = paste0(cond,"_vs_BASELINE")) %>%
#       rownames_to_column("gene")
#     return(out)
#   })
#   
#   test_cond %>%
#     bind_rows() %>%
#     mutate(cell_id = cell_id)
# })
# 
# # save a table
# df_test02_subset <- list_test02 %>%
#   bind_rows()
# 
# df_test02_subset |>
#   write_tsv("../../out/table/20_FindMarkers_cellID_treat_minpct5_logfcthr01_refBASELINE.tsv")

df_test02_subset <- read_tsv("../../out/table/20_FindMarkers_cellID_treat_minpct5_logfcthr01_refBASELINE.tsv")

# check the stats for the GOI
df_sig <- df_test02_subset %>%
  filter(gene %in% GOI) %>%
  mutate(col= case_when(p_val_adj<0.05 & abs(avg_log2FC)>0.5~"sig",T~"non-sig"))

# df_test02_subset %>%
df_sig %>%
  ggplot(aes(x = avg_log2FC,y = -log(p_val_adj))) +
  # geom_point(col="gray50",alpha=0.2,shape=1) +
  geom_point(aes(col = col)) +
  ggrepel::geom_text_repel(aes(label = gene),col="black") +
  facet_grid(cell_id~test) +
  theme_bw() +
  theme(strip.background = element_blank())+
  geom_vline(xintercept = c(0),col="gray",linetype = "dashed")+
  geom_hline(yintercept = -log(0.05),col="gray",linetype = "dashed")
# ggsave("../../out/image/26_FindMarkers_cellID_treat_minpct5_logfcthr01_refBASELINE.pdf",width = 17,height = 15)
