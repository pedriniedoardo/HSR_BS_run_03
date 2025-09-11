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
rownames(data.combined) %>% str_subset(pattern = "^MT-")
GOI <- c("MT-ND1","MT-ND2","MT-CO1","MT-CO2","MT-ATP8","MT-ATP6","MT-CO3","MT-ND3","MT-ND4L","MT-ND4","MT-ND5","MT-ND6","MT-CYB")
# set_02 <- c("ABCC1", "AKR1C3", "ALOX12", "ALOX15", "ALOX5", "CBR1", "CYP1A2", "CYP2C8", "CYP2C9", "CYP2D6", "CYP2E1", "CYP3A4", "CYP4F2", "CYP4F3", "DPEP1", "DPEP2", "DPEP3", "EPHX1", "EPHX2", "EPHX3", "GGT1", "GGT2", "GGT5", "GPX2", "GPX4", "GSTM4", "HPGDS", "LTA4H", "LTC4S", "PTGDS", "PTGES", "PTGES2", "PTGIS", "PTGS1", "PTGS2", "TBXAS1", "TXN")

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
ggsave("../../out/image/20_UMAPSeurat_annotationConfident_set01_count.pdf",width = 12,height = 10)

Idents(data.combined) <- "expertAnno.l1"
test <- DotPlot(data.combined,features = GOI)

df_test <- lapply(list(MT_genes = GOI),function(x){
  test$data %>% 
    filter(features.plot %in% x)
}) %>% 
  bind_rows(.id = "cell_type")

head(df_test)

df_test %>%
  # force the order
  # mutate(id = factor(id,levels = c(3,11,0,6,2,12,13,9,4,5,7,8,1,10))) %>% 
  # mutate(cell_type = factor(cell_type,levels = c("IMMUNE","OLIGO","OPC","ASTRO","VAS","NEURONS","EPENDYMA"))) %>% 
  ggplot(aes(x = features.plot,y = id)) +
  geom_point(aes(size = pct.exp, col = avg.exp.scaled))+
  scale_radius(range = c(0, 8)) +
  facet_grid(~cell_type,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_color_gradient(low = "lightgrey",high = "blue")
ggsave("../../out/image/20_DotPlot_MT_genes.pdf",width = 6,height = 6)

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
df_avg

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
# ggsave("../../out/image/20_avgExp_annotationConfident_set01.pdf",width = 12,height = 10)

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
ggsave("../../out/image/20_avgExpSplit_annotationConfident_set01.pdf",width = 12,height = 10)
