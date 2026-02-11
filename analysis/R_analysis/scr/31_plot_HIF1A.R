# AIM ---------------------------------------------------------------------
# explore the expression of HIF1A in the whole dataset

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
GOI <- c("HIF1A")

table(data.combined@meta.data$expertAnno.l1)

# generate the table for the plots ----------------------------------------
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

df_tot <- reduce(list(meta,UMAP1_df,df_exp),left_join, by="barcodes")
df_tot_avg <- df_tot %>% group_by(expertAnno.l1) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)

dim(df_tot)

head(data.combined@meta.data)

# plot the average expression per sample use the variable cell tyep per donor as grouping
# data.combined$group <- paste0(data.combined$harmonized_donor2,".",data.combined$cell_type2)
data.combined$group <- paste0(data.combined$treat_full,"-",data.combined$expertAnno.l1,"-",data.combined$harmonized_donor2)
# data.combined$group2 <- paste0(data.combined$harmonized_donor2,".",data.combined$treat,".",data.combined$cell_type2)
Idents(data.combined) <- "group"
DefaultAssay(data.combined) <- "RNA"

average_GOI <- AverageExpression(data.combined,features = GOI,group.by = c("group"))

# plot general UMAP -------------------------------------------------------
# build the plot using both info
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=expertAnno.l1),size=0.3) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = df_tot_avg,aes(x = UMAP_1,y = UMAP_2,label = expertAnno.l1),col="black")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_cowplot()

# no lab
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = UMAP_1,y = UMAP_2,col=expertAnno.l1),size=0.3) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  # ggrepel::geom_text_repel(data = df_tot_avg,aes(x = UMAP_1,y = UMAP_2,label = cell_type2),col="black")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_cowplot()

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
ggsave("../../out/image/31_UMAPggplot_annotationConfident_HIF1A_count.pdf",width = 13,height = 12)

# do the same using Seurat
FeaturePlot(data.combined,features = GOI,split.by = "treat_full",raster = T,order = T,ncol = 3)
# ggsave("../../out/image/06_UMAPSeurat_annotationConfident_HIF1A_count.pdf",width = 25,height = 3)

df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp)) + geom_point(alpha = 0.5,size = 0.2) +
  facet_wrap(~gene) +
  theme_cowplot() +
  scale_color_gradient(low = "gray",high = "blue") +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
ggsave("../../out/image/31_UMAPggplot_annotationConfident_HIF1A_count2.pdf",width = 6,height = 5)

# do the same using Seurat
FeaturePlot(data.combined,features = GOI,raster = T,order = T)

# try nebulosa option
plot_density(data.combined, GOI,reduction = "umap")
ggsave("../../out/image/31_UMAPSeurat_annotationConfident_HIF1A_countNebulosa.pdf",width = 6,height = 5)

df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp)) + geom_point(alpha = 0.5,size = 0.05) +
  facet_wrap(gene~treat_full) +
  theme_cowplot() +
  # scale_color_gradient(low = "gray",high = "blue") +
  scale_color_viridis_c(option = "turbo") +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
# ggsave("../../out/image/31_UMAPggplot_annotationConfident_HIF1A_count_alt.pdf",width = 13,height = 12)

# by min max normalized counts
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp_min_max)) + geom_point(alpha = 0.5,size = 0.2) +
  facet_wrap(gene~treat_full) +
  theme_cowplot() +
  scale_color_gradient(low = "gray",high = "blue") +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
# ggsave("../../out/image/31_UMAPggplot_annotationConfident_HIF1A_minmax.pdf",width = 13,height = 12)

# plot the category. being 0 or non zero per cell
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24","NMDA_36")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp_cat)) + geom_point(alpha = 0.5,size = 0.05) +
  # facet_wrap(gene~NMDA_time,nrow = 2) +
  # facet_rep_wrap(gene~treat,repeat.tick.labels = "all",nrow=3)+
  facet_wrap(gene~treat_full)+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_cowplot() +
  scale_color_manual(values = c("gray","blue")) +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
# ggsave("../../out/image/31_UMAPggplot_annotationConfident_HIF1A_proppos.pdf",width = 13,height = 12)

df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24","NMDA_36")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2,col = exp_cat)) + geom_point(alpha = 0.5,size = 0.05) +
  # facet_wrap(gene~NMDA_time,nrow = 2) +
  # facet_rep_wrap(gene~treat,repeat.tick.labels = "all",nrow=3)+
  facet_wrap(~gene)+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_cowplot() +
  scale_color_manual(values = c("gray","blue")) +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())

# ggsave("../../out/image/31_UMAPggplot_annotationConfident_HIF1A_proppos2.pdf",width = 6,height = 5)

# violin plot for GOI expression use macro categories
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24","NMDA_36")) %>%
  # this is the processing shown in the violinplot function
  # mutate(exp_fix = exp + rnorm(nrow(.))/100000) %>%
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=treat_full,y=exp_fix)) + 
  geom_violin(scale = "width")+
  geom_point(position=position_jitter(width = 0.2),alpha=0.01) +
  facet_wrap(~expertAnno.l1) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1,angle = 90)) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
# ggsave("../../out/image/31_violin_annotationConfident_HIF1A.pdf",width = 15,height = 10)

# try to depict the average expression there is roughly one sample per condition
df_avg <- average_GOI$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(gene = GOI) |> 
  pivot_longer(names_to = "group",values_to = "avg_exp",-gene) %>%
  # filter(!str_detect(group,pattern="doublet|unassigned")) |> 
  mutate(treat_full = str_extract(group,pattern = c("myelin|CSF.MS.24h|CSF.ctrl.24h|cytokine|CSF.MS.48h|TBHP|Fe|BASELINE"))) |> 
  mutate(donor = str_extract(group,pattern = c("donRR16|doublet|donRR25|unassigned|donRR24" ))) |> 
  mutate(expertAnno.l1 = str_extract(group,pattern = c("GLIA_IMM|OLIGO|NEU|PROG|MG|ASTRO|OPC")))

# confirm all the value have a label
df_avg %>%
  filter(is.na(expertAnno.l1))

# plot the average expression by cell annotation
# calculate the average per cluster
order_id <- df_avg %>% group_by(expertAnno.l1) %>% summarise(avg = median(avg_exp)) %>% arrange(desc(avg)) %>% pull(expertAnno.l1)

# df_avg %>%
#   filter(expertAnno.l1 == "NEU") %>%
#   group_by(treat_full) %>%
#   summarise()
# 
# df_avg %>%
#   filter(expertAnno.l1 == "NEU") %>%
#   group_by(donor) %>%
#   summarise(n = n())
# 
# df_avg %>%
#   filter(expertAnno.l1 == "NEU",donor == "donRR24")

df_avg |>
  mutate(expertAnno.l1 = factor(expertAnno.l1,levels = order_id)) %>%
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
  facet_wrap(~gene,scales = "free")+
  scale_y_continuous(trans = "log1p")
ggsave("../../out/image/31_dotplot_annotationConfident_HIF1A_expressionAvg_celltype_general.pdf",width = 6,height = 4)

# do the same as above but split by condition
df_avg |>
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=treat_full,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(~gene,scales = "free") +
  scale_y_continuous(trans = "log1p")
ggsave("../../out/image/31_dotplot_annotationConfident_HIF1A_expressionAvg_treat_general.pdf",width = 6,height = 4)

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
  facet_wrap(~expertAnno.l1,scales = "free")
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
# ggsave("../../out/image/06_dotplot_annotationConfident_HIF1A_expressionAvg_treatFull.pdf",width = 9,height = 9)

# plot splitting by treat
df_avg |>
  mutate(expertAnno.l1 = factor(expertAnno.l1,levels = order_id)) %>%
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=treat_full,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(~expertAnno.l1,scales = "free")
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
ggsave("../../out/image/31_dotplot_annotationConfident_HIF1A_expressionAvg_treat.pdf",width = 9,height = 9)

# try to keep the same scale
# calculate the median per annotation
df_avg_summary <- df_avg %>%
  group_by(expertAnno.l1) %>%
  summarise(med = median(avg_exp)) %>%
  ungroup() %>%
  mutate(expertAnno.l1 = fct_reorder(expertAnno.l1,desc(med)))

df_avg %>%
  mutate(expertAnno.l1 = factor(expertAnno.l1,levels = levels(df_avg_summary$expertAnno.l1))) %>%
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=treat_full,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  geom_hline(data = df_avg_summary,aes(yintercept = med),col="red",linetype="dashed") +
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(~expertAnno.l1,nrow=1)
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
ggsave("../../out/image/31_dotplot_annotationConfident_HIF1A_expressionAvg_treat_scale.pdf",width = 14,height = 3)


# -------------------------------------------------------------------------
# relative to Elif poster, martina suggested to the following modifications
# for the general trend of everage expression, remove the THBP treatment
# aggregate the MS conditions without splitting by time points

# show the groups
data.combined@meta.data %>%
  select(treat,expertAnno.l1,harmonized_donor2) %>%
  mutate(treat = as.factor(treat),
         expertAnno.l1 = as.factor(expertAnno.l1),
         harmonized_donor2 = as.factor(harmonized_donor2)) %>%
  summary()

# filter out the samples from `doublet` and `unassigned`
data.combined_sub <- subset(data.combined,subset = harmonized_donor2 %in% c("donRR16","donRR24","donRR25") & treat != "TBHP")

# generate the covariate for the aggregation
data.combined_sub$group2 <- paste0(data.combined_sub$treat,"-",data.combined_sub$expertAnno.l1,"-",data.combined_sub$harmonized_donor2)

# confirm the update
data.combined_sub@meta.data %>%
  select(harmonized_donor2,treat) %>%
  mutate(harmonized_donor2 = as.factor(harmonized_donor2),
         treat = as.factor(treat),) %>%
  summary()

# data.combined$group2 <- paste0(data.combined$harmonized_donor2,".",data.combined$treat,".",data.combined$cell_type2)
Idents(data.combined_sub) <- "group"
DefaultAssay(data.combined_sub) <- "RNA"

average_GOI2 <- AverageExpression(data.combined_sub,features = GOI,group.by = c("group2"))

df_avg2 <- average_GOI2$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(gene = GOI) |> 
  pivot_longer(names_to = "group2",values_to = "avg_exp",-gene) %>%
  # filter(!str_detect(group,pattern="doublet|unassigned")) |> 
  mutate(treat_full = str_extract(group2,pattern = c("myelin|CSF.ctrl|CSF.MS|cytokine|Fe|BASELINE"))) |> 
  mutate(donor = str_extract(group2,pattern = c("donRR16|donRR25|donRR24" ))) |> 
  mutate(expertAnno.l1 = str_extract(group2,pattern = c("GLIA_IMM|OLIGO|NEU|PROG|MG|ASTRO|OPC"))) %>%
  mutate(avg_exp_log1p = log1p(avg_exp))

# make the general plot
df_avg2 |>
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=treat_full,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(~gene,scales = "free") +
  scale_y_continuous(trans = "log1p")
ggsave("../../out/image/31_dotplot_annotationConfident_HIF1A_expressionAvg_treat_general2.pdf",width = 6,height = 4)

df_avg2 |>
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
  facet_wrap(~gene,scales = "free") +
  scale_y_continuous(trans = "log1p")
ggsave("../../out/image/31_dotplot_annotationConfident_HIF1A_expressionAvg_celltype_general2.pdf",width = 6,height = 4)

df_avg2 |>
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=treat_full,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(expertAnno.l1~gene,scales = "free") +
  scale_y_continuous(trans = "log1p")
ggsave("../../out/image/31_dotplot_annotationConfident_HIF1A_expressionAvg_celltypeSplitTreat_general2.pdf",width = 9,height = 9)

# generate the univariate analysis for the expression by controlling for the cell type and donor, use as reference the BASELINE
res_lm01 <- lm(data = df_avg2,formula = avg_exp~donor + expertAnno.l1 + treat_full)
summary(res_lm01)

res_lm02 <- lm(data = df_avg2,formula = avg_exp_log1p~donor + expertAnno.l1 + treat_full)
summary(res_lm02)

# try with random variable for the donor
dependent <- "avg_exp"
# explanatory <- c("expertAnno.l1","treat_full")
explanatory <- c("treat_full","expertAnno.l1")

# both are the same the same as:
# random_effect <- "(1 | donor)"
random_effect <- "(1 | donor)"

df_avg2 %>% 
  finalfit(dependent,
           explanatory, 
           random_effect = random_effect,
           metrics = TRUE)

df_avg2 %>%
  finalfit::coefficient_plot(
    dependent = dependent,
    explanatory = explanatory,
    random_effect = random_effect
  )

# try with random variable for the donor with log expression data
dependent0 <- "avg_exp_log1p"
# explanatory <- c("expertAnno.l1","treat_full")
explanatory0 <- c("treat_full","expertAnno.l1")

# both are the same the same as:
# random_effect <- "(1 | donor)"
random_effect0 <- "(1 | donor)"

df_avg2 %>% 
  finalfit(dependent0,
           explanatory0, 
           random_effect = random_effect0,
           metrics = TRUE)

df_avg2 %>%
  finalfit::coefficient_plot(
    dependent = dependent0,
    explanatory = explanatory0,
    random_effect = random_effect
  )

# try with all fixed effects
dependent2 <- "avg_exp"
# explanatory <- c("expertAnno.l1","treat_full")
explanatory2 <- c("donor","treat_full","expertAnno.l1")

df_avg2 %>% 
  finalfit(dependent2,
           explanatory2, 
           metrics = TRUE)

df_avg2 %>%
  finalfit::coefficient_plot(
    dependent = dependent2,
    explanatory = explanatory2
  )

# try with all fixed effects and log expression
dependent3 <- "avg_exp_log1p"
# explanatory <- c("expertAnno.l1","treat_full")
explanatory3 <- c("donor","treat_full","expertAnno.l1")

df_avg2 %>% 
  finalfit(dependent3,
           explanatory3, 
           metrics = TRUE)

pdf("../../out/image/31_model_HIF1A.pdf",width = 10,height = 5)
df_avg2 %>%
  finalfit::coefficient_plot(
    dependent = dependent3,
    explanatory = explanatory3
  )
dev.off()

# check the correlation between genes -------------------------------------
# martina also asked to check if the expression of HIF1A and DNMT3A

# filter out the samples from `doublet` and `unassigned`
data.combined_sub <- subset(data.combined,subset = harmonized_donor2 %in% c("donRR16","donRR24","donRR25") & treat != "TBHP")

# generate the covariate for the aggregation
data.combined_sub$group2 <- paste0(data.combined_sub$treat,"-",data.combined_sub$expertAnno.l1,"-",data.combined_sub$harmonized_donor2)

# confirm the update
data.combined_sub@meta.data %>%
  select(harmonized_donor2,treat) %>%
  mutate(harmonized_donor2 = as.factor(harmonized_donor2),
         treat = as.factor(treat),) %>%
  summary()

# data.combined$group2 <- paste0(data.combined$harmonized_donor2,".",data.combined$treat,".",data.combined$cell_type2)
Idents(data.combined_sub) <- "group"
DefaultAssay(data.combined_sub) <- "RNA"

average_GOI3 <- AverageExpression(data.combined_sub,features = c("HIF1A","DNMT3A"),group.by = c("group2"))

df_avg3 <- average_GOI3$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(names_to = "group2",values_to = "avg_exp",-gene) %>%
  # filter(!str_detect(group,pattern="doublet|unassigned")) |> 
  mutate(treat_full = str_extract(group2,pattern = c("myelin|CSF.ctrl|CSF.MS|cytokine|Fe|BASELINE"))) |> 
  mutate(donor = str_extract(group2,pattern = c("donRR16|donRR25|donRR24" ))) |> 
  mutate(expertAnno.l1 = str_extract(group2,pattern = c("GLIA_IMM|OLIGO|NEU|PROG|MG|ASTRO|OPC"))) %>%
  mutate(avg_exp_log1p = log1p(avg_exp))

# check the correlation in the dataset
df_avg3 %>%
  select(-avg_exp) %>%
  pivot_wider(names_from = gene,values_from = avg_exp_log1p) %>%
  ggplot(aes(x = HIF1A,y = DNMT3A)) + geom_point() + theme_bw() + geom_smooth(method = "lm")
ggsave("../../out/image/31_correlation_HIF_DNMT_global.pdf",width = 4,height = 4)

# run the stat
test <- df_avg3 %>%
  select(-avg_exp) %>%
  pivot_wider(names_from = gene,values_from = avg_exp_log1p)

cor.test(test$HIF1A,test$DNMT3A) %>%
  broom::tidy()
  
# try to split by cell type
df_avg3 %>%
  select(-avg_exp) %>%
  pivot_wider(names_from = gene,values_from = avg_exp_log1p) %>%
  ggplot(aes(x = HIF1A,y = DNMT3A)) + geom_point() + theme_bw() + geom_smooth(method = "lm") + facet_wrap(~expertAnno.l1,scales="free") + theme(strip.background = element_blank())
ggsave("../../out/image/31_correlation_HIF_DNMT_SplitCellType.pdf",width = 9,height = 9)

# run the stat
list_test <- df_avg3 %>%
  select(-avg_exp) %>%
  pivot_wider(names_from = gene,values_from = avg_exp_log1p) %>%
  split(f = .$expertAnno.l1)

lapply(list_test,function(test){
  cor.test(test$HIF1A,test$DNMT3A) %>%
    broom::tidy()
}) %>%
  bind_rows(.id = "cell_type")

# try to split by treatment
df_avg3 %>%
  select(-avg_exp) %>%
  pivot_wider(names_from = gene,values_from = avg_exp_log1p) %>%
  ggplot(aes(x = HIF1A,y = DNMT3A)) + geom_point() + theme_bw() + geom_smooth(method = "lm") + facet_wrap(~treat_full,scales="free") + theme(strip.background = element_blank())
ggsave("../../out/image/31_correlation_HIF_DNMT_SplitTreat.pdf",width = 9,height = 6)

# run the stat
list_test2 <- df_avg3 %>%
  select(-avg_exp) %>%
  pivot_wider(names_from = gene,values_from = avg_exp_log1p) %>%
  split(f = .$treat_full)

lapply(list_test2,function(test){
  cor.test(test$HIF1A,test$DNMT3A) %>%
    broom::tidy()
}) %>%
  bind_rows(.id = "treat_full")

df_avg3 %>%
  select(-avg_exp) %>%
  pivot_wider(names_from = gene,values_from = avg_exp_log1p) %>%
  ggplot(aes(x = HIF1A,y = DNMT3A)) + geom_point(aes(col=expertAnno.l1)) + theme_bw() + geom_smooth(method = "lm") + facet_wrap(~treat_full,scales="free") + theme(strip.background = element_blank())
ggsave("../../out/image/31_correlation_HIF_DNMT_SplitTreat2.pdf",width = 10,height = 6)

df_avg3 %>%
  ggplot(aes(x = expertAnno.l1,y = avg_exp_log1p)) + 
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~gene,ncol=1,scales="free") +
  geom_point(position = position_jitter(width = 0.2),shape = 1) + theme_bw() + theme(strip.background = element_blank(), axis.text.x = element_text(angle = 45,hjust = 1))
ggsave("../../out/image/31_boxplot_HIF_DNMT_SplitGene.pdf",width = 4,height = 6)
