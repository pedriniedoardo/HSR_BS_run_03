# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
# library(SeuratData)
# library(ggridges)
# library(ComplexHeatmap)

# read in the data --------------------------------------------------------
list_out1 <- readRDS("../../out/object/list_out_label_transfer_BSrun0102_SoupX_00200_06000_15.rds") %>%
  bind_rows(.id = "sample") %>% 
  mutate(filter_param = "00200_06000_15")
list_out2 <- readRDS("../../out/object/list_out_label_transfer_BSrun0102_SoupX_00600_06000_15.rds") %>% 
  bind_rows(.id = "sample") %>% 
  mutate(filter_param = "00600_06000_15")
list_out3 <- readRDS("../../out/object/list_out_label_transfer_BSrun0102_SoupX_01000_06000_15.rds") %>% 
  bind_rows(.id = "sample") %>% 
  mutate(filter_param = "01000_06000_15")

# read in the reference metadata for the run 01 and run 02
df_ref <- read_tsv(file = "../../out/table/ref_BS_meta.tsv")

# wrangling ---------------------------------------------------------------
# compare all the filtering options
df_all <- list(list_out1,
               list_out2,
               list_out3) %>%
  bind_rows()

# plot all the samples over the ref plot umap
plot <- ggplot(label= TRUE) +
  # reference layer
  geom_point(data = df_ref %>% dplyr::select(UMAP_1,UMAP_2),aes(x=UMAP_1,y=UMAP_2),col="gray",alpha=0.3,size=0.1)+
  # transfer layer
  geom_point(data = df_all,aes(x = refUMAP_1,y = refUMAP_2, col = predicted.id),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  theme_bw() +
  facet_grid(filter_param~ID)+
  guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))+theme(strip.background = element_blank())
ggsave(plot = plot,"../../out/image/UMAP_label_transfer_BS_panel.pdf",width = 40,height = 9)

# plot the ralative proportions across filtering stages
df_max <- df_all %>%
  group_by(ID,filter_param) %>%
  summarise(tot = n()) %>%
  ungroup() %>%
  filter(filter_param == "00200_06000_15") %>%
  dplyr::select(ID,tot)

df_summary <- df_all %>%
  group_by(ID,filter_param,predicted.id) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  left_join(df_max,by = "ID") %>%
  mutate(prop = n/tot) %>%
  group_by(predicted.id) %>%
  mutate(avg_prop = mean(prop)) %>%
  ungroup()

# check the proportions
df_summary %>%
  mutate(predicted.id = fct_reorder(predicted.id, avg_prop,.desc = T)) %>%
  ggplot(aes(x=predicted.id,y=prop,fill=filter_param))+
  geom_col(position = "dodge")+facet_wrap(~ID,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank())
ggsave("../../out/image/barplot_label_transfer_BS_panel.pdf",width = 15,height = 10)

# check if it is losing or gaining clusters in the comparison
df_summary2 <- df_summary %>% 
  # mutate(f = paste0(ID,"_",predicted.id)) %>% 
  split(f = paste0(.$ID,"_",.$predicted.id)) %>% 
  lapply(function(x){
    x %>%
      arrange(filter_param) %>% 
      mutate(FC_prop = prop/prop[1]) %>% 
      mutate(logFC_prop = log(FC_prop))
  }) %>% 
  bind_rows() %>% 
  ungroup() %>% 
  group_by(predicted.id) %>% 
  mutate(avg_logFC_prop = mean(logFC_prop)) %>% 
  ungroup()

# plot the general FC per cluster
df_summary2 %>%
  mutate(predicted.id = fct_reorder(predicted.id, avg_logFC_prop,.desc = T)) %>%
  ggplot(aes(x=predicted.id,y=logFC_prop,fill=filter_param))+
  geom_col(position = "dodge")+facet_wrap(~ID,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank())
ggsave("../../out/image/barplot_label_transfer_BS_panel2.pdf",width = 15,height = 10)

# plot a summary trend
df_summary2 %>%
  mutate(predicted.id = fct_reorder(predicted.id, avg_logFC_prop,.desc = T)) %>%
  ggplot(aes(x=predicted.id,y=logFC_prop,col=filter_param))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_dodge(width = 0.7))+
  theme_bw() +
  theme(strip.background = element_blank())
ggsave("../../out/image/barplot_label_transfer_BS_panel3.pdf",width = 6,height = 4)
