# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(ComplexHeatmap)
library(UpSetR)
library(gplots)

# plot genes dotplot ------------------------------------------------------
# save the filtered object
# read_tsv("../../out/table/DE_treatvsCSFctrl_pseudobulk_MG_shr.tsv") %>%
#   filter(conditionVsCSFctrl == "CSF.MS.24h_shr") %>%
#   filter(symbol %in% c("CSF1R","CTSA","CTSB","TREM2"))

data <- readRDS("../../out/object/ddsHTSeq_pseudobulk_MG_refCSFctrl.rds")

lut <- colData(data) %>%
  data.frame()

# GOI <- c(sig_B$symbol,sig_D$symbol)
GOI <- c("CSF1R","CTSA","CTSB","TREM2")
# GOI <- subset_genes

MR <- counts(data,normalized=T)%>%
  data.frame()%>%
  rownames_to_column("symbol") %>%
  pivot_longer(names_to = "sample",values_to = "exp",-symbol)%>%
  group_by(sample)%>%
  summarise(MR = sum(exp)/10^6)

# plot the data following the methods implemented in the plotCounts funciton from DESeq2
# Normalized counts plus a pseudocount of 0.5 are shown by default.
counts(data,normalized=T)%>%
  data.frame()%>%
  rownames_to_column("symbol") %>%
  filter(symbol %in% GOI) %>%
  pivot_longer(names_to = "sample",values_to = "count",-symbol) %>%
  # add the milion reads per sample
  left_join(MR,by = "sample") %>%
  left_join(lut,by = c("sample" = "samples")) %>%
  mutate(count_norm_adj = count + 0.5)%>%
  ggplot(aes(x=treat_full,y = count_norm_adj))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha=0.6)+facet_wrap(~symbol,scales = "free")+scale_y_log10()+ theme_bw()+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("../../out/image/boxplot_GOI.pdf",width = 6,height = 6)


counts(data,normalized=T)%>%
  data.frame()%>%
  rownames_to_column("symbol") %>%
  filter(symbol %in% GOI) %>%
  pivot_longer(names_to = "sample",values_to = "count",-symbol) %>%
  # add the milion reads per sample
  left_join(MR,by = "sample") %>%
  left_join(lut,by = "sample") %>%
  mutate(count_norm_adj = count + 0.5)%>%
  ggplot(aes(x=treat,y = count_norm_adj,label=clone))+
  geom_boxplot(outlier.shape = NA,alpha=0.5)+
  geom_point(position = position_jitter(width = 0.1),alpha=0.6)+facet_wrap(~symbol,scales = "free")+scale_y_log10()+ theme_bw()+
  geom_text_repel()+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
ggsave("../../out/image/boxplot_GOI_label.pdf",width = 6,height = 6)

# -------------------------------------------------------------------------
# are they significant in the comparison?
file <- dir("../../out/table/") %>%
  str_subset(pattern = "res_") %>%
  str_subset(pattern = ".txt") %>%
  str_subset(pattern = "shr",negate = F)
file 

# load the results 
results <- lapply(paste0("../../out/table/",file),function(x){
  read_tsv(x) 
}) %>%
  setNames(str_remove_all(file,pattern = ".txt")) %>% 
  bind_rows(.id = "comparison")

results %>% 
  filter(symbol %in% GOI)
