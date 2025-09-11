# # this script is run to explore the dataset individually and define a consensus threshold for the mito and reads counts
# # LIBRARIES ---------------------------------------------------------------
# library(scater)
# library(Seurat)
# library(tidyverse)
# library(robustbase)
# # library(SeuratData)
# library(patchwork)
# 
# # READ IN DATA ------------------------------------------------------------
# # file <- c("ECFC_P09")
# id_sample <- dir("../data/cellranger61/out/") %>%
#   str_subset(pattern = "^cr61",negate = T)
# 
# # load the LUT
# LUT <- read_csv("../data/LUT_samples.csv")
# # mutate(sample_id = str_pad(ID,width = 2,side = "left",pad = "0")) %>% 
# # mutate(sample_id = paste0(sample_id,"_cr_61"))
# 
# # do the preprocessing over all the dataset and save the objects
# # x <- "W8_18h_myelin_plus_untreated_multiplexed"
# list_datasc <- lapply(id_sample,function(x){
#   # to track the processing of the progress of the lapply
#   print(x)
#   data <- Read10X(data.dir = paste0("../data/cellranger61/out/",x,"/outs/filtered_feature_bc_matrix/"))
#   
#   datasc <- CreateSeuratObject(counts = data, project = LUT %>%
#                                  filter(sample == x) %>%
#                                  pull(sample), min.cells = 20, min.features = 200)
#   
#   # datasc <- CreateSeuratObject(counts = data, min.cells = 20, min.features = 200)
#   
#   # datasc@meta.data
#   datasc$percent.mt <- PercentageFeatureSet(datasc, pattern = "^MT-")
#   datasc$percent.ribo <- PercentageFeatureSet(datasc, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
#   
#   # label the cells based on the mt reads content
#   datasc$mt_bin <- datasc@meta.data %>%
#     mutate(test = case_when(percent.mt < 1~"low",
#                             percent.mt < 10~"mid",
#                             T ~ "high")) %>%
#     pull(test)
#   
#   datasc$treat <- LUT %>%
#     filter(sample == x) %>%
#     pull(treat)
#   
#   datasc$clone <- LUT %>%
#     filter(sample == x) %>%
#     pull(clone)
#   
#   datasc$doxy <- LUT %>%
#     filter(sample == x) %>%
#     pull(doxy)
#   
#   datasc$exposure <- LUT %>%
#     filter(sample == x) %>%
#     pull(exposure)
#   
#   # add the filtering variable based on the fixed threshold
#   datasc$test <- datasc@meta.data %>%
#     # mutate(test = percent.mt < 20 & nFeature_RNA > 700 & nFeature_RNA < 9000) %>%
#     # mutate(test = percent.mt > 1 & percent.mt < 10 & nFeature_RNA > 1000 & nFeature_RNA < 9000) %>%
#     mutate(test = percent.mt < 10 & nFeature_RNA > 600 & nFeature_RNA < 6000) %>% 
#     pull(test)
#   
#   # add the filtering variable based on the
#   stats <- cbind(log10(datasc@meta.data$nCount_RNA), log10(datasc@meta.data$nFeature_RNA),
#                  datasc@meta.data$percent.mt)
#   
#   # library(robustbase)
#   outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
#   #library(scater)
#   multi.outlier <- isOutlier(outlying, type = "higher")
#   # summary(multi.outlier)
#   
#   datasc$not_outlier <- !as.vector(multi.outlier)
#   
#   return(datasc)
# }) %>%
#   setNames(id_sample)
# 
# # plot QC -----------------------------------------------------------------
# # extract the metadata from each dataset
# meta_total <- lapply(list_datasc, function(x){
#   x@meta.data
# }) %>%
#   bind_rows(.id = "dataset") %>%
#   rownames_to_column("barcode")
# 
# meta_total %>%
#   write_tsv("../out/table/meta_datasc_test.tsv")
# 
# # how many cells are considered outliers
# meta_total %>%
#   dplyr::count(orig.ident,test)
# 
# meta_total %>%
#   dplyr::count(orig.ident,not_outlier)
# 
# # fixed threshold scatter nFeature vs percent.mt
# meta_total %>%
#   ggplot(aes(y = percent.mt,x = nFeature_RNA,col=test)) + geom_point(alpha=0.3) +
#   facet_wrap(~orig.ident) +
#   theme_bw() +
#   theme(strip.background = element_blank())
# # save the plot
# ggsave("../out/image/fixed_scatter_feature_mito_06_60.pdf",width = 12,height = 9)
# 
# meta_total %>%
#   ggplot(aes(y = percent.mt,x = nFeature_RNA,col=not_outlier)) +
#   geom_point(alpha=0.3) +
#   facet_wrap(~orig.ident) +
#   theme_bw() +
#   theme(strip.background = element_blank())
# # save the plot
# ggsave("../out/image/adaptive_scatter_feature_mito.pdf",width = 12,height = 9)
# 
# #
# meta_total %>%
#   gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt)) %>%
#   ggplot(aes(x=orig.ident,y=value)) +
#   geom_violin() +
#   geom_jitter(width = 0.2,alpha=0.01) +
#   facet_wrap(~var,scales = "free") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
#   theme(strip.background = element_blank())
# # save the plot
# ggsave("../out/image/fixed_boxplot_reads.pdf",width = 12,height = 4)
# 
# #
# meta_total %>%
#   gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt)) %>%
#   filter(var == "percent.mt") %>%
#   ggplot(aes(x=value))+geom_histogram(binwidth = 0.05) +
#   facet_wrap(orig.ident~var,scales = "free") +
#   theme_bw() +
#   scale_x_log10() +
#   geom_vline(xintercept = c(10),col="red",linetype="dashed") +
#   annotate("rect", xmin=0, xmax=10, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
#   theme(strip.background = element_blank())
# # save the plot
# ggsave("../out/image/fixed_histo_mito_10.pdf",width = 12,height = 9)
# 
# meta_total %>%
#   gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt)) %>%
#   filter(var == "nFeature_RNA") %>%
#   ggplot(aes(x=value)) + 
#   geom_histogram(binwidth = 0.05) +
#   facet_wrap(orig.ident~var,scales = "free") +
#   theme_bw() +
#   scale_x_log10() +
#   geom_vline(xintercept = c(1000,6000),col="red",linetype="dashed") +
#   annotate("rect", xmin=1000, xmax=6000, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
#   theme(strip.background = element_blank())
# # save the plot
# ggsave("../out/image/fixed_histo_features_10_60.pdf",width = 12,height = 9)
# 
# meta_total %>%
#   gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt)) %>%
#   filter(var == "nFeature_RNA") %>%
#   ggplot(aes(x=value)) + 
#   geom_histogram(binwidth = 0.05) +
#   facet_wrap(orig.ident~var,scales = "free") +
#   theme_bw() +
#   scale_x_log10() +
#   geom_vline(xintercept = c(600,6000),col="red",linetype="dashed") +
#   annotate("rect", xmin=600, xmax=6000, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
#   theme(strip.background = element_blank())
# # save the plot
# ggsave("../out/image/fixed_histo_features_06_60.pdf",width = 12,height = 9)
# 
# meta_total %>%
#   gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt)) %>%
#   filter(var == "nFeature_RNA") %>%
#   ggplot(aes(x=value)) + 
#   geom_histogram(binwidth = 0.05) +
#   facet_wrap(orig.ident~var,scales = "free") +
#   theme_bw() +
#   scale_x_log10() +
#   geom_vline(xintercept = c(200,6000),col="red",linetype="dashed") +
#   annotate("rect", xmin=200, xmax=6000, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
#   theme(strip.background = element_blank())
# # save the plot
# ggsave("../out/image/fixed_histo_features_02_60.pdf",width = 12,height = 9)
# 
# meta_total %>%
#   gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt)) %>%
#   filter(var == "nCount_RNA") %>%
#   ggplot(aes(x=value)) + 
#   geom_histogram(binwidth = 0.05) +
#   facet_wrap(orig.ident~var,scales = "free") +
#   theme_bw() +
#   scale_x_log10() +
#   # geom_vline(xintercept = c(500,5000),col="red",linetype="dashed") +
#   # annotate("rect", xmin=500, xmax=5000, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
#   theme(strip.background = element_blank())
# # save the plot
# ggsave("../out/image/fixed_histo_counts.pdf",width = 12,height = 9)
# 
# meta_total %>%
#   gather(var,value,c(nCount_RNA,nFeature_RNA,percent.mt,percent.ribo)) %>%
#   filter(var == "percent.ribo") %>%
#   ggplot(aes(x=value)) + 
#   geom_histogram(binwidth = 0.05) +
#   facet_wrap(orig.ident~var,scales = "free") +
#   theme_bw() +
#   scale_x_log10() +
#   # geom_vline(xintercept = c(500,5000),col="red",linetype="dashed") +
#   # annotate("rect", xmin=500, xmax=5000, ymin=0, ymax=Inf, alpha=0.1, fill="red")+
#   theme(strip.background = element_blank())
# # save the plot
# ggsave("../out/image/fixed_histo_ribo.pdf",width = 12,height = 9)
# 
# #
# meta_total %>%
#   dplyr::count(mt_bin,orig.ident)
# # we are going to select only the test T
# meta_total %>%
#   dplyr::count(orig.ident,mt_bin,test)
# 
# # color the bins for the amount of reads
# meta_total %>%
#   ggplot(aes(x = nCount_RNA,y = nFeature_RNA,col=mt_bin)) + geom_point(alpha=0.3) + facet_grid(orig.ident~mt_bin,scales = "free_y")+theme_bw() +
#   scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
#   theme(strip.background = element_blank())
# # save the plot
# # ggsave("out/image/fixed_threshold/fixed_scatter_mito.pdf",width = 10,height = 6)
