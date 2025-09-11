# # libraries ---------------------------------------------------------------
# library(tidyverse)
# library(ComplexHeatmap)
# 
# # read in the data --------------------------------------------------------
# # read in one run of genotyping
# test_01 <- read_tsv("../genetic_demux/results/02_vireo/W8_24h_CSF-controls_plus_untreated_multiplexed/GT_donors.vireo.vcf.gz",skip = 1)
# 
# # read in another run of genotyping
# test_02 <- read_tsv("../genetic_demux/results/02_vireo/W8_24h_cytokines/GT_donors.vireo.vcf.gz",skip = 1)
# 
# # are there any qual different from .
# table(test_01$QUAL)
# table(test_01$FILTER)
# 
# # set up a list for batch processing
# list_vcf <- list(test_01 = test_01,test_02 = test_02)
# 
# # x <- "test_01"
# # batch process to generate the matrix
# list_mat <- lapply(names(list_vcf),function(x){
#   
#   # keep track
#   print(x)
#   
#   # extract the GT field from every donor
#   test_sample <- list_vcf[[x]] %>% 
#     mutate(ref = paste0(`#CHROM`,"_",POS)) %>% 
#     select(ref,contains("donor")) %>% 
#     pivot_longer(names_to = "donor",values_to = "anno",-ref) %>% 
#     separate(anno,into = c("GT","AD","DP","PL"),sep = ":")
#   
#   # filter donors that are not likely to exist. meaning donors that have very few snp count
#   test_sample_stat <- test_sample %>% 
#     group_by(donor) %>% 
#     summarise(avg_snp_presence = mean(DP>0))
#   
#   # print the snp stats
#   print(test_sample_stat)
#   
#   id_donor <- test_sample_stat %>% 
#     filter(avg_snp_presence > 0.1) %>% 
#     pull(donor)
#   
#   test_sample_filter <- test_sample %>% 
#     dplyr::filter(donor %in% id_donor)
#   
#   # make a matrix
#   # make a matrix out of the GT
#   # mat_discrete <- test_sample_filter %>% 
#   #   dplyr::select(ref,donor,GT) %>% 
#   #   pivot_wider(names_from = donor,values_from = GT) %>% 
#   #   column_to_rownames("ref")
#   
#   # test per locus the number of phenotypes
#   summary_pheno <- test_sample_filter %>% 
#     group_by(ref) %>% 
#     summarise(test = length(unique(GT)))
#   
#   # filter only the locus where the number of phenotypes is > 2 to further reduce the pool
#   df_id1 <- summary_pheno %>% 
#     filter(test>2)
#   
#   # df_id1 <- summary_pheno %>% 
#   #   filter(test>2)
#   
#   # make a matrix out of the GT
#   mat_discrete2 <- test_sample_filter %>%
#     # filter only a subset of the ids
#     filter(ref %in% df_id1$ref) %>%
#     mutate(GT2 = case_when(GT=="1/1"~3,
#                            GT=="1/0"~2,
#                            GT=="0/0"~1)) %>% 
#     # show only the non concordant snps
#     dplyr::select(ref,donor,GT2) %>% 
#     pivot_wider(names_from = donor,values_from = GT2) %>% 
#     column_to_rownames("ref")
#   
#   return(mat_discrete2)
#   
# }) %>% 
#   setNames(names(list_vcf))
# 
# # attempt plot ------------------------------------------------------------
# test_sample <- list_mat[[1]]
# 
# # inner join the matrices
# # x <- "test_01"
# mat_join <- lapply(names(list_mat), function(x){
#   print(x)
#   obj <- list_mat[[x]]
#   colnames(obj) <- paste0(colnames(obj),".",x)
#   obj %>% 
#     rownames_to_column(var = "loc")
# }) %>% 
#   purrr::reduce(inner_join,by="loc") %>% 
#   column_to_rownames(var = "loc") %>% 
#   as.matrix()
# 
# # define the colors of the heatmap
# colors <- structure(1:3, names = c("1", "2", "3"))
# 
# # cluster the donors based on the type of GT annotation
# mat2 <- Heatmap(mat_join,
#                 name = "mat",
#                 col=colors,
#                 column_title = "test experiment",show_row_names = F,show_row_dend = F,show_heatmap_legend = F)
# 
# lgd <- Legend(labels = c("0/0", "1/0", "1/1"),
#               legend_gp = gpar(fill = 1:3),title = "GT")
# 
# draw(mat2, annotation_legend_list = list(lgd))
# 
# # plot PCA ----------------------------------------------------------------
# # calculate the variance for each gene 
# library(matrixStats)
# object <- as.matrix(mat_join)
# rv <- rowVars(object)
# 
# # select the ntop genes by variance 
# select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))] 
# 
# # perform a PCA on the data in assay(x) for the selected genes 
# pca <- prcomp(t(object[select,])) 
# 
# # the contribution to the total variance for each component 
# percentVar <- pca$sdev^2 / sum( pca$sdev^2 ) 
# 
# pca$x %>% 
#   data.frame() %>% 
#   rownames_to_column(var = "sample") %>% 
#   ggplot(aes(x=PC1,y=PC2,label=sample))+geom_point()+theme_bw()+ggrepel::geom_label_repel()
# 
