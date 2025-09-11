# libraries ---------------------------------------------------------------
library(tidyverse)
library(ComplexHeatmap)
library(matrixStats)

# read in the data --------------------------------------------------------
LUT_sample <- data.frame(id_sample = c("hBS_CTR4_MG","hBS_RR16_MG","hBS_RR24_MG","hBS_RR25_MG","W8_18h_myelin_plus_untreated_multiplexed","W8_24h_CSF-controls_plus_untreated_multiplexed","W8_24h_CSF-MS_plus_untreated_multiplexed","W8_24h_cytokines","W8_48h_CSF-MS_multiplexed","W8_48h_TBHP_and_TBHP-dasatinib_multiplexed","W8_6h_Fe_RSL3"),
                         id_sample_short = c("CTR4","RR16","RR24","RR25","myelin","CSFctrl","CSFms24h","cytokine","CSFms48h","TBHP","Fe"))

# read in one run of genotyping
list_vcf <- lapply(LUT_sample$id_sample,function(x){
  file <- paste0("../genetic_demux/results/02_vireo/",x,"/GT_donors.vireo.vcf.gz")
  read_tsv(file,skip = 1)
}) %>% 
  setNames(LUT_sample$id_sample_short)

# x <- "TBHP"
# batch process to generate the matrix
list_mat <- lapply(names(list_vcf),function(x){
  
  # keep track
  print(x)
  
  # extract the GT field from every donor
  test_sample <- list_vcf[[x]] %>% 
    # for some situations I also needed to add the REF to make sure the locus is unique
    mutate(ref = paste0(`#CHROM`,"_",POS,"_",REF)) %>% 
    select(ref,contains("donor")) %>% 
    pivot_longer(names_to = "donor",values_to = "anno",-ref) %>% 
    separate(anno,into = c("GT","AD","DP","PL"),sep = ":")
  
  # filter donors that are not likely to exist. meaning donors that have very few snp count
  test_sample_stat <- test_sample %>% 
    group_by(donor) %>% 
    summarise(avg_snp_presence = mean(DP>0))
  
  # print the snp stats
  print(test_sample_stat)
  
  id_donor <- test_sample_stat %>% 
    filter(avg_snp_presence > 0.5) %>% 
    pull(donor)
  
  test_sample_filter <- test_sample %>% 
    dplyr::filter(donor %in% id_donor)
  
  # following Elisa's suggestion, consider filtering the locus that have very few reads
  
  # make a matrix
  # make a matrix out of the GT
  # mat_discrete <- test_sample_filter %>% 
  #   dplyr::select(ref,donor,GT) %>% 
  #   pivot_wider(names_from = donor,values_from = GT) %>% 
  #   column_to_rownames("ref")
  
  # test per locus the number of phenotypes
  summary_pheno <- test_sample_filter %>% 
    group_by(ref) %>% 
    summarise(test = length(unique(GT)))
  
  # if possible, filter the number of phenotypes
  if(length(id_donor) == 1){
    # keep all the locus
    df_id1 <- summary_pheno %>% 
      filter(test > 0) 
  }else if(length(id_donor) == 2){
    # in case of two donors, keep only the locus with 2 levels
    df_id1 <- summary_pheno %>% 
      filter(test > 1) 
  } else {
    # for three or more donors, keep only the locus with all three possible combinations of alleles
    df_id1 <- summary_pheno %>% 
      filter(test > 2) 
  }
  
  # df_id1 <- summary_pheno %>% 
  #   filter(test>2)
  
  # make a matrix out of the GT
  mat_discrete2 <- test_sample_filter %>%
    # filter only a subset of the ids
    filter(ref %in% df_id1$ref) %>%
    # filter only for locus with 2 levels
    filter(str_detect(GT,pattern = "1/1|1/0|0/0")) %>% 
    mutate(GT2 = case_when(GT=="1/1"~2,
                           GT=="1/0"~1,
                           GT=="0/0"~0)) %>% 
    # show only the non concordant snps
    dplyr::select(ref,donor,GT2) %>% 
    pivot_wider(names_from = donor,values_from = GT2) %>% 
    column_to_rownames("ref")
  
  return(mat_discrete2)
  
}) %>% 
  setNames(names(list_vcf))

# attempt plot ------------------------------------------------------------
# test_sample <- list_mat[[1]]

# inner join the matrices
# x <- "test_01"
mat_join <- lapply(names(list_mat), function(x){
  print(x)
  obj <- list_mat[[x]]
  colnames(obj) <- paste0(colnames(obj),".",x)
  obj %>% 
    rownames_to_column(var = "loc")
}) %>% 
  purrr::reduce(inner_join,by="loc") %>% 
  column_to_rownames(var = "loc") %>% 
  as.matrix()

# define the colors of the heatmap
colors <- structure(1:3, names = c("0", "1", "2"))

# define the column annotation
LUT <- data.frame(colname = colnames(mat_join)) %>% 
  separate(col = colname,into = c("donor.random","experiment"),sep = "\\.",remove = F) %>% 
  mutate(experiment_factor = factor(experiment,levels = c("CSFctrl",
                                                          "CSFms24h",
                                                          "CSFms48h",
                                                          "cytokine",
                                                          "myelin",
                                                          "Fe",
                                                          "TBHP",
                                                          "CTR4",
                                                          "RR16",
                                                          "RR24",
                                                          "RR25")))

column_ha <- HeatmapAnnotation(experiment = LUT$experiment_factor,
                               col = list(experiment = c("CTR4" = "gray",
                                                         "RR16" = "gray50",
                                                         "RR24" = "gray20",
                                                         "RR25" = "black",
                                                         "myelin" = "pink",
                                                         "CSFctrl" = "yellow",
                                                         "CSFms24h" = "orange",
                                                         "cytokine" = "magenta",
                                                         "CSFms48h" = "red",
                                                         "TBHP" = "cyan",
                                                         "Fe" = "blue")))

# cluster the donors based on the type of GT annotation
mat2 <- Heatmap(mat_join,top_annotation = column_ha,
                name = "mat",
                col=colors,
                column_title = "BS 03 experiment",show_row_names = F,show_row_dend = F,show_heatmap_legend = F)

lgd <- Legend(labels = c("0/0", "1/0", "1/1"),
              legend_gp = gpar(fill = 1:3),title = "GT")

pdf("../../out/image/heatmap_donor_deconvolution.pdf",height = 10,width = 5)
draw(mat2, annotation_legend_list = list(lgd))
dev.off()

# plot PCA ----------------------------------------------------------------
# calculate the variance for each gene 
object <- as.matrix(mat_join)
rv <- rowVars(object)

# select the ntop genes by variance 
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))] 

# perform a PCA on the data in assay(x) for the selected genes 
# differently from DESeq2 approach the input is not scaled or centered. therefore we need to specify at lease the .scale argument as T
pca <- prcomp(t(object[select,]),scale. = T,center = T) 

# pull the varaince explained
summary(pca)

# the contribution to the total variance for each component 
percentVar <- pca$sdev^2 / sum( pca$sdev^2) 

pca$x %>% 
  data.frame() %>% 
  rownames_to_column(var = "sample") %>% 
  separate(col = sample,into = c("donor.random","experiment"),sep = "\\.",remove = F) %>% 
  mutate(experiment_factor = factor(experiment,levels = c("CSFctrl",
                                                          "CSFms24h",
                                                          "CSFms48h",
                                                          "cytokine",
                                                          "myelin",
                                                          "Fe",
                                                          "TBHP",
                                                          "CTR4",
                                                          "RR16",
                                                          "RR24",
                                                          "RR25"))) %>% 
  ggplot(aes(x=PC1,y=PC2,label=sample,col= experiment_factor))+geom_point()+theme_bw()+ggrepel::geom_label_repel(size = 5, # font size in the text labels
                                                                                                                 min.segment.length = 0, # draw all line segments
                                                                                                                 box.padding = unit(0.35, "lines"),
                                                                                                                 point.padding = unit(0.3, "lines"),max.overlaps = 25) +
  scale_color_manual(values = c("CSFctrl" = "gold","CSFms24h" = "orange","CSFms48h" = "red","cytokine" = "magenta","myelin" = "pink","Fe" = "blue","TBHP" = "cyan","CTR4" = "gray","RR16" = "gray50","RR24" = "gray20","RR25" = "black")) +
  xlab(paste0("PC1 ",round(percentVar[1]*100,digits = 1),"%"))+
  ylab(paste0("PC2 ",round(percentVar[2]*100,digits = 1),"%"))
ggsave(filename = "../../out/image/PCA_donor_deconvolution.pdf",width = 11,height = 10)
