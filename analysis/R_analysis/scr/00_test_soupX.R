# libraries ---------------------------------------------------------------
library(tidyverse)
library(SoupX)
library(DropletUtils)

# READ IN DATA ------------------------------------------------------------
id_sample <- dir("../../data/cellranger61/out/") %>%
  str_subset(pattern = "^cr61",negate = T)

# # load the LUT
# LUT <- read_csv("data/clinical_data.csv") %>% 
#   mutate(sample_id = str_pad(ID,width = 2,side = "left",pad = "0")) %>% 
#   mutate(sample_id = paste0(sample_id,"_cr_61"))

# do the preprocessing over all the dataset and save the output matrices
# x <- "01_cr_61"
# file <- paste0("../../raw_data/cellranger/",x,"/outs/")
# test <- load10X(dataDir = file)

lapply(id_sample,function(x){
  # track the progress
  print(x)
  # define the location of the output of cellranger
  file <- paste0("../../data/cellranger61/out/",x,"/outs/")
  # Load data and estimate soup profile
  sc <- load10X(file)
  # Estimate rho
  sc <- autoEstCont(sc)
  # Clean the data
  out <- adjustCounts(sc,roundToInt = T)
  # save the data
  DropletUtils:::write10xCounts(paste0("../../data/SoupX_default/",x), out)
})
