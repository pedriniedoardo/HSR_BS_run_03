# AIM ---------------------------------------------------------------------
# focus on the ferroptosis term and lysosome terms. pull the genes to build a custom signature object.

# libraries ---------------------------------------------------------------
library(tidyverse)
library(GSVA)
library(limma)
library(ComplexHeatmap)
library(AnnotationHub)
library(AnnotationDbi)
library(msigdbr)
library(ggrepel)
library(biomaRt)

# focus on the missing terms from Aletta ----------------------------------
# repriduce the query done in msigdbs to show which terms have been excluded from the analysis
# read in the list of terms provided by Aletta
TOI <- read_csv("../../data/GOterms_iron_myelin.csv")

# generate the signature file
gene_sets <- msigdbr(species = "Homo sapiens", category = "C5",subcategory = "GO:BP")
head(gene_sets)

# summarise the terms
test01_01_summary <- gene_sets %>%
  group_by(gs_exact_source,gs_name) %>%
  summarise(n = n())

#
test01_01_summary

# check that the terms provided are present in the current dataset
missing_go <- TOI %>%
  left_join(test01_01_summary,by = c("GO_id" = "gs_exact_source")) %>%
  mutate(test = !is.na(n)) %>%
  filter(is.na(n)) %>%
  pull(GO_id)

df_missing_go <- TOI %>%
  left_join(test01_01_summary,by = c("GO_id" = "gs_exact_source")) %>%
  mutate(test = !is.na(n))

TOI %>%
  left_join(test01_01_summary,by = c("GO_id" = "gs_exact_source")) %>%
  mutate(test = !is.na(n)) %>%
  print(n = 30)


# ferroptosis -------------------------------------------------------------
# read in the amigo output
df_test <- read_tsv("../../data/signatures/ferroptosis/amiGO_GO0110075.txt",col_names = F)

# how may unique genes are available in the table
df_test_summary <- df_test %>%
  group_by(X12,X19,X33,X36) %>%
  summarise() %>%
  arrange(X33)

# how many annotations are avaulable in the table
df_test$X19 %>%
  table()

# # compare this output with the one obtained from biomart.
# # the advantage of using biomart is that I can programmatically access the gene names
# 
# Select the Ensembl mart and the human dataset
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
# 
# # pull the anntaiton abailable in the database
# # keytypes(ensembl)
# # listAttributes(ensembl)
# 
# # pull the genes from individual terms
# list_test_go <- lapply(c("GO:0110075","GO:0110076","GO:0160020"), function(x){
#   genes <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id'),
#                  filters = 'go',
#                  values = x,
#                  mart = ensembl) %>%
#     mutate(GO_id = x)
#   
#   return(genes)
# })
# 
# list_test_go %>%
#   bind_rows() 
# 
# # make the gene unique
# list_test_go %>%
#   bind_rows() %>%
#   group_by(hgnc_symbol, GO_id) %>%
#   summarise()

# The GO term of interest
target_go_id_ferroptosis <- "GO:0110075"

# Use the 'go_parent_term' filter to get the term and all its children
genes_with_children_ferroptosis <- getBM(
  attributes = c('hgnc_symbol', 'go_id', 'name_1006'),
  filters = 'go_parent_term',
  values = target_go_id_ferroptosis,
  mart = ensembl
) %>%
  mutate(id = "ferroptosis_all")

# View the first few results
genes_with_children_ferroptosis %>%
  arrange(go_id,hgnc_symbol)

# lysosomes ---------------------------------------------------------------
# do the same for the lysososmes
# The GO term of interest
target_go_id_lysosome <- "GO:0097213"

# Use the 'go_parent_term' filter to get the term and all its children
genes_with_children_lysosome <- getBM(
  attributes = c('hgnc_symbol', 'go_id', 'name_1006'),
  filters = 'go_parent_term',
  values = target_go_id_lysosome,
  mart = ensembl
)

# View the first few results
genes_with_children_lysosome %>%
  arrange(go_id,hgnc_symbol)

# build siganture file ----------------------------------------------------

genes_with_children_ferroptosis %>%
  split(.$id) %>%
  lapply(function(x){
    x %>%
      dplyr::rename(Pathway = id, Genes = hgnc_symbol)
  }) %>%
  saveRDS("../../out/object/28_custom_signatures.rds")
  

