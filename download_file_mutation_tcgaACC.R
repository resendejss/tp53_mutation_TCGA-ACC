################################################################################
## Name: download_file_mutation_tcgaACC.R                                     ##
##                                                                            ##
## Author: Jean Resende (jean.s.s.resende@gmail.com)                          ##
## Date of last update: 12/18/2023                                            ##
##                                                                            ##
## Description: This script downloads the TCGA-ACC mutations available in the ##
##              maf file                                                      ##
################################################################################

# -- Required packages --

library(TCGAbiolinks)
library(maftools)
library(magrittr)

# -- Download MAF files -- 
query <- GDCquery(
  project = "TCGA-ACC",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")

GDCdownload(query)
maf <- GDCprepare(query)

maf_tools <- maf %>% read.maf()

save(maf, file = "maf.RData")
save(maf_tools, file = "maf_tools.RData")

# ok -- 

################################# -- end -- ####################################