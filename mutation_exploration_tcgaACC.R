################################################################################
## Name: mutation_exploration_tcgaACC.R                                       ##
##                                                                            ##
## Author: Jean Resende (jean.s.s.resende@gmail.com)                          ##
## Date of last update: 01/11/2023                                            ##
##                                                                            ##
## Description: This script downloads the TCGA-ACC mutations available in the ##
##              maf file                                                      ##
################################################################################

# -- Required packages --

# -- exploration --
load("maf.RData")
load("maf_tools.RData")
load("TCGA-ACC.RData")

length(unique(maf$Tumor_Sample_Barcode)) # number of samples (90)
length(unique(maf_tools@clinical.data$Tumor_Sample_Barcode)) # number of samples (90)
sort(table(as.factor(maf$Hugo_Symbol)), decreasing = T) # top of mutated genes (number of mutations)
length(sort(table(as.factor(maf$Hugo_Symbol)), decreasing = T)) # number of mutated genes (5510)
length(unique(maf$Hugo_Symbol)) # number of mutated genes (5510)

length(table(as.factor(maf$Hugo_Symbol))[table(as.factor(maf$Hugo_Symbol)) > 1])






