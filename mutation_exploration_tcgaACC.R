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
library(maftools)

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

# -- summary maf files --
metadata <- getSampleSummary(maf_tools) # types of mutations per sample
write.csv(metadata, file = "metadata.csv")

getGeneSummary(maf_tools)
getClinicalData(maf_tools)
getFields(maf_tools)

# -- TP53 and CTNNB1 --
barcode_tp53_mut <- unique(maf$Tumor_Sample_Barcode[maf$Hugo_Symbol == "TP53"])
barcode_ctnnb1_mut <- unique(maf$Tumor_Sample_Barcode[maf$Hugo_Symbol == "CTNNB1"])

qtd_samples_tp53_mut <- length(barcode_tp53_mut)
qtd_samples_ctnnb1_mut <- length(barcode_ctnnb1_mut)
qtd_all_samples_mut <- length(unique(maf$Tumor_Sample_Barcode))

qtd_all_samples_mut - qtd_samples_tp53_mut
qtd_all_samples_mut - qtd_samples_ctnnb1_mut

# -- samples RNA-Seq --
qtd_all_samples_rnaseq <- length(tcgaProject$barcode)
barcode_rnaseq <- tcgaProject$barcode

table(substr(barcode_tp53_mut, 1, 12) %in% substr(barcode_rnaseq, 1, 12))
table(substr(barcode_ctnnb1_mut, 1, 12) %in% substr(barcode_rnaseq, 1, 12))




