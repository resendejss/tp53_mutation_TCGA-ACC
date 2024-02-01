################################################################################
## Name: comparation_TP53_CTNNB1_tcgaACC.R                                    ##
##                                                                            ##
## Author: Jean Resende (jean.s.s.resende@gmail.com)                          ##
## Date of last update: 02/01/2024                                            ##
##                                                                            ##
## Description: This script downloads the TCGA-ACC mutations available in the ##
##              maf file                                                      ##
################################################################################

library(maftools)

#load("maf_tools.RData")
load("clinical.RData")
load("maf.RData")
load("TCGA-ACC.RData")

barc_mut_tp53 <-  maf$Tumor_Sample_Barcode[maf$Hugo_Symbol=="TP53"]
barc_mut_ctnnb1 <- maf$Tumor_Sample_Barcode[maf$Hugo_Symbol=="CTNNB1"]
barc_mut_tp53_ctnnb1 <- barc_mut_tp53[barc_mut_tp53 %in% barc_mut_ctnnb1]

expr_tp53_g1 <- tcgaProject@assays@data$tpm_unstrand[
  tcgaProject@rowRanges$gene_name == "TP53",
  substr(tcgaProject$barcode,1,12) %in% substr(barc_mut_tp53,1,12)]

expr_ctnnb1_g2 <- tcgaProject@assays@data$tpm_unstrand[
  tcgaProject@rowRanges$gene_name == "CTNNB1",
  substr(tcgaProject$barcode,1,12) %in% substr(barc_mut_ctnnb1,1,12)]

expr_tp53_g3 <- tcgaProject@assays@data$tpm_unstrand[
  tcgaProject@rowRanges$gene_name == "TP53",
  substr(tcgaProject$barcode,1,12) %in% substr(barc_mut_tp53,1,12) == FALSE]

expr_tp53_g4 <- tcgaProject@assays@data$tpm_unstrand[
  tcgaProject@rowRanges$gene_name == "CTNNB1",
  substr(tcgaProject$barcode,1,12) %in% substr(barc_mut_ctnnb1,1,12) == FALSE]


