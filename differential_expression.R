################################################################################
## Name: differential_expression.R                                            ##
##                                                                            ##
## Author: Jean Resende (jean.s.s.resende@gmail.com)                          ##
## Date of last update: 15/02/2024                                            ##
##                                                                            ##
## Description:                                                               ##
################################################################################

setwd("~/GitHub/tp53_mutation_TCGA-ACC")

load("TCGA-ACC.RData")
load("clinical.RData")
load("maf.RData")

# -- data construction --
barc_mut_tp53 <-  maf$Tumor_Sample_Barcode[maf$Hugo_Symbol=="TP53"]
barc_mut_ctnnb1 <- maf$Tumor_Sample_Barcode[maf$Hugo_Symbol=="CTNNB1"]


tcgaProject$mut_tp53 <- ifelse(
  substr(tcgaProject$barcode, 1, 12) %in% substr(barc_mut_tp53, 1, 12),
  TRUE,FALSE)

tcgaProject$mut_ctnnb1 <- ifelse(
  substr(tcgaProject$barcode, 1, 12) %in% substr(barc_mut_ctnnb1, 1, 12),
  TRUE,FALSE)


table(tcgaProject$mut_tp53)
table(tcgaProject$mut_ctnnb1)



