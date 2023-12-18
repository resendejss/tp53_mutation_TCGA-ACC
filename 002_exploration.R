################################################################################
## Name:  ---.R                                     ##
##                                                                            ##
## Author: Jean Resende (jean.s.s.resende@gmail.com)                          ##
## Date of last update: 12/18/2023                                            ##
##                                                                            ##
## Description: ----                                                      ##
################################################################################

load("maf_tools.RData")
load("maf.RData")
load("TCGA-ACC.RData")

# -- barcodes
length(unique(maf_tools@clinical.data$Tumor_Sample_Barcode))
length(unique(tcgaProject$barcode))

head(maf_tools@clinical.data$Tumor_Sample_Barcode)
head(tcgaProject$barcode)

# -- How many genes are mutated?
sort(table(as.factor(maf$Hugo_Symbol)), decreasing = T)
length(sort(table(as.factor(maf$Hugo_Symbol)), decreasing = T))

barc_tp53 <- maf$Tumor_Sample_Barcode[maf$Hugo_Symbol=="TP53"]
barc_tp53 <- unique(barc_tp53)

barc_ctnnb1 <- maf$Tumor_Sample_Barcode[maf$Hugo_Symbol=="CTNNB1"]
barc_ctnnb1 <- unique(barc_ctnnb1)

table(barc_ctnnb1 %in% barc_tp53)

table(substr(tcgaProject$barcode, 1, 12) %in%
        substr(maf_tools@clinical.data$Tumor_Sample_Barcode, 1, 12))

tp53_rnaseq <- unique(substr(
  maf$Tumor_Sample_Barcode[maf$Hugo_Symbol=="TP53"], 1, 12))

ctnnb1_rnaseq <- unique(substr(
  maf$Tumor_Sample_Barcode[maf$Hugo_Symbol=="CTNNB1"], 1, 12))

table(tp53_rnaseq %in% tcgaProject$patient)
table(ctnnb1_rnaseq %in% tcgaProject$patient)

table(ctnnb1_rnaseq %in% tp53_rnaseq)

table(tcgaProject$patient %in% tp53_rnaseq) 
table(tcgaProject$patient %in% ctnnb1_rnaseq)