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
load("clinical.RData")

length(unique(maf$Tumor_Sample_Barcode)) # number of samples (90)
length(unique(maf_tools@clinical.data$Tumor_Sample_Barcode)) # number of samples (90)
sort(table(as.factor(maf$Hugo_Symbol)), decreasing = T) # top of mutated genes (number of mutations)
length(sort(table(as.factor(maf$Hugo_Symbol)), decreasing = T)) # number of mutated genes (5510)
length(unique(maf$Hugo_Symbol)) # number of mutated genes (5510)

length(table(as.factor(maf$Hugo_Symbol))[table(as.factor(maf$Hugo_Symbol)) > 1])

# -- mutations
barcodes_pac_mutations <- unique(substr(maf$Tumor_Sample_Barcode, 1, 12))
all(barcodes_pac_mutations %in% clinical$submitter_id)
gender_mutations <- clinical$gender[clinical$submitter_id %in% 
                                      barcodes_pac_mutations]
table(gender_mutations)

age_mutations <- clinical$age_at_index[clinical$submitter_id %in% 
                                         barcodes_pac_mutations]
mean(age_mutations)
range(age_mutations)
  

# -- gene expression
barcodes_pac_expression <- tcgaProject$patient
all(barcodes_pac_expression %in% clinical$submitter_id)
gender_expression <- clinical$gender[clinical$submitter_id %in% 
                                       barcodes_pac_expression]
table(gender_expression)

age_expression <- clinical$age_at_index[clinical$submitter_id %in%
                                          barcodes_pac_expression]
mean(age_expression)
range(age_expression)

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

table(substr(barcode_tp53_mut, 1,12) %in% substr(tcgaProject$barcode, 1,12))
table(substr(barcode_ctnnb1_mut, 1,12) %in% substr(tcgaProject$barcode, 1,12)) 


# -- samples RNA-Seq --
qtd_all_samples_rnaseq <- length(tcgaProject$barcode)
barcode_rnaseq <- tcgaProject$barcode

table(substr(barcode_tp53_mut, 1, 12) %in% substr(barcode_rnaseq, 1, 12))
table(substr(barcode_ctnnb1_mut, 1, 12) %in% substr(barcode_rnaseq, 1, 12))

# -- qtd gene
length(unique(maf$Hugo_Symbol))

# -- variant classification
unique(maf$Variant_Classification)
table(maf$Variant_Classification)
length(maf$Variant_Classification)

length(unique(maf$Hugo_Symbol[maf$Variant_Classification == "Silent"]))
length(unique(maf$Hugo_Symbol[maf$Variant_Classification == "Missense_Mutation"]))
length(unique(maf$Hugo_Symbol[maf$Variant_Classification == "Splice_Site"]))
length(unique(maf$Hugo_Symbol[maf$Variant_Classification == "Frame_Shift_Del"]))
length(unique(maf$Hugo_Symbol[maf$Variant_Classification == "RNA"]))
length(unique(maf$Hugo_Symbol[maf$Variant_Classification == "In_Frame_Del"]))
length(unique(maf$Hugo_Symbol[maf$Variant_Classification == "Nonsense_Mutation"]))
length(unique(maf$Hugo_Symbol[maf$Variant_Classification == "Frame_Shift_Ins"]))
length(unique(maf$Hugo_Symbol[maf$Variant_Classification == "5'UTR"]))
length(unique(maf$Hugo_Symbol[maf$Variant_Classification == "Intron"]))
length(unique(maf$Hugo_Symbol[maf$Variant_Classification == "Splice_Region"]))
length(unique(maf$Hugo_Symbol[maf$Variant_Classification == "3'Flank"]))
length(unique(maf$Hugo_Symbol[maf$Variant_Classification == "Translation_Start_Site"]))
length(unique(maf$Hugo_Symbol[maf$Variant_Classification == "Nonstop_Mutation"]))
length(unique(maf$Hugo_Symbol[maf$Variant_Classification == "3'UTR"]))
length(unique(maf$Hugo_Symbol[maf$Variant_Classification == "5'Flank"]))
length(unique(maf$Hugo_Symbol[maf$Variant_Classification == "In_Frame_Ins"]))
