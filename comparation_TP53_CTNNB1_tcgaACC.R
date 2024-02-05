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

# -- data construction --
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


data <- data.frame(barcode = tcgaProject$barcode,
                   gender = tcgaProject$gender,
                   stage = tcgaProject$ajcc_pathologic_stage,
                   age = tcgaProject$age_at_index,
                   vital_status = tcgaProject$vital_status,
                   tp53 = tcgaProject@assays@data$tpm_unstrand[
                     tcgaProject@rowRanges$gene_name == "TP53",],
                   ctnnb1 = tcgaProject@assays@data$tpm_unstrand[
                     tcgaProject@rowRanges$gene_name == "CTNNB1",])

data$tp53_mut_g1 <- ifelse(substr(data$barcode,1,12) %in% 
                             substr(barc_mut_tp53,1,12), "yes","no")

data$ctnnb1_mut_g2 <- ifelse(substr(data$barcode,1,12) %in% 
                             substr(barc_mut_ctnnb1,1,12), "yes","no")

data$tp53_g3 <- ifelse(substr(data$barcode,1,12) %in% 
                               substr(barc_mut_tp53,1,12), "no","yes")

data$ctnnb1_g4 <- ifelse(substr(data$barcode,1,12) %in% 
                         substr(barc_mut_ctnnb1,1,12), "no","yes")

# -- normality test --
shapiro.test(data$tp53) # p-value = 1.671e-05
nortest::ad.test(data$tp53) # p-value = 0.0004598

shapiro.test(data$ctnnb1) # p-value = 9.876e-12
nortest::ad.test(data$ctnnb1) # p-value < 2.2e-16

hist(data$tp53)
hist(data$ctnnb1)

# -- comparation --
library(ggplot2)
library(ggpubr)

ggplot(data, aes(x=tp53_mut_g1, y=log2(tp53), fill=tp53_mut_g1))+
  ggtitle("Mutation - TP53")+
  geom_boxplot(width=0.4, lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("Expression (log2(TPM))")+
  scale_x_discrete("")+
  scale_fill_manual(values=c("#8dd3c7","#ffffb3"))+
  geom_jitter(width = 0.1, alpha = 0.2)+
  stat_compare_means(method="wilcox.test", label="p.format")+
  theme(legend.position = "bottom", 
        legend.title = element_text(face="bold", size=10),
        legend.key=element_rect(size=10, color=NA), 
        legend.key.size=unit(8,"mm"),
        legend.text=element_text(size=10), 
        legend.direction = "horizontal",
        legend.box = "horizontal")

pdf(file = "fig_mutation_ctnnb1.pdf", width = 7, height = 5)
ggplot(data, aes(x=ctnnb1_mut_g2, y=log2(ctnnb1), fill=ctnnb1_mut_g2))+
  ggtitle("Mutation - CTNNB1")+
  geom_boxplot(width=0.4, lwd=0.2, outlier.size = 0.5)+
  scale_y_continuous("Expression (log2(TPM))")+
  scale_x_discrete("")+
  scale_fill_manual(values=c("#bebada","#fb8072"))+
  #geom_jitter(width = 0.1, alpha = 0.2)+
  geom_rug(data = data[data$ctnnb1_mut_g2 == "no",],
           aes(x=NULL),
           sides = "l",
           colour = "#bebada")+
  geom_rug(data = data[data$ctnnb1_mut_g2 == "yes",],
           aes(x=NULL),
           sides = "r",
           colour = "#fb8072")+
  stat_compare_means(method="wilcox.test",
                     label="p.format",
                     label.x = c(1.3))+
  theme(legend.position = "bottom", 
        legend.title = element_text(face="bold", size=12),
        legend.key=element_rect(size=12, color=NA), 
        legend.key.size=unit(8,"mm"),
        legend.text=element_text(size=12), 
        legend.direction = "horizontal",
        legend.box = "horizontal",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14))
dev.off()

