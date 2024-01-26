################################################################################
## Name: figure_1.R                                                           ##
##                                                                            ##
## Author: Jean Resende (jean.s.s.resende@gmail.com)                          ##
## Date of last update: 01/24/2023                                            ##
##                                                                            ##
## Description: This script produces the images in figure 1                   ##
################################################################################

library(maftools)

load("maf_tools.RData")
load("clinical.RData")
load("maf.RData")

Tumor_Sample_Barcode <- unique(maf$Tumor_Sample_Barcode)
#Tumor_Sample_Barcode_clinical <- Tumor_Sample_Barcode[sub]

idx <- match(clinical$submitter_id, substr(Tumor_Sample_Barcode, 1,12))
idx <- idx[!is.na(idx)]

clinical$Tumor_Sample_Barcode <- NA
clinical$Tumor_Sample_Barcode[idx] <- Tumor_Sample_Barcode

acc <- read.maf(maf=maf, clinicalData = clinical)

# -- Figure 1A
col = RColorBrewer::brewer.pal(n = 7, name = 'Accent')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 
               'Multi_Hit', 'Frame_Shift_Ins', 'Splice_Site',
               'In_Frame_Del')

#gendercolors <- c("#bdbdbd","#000000")
#names(gendercolors) <- c("female","male")
#gendercolors <- list(gender = gendercolors)

#stagecolors <- RColorBrewer::brewer.pal(n=4,name = "Greys")
#names(stagecolors) <- c("stage I", "Stage II", "Stage III", "Stage IV")
#stagecolors <- list(ajcc_pathologic_stage = stagecolors)

#oncoplot(maf = acc, colors = col,
#         clinicalFeatures = c('gender','ajcc_pathologic_stage'), top = 10, 
#         removeNonMutated = T, annotationColor = c(gendercolors, stagecolors))

oncoplot(maf = acc, colors = col, removeNonMutated = T, top = 10)

# -- Figure 1B
lollipopPlot(acc, gene = 'TP53', AACol = 'HGVSp')
lollipopPlot(acc, gene = 'CTNNB1', AACol = 'HGVSp')

lollipopPlot(maf = acc, gene = 'TP53', AACol = 'HGVSp',
             showMutationRate = TRUE, domainLabelSize = 1, defaultYaxis = FALSE)

lollipopPlot(maf = acc, gene = 'CTNNB1', AACol = 'HGVSp',
             showMutationRate = TRUE, domainLabelSize = 1, defaultYaxis = FALSE)

# -- survival
laml.mutload = tcgaCompare(maf = acc, cohortName = 'TCGA-ACC')
acc.sig <- oncodrive(maf = acc, AACol = "HGVSp", minMut = 5)
plotOncodrive(res = acc.sig, useFraction = F)

clinical_survival = acc@clinical.data
clinical_survival$Overall_Survival_Status <- ifelse(clinical_teste$vital_status == "Dead", 1,0)

mafSurvival(maf = acc, clinicalData = clinical_teste, genes = 'TP53', 
            time = 'days_to_last_follow_up', Status = 'Overall_Survival_Status', 
            showConfInt = TRUE)

mafSurvival(maf = acc, clinicalData = clinical_teste, genes = 'CTNNB1', 
            time = 'days_to_last_follow_up', Status = 'Overall_Survival_Status', 
            showConfInt = TRUE)








