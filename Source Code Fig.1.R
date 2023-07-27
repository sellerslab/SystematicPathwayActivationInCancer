######## Source code Figure1

#load code
library(readr)
library(tidyverse)
library(ggrepel)
library(readxl)
library(tidyverse)
library(select)
library(miceadds)
require(ggrepel)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(magrittr)
library(reshape2)
library(circlize)
library(ComplexHeatmap)
`%ni%` <- Negate(`%in%`)

##setwd to "fig1_plottingfiles"
####Figure 1a GENIE analysis

#Load Mutation Datasets
#read Foundation1
#data accessed from https://gdc.cancer.gov/about-gdc/contributed-genomic-data-cancer-research/foundation-medicine/foundation-medicine
F1_Maf <- read_rds("Foundation1_combinedMAF.rds") 
F1_clinical <- read_rds(("Foundation1_combinedsampleID.rds")) 

#read COSMIC hotspot mutation for filtering
#accessed from https://cancer.sanger.ac.uk/cosmic/file_download_info?data=GRCh38%2Fcosmic%2Fv89%2FCosmicMutantExportCensus.tsv.gz
CosmicMutant <- read_table2("CosmicMutantExportCensus.tsv") 
COSMIC_codon <- paste0(CosmicMutant$Gene,CosmicMutant$`1`)
COSMIC_AA <- paste0(CosmicMutant$Gene,CosmicMutant$Site_1)
F1_mut_codon <- paste0(F1_Maf$Hugo_Symbol,F1_Maf$HGVSc)
F1_Maf_mut_AA <- paste0(F1_Maf$Hugo_Symbol,F1_Maf$HGVSp_Short)
F1.mut.hot <- F1_Maf[((F1_mut_codon %in% COSMIC_codon) | (F1_Maf_mut_AA %in% COSMIC_AA)),] 

#read GENIE 7.0
#accessed from https://www.synapse.org/#!Synapse:syn7222066/wiki/410924
Genie7.clinical <- read.delim("GENIE7public/data_clinical_sample_7.0-public.txt", comment.char="#")
genie7.mut <- read_delim("GENIE7public/data_mutations_extended_7.0-public.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE) ##mutation data
genie_centers <- unique(genie7.mut$Center) 
genie_mut_codon <- paste0(genie7.mut$Hugo_Symbol,genie7.mut$HGVSc)
genie_mut_AA <- paste0(genie7.mut$Hugo_Symbol,genie7.mut$HGVSp_Short)
genie.mut.hot <- genie7.mut[((genie_mut_AA %in% COSMIC_AA)),] 

#read Pan-can TCGA
#accessed from https://xenabrowser.net/datapages/?cohort=TCGA%20Pan-Cancer%20(PANCAN)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
TCGA_mutation <- read_delim("mc3.v0.2.8.PUBLIC.xena", "\t", escape_double = FALSE, trim_ws = TRUE)
saveRDS(TCGA_mutation, file = "TCGA_mutation.rds",compress = TRUE)
pancanTCGA_clinical <- read_delim("Pan Can TCGA/TCGA_phenotype_denseDataOnlyDownload (1).tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
TCGA_Maf_mut_AA <- paste0(TCGA_mutation$gene,TCGA_mutation$Amino_Acid_Change)
TCGA.mut.hot <- TCGA_mutation[(TCGA_Maf_mut_AA %in% COSMIC_AA),] 


#Select Analysis in specific cohorts of cancer by combining datasets first

#lung
#subset lung cancer samples and case list from TCGA
#sample subset
table(pancanTCGA_clinical$`_primary_disease`)
LUAD_pancanTCGA_clinical <- pancanTCGA_clinical %>% filter(`_primary_disease`=="lung adenocarcinoma")
panlung_pancanTCGA_clinical <- pancanTCGA_clinical %>% filter(`_primary_disease` %in% c("lung adenocarcinoma","lung squamous cell carcinoma"))

Lung_F1_clinical <- F1_clinical %>% filter(`cases.primary_site`=="Bronchus And Lung")
M1.f1_hot_lung <- table(((F1.mut.hot %>% dplyr::select(case_id,Hugo_Symbol))))
AL1.f1_hot_lung <- new.AL(M1.f1_hot_lung )
M2.f1_hot_lung <- as.character(Lung_F1_clinical$diagnoses.primary_diagnosis[match(rownames(M1.f1_hot_lung), F1_clinical$case_id)])
names(M2.f1_hot_lung)<- rownames(M1.f1_hot_lung)
AL1.f1_hot_lung$samples$tumor.subtype <- M2.f1_hot_lung
M3.f1_hot_lung <- rep("MUTATION", ncol(M1.f1_hot_lung))
names(M3.f1_hot_lung) <- colnames(M1.f1_hot_lung)
AL1.f1_hot_lung$alterations$alteration.type <- M3.f1_hot_lung
result.F1.allhot = select::select(M = AL1.13$am, sample.class = AL1.13$samples$tumor.subtype, alteration.class = AL1.13$alterations$alteration.type,folder = "F1 hot_lung", save.intermediate.files = T, verbose = T, remove.0.samples = TRUE)

result.F1.hot.lung <- result.F1.allhot

#subset lung cancer samples and case list from DFCI and MSKCC
DM.hotmut <- genie.mut.hot %>% filter(Center %in% c("DFCI","MSK")) #MSKCC and DFCI only

#clinical processing
lung.Genie7.clinical <- Genie7.clinical %>% filter(CANCER_TYPE == "Non-Small Cell Lung Cancer")
Lung_F1_clinical <- F1_clinical %>% filter(`cases.primary_site`=="Bronchus And Lung")
panlung_pancanTCGA_clinical <- pancanTCGA_clinical %>% filter(`_primary_disease` %in% c("lung adenocarcinoma","lung squamous cell carcinoma"))

M2.f1_hot_lung <- as.character(Lung_F1_clinical$diagnoses.primary_diagnosis[match(rownames(M1.f1_hot_lung), F1_clinical$case_id)])


combined.clinical.lung <- lung.Genie7.clinical$PATIENT_ID[lung.Genie7.clinical$PATIENT_ID %in% DM.hotmut$Tumor_Sample_Barcode]
#combine samples
T.hotmut.lung <- TCGA.mut.hot%>% dplyr::select(gene,sample) %>% filter(sample %in% panlung_pancanTCGA_clinical$sample)
F.hotmut.lung <- F1.mut.hot %>% dplyr::select(gene=Hugo_Symbol,sample = case_id) %>% filter(sample %in% Lung_F1_clinical$case_id)
DM.hotmut.lung <- DM.hotmut %>% dplyr::select(gene=Hugo_Symbol,sample = Tumor_Sample_Barcode) %>%filter(sample %in% lung.Genie7.clinical$SAMPLE_ID)
commongene.lung <- Reduce(intersect, list(unique(T.hotmut.lung$gene),unique(F.hotmut.lung$gene),unique(DM.hotmut.lung$gene)))
combined.hotmut.lung <- bind_rows(T.hotmut.lung,F.hotmut.lung,DM.hotmut.lung) %>% filter(gene %in%commongene.lung) #combined filtered sample

#run GENIE
M1.COMBINED_hot_lung <- table(combined.hotmut.lung)
AL1.COMBINED_hot_lung <- new.AL(M1.COMBINED_hot_lung,samples.on.row = FALSE)
?new.AL
M2.COMBINED_hot_lung <- rep("Lung Cancer", ncol(M1.COMBINED_hot_lung))
names(M2.COMBINED_hot_lung)<- colnames(M1.COMBINED_hot_lung)
AL1.COMBINED_hot_lung$samples$tumor.subtype <- M2.COMBINED_hot_lung

M3.COMBINED_hot_lung <- rep("MUTATION", nrow(M1.COMBINED_hot_lung))
names(M3.COMBINED_hot_lung) <- rownames(M1.COMBINED_hot_lung)
AL1.COMBINED_hot_lung$alterations$alteration.type <- M3.COMBINED_hot_lung
setwd("~/Sellers Lab Dropbox/Liang Chang/Bioinfo analysis/20191219 Comprehensive Bioinformatics of ME/Results")
result.COMBINED.hot_lung = select::select(M = AL1.COMBINED_hot_lung$am, sample.class = AL1.COMBINED_hot_lung$samples$tumor.subtype, alteration.class = AL1.COMBINED_hot_lung$alterations$alteration.type,folder = "COMBINED hot_lung", save.intermediate.files = T, verbose = T, remove.0.samples = TRUE)


#Colon
colon.Genie7.clinical <- Genie7.clinical %>% filter(CANCER_TYPE == "Colorectal Cancer")
table(pancanTCGA_clinical$`_primary_disease`)
Colon_F1_clinical <- F1_clinical %>% filter(`cases.primary_site`%in% c("Colon","Rectum"))
colon_pancanTCGA_clinical <- pancanTCGA_clinical %>% filter(`_primary_disease` %in% c("rectum adenocarcinoma","colon adenocarcinoma"))

#combine samples
T.hotmut.colon <- TCGA.mut.hot%>% dplyr::select(sample,gene,) %>% filter(sample %in% colon_pancanTCGA_clinical$sample)
F.hotmut.colon <- F1.mut.hot %>% dplyr::select(sample = case_id,gene=Hugo_Symbol) %>% filter(sample %in% Colon_F1_clinical$case_id)
DM.hotmut.colon <- DM.hotmut %>% dplyr::select(sample = Tumor_Sample_Barcode,gene=Hugo_Symbol) %>%filter(sample %in% colon.Genie7.clinical$SAMPLE_ID)
commongene.colon <- Reduce(intersect, list(unique(T.hotmut.colon$gene),unique(F.hotmut.colon$gene),unique(DM.hotmut.colon$gene)))
combined.hotmut.colon <- bind_rows(T.hotmut.colon,F.hotmut.colon,DM.hotmut.colon) %>% filter(gene %in%commongene.colon) #combined filtered sample
#run GENIE
M1.COMBINED_hot_colon <- table(combined.hotmut.colon)
AL1.COMBINED_hot_colon <- new.AL(M1.COMBINED_hot_colon)

M2.COMBINED_hot_colon <- rep("Colon Cancer", nrow(M1.COMBINED_hot_colon))
names(M2.COMBINED_hot_colon)<- rownames(M1.COMBINED_hot_colon)
AL1.COMBINED_hot_colon$samples$tumor.subtype <- M2.COMBINED_hot_colon

M3.COMBINED_hot_colon <- rep("MUTATION", ncol(M1.COMBINED_hot_colon))
names(M3.COMBINED_hot_colon) <- colnames(M1.COMBINED_hot_colon)
AL1.COMBINED_hot_colon$alterations$alteration.type <- M3.COMBINED_hot_colon
setwd("~/Sellers Lab Dropbox/Liang Chang/Bioinfo analysis/20191219 Comprehensive Bioinformatics of ME/Results")
result.COMBINED.hot_colon = select::select(M = AL1.COMBINED_hot_colon$am, sample.class = AL1.COMBINED_hot_colon$samples$tumor.subtype, alteration.class = AL1.COMBINED_hot_colon$alterations$alteration.type,folder = "COMBINED hot_colon", save.intermediate.files = T, verbose = T, remove.0.samples = TRUE)


#breast
breast.Genie7.clinical <- Genie7.clinical %>% filter(CANCER_TYPE == "Breast Cancer")
Breast_F1_clinical <- F1_clinical %>% filter(`cases.primary_site`%in% c("Breast"))
breast_pancanTCGA_clinical <- pancanTCGA_clinical %>% filter(`_primary_disease` %in% c("breast invasive carcinoma"))

#combine samples
T.hotmut.breast <- TCGA.mut.hot%>% dplyr::select(gene,sample) %>% filter(sample %in% breast_pancanTCGA_clinical$sample)
F.hotmut.breast <- F1.mut.hot %>% dplyr::select(gene=Hugo_Symbol,sample = case_id) %>% filter(sample %in% Breast_F1_clinical$case_id)
DM.hotmut.breast <- DM.hotmut %>% dplyr::select(gene=Hugo_Symbol,sample = Tumor_Sample_Barcode) %>%filter(sample %in% breast.Genie7.clinical$SAMPLE_ID)
commongene.breast <- Reduce(intersect, list(unique(T.hotmut.breast$gene),unique(F.hotmut.breast$gene),unique(DM.hotmut.breast$gene)))
combined.hotmut.breast <- bind_rows(T.hotmut.breast,F.hotmut.breast,DM.hotmut.breast) %>% filter(gene %in%commongene.breast) #combined filtered sample

#run GENIE
M1.COMBINED_hot_breast <- table(combined.hotmut.breast %>% select(sample,gene))
AL1.COMBINED_hot_breast <- new.AL(M1.COMBINED_hot_breast)

M2.COMBINED_hot_breast <- rep("Breast Cancer", nrow(M1.COMBINED_hot_breast))
names(M2.COMBINED_hot_breast)<- rownames(M1.COMBINED_hot_breast)
AL1.COMBINED_hot_breast$samples$tumor.subtype <- M2.COMBINED_hot_breast

M3.COMBINED_hot_breast <- rep("MUTATION", ncol(M1.COMBINED_hot_breast))
names(M3.COMBINED_hot_breast) <- colnames(M1.COMBINED_hot_breast)
AL1.COMBINED_hot_breast$alterations$alteration.type <- M3.COMBINED_hot_breast
setwd("~/Sellers Lab Dropbox/Liang Chang/Bioinfo analysis/20191219 Comprehensive Bioinformatics of ME/Results")
result.COMBINED.hot_breast = select::select(M = AL1.COMBINED_hot_breast$am, sample.class = AL1.COMBINED_hot_breast$samples$tumor.subtype, alteration.class = AL1.COMBINED_hot_breast$alterations$alteration.type,folder = "COMBINED hot_breast", save.intermediate.files = T, verbose = T, remove.0.samples = TRUE)

#brain
GBM.Genie7.clinical <- Genie7.clinical %>% filter(CANCER_TYPE == "Glioma")
GBM_F1_clinical <- F1_clinical %>% filter(`cases.disease_type`%in% c("Gliomas"))
GBM_pancanTCGA_clinical <- pancanTCGA_clinical %>% filter(`_primary_disease` %in% c("glioblastoma multiforme","brain lower grade glioma"))

#combine samples
T.hotmut.GBM <- TCGA.mut.hot%>% dplyr::select(gene,sample) %>% filter(sample %in% GBM_pancanTCGA_clinical$sample)
F.hotmut.GBM <- F1.mut.hot %>% dplyr::select(gene=Hugo_Symbol,sample = case_id) %>% filter(sample %in% GBM_F1_clinical$case_id)
DM.hotmut.GBM <- DM.hotmut %>% dplyr::select(gene=Hugo_Symbol,sample = Tumor_Sample_Barcode) %>%filter(sample %in% GBM.Genie7.clinical$SAMPLE_ID)
commongene.GBM <- Reduce(intersect, list(unique(T.hotmut.GBM$gene),unique(F.hotmut.GBM$gene),unique(DM.hotmut.GBM$gene)))
combined.hotmut.GBM <- bind_rows(T.hotmut.GBM,F.hotmut.GBM,DM.hotmut.GBM) %>% filter(gene %in%commongene.GBM) #combined filtered sample

#run GENIE
M1.COMBINED_hot_GBM <- table(combined.hotmut.GBM %>% select(sample,gene))
AL1.COMBINED_hot_GBM <- new.AL(M1.COMBINED_hot_GBM)

M2.COMBINED_hot_GBM <- rep("GBM Cancer", nrow(M1.COMBINED_hot_GBM))
names(M2.COMBINED_hot_GBM)<- rownames(M1.COMBINED_hot_GBM)
AL1.COMBINED_hot_GBM$samples$tumor.subtype <- M2.COMBINED_hot_GBM

M3.COMBINED_hot_GBM <- rep("MUTATION", ncol(M1.COMBINED_hot_GBM))
names(M3.COMBINED_hot_GBM) <- colnames(M1.COMBINED_hot_GBM)
AL1.COMBINED_hot_GBM$alterations$alteration.type <- M3.COMBINED_hot_GBM
setwd("~/Sellers Lab Dropbox/Liang Chang/Bioinfo analysis/20191219 Comprehensive Bioinformatics of ME/Results")
result.COMBINED.hot_GBM = select::select(M = AL1.COMBINED_hot_GBM$am, sample.class = AL1.COMBINED_hot_GBM$samples$tumor.subtype, alteration.class = AL1.COMBINED_hot_GBM$alterations$alteration.type,folder = "COMBINED hot_GBM", save.intermediate.files = T, verbose = T, remove.0.samples = TRUE)

#skin melanoma
SKCM.Genie7.clinical <- Genie7.clinical %>% filter(ONCOTREE_CODE%in% c("SKCM","MEL","MUP","ACRM"))#remove uveal
SKCM_F1_clinical <- F1_clinical %>% filter(`cases.disease_type`%in% c("Nevi and Melanomas")) %>% filter(cases.primary_site != "Eye And Adnexa")#remove uveal
SKCM_pancanTCGA_clinical <- pancanTCGA_clinical %>% filter(`_primary_disease` %in% c("skin cutaneous melanoma"))

#combine samples
T.hotmut.SKCM <- TCGA.mut.hot%>% dplyr::select(gene,sample) %>% filter(sample %in% SKCM_pancanTCGA_clinical$sample)
F.hotmut.SKCM <- F1.mut.hot %>% dplyr::select(gene=Hugo_Symbol,sample = case_id) %>% filter(sample %in% SKCM_F1_clinical$case_id)
DM.hotmut.SKCM <- DM.hotmut %>% dplyr::select(gene=Hugo_Symbol,sample = Tumor_Sample_Barcode) %>%filter(sample %in% SKCM.Genie7.clinical$SAMPLE_ID)
commongene.SKCM <- Reduce(intersect, list(unique(T.hotmut.SKCM$gene),unique(F.hotmut.SKCM$gene),unique(DM.hotmut.SKCM$gene)))
combined.hotmut.SKCM <- bind_rows(T.hotmut.SKCM,F.hotmut.SKCM,DM.hotmut.SKCM) %>% filter(gene %in%commongene.SKCM) #combined filtered sample

#run GENIE
M1.COMBINED_hot_SKCM <- table(combined.hotmut.SKCM %>% select(sample,gene))
AL1.COMBINED_hot_SKCM <- new.AL(M1.COMBINED_hot_SKCM)

M2.COMBINED_hot_SKCM <- rep("SKCM Cancer", nrow(M1.COMBINED_hot_SKCM))
names(M2.COMBINED_hot_SKCM)<- rownames(M1.COMBINED_hot_SKCM)
AL1.COMBINED_hot_SKCM$samples$tumor.subtype <- M2.COMBINED_hot_SKCM

M3.COMBINED_hot_SKCM <- rep("MUTATION", ncol(M1.COMBINED_hot_SKCM))
names(M3.COMBINED_hot_SKCM) <- colnames(M1.COMBINED_hot_SKCM)
AL1.COMBINED_hot_SKCM$alterations$alteration.type <- M3.COMBINED_hot_SKCM
setwd("~/Sellers Lab Dropbox/Liang Chang/Bioinfo analysis/20191219 Comprehensive Bioinformatics of ME/Results")
result.COMBINED.hot_SKCM = select::select(M = AL1.COMBINED_hot_SKCM$am, sample.class = AL1.COMBINED_hot_SKCM$samples$tumor.subtype, alteration.class = AL1.COMBINED_hot_SKCM$alterations$alteration.type,folder = "COMBINED hot_SKCM", save.intermediate.files = T, verbose = F, remove.0.samples = TRUE)


##pancreatic
PANCREATIC.Genie7.clinical <- Genie7.clinical %>% filter(CANCER_TYPE == "Pancreatic Cancer")
PANCREATIC_F1_clinical <- F1_clinical %>% filter(`cases.primary_site`%in% c("Pancreas"))
PANCREATIC_pancanTCGA_clinical <- pancanTCGA_clinical %>% filter(`_primary_disease` %in% c("pancreatic adenocarcinoma"))

#combine samples
T.hotmut.PANCREATIC <- TCGA.mut.hot%>% dplyr::select(gene,sample) %>% filter(sample %in% PANCREATIC_pancanTCGA_clinical$sample)
F.hotmut.PANCREATIC <- F1.mut.hot %>% dplyr::select(gene=Hugo_Symbol,sample = case_id) %>% filter(sample %in% PANCREATIC_F1_clinical$case_id)
DM.hotmut.PANCREATIC <- DM.hotmut %>% dplyr::select(gene=Hugo_Symbol,sample = Tumor_Sample_Barcode) %>%filter(sample %in% PANCREATIC.Genie7.clinical$SAMPLE_ID)
commongene.PANCREATIC <- Reduce(intersect, list(unique(T.hotmut.PANCREATIC$gene),unique(F.hotmut.PANCREATIC$gene),unique(DM.hotmut.PANCREATIC$gene)))
combined.hotmut.PANCREATIC <- bind_rows(T.hotmut.PANCREATIC,F.hotmut.PANCREATIC,DM.hotmut.PANCREATIC) %>% filter(gene %in%commongene.PANCREATIC) #combined filtered sample

#run GENIE
M1.COMBINED_hot_PANCREATIC <- table(combined.hotmut.PANCREATIC %>% select(sample,gene))
AL1.COMBINED_hot_PANCREATIC <- new.AL(M1.COMBINED_hot_PANCREATIC)

M2.COMBINED_hot_PANCREATIC <- rep("PANCREATIC Cancer", nrow(M1.COMBINED_hot_PANCREATIC))
names(M2.COMBINED_hot_PANCREATIC)<- rownames(M1.COMBINED_hot_PANCREATIC)
AL1.COMBINED_hot_PANCREATIC$samples$tumor.subtype <- M2.COMBINED_hot_PANCREATIC

M3.COMBINED_hot_PANCREATIC <- rep("MUTATION", ncol(M1.COMBINED_hot_PANCREATIC))
names(M3.COMBINED_hot_PANCREATIC) <- colnames(M1.COMBINED_hot_PANCREATIC)
AL1.COMBINED_hot_PANCREATIC$alterations$alteration.type <- M3.COMBINED_hot_PANCREATIC
setwd("~/Sellers Lab Dropbox/Liang Chang/Bioinfo analysis/20191219 Comprehensive Bioinformatics of ME/Results")
result.COMBINED.hot_PANCREATIC = select::select(M = AL1.COMBINED_hot_PANCREATIC$am, sample.class = AL1.COMBINED_hot_PANCREATIC$samples$tumor.subtype, alteration.class = AL1.COMBINED_hot_PANCREATIC$alterations$alteration.type,folder = "COMBINED hot_PANCREATIC", save.intermediate.files = T, verbose = T, remove.0.samples = TRUE)


##prostate
PROSTATE.Genie7.clinical <- Genie7.clinical %>% filter(CANCER_TYPE == "Prostate Cancer")
PROSTATE_F1_clinical <- F1_clinical %>% filter(`cases.primary_site`%in% c("Prostate Gland"))
PROSTATE_pancanTCGA_clinical <- pancanTCGA_clinical %>% filter(`_primary_disease` %in% c("prostate adenocarcinoma"))

#combine samples
T.hotmut.PROSTATE <- TCGA.mut.hot%>% dplyr::select(gene,sample) %>% filter(sample %in% PROSTATE_pancanTCGA_clinical$sample)
F.hotmut.PROSTATE <- F1.mut.hot %>% dplyr::select(gene=Hugo_Symbol,sample = case_id) %>% filter(sample %in% PROSTATE_F1_clinical$case_id)
DM.hotmut.PROSTATE <- DM.hotmut %>% dplyr::select(gene=Hugo_Symbol,sample = Tumor_Sample_Barcode) %>%filter(sample %in% PROSTATE.Genie7.clinical$SAMPLE_ID)
commongene.PROSTATE <- Reduce(intersect, list(unique(T.hotmut.PROSTATE$gene),unique(F.hotmut.PROSTATE$gene),unique(DM.hotmut.PROSTATE$gene)))
combined.hotmut.PROSTATE <- bind_rows(T.hotmut.PROSTATE,F.hotmut.PROSTATE,DM.hotmut.PROSTATE) %>% filter(gene %in%commongene.PROSTATE) #combined filtered sample

#run GENIE
M1.COMBINED_hot_PROSTATE <- table(combined.hotmut.PROSTATE %>% select(sample,gene))
AL1.COMBINED_hot_PROSTATE <- new.AL(M1.COMBINED_hot_PROSTATE)

M2.COMBINED_hot_PROSTATE <- rep("PROSTATE Cancer", nrow(M1.COMBINED_hot_PROSTATE))
names(M2.COMBINED_hot_PROSTATE)<- rownames(M1.COMBINED_hot_PROSTATE)
AL1.COMBINED_hot_PROSTATE$samples$tumor.subtype <- M2.COMBINED_hot_PROSTATE

M3.COMBINED_hot_PROSTATE <- rep("MUTATION", ncol(M1.COMBINED_hot_PROSTATE))
names(M3.COMBINED_hot_PROSTATE) <- colnames(M1.COMBINED_hot_PROSTATE)
AL1.COMBINED_hot_PROSTATE$alterations$alteration.type <- M3.COMBINED_hot_PROSTATE
setwd("~/Sellers Lab Dropbox/Liang Chang/Bioinfo analysis/20191219 Comprehensive Bioinformatics of ME/Results")
result.COMBINED.hot_PROSTATE = select::select(M = AL1.COMBINED_hot_PROSTATE$am, sample.class = AL1.COMBINED_hot_PROSTATE$samples$tumor.subtype, alteration.class = AL1.COMBINED_hot_PROSTATE$alterations$alteration.type,folder = "COMBINED hot_PROSTATE", save.intermediate.files = T, verbose = T, remove.0.samples = TRUE)

#ovarian
OVARIAN.Genie7.clinical <- Genie7.clinical %>% filter(CANCER_TYPE == "Ovarian Cancer")
OVARIAN_F1_clinical <- F1_clinical %>% filter(`cases.primary_site`%in% c("Ovary"))
OVARIAN_pancanTCGA_clinical <- pancanTCGA_clinical %>% filter(`_primary_disease` %in% c("ovarian serous cystadenocarcinoma"))

#combine samples
T.hotmut.OVARIAN <- TCGA.mut.hot%>% dplyr::select(gene,sample) %>% filter(sample %in% OVARIAN_pancanTCGA_clinical$sample)
F.hotmut.OVARIAN <- F1.mut.hot %>% dplyr::select(gene=Hugo_Symbol,sample = case_id) %>% filter(sample %in% OVARIAN_F1_clinical$case_id)
DM.hotmut.OVARIAN <- DM.hotmut %>% dplyr::select(gene=Hugo_Symbol,sample = Tumor_Sample_Barcode) %>%filter(sample %in% OVARIAN.Genie7.clinical$SAMPLE_ID)
commongene.OVARIAN <- Reduce(intersect, list(unique(T.hotmut.OVARIAN$gene),unique(F.hotmut.OVARIAN$gene),unique(DM.hotmut.OVARIAN$gene)))
combined.hotmut.OVARIAN <- bind_rows(T.hotmut.OVARIAN,F.hotmut.OVARIAN,DM.hotmut.OVARIAN) %>% filter(gene %in%commongene.OVARIAN) #combined filtered sample

#run GENIE
M1.COMBINED_hot_OVARIAN <- table(combined.hotmut.OVARIAN %>% select(sample,gene))
AL1.COMBINED_hot_OVARIAN <- new.AL(M1.COMBINED_hot_OVARIAN)

M2.COMBINED_hot_OVARIAN <- rep("OVARIAN Cancer", nrow(M1.COMBINED_hot_OVARIAN))
names(M2.COMBINED_hot_OVARIAN)<- rownames(M1.COMBINED_hot_OVARIAN)
AL1.COMBINED_hot_OVARIAN$samples$tumor.subtype <- M2.COMBINED_hot_OVARIAN

M3.COMBINED_hot_OVARIAN <- rep("MUTATION", ncol(M1.COMBINED_hot_OVARIAN))
names(M3.COMBINED_hot_OVARIAN) <- colnames(M1.COMBINED_hot_OVARIAN)
AL1.COMBINED_hot_OVARIAN$alterations$alteration.type <- M3.COMBINED_hot_OVARIAN
setwd("~/Sellers Lab Dropbox/Liang Chang/Bioinfo analysis/20191219 Comprehensive Bioinformatics of ME/Results")
result.COMBINED.hot_OVARIAN = select::select(M = AL1.COMBINED_hot_OVARIAN$am, sample.class = AL1.COMBINED_hot_OVARIAN$samples$tumor.subtype, alteration.class = AL1.COMBINED_hot_OVARIAN$alterations$alteration.type,folder = "COMBINED hot_OVARIAN", save.intermediate.files = T, verbose = T, remove.0.samples = TRUE)

##bladder
BLADDER.Genie7.clinical <- Genie7.clinical %>% filter(CANCER_TYPE == "Bladder Cancer")
BLADDER_F1_clinical <- F1_clinical %>% filter(`cases.primary_site`%in% c("Bladder"))
BLADDER_pancanTCGA_clinical <- pancanTCGA_clinical %>% filter(`_primary_disease` %in% c("bladder urothelial carcinoma"))

#combine samples
T.hotmut.BLADDER <- TCGA.mut.hot%>% dplyr::select(gene,sample) %>% filter(sample %in% BLADDER_pancanTCGA_clinical$sample)
F.hotmut.BLADDER <- F1.mut.hot %>% dplyr::select(gene=Hugo_Symbol,sample = case_id) %>% filter(sample %in% BLADDER_F1_clinical$case_id)
DM.hotmut.BLADDER <- DM.hotmut %>% dplyr::select(gene=Hugo_Symbol,sample = Tumor_Sample_Barcode) %>%filter(sample %in% BLADDER.Genie7.clinical$SAMPLE_ID)
commongene.BLADDER <- Reduce(intersect, list(unique(T.hotmut.BLADDER$gene),unique(F.hotmut.BLADDER$gene),unique(DM.hotmut.BLADDER$gene)))
combined.hotmut.BLADDER <- bind_rows(T.hotmut.BLADDER,F.hotmut.BLADDER,DM.hotmut.BLADDER) %>% filter(gene %in%commongene.BLADDER) #combined filtered sample

#run GENIE
M1.COMBINED_hot_BLADDER <- table(combined.hotmut.BLADDER %>% select(sample,gene))
AL1.COMBINED_hot_BLADDER <- new.AL(M1.COMBINED_hot_BLADDER)

M2.COMBINED_hot_BLADDER <- rep("BLADDER Cancer", nrow(M1.COMBINED_hot_BLADDER))
names(M2.COMBINED_hot_BLADDER)<- rownames(M1.COMBINED_hot_BLADDER)
AL1.COMBINED_hot_BLADDER$samples$tumor.subtype <- M2.COMBINED_hot_BLADDER

M3.COMBINED_hot_BLADDER <- rep("MUTATION", ncol(M1.COMBINED_hot_BLADDER))
names(M3.COMBINED_hot_BLADDER) <- colnames(M1.COMBINED_hot_BLADDER)
AL1.COMBINED_hot_BLADDER$alterations$alteration.type <- M3.COMBINED_hot_BLADDER
setwd("~/Sellers Lab Dropbox/Liang Chang/Bioinfo analysis/20191219 Comprehensive Bioinformatics of ME/Results")
result.COMBINED.hot_BLADDER = select::select(M = AL1.COMBINED_hot_BLADDER$am, sample.class = AL1.COMBINED_hot_BLADDER$samples$tumor.subtype, alteration.class = AL1.COMBINED_hot_BLADDER$alterations$alteration.type,folder = "COMBINED hot_BLADDER", save.intermediate.files = T, verbose = T, remove.0.samples = TRUE)

##uterine
UTERINE.Genie7.clinical <- Genie7.clinical %>% filter(CANCER_TYPE == "Endometrial Cancer")
UTERINE_F1_clinical <- F1_clinical %>% filter(`cases.primary_site`%in% c("Uterus, NOS"))
UTERINE_pancanTCGA_clinical <- pancanTCGA_clinical %>% filter(`_primary_disease` %in% c("uterine corpus endometrioid carcinoma"))

#combine samples
T.hotmut.UTERINE <- TCGA.mut.hot%>% dplyr::select(gene,sample) %>% filter(sample %in% UTERINE_pancanTCGA_clinical$sample)
F.hotmut.UTERINE <- F1.mut.hot %>% dplyr::select(gene=Hugo_Symbol,sample = case_id) %>% filter(sample %in% UTERINE_F1_clinical$case_id)
DM.hotmut.UTERINE <- DM.hotmut %>% dplyr::select(gene=Hugo_Symbol,sample = Tumor_Sample_Barcode) %>%filter(sample %in% UTERINE.Genie7.clinical$SAMPLE_ID)
commongene.UTERINE <- Reduce(intersect, list(unique(T.hotmut.UTERINE$gene),unique(F.hotmut.UTERINE$gene),unique(DM.hotmut.UTERINE$gene)))
combined.hotmut.UTERINE <- bind_rows(T.hotmut.UTERINE,F.hotmut.UTERINE,DM.hotmut.UTERINE) %>% filter(gene %in%commongene.UTERINE) #combined filtered sample

#run GENIE
M1.COMBINED_hot_UTERINE <- table(combined.hotmut.UTERINE %>% select(sample,gene))
AL1.COMBINED_hot_UTERINE <- new.AL(M1.COMBINED_hot_UTERINE)

M2.COMBINED_hot_UTERINE <- rep("UTERINE Cancer", nrow(M1.COMBINED_hot_UTERINE))
names(M2.COMBINED_hot_UTERINE)<- rownames(M1.COMBINED_hot_UTERINE)
AL1.COMBINED_hot_UTERINE$samples$tumor.subtype <- M2.COMBINED_hot_UTERINE

M3.COMBINED_hot_UTERINE <- rep("MUTATION", ncol(M1.COMBINED_hot_UTERINE))
names(M3.COMBINED_hot_UTERINE) <- colnames(M1.COMBINED_hot_UTERINE)
AL1.COMBINED_hot_UTERINE$alterations$alteration.type <- M3.COMBINED_hot_UTERINE
setwd("~/Sellers Lab Dropbox/Liang Chang/Bioinfo analysis/20191219 Comprehensive Bioinformatics of ME/Results")
result.COMBINED.hot_UTERINE = select::select(M = AL1.COMBINED_hot_UTERINE$am, sample.class = AL1.COMBINED_hot_UTERINE$samples$tumor.subtype, alteration.class = AL1.COMBINED_hot_UTERINE$alterations$alteration.type,folder = "COMBINED hot_UTERINE", save.intermediate.files = T, verbose = T, remove.0.samples = TRUE)

##Gastric

`%notin%` <- Negate(`%in%`)
STOMACHE.Genie7.clinical <- Genie7.clinical %>% filter(CANCER_TYPE == "Esophagogastric Cancer") %>% filter(ONCOTREE_CODE %notin% c('ESCA', 'GEJ', 'EGC', 'ESCC', 'EPDCA'))
STOMACHE_F1_clinical <- F1_clinical %>% filter(`cases.primary_site`%in% c("Stomach"))
STOMACHE_pancanTCGA_clinical <- pancanTCGA_clinical %>% filter(`_primary_disease` %in% c("stomach adenocarcinoma"))

#combine samples
T.hotmut.STOMACHE <- TCGA.mut.hot%>% dplyr::select(gene,sample) %>% filter(sample %in% STOMACHE_pancanTCGA_clinical$sample)
F.hotmut.STOMACHE <- F1.mut.hot %>% dplyr::select(gene=Hugo_Symbol,sample = case_id) %>% filter(sample %in% STOMACHE_F1_clinical$case_id)
DM.hotmut.STOMACHE <- DM.hotmut %>% dplyr::select(gene=Hugo_Symbol,sample = Tumor_Sample_Barcode) %>%filter(sample %in% STOMACHE.Genie7.clinical$SAMPLE_ID)
commongene.STOMACHE <- Reduce(intersect, list(unique(T.hotmut.STOMACHE$gene),unique(F.hotmut.STOMACHE$gene),unique(DM.hotmut.STOMACHE$gene)))
combined.hotmut.STOMACHE <- bind_rows(T.hotmut.STOMACHE,F.hotmut.STOMACHE,DM.hotmut.STOMACHE) %>% filter(gene %in%commongene.STOMACHE) #combined filtered sample

#run GENIE
M1.COMBINED_hot_STOMACHE <- table(combined.hotmut.STOMACHE %>% select(sample,gene))
AL1.COMBINED_hot_STOMACHE <- new.AL(M1.COMBINED_hot_STOMACHE)

M2.COMBINED_hot_STOMACHE <- rep("STOMACHE Cancer", nrow(M1.COMBINED_hot_STOMACHE))
names(M2.COMBINED_hot_STOMACHE)<- rownames(M1.COMBINED_hot_STOMACHE)
AL1.COMBINED_hot_STOMACHE$samples$tumor.subtype <- M2.COMBINED_hot_STOMACHE

M3.COMBINED_hot_STOMACHE <- rep("MUTATION", ncol(M1.COMBINED_hot_STOMACHE))
names(M3.COMBINED_hot_STOMACHE) <- colnames(M1.COMBINED_hot_STOMACHE)
AL1.COMBINED_hot_STOMACHE$alterations$alteration.type <- M3.COMBINED_hot_STOMACHE
setwd("~/Sellers Lab Dropbox/Liang Chang/Bioinfo analysis/20191219 Comprehensive Bioinformatics of ME/Results")
result.COMBINED.hot_STOMACHE = select::select(M = AL1.COMBINED_hot_STOMACHE$am, sample.class = AL1.COMBINED_hot_STOMACHE$samples$tumor.subtype, alteration.class = AL1.COMBINED_hot_STOMACHE$alterations$alteration.type,folder = "COMBINED hot_STOMACHE", save.intermediate.files = T, verbose = T, remove.0.samples = TRUE)

##esophageal
ESOPHAGEAL.Genie7.clinical <- Genie7.clinical %>% filter(CANCER_TYPE == "Esophagogastric Cancer") %>% filter(ONCOTREE_CODE %in% c('ESCA', 'GEJ', 'EGC', 'ESCC', 'EPDCA'))
ESOPHAGEAL_F1_clinical <- F1_clinical %>% filter(`cases.primary_site`%in% c("Stomach"))
ESOPHAGEAL_pancanTCGA_clinical <- pancanTCGA_clinical %>% filter(`_primary_disease` %in% c("stomach adenocarcinoma"))

#combine samples
T.hotmut.ESOPHAGEAL <- TCGA.mut.hot%>% dplyr::select(gene,sample) %>% filter(sample %in% ESOPHAGEAL_pancanTCGA_clinical$sample)
F.hotmut.ESOPHAGEAL <- F1.mut.hot %>% dplyr::select(gene=Hugo_Symbol,sample = case_id) %>% filter(sample %in% ESOPHAGEAL_F1_clinical$case_id)
DM.hotmut.ESOPHAGEAL <- DM.hotmut %>% dplyr::select(gene=Hugo_Symbol,sample = Tumor_Sample_Barcode) %>%filter(sample %in% ESOPHAGEAL.Genie7.clinical$SAMPLE_ID)
commongene.ESOPHAGEAL <- Reduce(intersect, list(unique(T.hotmut.ESOPHAGEAL$gene),unique(F.hotmut.ESOPHAGEAL$gene),unique(DM.hotmut.ESOPHAGEAL$gene)))
combined.hotmut.ESOPHAGEAL <- bind_rows(T.hotmut.ESOPHAGEAL,F.hotmut.ESOPHAGEAL,DM.hotmut.ESOPHAGEAL) %>% filter(gene %in%commongene.ESOPHAGEAL) #combined filtered sample

#run GENIE
M1.COMBINED_hot_ESOPHAGEAL <- table(combined.hotmut.ESOPHAGEAL %>% select(sample,gene))
AL1.COMBINED_hot_ESOPHAGEAL <- new.AL(M1.COMBINED_hot_ESOPHAGEAL)

M2.COMBINED_hot_ESOPHAGEAL <- rep("ESOPHAGEAL Cancer", nrow(M1.COMBINED_hot_ESOPHAGEAL))
names(M2.COMBINED_hot_ESOPHAGEAL)<- rownames(M1.COMBINED_hot_ESOPHAGEAL)
AL1.COMBINED_hot_ESOPHAGEAL$samples$tumor.subtype <- M2.COMBINED_hot_ESOPHAGEAL

M3.COMBINED_hot_ESOPHAGEAL <- rep("MUTATION", ncol(M1.COMBINED_hot_ESOPHAGEAL))
names(M3.COMBINED_hot_ESOPHAGEAL) <- colnames(M1.COMBINED_hot_ESOPHAGEAL)
AL1.COMBINED_hot_ESOPHAGEAL$alterations$alteration.type <- M3.COMBINED_hot_ESOPHAGEAL
setwd("~/Sellers Lab Dropbox/Liang Chang/Bioinfo analysis/20191219 Comprehensive Bioinformatics of ME/Results")
result.COMBINED.hot_ESOPHAGEAL = select::select(M = AL1.COMBINED_hot_ESOPHAGEAL$am, sample.class = AL1.COMBINED_hot_ESOPHAGEAL$samples$tumor.subtype, alteration.class = AL1.COMBINED_hot_ESOPHAGEAL$alterations$alteration.type,folder = "COMBINED hot_ESOPHAGEAL", save.intermediate.files = T, verbose = T, remove.0.samples = TRUE)



###Load files after SELECT Analysis (from previous section)
group.colors <- c("Different pathways" = "gray23","Wnt/B-catenin signaling" = "blue", "MAPK signaling" = "red","PI3K signaling" = "green4")

filenames <- list.files(path=getwd(),full.names = TRUE,pattern = "\\.RData$")
filenames_short <- list.files(path=getwd())
result <- list()
for (i in 1:length(filenames)){
  load(filenames[i])
  result[[i]] <- alpi
}

for (i in 1:12){
  result[[i]] <- result[[i]] %>% add_column(type = rep(str_sub(filenames_short[i],6,-10),nrow(result[[i]]))) 
}

all_hotspot <- data.frame()
for (i in 1:12){
  all_hotspot <- bind_rows(all_hotspot,result[[i]]) 
}


hotspot_summary_noFDRfilter <- all_hotspot %>% mutate(support1_coverage = support_1/((support_1+support_2-overlap)/freq_coverage),support2_coverage = support_2/((support_1+support_2-overlap)/freq_coverage)) %>% filter(support2_coverage >0.03 & support1_coverage>0.03) %>% mutate(fold_difference = (r_overlap+1)/(overlap+1)) %>% mutate(label = paste(name,type,sep = " ")) %>% mutate (case_coverage = (support_1+support_2-overlap)) %>% mutate(cohort_freq = (support_1+support_2-overlap)/freq_coverage)

pathway_info <- read_csv("Final_ORFlist_Bioinfo.csv") %>% dplyr::select(gene, consensus_symbol) %>% add_row(gene ="PTEN",consensus_symbol="PI3K signaling") %>% add_row(gene ="APC",consensus_symbol="Wnt/B-catenin signaling") %>% add_row(gene ="PIK3R2",consensus_symbol="PI3K signaling") 

pathway_1_nofilter <- pathway_info$consensus_symbol[match(hotspot_summary_noFDRfilter$SFE_1,pathway_info$gene)]%>% recode('RTK signaling' = 'MAPK signaling')
pathway_2_nofilter <- pathway_info$consensus_symbol[match(hotspot_summary_noFDRfilter$SFE_2,pathway_info$gene)]%>% recode('RTK signaling' = 'MAPK signaling')

same_pathway_nofilter <- (pathway_1_nofilter == pathway_2_nofilter) %>% replace_na(FALSE)

what_p_nofilter <- c()
for(i in 1:length(same_pathway_nofilter)){
  if (same_pathway_nofilter [i] == TRUE){
    what_p_nofilter[i] = pathway_1_nofilter[i]
  }
  else what_p_nofilter[i] = "Different pathways"
}

hotspot_summary_ME_noFDRfilter <- hotspot_summary_noFDRfilter %>% 
  add_column(pathway_1=pathway_1_nofilter,pathway_2=pathway_2_nofilter,
             same_pathway =same_pathway_nofilter,Gene_Pair_Same_Pathway=what_p_nofilter)%>% 
  filter(direction == "ME",type %ni% c("Pancreatic","Uterine")) %>% 
  mutate(Gene_Pair_Same_Pathway = ifelse((Gene_Pair_Same_Pathway != "Different pathways"&
                                            FDR== TRUE&APC_good == TRUE)|(fold_difference>6&
                                                                            FDR== TRUE&APC_good == TRUE),Gene_Pair_Same_Pathway,"")) %>%
  mutate(LabelName = ifelse(Gene_Pair_Same_Pathway!="",name,""))


hotspot_summary_CO_noFDRfilter <- hotspot_summary_noFDRfilter %>% 
  add_column(pathway_1=pathway_1_nofilter,pathway_2=pathway_2_nofilter,
             same_pathway =same_pathway_nofilter,Gene_Pair_Same_Pathway=what_p_nofilter)%>% 
  filter(direction == "CO",type %ni% c("Pancreatic","Uterine")) 
#final dataset for visualization



#4 visualization
group.colors <- c("Different pathways" = "gray23","Wnt/B-catenin signaling" = "blue", "MAPK signaling" = "red","PI3K signaling" = "green4")

pos <- position_jitter(width = 0.7, seed = 50)

Fig1a <- ggplot(data = hotspot_summary_ME_noFDRfilter, aes(y = fold_difference, x = 1)) +
  facet_wrap(. ~factor(type, levels=c('Bladder','Brain','Breast','Colorectal','Esophageal','Gastric','Lung',
                                      'Prostate','Melanoma')), nrow = 1) + 
  geom_jitter(data = (hotspot_summary_ME_noFDRfilter %>% filter(FDR!= TRUE, APC_good != TRUE)),
              alpha = 0.1,color = "grey",size = 4,position = pos)+
  geom_jitter(data = hotspot_summary_ME_noFDRfilter %>% filter(FDR== TRUE, APC_good == TRUE),
              alpha = 1,aes(color = Gene_Pair_Same_Pathway),size = 4,position = pos)+
  scale_color_manual(values=group.colors)+
  xlab("Cancer_type")+
  ylab("Co-occurance Fold Differences (Expected/Observed)")+
  theme_classic()+
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        legend.key.size = unit(2, 'cm'),
        legend.title = element_blank(),
        legend.text = element_blank(),
        plot.background = element_blank(), 
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA,size = 2),
        axis.text = element_text(face="bold"),
        axis.text.y = element_text(size = 19),
        axis.title.y = element_text(size = 19, colour = "black"),
        text = element_text(size=13,face="bold"),
        strip.text.x = element_text(size = 16, colour = "black", angle = 0))

ggsave(Fig1a, filename = paste0("Fig1a_new4_SELECT_Effect_",Sys.Date(),".pdf"), width = 15, height = 7)



#Figure 1d 
add_cor_labels = FALSE 
day11_orf_metadata <- read.csv("Fig1e_Biorep_metadata.csv",
                               header = FALSE)
day11_tech_rep <- readRDS("Fig1e_Biorep_rawdata.rds") %>%
  dplyr::select(CCLE_name,starts_with("11_")) %>% column_to_rownames(var = "CCLE_name") %>%
  select(day11_orf_metadata$V1)
colnames(day11_tech_rep) [match(colnames(day11_tech_rep) ,day11_orf_metadata$V1)] <- day11_orf_metadata$V2

cormat <- round(cor(day11_tech_rep),2)
melted_cormat <- melt(cormat)

library(RColorBrewer)

names_ann <- c("","BRAF_V600E","","","p16","" ,"","KRAS_G12V","" ,"","GFP","","","p53","","","NIC","","","MEK1-DD","","" ,"EGFR_exon19del" ,"",
               "","ERK2","","" ,"AKT_E17K","","","CTNNB1","","" ,"PTEN", "")

pdf(paste0('Fig1d_heatmap.pdf'), width = 13,height = 13)
# col = colorRamp2(c(0, 1), c("white","#A63446"))

h = Heatmap(cormat,column_names_gp = grid::gpar(fontsize = 20),
            row_names_gp = grid::gpar(fontsize = 20),show_column_names = T,
            show_row_names = T,column_names_rot = 45,
            row_order = colnames(cormat), 
            column_order = colnames(cormat),column_labels = names_ann,row_labels = names_ann,
            col = colorRamp2(c(0, 0.5, 1), c("white","orange","#A63446")),name = "Pearson\nCorrelation",
            column_names_centered =  T,
            heatmap_legend_param = list(
              legend_height = unit(6, "cm"),
              grid_width = unit(1, "cm"),
              border = "black",
              title = "Pearson\nCorrelation",title_gp = gpar(fontsize = 14,fontface = "bold"),
              labels_gp = gpar(fontsize = 17)))
draw(h, heatmap_legend_side = "left")
dev.off()



