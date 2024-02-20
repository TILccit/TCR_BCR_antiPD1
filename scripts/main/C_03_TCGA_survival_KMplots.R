#! /usr/bin/env Rscript
## Clear R-workspace
# rm(list=ls(all=TRUE))

## Close all graphic devices
graphics.off()
#####################
### Load packages ###
#####################
library(pacman)

pacman::p_load(extrafont,DT,dplyr,tibble,survminer,survival, ggplot2, hrbrthemes, grid, gridExtra, cowplot,Amelia,data.table)

extrafont::loadfonts(device="win")
windowsFonts(sans="Palatino Linotype")
loadfonts(device="win")
loadfonts(device="postscript")


###########################
### Main analysis paths ###
###########################
# Main
scriptsPath <- paste0("scripts/")
scriptsMainPath<- paste0(scriptsPath, "main/")
scriptsFunctionsPath <- paste0(scriptsPath,"functions/")
projectDataPath <- paste0("data/")
# Input
tcgaInputData <- paste0(projectDataPath,"TCGA/")
antiPDL1publicInputData <- paste0(projectDataPath,"antiPD(L)1_public/")
otherInputData <- paste0(projectDataPath,"other/")
referenceInputData <- paste0(projectDataPath,"reference/")
# Output
projectOutDataPath <- paste0("output/data_files/")
tcgaIntermediateData <- paste0(projectOutDataPath,"TCGA/")
# antiPDL1publicIntermediateData <- paste0(projectOutDataPath,"antiPD(L)1_public/")

# Session/dependencies
sessionInfoPath <- paste0("session_info/")

######################################
## Create intermediate output paths ##
######################################
if (!dir.exists(paste0(projectOutDataPath,"survival/"))) {dir.create(paste0(projectOutDataPath,"survival/"))}


##################################
### LOAD SOURCE FUNCTIONS FILE ###
##################################
source(paste0(scriptsFunctionsPath,"TCGA_survival_KMplots_functions.R"), local = T)


#####################
### File suffixes ###
#####################
Rdata.suffix <- ".RData"

# Data loading
tcgaData.tp <- loadRData(paste0(tcgaIntermediateData, "tcga.tp.Filt.final",Rdata.suffix)) %>% as.data.frame()
# tcga.drug.notICI <-loadRData(paste0(tcgaIntermediateData, "tcga.tp.nonFilt.final_ICI4Tissues.RData"))
tcga.drug.ici.To.Exclude2 <- loadRData(paste0(tcgaIntermediateData, "tcga_patientsNOTici_toEXCLUDE_4Tissues.RData"))

# Filter data to exclude immunotherapy
table(tcgaData.tp$project)
# table(tcga.drug.notICI$project)

# See which patients we need to remove
tcgaData.tp %>% dplyr::filter(bcr_patient_barcode %in% tcga.drug.ici.To.Exclude2) %>% droplevels() %>% distinct()
# Remove patients that received aPD1 immunotherapy.
tcgaData.tp.clean <- tcgaData.tp %>% dplyr::filter(bcr_patient_barcode %notin% tcga.drug.ici.To.Exclude2) %>% 
  droplevels() %>% 
  distinct()
dim(tcgaData.tp)
dim(tcgaData.tp.clean)

# Define features
featuresSelected <- c("bcr_patient_barcode", "project","PD_L1","logTMB","BCR_Richness","TCR_Richness",
                      "B_cells.CIBER", "T_Cells_CD8.CIBER", "CD4cells.CIBER","TIS_sigscore", "TumorPurity","OS","OS.time","PFI","PFI.time","age_at_initial_pathologic_diagnosis","gender")

extraFeatures <- c("StromalScore", "ImmuneScore")

featuresIn <- c("TCR_Richness","BCR_Richness","logTMB","PD_L1","TIS_sigscore",
                "CD8cells", "CD4cells","B_cells",  "TumorPurity")
featuresInLabels <- c("TCR Richness","BCR Richness","logTMB","PD-L1 expression","T-cell inflamed GEP", "CD8+ T cells infiltration","CD4+ T cells infiltration","B cells infiltration",
                      "Tumor Purity")
## FOR TP data
tcgaData.tp.sub <- tcgaData.tp.clean %>% dplyr::select(all_of(featuresSelected)) %>% droplevels() %>% distinct()
tcgaData.tp.sub.extra <- tcgaData.tp.clean %>% dplyr::select(all_of(c(featuresSelected,extraFeatures))) %>% droplevels() %>% distinct()

colnames(tcgaData.tp.sub)[7:9] <- gsub(".CIBER","",colnames(tcgaData.tp.sub)[7:9])
colnames(tcgaData.tp.sub)[8] <- "CD8cells"

# Keep rows where at least on of the features of interested is not NA
tcgaData.tp.sub <-tcgaData.tp.sub[rowSums(is.na(tcgaData.tp.sub[,3:11])) != 9,]
tcgaData.tp.sub.extra <-tcgaData.tp.sub.extra[rowSums(is.na(tcgaData.tp.sub.extra[,3:11])) != 9,]
dim(tcgaData.tp.sub)
dim(tcgaData.tp.sub.extra)

# Let's summarize the TCGA biospecimen data we just loaded
# Number of patients per project
tcgaData.tp.patient.sub.table <-tcgaData.tp.sub %>% dplyr::select(project,bcr_patient_barcode) %>% group_by(project) %>% distinct() %>% dplyr::summarize(n())


# Define cohorts
mycancertypes <- levels(as.factor(tcgaData.tp.sub$project))
# mycancertypes <- factor(mycancertypes,levels = mycancertypes)
## Exclude THYM and DLBC
mycancertypes.reduced <- mycancertypes[-which(mycancertypes %in% c("THYM", "DLBC","COAD_MSI_low","BRCA_Her2","BRCA_LuminalA", "BRCA_LuminalB","BRCA","BRCA_Normal","COAD","UCEC","LAML"))]

## TCGA histologies aligned with aPD1 data
mycancertypes.final <- c("SKCM", "KIRC","BLCA","STAD")
# mycancertypes.final <- factor(mycancertypes.final, levels = mycancertypes.final)
## Subset data to reduced types
tcgaData.tp.sub.filt <- tcgaData.tp.sub %>% dplyr::filter(project %in% mycancertypes.reduced) %>% droplevels() %>% as.data.frame()
tcgaData.tp.sub.extra.filt<- tcgaData.tp.sub.extra %>% dplyr::filter(project %in% mycancertypes.reduced) %>% droplevels() %>% as.data.frame()

# Make project a factor with levels
tcgaData.tp.sub.filt$project <- factor(tcgaData.tp.sub.filt$project, levels=mycancertypes.reduced)
tcgaData.tp.sub.extra.filt$project <- factor(tcgaData.tp.sub.extra.filt$project, levels=mycancertypes.reduced)

#CHECK DUPLICATED PATIENTS
which(duplicated(tcgaData.tp.sub.filt$bcr_patient_barcode))
tcgaData.tp.sub.filt[tcgaData.tp.sub.filt$bcr_patient_barcode==tcgaData.tp.sub.filt$bcr_patient_barcode[which(duplicated(tcgaData.tp.sub.filt$bcr_patient_barcode))],]

# OS units
tcgaData.tp.sub.filt$OS.time <-tcgaData.tp.sub.filt$OS.time/(365.24/12)

# Run univariate analyses
tcga.SURV<-sapply(1:length(featuresIn), function(x) runMultiSurvival_TCGA(DF = tcgaData.tp.sub.filt,var = featuresIn[x],varSelection =mycancertypes.final,varSelected = "project",survType = "OS", sampleID = "bcr_patient_barcode",varPlotName = featuresInLabels[x], printP=TRUE,saveP = TRUE,savePath =paste0(projectOutDataPath,"survival/") ,onlySign=FALSE,selectFacets=mycancertypes.final), simplify = FALSE)

tcga.SURV.df<-rbindlist(tcga.SURV)
####################
### Session info ###
####################

session_info <- sessionInfo()
writeLines(capture.output(session_info), paste0(sessionInfoPath,"C_03_TCGA_survival_KMplots.txt"))

sessionInfo()