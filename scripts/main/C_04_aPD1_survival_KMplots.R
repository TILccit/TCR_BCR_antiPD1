#! /usr/bin/env Rscript
## Clear R-workspace
# rm(list=ls(all=TRUE))

## Close all graphic devices
graphics.off()
#####################
### Load packages ###
#####################
library(pacman)

pacman::p_load(DT,dplyr,tibble,tidyverse,ggplot2,hrbrthemes,
               ggpubr, reshape2,cowplot,ggrepel,cowplot,grid,
               RColorBrewer, pander,gridExtra,
               reconPlots,flextable, extrafont,Amelia, estimate,strex,plyr,survival,survminer, data.table)

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
# tcgaIntermediateData <- paste0(projectOutDataPath,"TCGA/")
antiPDL1publicIntermediateData <- paste0(projectOutDataPath,"antiPD(L)1_public/")

# Session/dependencies
sessionInfoPath <- paste0("session_info/")

######################################
## Create intermediate output paths ##
######################################
if (!dir.exists(paste0(projectOutDataPath,"survival/"))) {dir.create(paste0(projectOutDataPath,"survival/"))}


##################################
### LOAD SOURCE FUNCTIONS FILE ###
##################################
source(paste0(scriptsFunctionsPath,"aPD1_survival_KMplots_functions.R"), local=T)

#####################
### File suffixes ###
#####################
Rdata.suffix <- ".RData"

######################
### Model variables ###
######################

modelTypeA <- "_TCRArichnessIGH.ds"
modelTypeB <- "_noLIU"
modelTypeC <- "_mintcr4bcr10"

# Load data
immdata.TcrBcr.full.red <- loadRData(paste0(antiPDL1publicIntermediateData,"immdata.TcrBcr.full.red",modelTypeA,modelTypeB,modelTypeC,".RData"))
immdata.TcrBcr.full.red$response <- ifelse(immdata.TcrBcr.full.red$response=="CR+PR","CRPR",immdata.TcrBcr.full.red$response)

## EXTRACT RESPONSE, NOMINAL, MEASURE VARIABLES FROM DATA
response.var <- "response"
nominal.vars <- c("dataset","tissue","gender")
measure.vars <-colnames(immdata.TcrBcr.full.red)[c(2:10,12,32:43)]
measure.vars <- measure.vars[c(6,22,1:3,20,21,4:5,10,7:9,11:19)]


## CHECK INFINITE LOG TMB VALUES
immdata.TcrBcr.full.red[which(is.infinite(immdata.TcrBcr.full.red$logTMB)),]

## Extract tissues in data so far
tissues <- unique(immdata.TcrBcr.full.red$tissue)

# Which are the covariates I am considering
covariatesIn <- measure.vars[c(1:2,4:7,14:16)]
responseVar <- "response"

mycancertypes.final <-unique(immdata.TcrBcr.full.red$tissue)
mycancertypes.final <- mycancertypes.final[c(2,4,1,3)]

# Extract features
colnames(immdata.TcrBcr.full.red)[c(7,43,32,4)] <- c("TCR_Richness","BCR_Richness","B_cells","PD_L1")
featuresIn <- colnames(immdata.TcrBcr.full.red)[c(4,42,43,7,32,33,34,3,41)]
featuresInLabels <- c("PD-L1 expression","logTMB","BCR Richness","TCR Richness","B cells infiltration","CD8+ T cells infiltration","CD4+ T cells infiltration",
                      "T-cell inflamed GEP", "Tumor Purity")

tissues<- unique(immdata.TcrBcr.full.red$tissue)
tissues <- tissues[c(2,4,1,3)]

immdata.TcrBcr.full.red$tissue <- ifelse(immdata.TcrBcr.full.red$tissue=="melanoma", "SKCM anti-PD1/L1",
                                         ifelse(immdata.TcrBcr.full.red$tissue=="RCC","KIRC anti-PD1/L1",
                                                ifelse(immdata.TcrBcr.full.red$tissue=="bladder","BLCA anti-PD1/L1",                                             ifelse(immdata.TcrBcr.full.red$tissue=="gastric","STAD anti-PD1/L1",immdata.TcrBcr.full.red$tissue))))

immdata.TcrBcr.full.red$tissue <- factor(immdata.TcrBcr.full.red$tissue, levels=c("SKCM anti-PD1/L1","KIRC anti-PD1/L1","BLCA anti-PD1/L1","STAD anti-PD1/L1"))
tissues<- levels(immdata.TcrBcr.full.red$tissue)

aPD1.SURV<-sapply(1:length(featuresIn), function(x) runMultiSurvival_aPD1(DF =immdata.TcrBcr.full.red,var = featuresIn[x],varSelection =tissues,varSelected = "tissue",survType = "OS", sampleID = "run_accession",varPlotName = featuresInLabels[x], 
                                                                          printP=TRUE,saveP = TRUE,savePath=paste0(projectOutDataPath,"survival/"),onlySign=FALSE,selectFacets=tissues), simplify = FALSE)

aPD1.SURV.df<-rbindlist(aPD1.SURV)

####################
### Session info ###
####################

session_info <- sessionInfo()
writeLines(capture.output(session_info), paste0(sessionInfoPath,"C_04_aPD1_survival_KMplots.txt"))

sessionInfo()