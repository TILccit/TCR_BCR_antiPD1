#!/services/tools/R/4.1.0-GCC-MKL/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied: model algorithm", call.=FALSE)
} 
## HOW TO RUN SCRIPT FROM COMMAND LINE IN COMPUTEROME #####
# module load intel/perflibs/2020_update4 gcc/11.1.0 R/4.1.0
# Rscript ./E_01_aPD1_RFE_multivariate_Scheme_1.R NB

# pkgFile <- "https://cran.r-project.org/src/contrib/Archive/OutlierDetection/OutlierDetection_0.1.1.tar.gz" 
#https://cran.r-project.org/src/contrib/Archive/OutlierDetection/OutlierDetection_0.1.1.tar.gz
# devtools::install_github('jhmadsen/DDoutlier')
# install.packages(pkgs="https://cran.r-project.org/src/contrib/circular_0.4-93.tar.gz", type="source", repos=NULL)
# install.packages(pkgs="https://cran.r-project.org/src/contrib/depth_2.1-1.1.tar.gz", type="source", repos=NULL)
# install.packages(pkgs="https://cran.r-project.org/src/contrib/depthTools_0.7.tar.gz", type="source", repos=NULL)
# install.packages(pkgs="https://cran.r-project.org/src/contrib/ldbod_0.1.2.tar.gz", type="source", repos=NULL)
# install.packages(pkgs=pkgFile, type="source", repos=NULL)
# remotes::install_github("rubak/spatstat.revdep", subdir = "OutlierDetection")

# library(devtools)
# devtools::install_github("dalpozz/unbalanced")
# require("devtools")
# install_github("tidymodels/themis")
# install.packages("https://cran.r-project.org/src/contrib/Archive/DMwR/DMwR_0.1.0.tar.gz", repos = NULL, type="source")

# library(remotes)
# install_version("OutlierDetection", "0.1.1")
# packageurl <- "https://cran.r-project.org/src/contrib/Archive/caret/caret_6.0-89.tar.gz"
# # packageurl <- "http://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_0.9.1.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")
# install_version("caret", '6.0.88')
# install.packages("caret", version='6.0.88')
#####################
### Load packages ####
#####################
library(pacman)
pacman::p_load(DT,stringr,dplyr,plyr,tibble,tidyverse,data.table, strex,Hmisc,
               ggpubr,mlbench,Amelia,caret, DMwR, DataExplorer,FactoMineR,factoextra, parallel,doParallel,scales,precrec, rlist,RColorBrewer
                )

###########################
### Main analysis paths ####
###########################
projectDir <- 'D:/BIOINF/PROJECTS/2_ANALYSES/18_GIT_PAPERS/TCR_BCR_antiPD1/'
setwd(projectDir)
# Main
scriptsPath <- paste0("scripts/")
scriptsFunctionsPath <- paste0(scriptsPath,"functions/")
projectDataPath <- paste0("data/")
# Input
# tcgaInputData <- paste0(projectDataPath,"TCGA/")
# antiPDL1publicInputData <- paste0(projectDataPath,"antiPD(L)1_public/")
# otherInputData <- paste0(projectDataPath,"other/")
# referenceInputData <- paste0(projectDataPath,"reference/")
# Output - Intermediate output
projectOutDataPath <- paste0("output/data_files/")
tcgaIntermediateData <- paste0(projectOutDataPath,"TCGA/")
antiPDL1publicIntermediateData <- paste0(projectOutDataPath,"antiPD(L)1_public/")
#
figuresPath <- paste0("output/figures/")
tablesPath <- paste0("output/tables/")
modelsPath <- paste0("output/models/")

if (!dir.exists(paste0(projectDir,modelsPath,"paramOptim/"))) {dir.create(paste0(projectDir,modelsPath,"paramOptim/"))}

##########################
### FOLDER/FILES NAMES ####
##########################
# Create list with folders names for the different datasets.
datasets <- c("Hugo","Riaz","Gide","Liu","Rod","Miao","Immotion150","Kim","Mari", "Powles","Braun")
data.foldersList <- as.list(c("Hugo_MEL","Riaz_MEL","Gide_MEL", "Liu_MEL","Rodriguez_MEL","Miao_RCC","Immo150_RCC","Kim_GAS","Mari_BLA","Powles_BLA","Braun_RCC"))
names(data.foldersList)<- c("Hugo","Riaz","Gide", "Liu","Rod", "Miao","Immotion150","Kim","Mari","Powles","Braun")

##################################
### LOAD SOURCE FUNCTIONS FILE ####
##################################
source(paste0(scriptsFunctionsPath,"aPD1_RFE_multivariate_functions.R"))

#####################
### File suffixes ###
#####################
Rdata.suffix <- ".RData"

######################
# Model variables ####
######################

modelTypeA <- "_TCRArichnessIGH.ds"
modelTypeB <- "_noLIU"
modelTypeC <- "_mintcr4bcr10"

#################
### Load data ####
#################

if (file.exists(paste0(antiPDL1publicIntermediateData,"samples.aPD1.f.tpm.dataset","_ALLresp", Rdata.suffix))){
  samples_datasets <-loadRData(paste0(antiPDL1publicIntermediateData,"samples.aPD1.f.tpm.dataset","_ALLresp", Rdata.suffix))
  samples_tissues <- loadRData(paste0(antiPDL1publicIntermediateData,"samples.aPD1.tpm.batch.TISSUE","_ALLresp",Rdata.suffix))
  samples_all.merged <- bind_rows(samples_tissues)
}else{
  print("Run script B_01_aPD1_data_count_deconvolution_processing.R to create table with sample metadata")
}


if (file.exists(paste0(antiPDL1publicIntermediateData,"samples.aPD1.estimate.CIBER.ABS",Rdata.suffix)) && file.exists(paste0(antiPDL1publicIntermediateData,"samples.aPD1.estimate.CIBER.REL",Rdata.suffix))){
  # Load also metadata with estimate and cibersort
  samples_meta.Abs <- loadRData(paste0(antiPDL1publicIntermediateData,"samples.aPD1.estimate.CIBER.ABS",Rdata.suffix))
  samples_meta.Rel <- loadRData(paste0(antiPDL1publicIntermediateData,"samples.aPD1.estimate.CIBER.REL",Rdata.suffix))
  
  samples_meta.Abs.compact <- loadRData(paste0(antiPDL1publicIntermediateData,"samples.aPD1.estimate.CIBER.ABS.compact",Rdata.suffix))
}else{
  print("Run script B_01_aPD1_data_count_deconvolution_processing.RMD to get CIBERSORT results")
}


# Fix all samples metadata tables
samples_all.merged.fix <-fix_immoRunAcc(samples_all.merged,datasets[1:10],7)

samples_meta.Abs.fix <- fix_immoRunAcc(samples_meta.Abs,datasets[1:10],7)
samples_meta.Rel.fix <- fix_immoRunAcc(samples_meta.Rel,datasets[1:10],7)
samples_meta.Abs.compact.fix <- fix_immoRunAcc(samples_meta.Abs.compact,datasets[1:10],7)

if (file.exists(paste0(antiPDL1publicIntermediateData,"tcrAclones_downSampled.divEst.aPD1.FINALmin4.v1",Rdata.suffix))){
  # Load data
  immdata.TcrA <- loadRData(paste0(antiPDL1publicIntermediateData,"tcrAclones_downSampled.divEst.aPD1.FINALmin4.v1",Rdata.suffix))
  # TCRB
  immdata.TcrB <- loadRData(paste0(antiPDL1publicIntermediateData,"tcrBclones_downSampled.divEst.aPD1.FINALmin4.v1",Rdata.suffix))

  # BCR-IGH
  immdata.IgH <- loadRData(paste0(antiPDL1publicIntermediateData,"igHclones_downSampled.divEst.aPD1.FINALmin10.v1",Rdata.suffix))
}else{
  print("Run script B_02_aPD1_TCR_BCR_clones_analysis.RMD to process mixcr data and get diversiy measure tables")
}


###############################
### Processing data tables ####
###############################
# FOR EACH DATASET The diversityEst object is a list, so we need to unlist that and have it as a dataframe, then
# we need to merge that with the metadata object which is already a dataframe.
# After that we will merge all the datasets, concatenate their dataframes.

# TCRA
immdata.TcrA.clean <-sapply(datasets[1:10], function(x) extractReduce_diversityDF(x, immdata.TcrA[[x]]), simplify = FALSE)

# TCRB
immdata.TcrB.clean <-sapply(datasets[1:10], function(x) extractReduce_diversityDF(x, immdata.TcrB[[x]]), simplify = FALSE)

# IGH
# TCRB
immdata.IgH.clean <-sapply(datasets[1:10], function(x) extractReduce_diversityDF(x, immdata.IgH[[x]]), simplify = FALSE)


## TRANSFORM LIST OF DATASETS TO DATAFRAME
immdata.TcrA <-ldply(immdata.TcrA.clean,.id = NULL)

immdata.TcrB <-ldply(immdata.TcrB.clean,.id = NULL)

immdata.IgH <-ldply(immdata.IgH.clean,.id = NULL)

## Filter for q=1
immdata.TcrA <-immdata.TcrA %>% filter(Q.hill ==1) %>% droplevels()

immdata.TcrB <-immdata.TcrB %>% filter(Q.hill ==1) %>% droplevels()

immdata.IgH <-immdata.IgH %>% filter(Q.hill ==1) %>% droplevels()

## Manual clonality measures, filter for clone columns, survival, tmb
tcrSamTableMerged.red <- tcrSamTableMerged %>% select(-cloneId,-allVHitsWithScore,-nSeqCDR3, -aaSeqCDR3,- chain,-cloneCount,-cloneFraction,-dominant,-dominantMetric,-pcTotTCR,-pcTotTCRchain,-pcTotReads,-age, -gender, -disease_status, -biopsy_site,
                                                      -primary_tumor, -treatment,
                                                      -prior_treatment, -response, -OS.time, -OS,
                                                      -PFS.time, -PFS, -pre_on_treatment,
                                                      -total_muts, -nonsyn_muts, -clonal_muts, -subclonal_muts, -TMB, -tissue, -dataset,-title) %>% distinct()

tcrSamTableMerged.full <- merge(tcrSamTableMerged.red, samples_meta.Abs.compact.fix, by="run_accession",all=TRUE)

### Merge DFs with the metadata object (.full)
## TCRA
immdata.TcrA.full <- merge(immdata.TcrA %>% select(run_accession, TuTACK_sigscore, TIS_sigscore, PDL1, CD8A, CD8B,
                                                   Estimator.chao, SD.chao, Conf.95.lo.chao, Conf.95.hi.chao, Q.hill, Value.hill, Clones.d50, Percentage.d50,
                                                   Value.trueDiv, Clones.dxx, Percentage.dxx, modEntr.ds, TCRArichness.ds),
                           samples_meta.Abs.compact.fix, by="run_accession",all.x=TRUE) #  THIS IS WHY I HAVE DATA THAT ARE NA, EXRTA SAMPLES FROM THE ONES DOWNSAMPLED., setting to all.x keeps only the downsampled information



## TCRB
immdata.TcrB.full <- merge(immdata.TcrB %>% select(run_accession, TuTACK_sigscore, TIS_sigscore, PDL1, CD8A, CD8B,
                                                   Estimator.chao, SD.chao, Conf.95.lo.chao, Conf.95.hi.chao, Q.hill, Value.hill, Clones.d50, Percentage.d50,
                                                   Value.trueDiv, Clones.dxx, Percentage.dxx, modEntr.ds, TCRBrichness.ds),
                           samples_meta.Abs.compact.fix, by="run_accession",all.x=TRUE)


## IGH
immdata.IgH.full <- merge(immdata.IgH %>% select(run_accession, TuTACK_sigscore, TIS_sigscore, PDL1, CD8A, CD8B,
                                                 Estimator.chao, SD.chao, Conf.95.lo.chao, Conf.95.hi.chao, Q.hill, Value.hill, Clones.d50, Percentage.d50,
                                                 Value.trueDiv, Clones.dxx, Percentage.dxx, modEntr.ds, BCR.IGHrichness.ds),
                          samples_meta.Abs.compact.fix, by="run_accession",all.x=TRUE)


## CALCULATE TUMOR PURITY
# Manual clonality measures
tcrSamTableMerged.full$TumorPurity <- cos(0.6049872018+0.0001467884*tcrSamTableMerged.full$ESTIMATEScore)
# Immunarch diversity measures
immdata.TcrA.full$TumorPurity <-cos(0.6049872018+0.0001467884*immdata.TcrA.full$ESTIMATEScore)
# immdata.TcrA.noDS.full$TumorPurity <-cos(0.6049872018+0.0001467884*immdata.TcrA.noDS.full$ESTIMATEScore)


immdata.TcrB.full$TumorPurity <-cos(0.6049872018+0.0001467884*immdata.TcrB.full$ESTIMATEScore)
# immdata.TcrB.noDS.full$TumorPurity <-cos(0.6049872018+0.0001467884*immdata.TcrB.noDS.full$ESTIMATEScore)

immdata.IgH.full$TumorPurity <-cos(0.6049872018+0.0001467884*immdata.IgH.full$ESTIMATEScore)
# immdata.IgH.noDS.full$TumorPurity <-cos(0.6049872018+0.0001467884*immdata.IgH.noDS.full$ESTIMATEScore)

## Correct TMB INF values
immdata.TcrA.full$TMB[ which(is.infinite(log(immdata.TcrA.full$TMB)))] <- min(immdata.TcrA.full$TMB[!immdata.TcrA.full$TMB==0], na.rm = TRUE)

immdata.TcrB.full$TMB[ which(is.infinite(log(immdata.TcrB.full$TMB)))] <- min(immdata.TcrB.full$TMB[!immdata.TcrB.full$TMB==0], na.rm = TRUE)

immdata.IgH.full$TMB[ which(is.infinite(log(immdata.IgH.full$TMB)))] <- min(immdata.IgH.full$TMB[!immdata.IgH.full$TMB==0], na.rm = TRUE)

## CALCULATE LOG TMB
# Manual clonality measures
tcrSamTableMerged.full$logTMB <-log(tcrSamTableMerged.full$TMB)

# Immunarch diversity measures
immdata.TcrA.full$logTMB <-log(immdata.TcrA.full$TMB)

immdata.TcrB.full$logTMB <-log(immdata.TcrB.full$TMB)

immdata.IgH.full$logTMB <-log(immdata.IgH.full$TMB)

## FIND COMMON SAMPLES ACROSS DIFF DIVERSITY MEASURES I AM CONSIDERING
samples.common <- intersect(immdata.TcrA.full$run_accession, immdata.IgH.full$run_accession)
length(samples.common)

## FILTER DATA FOR COMMON SAMPLES (.red)
immdata.TcrA.full.red <- immdata.TcrA.full %>% filter(run_accession %in% samples.common) %>% droplevels()
immdata.IgH.full.red <- immdata.IgH.full %>% filter(run_accession %in% samples.common) %>% droplevels()

## MERGE DATA WITH ALL DIVERSITY MEASURES AND COVARIATES (.TcrBcr.full.red)
immdata.TcrBcr.full.red <- merge( immdata.TcrA.full.red %>% select(run_accession, TuTACK_sigscore, TIS_sigscore, PDL1,
                                                                   CD8A,CD8B, TCRArichness.ds, StromalScore, ImmuneScore,ESTIMATEScore,
                                                                   title,age,gender, disease_status,biopsy_site, primary_tumor,
                                                                   treatment, prior_treatment, response, OS.time, OS, PFS.time, PFS,
                                                                   pre_on_treatment, total_muts, nonsyn_muts, clonal_muts, subclonal_muts,TMB,
                                                                   tissue,dataset, Bcells, CD8cells, CD4cells, Tregs, NKcells, DCcells, Monocytes,Granulocytes, Macrophages, TumorPurity, logTMB),
                                  immdata.IgH.full.red %>% select(run_accession, BCR.IGHrichness.ds), by = "run_accession",all=TRUE)

## FILTER DATA FOR PD, CR, PR (.sub)
immdata.TcrBcr.full.red.sub <- immdata.TcrBcr.full.red %>% filter(response %in% c("PD","CR+PR")) %>% droplevels()
# FIX RESPONSE INFORMATION
immdata.TcrBcr.full.red.sub$response <- ifelse(immdata.TcrBcr.full.red.sub$response=="CR+PR","CRPR",immdata.TcrBcr.full.red.sub$response)


## EXTRACT RESPONSE, NOMINAL, MEASURE VARIABLES FROM DATA
response.var <- "response"
nominal.vars <- c("dataset","tissue","gender")
measure.vars <-colnames(immdata.TcrBcr.full.red.sub)[c(2:10,12,32:43)]
measure.vars <- measure.vars[c(6,22,1:3,20,21,4:5,10,7:9,11:19)]


## CHECK INFINITE LOG TMB VALUES
immdata.TcrBcr.full.red.sub[which(is.infinite(immdata.TcrBcr.full.red.sub$logTMB)),]

## Extract tissues in data so far
tissues <- unique(immdata.TcrBcr.full.red.sub$tissue)

# Keep samples with TMB (.tmb)
immdata.TcrBcr.full.red.sub.tmb <- immdata.TcrBcr.full.red.sub %>% filter(complete.cases(logTMB)) %>% droplevels();nrow(immdata.TcrBcr.full.red.sub.tmb)
# ONLY MELANOMA, WITH TMB (.tmb.MEL)
immdata.TcrBcr.full.red.sub.tmb.MEL <- immdata.TcrBcr.full.red.sub.tmb %>% filter(tissue=="melanoma") %>% droplevels();nrow(immdata.TcrBcr.full.red.sub.tmb.MEL)
## ALL MELANOMA AND BLADDER, WITH TMB (.tmb.ALL)
immdata.TcrBcr.full.red.sub.tmb.ALL <- immdata.TcrBcr.full.red.sub.tmb %>% filter(tissue %in% c("melanoma","bladder")) %>% droplevels();nrow(immdata.TcrBcr.full.red.sub.tmb.ALL)
# ALL DATA, MELANOMA AND BLADDER, BUT LEAVE MISSING TMB (.ALL)
immdata.TcrBcr.full.red.sub.ALL <- immdata.TcrBcr.full.red.sub %>% filter(tissue %in% c("melanoma","bladder")) %>% droplevels();nrow(immdata.TcrBcr.full.red.sub.ALL)


# covariates included
covariatesIn <- measure.vars[c(1:2,4:7,14:16)]
responseVar <- "response"


## DATA TO WORK WITH
#####################################
## (.TcrBcr.full.red) : Merged dataframes of all different diversity metrics, with other data
#=======================#
# table(immdata.TcrBcr.full.red$tissue)
#
# bladder  gastric melanoma      RCC
# 345       43      223       95

## (.TcrBcr.full.red.sub)
#=======================#
# table(immdata.TcrBcr.full.red.sub$tissue) : Merged dataframes of all different diversity metrics, with other data, only PD, CR, PR
#
# bladder  gastric melanoma      RCC
# 249       29      193       56

## (.TcrBcr.full.red.sub.tmb) : All tissues, only PD, CR, PR, ONLY THOSE WITH TMB
#=======================#
# table(immdata.TcrBcr.full.red.sub.tmb$tissue)
#
# bladder melanoma      RCC
# 216      130       34

## (.TcrBcr.full.red.sub.tmb.ALL): All melanoma & bladder, PD, CR, PR, with TMB
#=======================#
# table(immdata.TcrBcr.full.red.sub.tmb.ALL$tissue)
#
# bladder melanoma
# 216      130

## (.TcrBcr.full.red.sub.ALL): All melanoma & bladder, PD, CR, PR (with missing TMB)
# table(immdata.TcrBcr.full.red.sub.ALL$tissue)
#
# bladder melanoma
# 249      193
#####################################
#####################################

## INSTEAD OF FILTERINF FOR NO MISSING TMB, WE DO THAT FOR ALL COVARIATES OF INTEREST
# Remove rows that have outliers/missing values since I get errors with Knn imputation
immdata.TcrBcr.full.red.sub.ALL.filt <- immdata.TcrBcr.full.red.sub.ALL %>% filter_at(vars(covariatesIn), all_vars(complete.cases(.))) %>% droplevels()
# DOING THE SAME FOR THE ORIGINAL DATA - USE THIS WHEN TESTING
# these data include other tissues apart from melanoma and bladder
immdata.TcrBcr.full.red.sub.filt <- immdata.TcrBcr.full.red.sub %>% filter_at(vars(covariatesIn), all_vars(complete.cases(.))) %>% droplevels()

# FINAL DATA FOR MODELLING ####
# Remove the Liu dataset for the model testing
immdata.TcrBcr.BLAD.MEL.ready <-immdata.TcrBcr.full.red.sub.ALL.filt %>% filter(dataset %notin% c("Liu") ) %>% droplevels() # All melanoma & bladder, PD, CR, PR, complete data
immdata.TcrBcr.ALL.ready <- immdata.TcrBcr.full.red.sub.filt %>% filter(dataset %notin% c("Liu") ) %>% droplevels() # All tissues, PD, CR, PR, complete data



###########################################

### SETTING THE GRID ####
sample.params1 <- c("smote","down","up")
sample.params2 <- c("in")
kfold.params <- c(5,10)
repeat.params <-c(30)#10,
split.params <-c(0.66,0.7,0.75,0.8)
outlier.remove <-c(FALSE)
outlier.params <-c("")
scale.params<-c(TRUE)

# The grid
myparamGrid1 <- expand.grid(sample.params1,sample.params2, kfold.params,repeat.params,split.params,
                            outlier.remove,outlier.params,
                            scale.params,stringsAsFactors=TRUE)
dim(myparamGrid1)
myparamGrid1 <- myparamGrid1 %>% distinct()
dim(myparamGrid1)


####################################################
## RUN THE MODELS: A  - allow parallel in model ####
####################################################
# NB #####
if (args[1]=="NB") {
  ## NAIVE BAYES ##
  start_time <- Sys.time()
  # Apply model parameter optimization
  ref_rf_res<- apply(myparamGrid1 , 1,function(x) fullOptimalCutoffModelResult.RFE(df =  immdata.TcrBcr.BLAD.MEL.ready ,
                                                                                  df.extra = immdata.TcrBcr.ALL.ready ,
                                                                                  selectedSplit = as.numeric(x[5]),ClassLevels = c("CRPR","PD"),
                                                                                  PositiveClass = "CRPR",r = as.numeric(x[4]),f = as.numeric(x[3]),
                                                                                  METRIC = "ROC",
                                                                                  subSamplingParams1 = as.character(x[1]),subSamplingParams2 =as.character(x[2]),
                                                                                  cD = FALSE,rO =  as.logical(x[6]),
                                                                                  quant = .99,outMethod = as.character(x[7]),pCl = FALSE,
                                                                                  filt = TRUE,filtVar = "dataset",filtSel = c("Powles","Gide","Mari","Rod","Hugo","Riaz"),
                                                                                  cr = FALSE,crt = 0.85,tissue="none",
                                                                                  doScaling=TRUE,modelName='',
                                                                                  applyStratify = FALSE,
                                                                                  mainPredictors = covariatesIn,
                                                                                  rfeSizeVector = c(1:9),
                                                                                  selectedRFEModel = "NB",
                                                                                  yVar = "response",
                                                                                  sampleID = "run_accession"))

  end_time <- Sys.time()
  end_time - start_time
  
  ref_rf_res<-as.data.frame(t(ref_rf_res))
  ref_rf_res <-ref_rf_res %>%
    mutate_if(is.factor,as.character)
  ref_rf_res <- ref_rf_res %>% mutate_at(colnames(ref_rf_res)[c(2:5,8,13,16:18)], as.numeric)
  
  save(ref_rf_res,file=paste0(modelsPath,"optimizationParamsModels","_RFE","_NB",modelTypeA,modelTypeB, modelTypeC,Rdata.suffix))
  
# LR ##### 
}else if(args[1]=="LR"){
  ## LINEAR REGRESSION ##
  start_time <- Sys.time()
  # Apply model parameter optimization
  ref_rf_res<- apply(myparamGrid1, 1,function(x) fullOptimalCutoffModelResult.RFE(df =  immdata.TcrBcr.BLAD.MEL.ready ,
                                                                                  df.extra = immdata.TcrBcr.ALL.ready ,
                                                                                  selectedSplit = as.numeric(x[5]),ClassLevels = c("CRPR","PD"),
                                                                                  PositiveClass = "CRPR",r = as.numeric(x[4]),f = as.numeric(x[3]),
                                                                                  METRIC = "ROC",
                                                                                  subSamplingParams1 = as.character(x[1]),subSamplingParams2 =as.character(x[2]),
                                                                                  cD = FALSE,rO =  as.logical(x[6]),
                                                                                  quant = .99,outMethod = as.character(x[7]),pCl = FALSE,
                                                                                  filt = TRUE,filtVar = "dataset",filtSel = c("Powles","Gide","Mari","Rod","Hugo","Riaz"),
                                                                                  cr = FALSE,crt = 0.85,tissue="none",
                                                                                  doScaling=TRUE,modelName='',
                                                                                  applyStratify = FALSE,
                                                                                  mainPredictors = covariatesIn,
                                                                                  rfeSizeVector = c(1:9),
                                                                                  selectedRFEModel = "LR",
                                                                                  yVar = "response",
                                                                                  sampleID = "run_accession"))

  end_time <- Sys.time()
  end_time - start_time

  ref_rf_res<-as.data.frame(t(ref_rf_res))
  ref_rf_res <-ref_rf_res %>%
    mutate_if(is.factor,as.character)
  ref_rf_res <- ref_rf_res %>% mutate_at(colnames(ref_rf_res)[c(2:5,8,13,16:18)], as.numeric)



  save(ref_rf_res,file=paste0(modelsPath,"optimizationParamsModels","_RFE","_LR",modelTypeA,modelTypeB, modelTypeC,Rdata.suffix))
# RF #####
}else if(args[1]=="RF"){
  ## RANDOM FOREST ##
  start_time <- Sys.time()
  # Apply model parameter optimization
  ref_rf_res<- apply(myparamGrid1, 1,function(x) fullOptimalCutoffModelResult.RFE(df =  immdata.TcrBcr.BLAD.MEL.ready ,
                                                                                  df.extra = immdata.TcrBcr.ALL.ready ,
                                                                                  selectedSplit = as.numeric(x[5]),ClassLevels = c("CRPR","PD"),
                                                                                  PositiveClass = "CRPR",r = as.numeric(x[4]),f = as.numeric(x[3]),
                                                                                  METRIC = "ROC",
                                                                                  subSamplingParams1 = as.character(x[1]),subSamplingParams2 =as.character(x[2]),
                                                                                  cD = FALSE,rO =  as.logical(x[6]),
                                                                                  quant = .99,outMethod = as.character(x[7]),pCl = FALSE,
                                                                                  filt = TRUE,filtVar = "dataset",filtSel = c("Powles","Gide","Mari","Rod","Hugo","Riaz"),
                                                                                  cr = FALSE,crt = 0.85,tissue="none",
                                                                                  doScaling=TRUE,modelName='',
                                                                                  applyStratify = FALSE,
                                                                                  mainPredictors = covariatesIn,
                                                                                  rfeSizeVector = c(1:9),
                                                                                  selectedRFEModel = "RF",
                                                                                  yVar = "response",
                                                                                  sampleID = "run_accession"))

  end_time <- Sys.time()
  end_time - start_time

  ref_rf_res<-as.data.frame(t(ref_rf_res))
  ref_rf_res <-ref_rf_res %>%
    mutate_if(is.factor,as.character)
  ref_rf_res <- ref_rf_res %>% mutate_at(colnames(ref_rf_res)[c(2:5,8,13,16:18)], as.numeric)


  save(ref_rf_res,file=paste0(modelsPath,"optimizationParamsModels","_RFE","_RF",modelTypeA,modelTypeB, modelTypeC,Rdata.suffix))
# SVM radial ####
}else if(args[1]=="SVMr"){
  ## SVM radial ##
  start_time <- Sys.time()
  # Apply model parameter optimization
  ref_rf_res<- apply(myparamGrid1, 1,function(x) fullOptimalCutoffModelResult.RFE(df =  immdata.TcrBcr.BLAD.MEL.ready ,
                                                                                  df.extra = immdata.TcrBcr.ALL.ready ,
                                                                                  selectedSplit = as.numeric(x[5]),ClassLevels = c("CRPR","PD"),
                                                                                  PositiveClass = "CRPR",r = as.numeric(x[4]),f = as.numeric(x[3]),
                                                                                  METRIC = "ROC",
                                                                                  subSamplingParams1 = as.character(x[1]),subSamplingParams2 =as.character(x[2]),
                                                                                  cD = FALSE,rO =  as.logical(x[6]),
                                                                                  quant = .99,outMethod = as.character(x[7]),pCl = FALSE,
                                                                                  filt = TRUE,filtVar = "dataset",filtSel = c("Powles","Gide","Mari","Rod","Hugo","Riaz"),
                                                                                  cr = FALSE,crt = 0.85,tissue="none",
                                                                                  doScaling=TRUE,modelName='',
                                                                                  applyStratify = FALSE,
                                                                                  mainPredictors = covariatesIn,
                                                                                  rfeSizeVector = c(1:9),
                                                                                  selectedRFEModel = "SVMr",
                                                                                  yVar = "response",
                                                                                  sampleID = "run_accession"))

  end_time <- Sys.time()
  end_time - start_time

  ref_rf_res<-as.data.frame(t(ref_rf_res))
  ref_rf_res <-ref_rf_res %>%
    mutate_if(is.factor,as.character)
  ref_rf_res <- ref_rf_res %>% mutate_at(colnames(ref_rf_res)[c(2:5,8,13,16:18)], as.numeric)



  save(ref_rf_res,file=paste0(modelsPath,"optimizationParamsModels","_RFE","_SVMr",modelTypeA,modelTypeB, modelTypeC,Rdata.suffix))
# SVM linear #####
}else if(args[1]=="SVMl"){
  ## SVMl ##
  start_time <- Sys.time()
  # Apply model parameter optimization
  ref_rf_res<- apply(myparamGrid1, 1,function(x) fullOptimalCutoffModelResult.RFE(df =  immdata.TcrBcr.BLAD.MEL.ready ,
                                                                                  df.extra = immdata.TcrBcr.ALL.ready ,
                                                                                  selectedSplit = as.numeric(x[5]),ClassLevels = c("CRPR","PD"),
                                                                                  PositiveClass = "CRPR",r = as.numeric(x[4]),f = as.numeric(x[3]),
                                                                                  METRIC = "ROC",
                                                                                  subSamplingParams1 = as.character(x[1]),subSamplingParams2 =as.character(x[2]),
                                                                                  cD = FALSE,rO =  as.logical(x[6]),
                                                                                  quant = .99,outMethod = as.character(x[7]),pCl = FALSE,
                                                                                  filt = TRUE,filtVar = "dataset",filtSel = c("Powles","Gide","Mari","Rod","Hugo","Riaz"),
                                                                                  cr = FALSE,crt = 0.85,tissue="none",
                                                                                  doScaling=TRUE,modelName='',
                                                                                  applyStratify = FALSE,
                                                                                  mainPredictors = covariatesIn,
                                                                                  rfeSizeVector = c(1:9),
                                                                                  selectedRFEModel = "SVMl",
                                                                                  yVar = "response",
                                                                                  sampleID = "run_accession"))

  end_time <- Sys.time()
  end_time - start_time

  ref_rf_res<-as.data.frame(t(ref_rf_res))
  ref_rf_res <-ref_rf_res %>%
    mutate_if(is.factor,as.character)
  ref_rf_res <- ref_rf_res %>% mutate_at(colnames(ref_rf_res)[c(2:5,8,13,16:18)], as.numeric)




  save(ref_rf_res,file=paste0(modelsPath,"optimizationParamsModels","_RFE","_SVMl",modelTypeA,modelTypeB, modelTypeC,Rdata.suffix))
# KNN ####
}else if(args[1]=="KNN"){
  ## KNN ##
  start_time <- Sys.time()
  # Apply model parameter optimization
  ref_rf_res<- apply(myparamGrid1, 1,function(x) fullOptimalCutoffModelResult.RFE(df =  immdata.TcrBcr.BLAD.MEL.ready ,
                                                                                  df.extra = immdata.TcrBcr.ALL.ready ,
                                                                                  selectedSplit = as.numeric(x[5]),ClassLevels = c("CRPR","PD"),
                                                                                  PositiveClass = "CRPR",r = as.numeric(x[4]),f = as.numeric(x[3]),
                                                                                  METRIC = "ROC",
                                                                                  subSamplingParams1 = as.character(x[1]),subSamplingParams2 =as.character(x[2]),
                                                                                  cD = FALSE,rO =  as.logical(x[6]),
                                                                                  quant = .99,outMethod = as.character(x[7]),pCl = FALSE,
                                                                                  filt = TRUE,filtVar = "dataset",filtSel = c("Powles","Gide","Mari","Rod","Hugo","Riaz"),
                                                                                  cr = FALSE,crt = 0.85,tissue="none",
                                                                                  doScaling=TRUE,modelName='',
                                                                                  applyStratify = FALSE,
                                                                                  mainPredictors = covariatesIn,
                                                                                  rfeSizeVector = c(1:9),
                                                                                  selectedRFEModel = "KNN",
                                                                                  yVar = "response",
                                                                                  sampleID = "run_accession"))

  end_time <- Sys.time()
  end_time - start_time

  ref_rf_res<-as.data.frame(t(ref_rf_res))
  ref_rf_res <-ref_rf_res %>%
    mutate_if(is.factor,as.character)
  ref_rf_res <- ref_rf_res %>% mutate_at(colnames(ref_rf_res)[c(2:5,8,13,16:18)], as.numeric)




  save(ref_rf_res,file=paste0(modelsPath,"optimizationParamsModels","_RFE","_KNN",modelTypeA,modelTypeB, modelTypeC,Rdata.suffix))
# SINGLE models- GET OPTIMIZATION PARAMETERS #####
}else if(args[1]=="single"){
  #######################
  ## RUN SINGLE MODELS ##
  #######################
  rfe_rf <- loadRData(paste0(modelsPath,"optimizationParamsModels","_RFE","_RF",modelTypeA,modelTypeB, modelTypeC,Rdata.suffix))
  # rfe_rf[which(rfe_rf$ROC==max(rfe_rf$ROC)),]
  # hist(rfe_rf$ROC)
  # rfe_rf[order(rfe_rf$ROC, decreasing = TRUE),] %>% head()
  rfe_rf.opt <-rfe_rf[which(rfe_rf$ROC==max(rfe_rf$ROC)),]
  # #
  # #
  rfe_lr <- loadRData(paste0(modelsPath,"optimizationParamsModels","_RFE","_LR",modelTypeA,modelTypeB, modelTypeC,Rdata.suffix))
  # rfe_lr[which(rfe_lr$ROC==max(rfe_lr$ROC)),]
  # hist(rfe_lr$ROC)
  # rfe_lr[order(rfe_lr$ROC, decreasing = TRUE),] %>% head()
  rfe_lr.opt <-rfe_lr[which(rfe_lr$ROC==max(rfe_lr$ROC)),]
  # #
  # #
  # #
  rfe_nb <- loadRData(paste0(modelsPath,"optimizationParamsModels","_RFE","_NB",modelTypeA,modelTypeB, modelTypeC,Rdata.suffix))
  # rfe_nb[which(rfe_nb$ROC==max(rfe_nb$ROC)),]
  # hist(rfe_nb$ROC)
  # rfe_nb[order(rfe_nb$ROC, decreasing = TRUE),] %>% head()
  rfe_nb.opt <-rfe_nb[which(rfe_nb$ROC==max(rfe_nb$ROC)),]

 
  rfe_svmr <- loadRData(paste0(modelsPath,"optimizationParamsModels","_RFE","_SVMr",modelTypeA,modelTypeB, modelTypeC,Rdata.suffix))
  # rfe_svmr[which(rfe_svmr$ROC==max(rfe_svmr$ROC)),]
  # hist(rfe_svmr$ROC)
  # rfe_svmr[order(rfe_svmr$ROC, decreasing = TRUE),] %>% head()
  rfe_svmr.opt <-rfe_svmr[which(rfe_svmr$ROC==max(rfe_svmr$ROC)),]
  # #
  rfe_svml <- loadRData(paste0(modelsPath,"optimizationParamsModels","_RFE","_SVMl",modelTypeA,modelTypeB, modelTypeC,Rdata.suffix))
  # rfe_svml[which(rfe_svml$ROC==max(rfe_svml$ROC)),]
  # hist(rfe_svml$ROC)
  # rfe_svml[order(rfe_svml$ROC, decreasing = TRUE),] %>% head()
  rfe_svml.opt <-rfe_svml[which(rfe_svml$ROC==max(rfe_svml$ROC)),]
  # #
  rfe_knn <- loadRData(paste0(modelsPath,"optimizationParamsModels","_RFE","_KNN",modelTypeA,modelTypeB, modelTypeC,Rdata.suffix))
  # rfe_knn[which(rfe_knn$ROC==max(rfe_knn$ROC)),]
  # hist(rfe_knn$ROC)
  # rfe_knn[order(rfe_knn$ROC, decreasing = TRUE),] %>% head()
  rfe_knn.opt <-rfe_knn[which(rfe_knn$ROC==max(rfe_knn$ROC)),]
  # #
  # #
  rfe.models.opt <- list(`Random Forest` = rfe_rf.opt,
                         `Linear Logistic Regression` = rfe_lr.opt,
                         `Naive Bayes` = rfe_nb.opt,
                         `SVM Radial` = rfe_svmr.opt,
                         `SVM Linear` = rfe_svml.opt,
                         `K-Nearest Neighbours` = rfe_knn.opt)

  rfe.models.opt.DF <-bind_rows(rfe.models.opt, .id = "model")

  # Find optimal parameters and model
  # rfe.models.opt.DF2 <- rbind(rfe.models.opt.DF,rfe.models.opt.DF)
  rfe.models.opt.final <-rfe.models.opt.DF %>% filter(ROC == max(ROC)) %>% slice(1)
  
  optimal.model <- rfe.models.opt.final$model
  print(paste0(">>>>Model giving best results is: ", optimal.model,"<<<<<<<"))
  
  kFold.opt = rfe.models.opt.final$kfolds
  Sampling.opt =rfe.models.opt.final$subsampling
  split.opt =rfe.models.opt.final$split
  optList = list(kfold = kFold.opt,
                 sampling=Sampling.opt,
                 split = split.opt)
  model_rf <-trainRUN.RFE(DF = immdata.TcrBcr.BLAD.MEL.ready,
                          split = split.opt, ClassLevels = c("CRPR","PD"), PositiveClass = "CRPR",LevelIndex = 1,
                          seed.partition = 3456,seed.model = 5627,setCVseeds=NULL,
                          dummy=FALSE,dummyVars=c("tissue"),categVars=c("tissue","dataset"),
                          obsColumn = "run_accession",responseVar = "response",
                          REPEATS=30,KFOLDS=kFold.opt,
                          metricOpt = "ROC",subsampMethod=Sampling.opt,inoutSub="in",
                          modelDesign=paste(covariatesIn,collapse = " + "),controlDum=FALSE,
                          model.name = "RF",mainFeats=covariatesIn,
                          removeOutliers=FALSE,quantileClass1=.95,quantileClass2=.95,
                          selectOutlierMethod= "meanDist",
                          stratify=FALSE,fullResults=TRUE,
                          positiveCoefs=FALSE,
                          filter=TRUE,filterVar="dataset",filterSelected=c("Powles","Gide","Mari","Rod","Hugo","Riaz"),
                          CorrRemove=FALSE,corrThres=.85,
                          selectedScale=TRUE,
                          sizeVector = c(1:9),
                          rfeModel = "RF")
  save(model_rf,file=paste0(modelsPath,"rf_optParams",modelTypeA,modelTypeB,modelTypeC,Rdata.suffix))

  model_lr <-trainRUN.RFE(DF = immdata.TcrBcr.BLAD.MEL.ready,
                          split = split.opt, ClassLevels = c("CRPR","PD"), PositiveClass = "CRPR",LevelIndex = 1,
                          seed.partition = 3456,seed.model = 5627,setCVseeds=NULL,
                          dummy=FALSE,dummyVars=c("tissue"),categVars=c("tissue","dataset"),
                          obsColumn = "run_accession",responseVar = "response",
                          REPEATS=30,KFOLDS=kFold.opt,
                          metricOpt = "ROC",subsampMethod=Sampling.opt,inoutSub="in",
                          modelDesign=paste(covariatesIn,collapse = " + "),controlDum=FALSE,
                          model.name = "LR",mainFeats=covariatesIn,
                          removeOutliers=FALSE,quantileClass1=.95,quantileClass2=.95,
                          selectOutlierMethod= "meanDist",
                          stratify=FALSE,fullResults=TRUE,
                          positiveCoefs=FALSE,
                          filter=TRUE,filterVar="dataset",filterSelected=c("Powles","Gide","Mari","Rod","Hugo","Riaz"),
                          CorrRemove=FALSE,corrThres=.85,
                          selectedScale=TRUE,
                          sizeVector = c(1:9),
                          rfeModel = "LR")
  save(model_lr,file=paste0(modelsPath,"lr_optParams",modelTypeA,modelTypeB,modelTypeC,Rdata.suffix))

  model_nb <-trainRUN.RFE(DF = immdata.TcrBcr.BLAD.MEL.ready,
                          split = split.opt, ClassLevels = c("CRPR","PD"), PositiveClass = "CRPR",LevelIndex = 1,
                          seed.partition = 3456,seed.model = 5627,setCVseeds=NULL,
                          dummy=FALSE,dummyVars=c("tissue"),categVars=c("tissue","dataset"),
                          obsColumn = "run_accession",responseVar = "response",
                          REPEATS=30,KFOLDS=kFold.opt,
                          metricOpt = "ROC",subsampMethod=Sampling.opt,inoutSub="in",
                          modelDesign=paste(covariatesIn,collapse = " + "),controlDum=FALSE,
                          model.name = "NB",mainFeats=covariatesIn,
                          removeOutliers=FALSE,quantileClass1=.95,quantileClass2=.95,
                          selectOutlierMethod= "meanDist",
                          stratify=FALSE,fullResults=TRUE,
                          positiveCoefs=FALSE,
                          filter=TRUE,filterVar="dataset",filterSelected=c("Powles","Gide","Mari","Rod","Hugo","Riaz"),
                          CorrRemove=FALSE,corrThres=.85,
                          selectedScale=TRUE,
                          sizeVector = c(1:9),
                          rfeModel = "NB")
  save(model_nb,file=paste0(modelsPath,"nb_optParams",modelTypeA,modelTypeB,modelTypeC,Rdata.suffix))

  model_svmr <-trainRUN.RFE(DF = immdata.TcrBcr.BLAD.MEL.ready,
                            split = split.opt, ClassLevels = c("CRPR","PD"), PositiveClass = "CRPR",LevelIndex = 1,
                            seed.partition = 3456,seed.model = 5627,setCVseeds=NULL,
                            dummy=FALSE,dummyVars=c("tissue"),categVars=c("tissue","dataset"),
                            obsColumn = "run_accession",responseVar = "response",
                            REPEATS=30,KFOLDS=kFold.opt,
                            metricOpt = "ROC",subsampMethod=Sampling.opt,inoutSub="in",
                            modelDesign=paste(covariatesIn,collapse = " + "),controlDum=FALSE,
                            model.name = "SVMr",mainFeats=covariatesIn,
                            removeOutliers=FALSE,quantileClass1=.95,quantileClass2=.95,
                            selectOutlierMethod= "meanDist",
                            stratify=FALSE,fullResults=TRUE,
                            positiveCoefs=FALSE,
                            filter=TRUE,filterVar="dataset",filterSelected=c("Powles","Gide","Mari","Rod","Hugo","Riaz"),
                            CorrRemove=FALSE,corrThres=.85,
                            selectedScale=TRUE,
                            sizeVector = c(1:9),
                            rfeModel = "SVMr")
  save(model_svmr,file=paste0(modelsPath,"svmr_optParams",modelTypeA,modelTypeB,modelTypeC,Rdata.suffix))

  model_svml <-trainRUN.RFE(DF = immdata.TcrBcr.BLAD.MEL.ready,
                            split = split.opt, ClassLevels = c("CRPR","PD"), PositiveClass = "CRPR",LevelIndex = 1,
                            seed.partition = 3456,seed.model = 5627,setCVseeds=NULL,
                            dummy=FALSE,dummyVars=c("tissue"),categVars=c("tissue","dataset"),
                            obsColumn = "run_accession",responseVar = "response",
                            REPEATS=30,KFOLDS=kFold.opt,
                            metricOpt = "ROC",subsampMethod=Sampling.opt,inoutSub="in",
                            modelDesign=paste(covariatesIn,collapse = " + "),controlDum=FALSE,
                            model.name = "SVMl",mainFeats=covariatesIn,
                            removeOutliers=FALSE,quantileClass1=.95,quantileClass2=.95,
                            selectOutlierMethod= "meanDist",
                            stratify=FALSE,fullResults=TRUE,
                            positiveCoefs=FALSE,
                            filter=TRUE,filterVar="dataset",filterSelected=c("Powles","Gide","Mari","Rod","Hugo","Riaz"),
                            CorrRemove=FALSE,corrThres=.85,
                            selectedScale=TRUE,
                            sizeVector = c(1:9),
                            rfeModel = "SVMl")
  save(model_svml,file=paste0(modelsPath,"svml_optParams",modelTypeA,modelTypeB,modelTypeC,Rdata.suffix))

  model_knn <-trainRUN.RFE(DF = immdata.TcrBcr.BLAD.MEL.ready,
                           split = split.opt, ClassLevels = c("CRPR","PD"), PositiveClass = "CRPR",LevelIndex = 1,
                           seed.partition = 3456,seed.model = 5627,setCVseeds=NULL,
                           dummy=FALSE,dummyVars=c("tissue"),categVars=c("tissue","dataset"),
                           obsColumn = "run_accession",responseVar = "response",
                           REPEATS=30,KFOLDS=kFold.opt,
                           metricOpt = "ROC",subsampMethod=Sampling.opt,inoutSub="in",
                           modelDesign=paste(covariatesIn,collapse = " + "),controlDum=FALSE,
                           model.name = "KNN",mainFeats=covariatesIn,
                           removeOutliers=FALSE,quantileClass1=.95,quantileClass2=.95,
                           selectOutlierMethod= "meanDist",
                           stratify=FALSE,fullResults=TRUE,
                           positiveCoefs=FALSE,
                           filter=TRUE,filterVar="dataset",filterSelected=c("Powles","Gide","Mari","Rod","Hugo","Riaz"),
                           CorrRemove=FALSE,corrThres=.85,
                           selectedScale=TRUE,
                           sizeVector = c(1:9),
                           rfeModel = "KNN")
  save(model_knn,file=paste0(modelsPath,"knn_optParams",modelTypeA,modelTypeB,modelTypeC,Rdata.suffix))

}

######################################################
## RUN THE MODELS: B - run each model in separate core
######################################################

# library(parallel)
# library(doParallel)
# cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
# registerDoParallel(cluster)
# stopCluster(cluster)
