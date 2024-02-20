projectDir <- 'D:/BIOINF/PROJECTS/2_ANALYSES/18_GIT_PAPERS/TCR_BCR_antiPD1/'
setwd(projectDir)
#####################
### Load packages ###
#####################
library(pacman)
pacman::p_load(DT,dplyr,tibble,tidyverse,ggplot2,hrbrthemes,
               ggpubr, reshape2,cowplot,ggrepel,cowplot,grid,
               RColorBrewer, pander,gridExtra,
               reconPlots,flextable, extrafont,Amelia, estimate,strex,plyr)



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
antiPDL1publicIntermediateData <- paste0(projectOutDataPath,"antiPD(L)1_public/")

# Session/dependencies
sessionInfoPath <- paste0("session_info/")


##########################
### FOLDER/FILES NAMES ###
##########################
# Create list with folders names for the different datasets.
datasets <- c("Hugo","Riaz","Gide","Liu","Rod","Miao","Immotion150","Kim","Mari", "Powles","Braun")
data.foldersList <- as.list(c("Hugo_MEL","Riaz_MEL","Gide_MEL", "Liu_MEL","Rodriguez_MEL","Miao_RCC","Immo150_RCC","Kim_GAS","Mari_BLA","Powles_BLA","Braun_RCC"))
names(data.foldersList)<- c("Hugo","Riaz","Gide", "Liu","Rod", "Miao","Immotion150","Kim","Mari","Powles","Braun")


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

#################
### Functions ###
#################
`%notin%` <- Negate(`%in%`)


##########################
## Loading RData object ##
##########################


loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}


############################
## Uppercase first letter ##
############################
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

extractReduce_diversityDF <- function(datasetName, immDF){
  # datasetName <- datasets[1:10][2]
  # immDF <- immdata.apd1.proc.TcrA.divEst[[datasetName]]
  
  # print(datasetName)
  # Get metadata df
  metaDF <- immDF$meta
  
  # Get diversity DF
  # First need to process the diversityEst object, to have the samples ids as rows, or with same sample column
  immDF$diversityEst$chao <- immDF$diversityEst$chao %>% as.data.frame() %>% rownames_to_column(var="Sample")
  colnames(immDF$diversityEst$chao)[2:5] <- paste0(colnames(immDF$diversityEst$chao)[2:5],".chao")
  # Hill
  colnames(immDF$diversityEst$hill)[2:3] <- paste0(colnames(immDF$diversityEst$hill)[2:3],".hill")
  #d50
  immDF$diversityEst$d50 <- immDF$diversityEst$d50 %>% as.data.frame() %>% rownames_to_column(var="Sample")
  colnames(immDF$diversityEst$d50)[2:3] <- paste0(colnames(immDF$diversityEst$d50)[2:3],".d50")
  # trueDiv
  colnames(immDF$diversityEst$trueDiv)[2] <- paste0(colnames(immDF$diversityEst$trueDiv)[2],".trueDiv")
  #dxx
  immDF$diversityEst$dxx <- immDF$diversityEst$dxx %>% as.data.frame() %>% rownames_to_column(var="Sample")
  colnames(immDF$diversityEst$dxx)[2:3] <- paste0(colnames(immDF$diversityEst$dxx)[2:3],".dxx")
  # modeledEntropy
  colnames(immDF$diversityEst[[6]])[1] <- "Sample"
  # Sampleid
  colnames(immDF$diversityEst[[7]])[1] <- "Sample"
  
  
  divDF <-Reduce(function(x, y) merge(x, y, all=TRUE), immDF$diversityEst) #  KEEP IN MIND HOW WE HAVE FOR HILL ESTIMATE DIFFERENT Qs
  
  # Merge metadata with diversity 
  immDF.new <- merge(metaDF, divDF, by.x = "run_accession", by.y = "Sample", all = TRUE)
  immDF.new
}

fix_immoRunAcc <- function(DF,datasetsVec,datasetNum){
  # DF <- samples_all.merged
  # datasetsVec <- datasets
  # datasetNum <- 7
  immoRunAcc <- DF %>% dplyr::filter(dataset==datasetsVec[datasetNum]) %>% droplevels()
  immoRunAcc<- as.character(immoRunAcc$run_accession)
  
  immoRunAcc.fix<-str_after_nth(immoRunAcc, "_", 7)
  immoRunAcc.map <- as.data.frame(cbind(immoRunAcc,immoRunAcc.fix))
  colnames(immoRunAcc.map) <- c("run_accession_full","run_accession_fixed")
  immoRunAcc.map$run_accession_full <- as.character(immoRunAcc.map$run_accession_full)
  immoRunAcc.map$run_accession_fixed<- as.character(immoRunAcc.map$run_accession_fixed)
  
  DF.new <- DF %>% dplyr::filter(dataset %in% datasetsVec[-datasetNum]) %>% droplevels()
  
  DF.Immo <- DF %>% dplyr::filter(dataset %in% datasetsVec[datasetNum]) %>% droplevels()
  
  DF.Immo$run_accession <- immoRunAcc.map$run_accession_fixed[match(DF.Immo$run_accession,immoRunAcc.map$run_accession_full)]
  # Bind fixed
  DF.fix <- rbind(DF.new,DF.Immo)
  DF.fix
}
######################
##### Load data #####

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
  # immdata.TcrA.noDS <- loadRData(paste0(projectRDataPath,"tcrAclones_downSampled.divEst.aPD1.FINALmin1",Rdata.suffix))
  # TCRB
  immdata.TcrB <- loadRData(paste0(antiPDL1publicIntermediateData,"tcrBclones_downSampled.divEst.aPD1.FINALmin4.v1",Rdata.suffix))
  # immdata.TcrB.noDS <- loadRData(paste0(projectRDataPath,"tcrBclones_downSampled.divEst.aPD1.FINALmin1",Rdata.suffix))
  
  # BCR-IGH
  immdata.IgH <- loadRData(paste0(antiPDL1publicIntermediateData,"igHclones_downSampled.divEst.aPD1.FINALmin10.v1",Rdata.suffix))
  # immdata.IgH.noDS <- loadRData(paste0(projectRDataPath,"igHclones_downSampled.divEst.aPD1.FINALmin1",Rdata.suffix))
}else{
  print("Run script B_02_aPD1_TCR_BCR_clones_analysis.RMD to process mixcr data and get diversiy measure tables")
}


# FOR EACH DATASET The diversityEst object is a list, so we need to unlist that and have it as a dataframe, then
# we need to merge that with the metadata object which is already a dataframe.
# After that we will merge all the datasets, concatenate their dataframes.

# TCRA
immdata.TcrA.clean <-sapply(datasets[1:10], function(x) extractReduce_diversityDF(x, immdata.TcrA[[x]]), simplify = FALSE)
# immdata.TcrA.noDS.clean <-sapply(datasets[1:10], function(x) extractReduce_diversityDF(x, immdata.TcrA.noDS[[x]]), simplify = FALSE)

# TCRB
immdata.TcrB.clean <-sapply(datasets[1:10], function(x) extractReduce_diversityDF(x, immdata.TcrB[[x]]), simplify = FALSE)
# immdata.TcrB.noDS.clean <-sapply(datasets[1:10], function(x) extractReduce_diversityDF(x, immdata.TcrB.noDS[[x]]), simplify = FALSE)

# IGH
# TCRB
immdata.IgH.clean <-sapply(datasets[1:10], function(x) extractReduce_diversityDF(x, immdata.IgH[[x]]), simplify = FALSE)
# immdata.IgH.noDS.clean <-sapply(datasets[1:10], function(x) extractReduce_diversityDF(x, immdata.IgH.noDS[[x]]), simplify = FALSE)


## TRANSFORM LIST OF DATASETS TO DATAFRAME
immdata.TcrA <-ldply(immdata.TcrA.clean,.id = NULL)
# immdata.TcrA.noDS <-ldply(immdata.TcrA.noDS.clean,.id = NULL)

immdata.TcrB <-ldply(immdata.TcrB.clean,.id = NULL)
# immdata.TcrB.noDS <-ldply(immdata.TcrB.noDS.clean,.id = NULL)

immdata.IgH <-ldply(immdata.IgH.clean,.id = NULL)
# immdata.IgH.noDS <-ldply(immdata.IgH.noDS.clean,.id = NULL)


## Filter for q=1
immdata.TcrA <-immdata.TcrA %>% dplyr::filter(Q.hill ==1) %>% droplevels()
# immdata.TcrA.noDS <-immdata.TcrA.noDS %>% dplyr::filter(Q.hill ==1) %>% droplevels()

immdata.TcrB <-immdata.TcrB %>% dplyr::filter(Q.hill ==1) %>% droplevels()
# immdata.TcrB.noDS <-immdata.TcrB.noDS %>% dplyr::filter(Q.hill ==1) %>% droplevels()

immdata.IgH <-immdata.IgH %>% dplyr::filter(Q.hill ==1) %>% droplevels()
# immdata.IgH.noDS <-immdata.IgH.noDS %>% dplyr::filter(Q.hill ==1) %>% droplevels()


## Now add Cibersort information:
### Merge DFs with the metadata object (.full)
## TCRA
immdata.TcrA.full <- merge(immdata.TcrA %>% dplyr::select(run_accession, TuTACK_sigscore, TIS_sigscore, PDL1, CD8A, CD8B,  
                                                   Estimator.chao, SD.chao, Conf.95.lo.chao, Conf.95.hi.chao, Q.hill, Value.hill, Clones.d50, Percentage.d50, 
                                                   Value.trueDiv, Clones.dxx, Percentage.dxx, modEntr.ds, TCRArichness.ds), 
                           samples_meta.Abs.compact.fix, by="run_accession",all.x=TRUE) #  THIS IS WHY I HAVE DATA THAT ARE NA, EXRTA SAMPLES FROM THE ONES DOWNSAMPLED., setting to all.x keeps only the downsampled information


## TCRB
immdata.TcrB.full <- merge(immdata.TcrB %>% dplyr::select(run_accession, TuTACK_sigscore, TIS_sigscore, PDL1,CD8A, CD8B,
                                                   Estimator.chao, SD.chao, Conf.95.lo.chao, Conf.95.hi.chao, Q.hill, Value.hill, Clones.d50, Percentage.d50, 
                                                   Value.trueDiv, Clones.dxx, Percentage.dxx, modEntr.ds, TCRBrichness.ds), 
                           samples_meta.Abs.compact.fix, by="run_accession",all.x=TRUE)



## IGH
immdata.IgH.full <- merge(immdata.IgH %>% dplyr::select(run_accession, TuTACK_sigscore, TIS_sigscore, PDL1, CD8A, CD8B, 
                                                 Estimator.chao, SD.chao, Conf.95.lo.chao, Conf.95.hi.chao, Q.hill, Value.hill, Clones.d50, Percentage.d50, 
                                                 Value.trueDiv, Clones.dxx, Percentage.dxx, modEntr.ds, BCR.IGHrichness.ds), 
                          samples_meta.Abs.compact.fix, by="run_accession",all.x=TRUE)


## CALCULATE TUMOR PURITY
# # Manual clonality measures
# tcrSamTableMerged.full$TumorPurity <- cos(0.6049872018+0.0001467884*tcrSamTableMerged.full$ESTIMATEScore)
# Immunarch diversity measures
immdata.TcrA.full$TumorPurity <-cos(0.6049872018+0.0001467884*immdata.TcrA.full$ESTIMATEScore)
# immdata.TcrA.noDS.full$TumorPurity <-cos(0.6049872018+0.0001467884*immdata.TcrA.noDS.full$ESTIMATEScore)


immdata.TcrB.full$TumorPurity <-cos(0.6049872018+0.0001467884*immdata.TcrB.full$ESTIMATEScore)
# immdata.TcrB.noDS.full$TumorPurity <-cos(0.6049872018+0.0001467884*immdata.TcrB.noDS.full$ESTIMATEScore)

immdata.IgH.full$TumorPurity <-cos(0.6049872018+0.0001467884*immdata.IgH.full$ESTIMATEScore)
# immdata.IgH.noDS.full$TumorPurity <-cos(0.6049872018+0.0001467884*immdata.IgH.noDS.full$ESTIMATEScore)

## IMPORTANT NOTE
# I have an issue with one sample from the Mariathasan dataset, where TMB(TMB per MB) from the clinical metadata file (Foundation one, gene panel) is assigned a value of zero, so when I try to log it I get Infinite value, which forces me to use a pseudocount, which I don't want to. After discussing  with Zsofi, we decide to assign for this particular sample the smallest TMB found in our data, which is very close to zero, so that we avoid this issue. In the Mariathasan paper, they considered patients with zero TMB as outliers and removed them from analysis. I do not want to lose patients at this moment, so I am not excluding this patient.

## Correct TMB INF values
# immdata.TcrA.full[ which(is.infinite(log(immdata.TcrA.full$TMB))),]

immdata.TcrA.full$TMB[ which(is.infinite(log(immdata.TcrA.full$TMB)))] <- min(immdata.TcrA.full$TMB[!immdata.TcrA.full$TMB==0], na.rm = TRUE)


# immdata.TcrB.full[ which(is.infinite(log(immdata.TcrB.full$TMB))),]

immdata.TcrB.full$TMB[ which(is.infinite(log(immdata.TcrB.full$TMB)))] <- min(immdata.TcrB.full$TMB[!immdata.TcrB.full$TMB==0], na.rm = TRUE)


# immdata.IgH.full[ which(is.infinite(log(immdata.IgH.full$TMB))),]

immdata.IgH.full$TMB[ which(is.infinite(log(immdata.IgH.full$TMB)))] <- min(immdata.IgH.full$TMB[!immdata.IgH.full$TMB==0], na.rm = TRUE)


## CALCULATE LOG TMB
# # Manual clonality measures
# tcrSamTableMerged.full$logTMB <-log(tcrSamTableMerged.full$TMB)

# Immunarch diversity measures
immdata.TcrA.full$logTMB <-log(immdata.TcrA.full$TMB)
# immdata.TcrA.noDS.full$logTMB <-log(immdata.TcrA.noDS.full$TMB)

immdata.TcrB.full$logTMB <-log(immdata.TcrB.full$TMB)
# immdata.TcrB.noDS.full$logTMB <-log(immdata.TcrB.noDS.full$TMB)

immdata.IgH.full$logTMB <-log(immdata.IgH.full$TMB)
# immdata.IgH.noDS.full$logTMB <-log(immdata.IgH.noDS.full$TMB)

## keep samples that we have info on TCRA and IGH.

## FIND COMMON SAMPLES ACROSS DIFF DIVERSITY MEASURES I AM CONSIDERING
samples.common <- intersect(immdata.TcrA.full$run_accession, immdata.IgH.full$run_accession)
length(samples.common)

## FILTER DATA FOR COMMON SAMPLES (.red)
immdata.TcrA.full.red <- immdata.TcrA.full %>% dplyr::filter(run_accession %in% samples.common) %>% droplevels()
immdata.IgH.full.red <- immdata.IgH.full %>% dplyr::filter(run_accession %in% samples.common) %>% droplevels()

## MERGE DATA WITH ALL DIVERSITY MEASURES AND COVARIATES (.TcrBcr.full.red)
immdata.TcrBcr.full.red <- merge(immdata.TcrA.full.red %>% dplyr::select(run_accession, TuTACK_sigscore, TIS_sigscore, PDL1,
                                                                   CD8A,CD8B, TCRArichness.ds, StromalScore, ImmuneScore,ESTIMATEScore,
                                                                   title,age,gender, disease_status,biopsy_site, primary_tumor,
                                                                   treatment, prior_treatment, response, OS.time, OS, PFS.time, PFS,
                                                                   pre_on_treatment, total_muts, nonsyn_muts, clonal_muts, subclonal_muts,TMB,
                                                                   tissue,dataset, Bcells, CD8cells, CD4cells, Tregs, NKcells, DCcells, Monocytes,
                                                                   Granulocytes, Macrophages, TumorPurity, logTMB),
                                  immdata.IgH.full.red %>% dplyr::select(run_accession, BCR.IGHrichness.ds), by = "run_accession",all=TRUE) %>% distinct() %>% droplevels()


#### SAVE DATA ### 
save(immdata.TcrBcr.full.red, file = paste0(antiPDL1publicIntermediateData,"immdata.TcrBcr.full.red",modelTypeA,modelTypeB,modelTypeC,".RData" ))

####################
### Session info ###
####################

session_info <- sessionInfo()
writeLines(capture.output(session_info), paste0(sessionInfoPath,"B_03_aPD1_merge_ready_data.txt"))

sessionInfo()

