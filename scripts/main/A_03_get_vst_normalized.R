#! /usr/bin/env Rscript

#####################
### Load packages ###
#####################
# Run below for installing RTCGA first
if (!require(devtools)) {
  install.packages("devtools")
  require(devtools)
}

if (!require(TCGAbiolinks)) {
  BiocManager::install("BioinformaticsFMRP/TCGAbiolinks" )
  require(TCGAbiolinks)
}


# install_version("TCGAbiolinks", version = "2.15.3")


library(pacman)
pacman::p_load(DESeq2,TCGAbiolinks)
###########################
### Main analysis paths ###
###########################
projectDir <- 'D:/BIOINF/PROJECTS/2_ANALYSES/18_GIT_PAPERS/TCR_BCR_antiPD1/'
setwd(projectDir)

projectDataPath <- paste0("data/")
tcgaInputData <- paste0(projectDir, projectDataPath,"TCGA/")
htseqCounts_Path <-paste0(tcgaInputData, "htseq_vst/")

# Session/dependencies
sessionInfoPath <- paste0("session_info/")

# setwd(tcgaInputData)

#################
### Functions ###
#################
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

vst.htseq <- function(countdata,conditionLabel){
  condition <- factor(rep(conditionLabel, ncol(countdata)))
  dds <- DESeqDataSetFromMatrix(countdata, DataFrame(condition), ~ 1)
  vstdata <- DESeq2::vst( dds,blind = TRUE)
  assay(vstdata)
}

get_vst_tcga<- function(selectedcancer,condition, outputPath){
  print(paste0(selectedcancer,"_RNASeqCounts.RData"))
  fname <- file.path(outputPath,
                     paste0(selectedcancer,"_RNASeqCounts.RData"))
  if (file.exists(fname)){
    load(paste0(outputPath,selectedcancer,"_RNASeqCounts.RData"))
    data<- df_all
    if(length(colnames(data)) <=1){
      print("One or less samples available, cannot normalize!Skipping...")
    }else{
      vst.DF <-vst.htseq(data,condition)
      save(vst.DF,file=paste0(outputPath,selectedcancer,"_RNASeqCountsVST.RData"))
    }
    
  }else{
    print(paste0("No htseq count data for ", selectedcancer))
  }
}

################################
### Calculate vst normalized ###
################################

tcga_projects <-TCGAbiolinks:::getGDCprojects()$project_id
mycancertypes <- tcga_projects[grepl("^TCGA-",tcga_projects)] # total 33

save.vst.data <- sapply(mycancertypes, function(x) get_vst_tcga(x,"tumor", outputPath = htseqCounts_Path))


####################
### Session info ###
####################

session_info <- sessionInfo()
writeLines(capture.output(session_info), paste0(projectDir,sessionInfoPath,"A_03_get_vst_normalized.txt"))

sessionInfo()

