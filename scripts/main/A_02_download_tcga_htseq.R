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
pacman::p_load(TCGAbiolinks,SummarizedExperiment, parallel)
###########################
### Main analysis paths ###
###########################
projectDir <- 'D:/BIOINF/PROJECTS/2_ANALYSES/18_GIT_PAPERS/TCR_BCR_antiPD1/'
setwd(projectDir)

projectDataPath <- paste0("data/")
tcgaInputData <- paste0(projectDir, projectDataPath,"TCGA/")
htseqFPKM_UQ_Path <-paste0(tcgaInputData, "htseq_fpkm_uq")
htseqCounts_Path <-paste0(tcgaInputData, "htseq_vst")

# Session/dependencies
sessionInfoPath <- paste0("session_info/")

# setwd(tcgaInputData)
######################################
## Create intermediate output paths ##
######################################
if (!dir.exists(htseqFPKM_UQ_Path)) {dir.create(htseqFPKM_UQ_Path)}
if (!dir.exists(htseqCounts_Path)) {dir.create(htseqCounts_Path)}


#################
### Functions ###
#################
download_htseq_data <- function(selectedcancer,type="HTSeq - Counts",name="Counts"){ 
  # type-- "HTSeq - Counts", "HTSeq - FPKM", "HTSeq - FPKM-UQ" | name---"Counts", "FPKM", "FPKM_UQ"
  
  # selectedcancer <- mycancertypes[1]
  # type="HTSeq - Counts"
  # name="Counts"
  # print(selectedcancer)
  
  query <- GDCquery(#project = paste0('TCGA-',selectedcancer),
    project = selectedcancer,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    sample.type = "Primary Tumor",
    workflow.type = type)
  
  GDCdownload(query)
  data <- GDCprepare(query)
  clin <- colData(data)
  df_all <- assay(data)
  
  save(df_all,clin,file=paste0(tcga_dir, selectedcancer,"_RNASeq",name,".RData"))
}

###################
### DOWNLOADING ###
###################
# Get TCGA cohort names
tcga_projects <-TCGAbiolinks:::getGDCprojects()$project_id
mycancertypes <- tcga_projects[grepl("^TCGA-",tcga_projects)] # total 33

#####################
### HTSEQ FPKM UQ ###
#####################
# Set dir
tcga_dir <- htseqCounts_Path

# Set up parallel cores
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
# all_var<-ls()
# clusterExport(cl, all_var)
clusterEvalQ(cl, c(library(TCGAbiolinks),library(SummarizedExperiment)))
# Export the play() function to the cluster
clusterExport(cl,c("download_htseq_data","tcga_dir"))

download_all <- parSapply(cl,mycancertypes,function(x) download_htseq_data(selectedcancer=x,type="HTSeq - FPKM-UQ",name="FPKM_UQ"))
stopCluster(cl)

#####################
### HTSEQ Counts ###
#####################
# Set dir
tcga_dir <- htseqFPKM_UQ_Path

# Set up parallel cores
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
# all_var<-ls()
# clusterExport(cl, all_var)
clusterEvalQ(cl, c(library(TCGAbiolinks),library(SummarizedExperiment)))
# Export the play() function to the cluster
clusterExport(cl,c("download_htseq_data","tcga_dir"))

download_all <- parSapply(cl,mycancertypes,function(x) download_htseq_data(selectedcancer=x,type="HTSeq - Counts",name="Counts"))
stopCluster(cl)

####################
### Session info ###
####################

session_info <- sessionInfo()
writeLines(capture.output(session_info), paste0(projectDir,sessionInfoPath,"A_02_download_tcga_htseq.txt"))

sessionInfo()

