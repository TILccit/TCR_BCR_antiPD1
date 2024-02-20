#! /usr/bin/env Rscript

# Using RTCGA package to download RNAseq data that are included in RTCGA.rnaseq package
# Rna-seq data format is explained here https://wiki.nci.nih.gov/display/TCGA/RNASeq+Version+2

#####################
### Load packages ###
#####################
# Run below for installing RTCGA first
if (!require(devtools)) {
  install.packages("devtools")
  require(devtools)
}

if (!require(RTCGA)) {
  install_github("RTCGA/RTCGA", build_vignettes = TRUE)
  require(RTCGA)
}


library(pacman)
pacman::p_load(RTCGA, dplyr,stringr, tidyr, splitstackshape)#RTCGA.rnaseq,
###########################
### Main analysis paths ###
###########################
projectDir <- 'D:/BIOINF/PROJECTS/2_ANALYSES/18_GIT_PAPERS/TCR_BCR_antiPD1/'
setwd(projectDir)

projectDataPath <- paste0("data/")
tcgaInputData <- paste0(projectDir, projectDataPath,"TCGA/")
rsemPath <-paste0(tcgaInputData, "rsem")

# Session/dependencies
sessionInfoPath <- paste0("session_info/")

setwd(tcgaInputData)
######################################
## Create intermediate output paths ##
######################################
if (!dir.exists(rsemPath)) {dir.create(rsemPath)}


#################
### Functions ###
#################
get_rsem_table<-function(rsemObj){
  rsem_df <-as.data.frame(t(rsemObj))
  rsem_df[-1,]
}

###################
### DOWNLOADING ###
###################

# All cohort names can be checked using:
  
(cohorts <- infoTCGA() %>% 
     rownames() %>% 
     sub("-counts", "", x=.))
#For all cohorts the following code downloads the RNAseq data.

#Downloading RNAseq files
releaseDate <- "2015-11-01"
sapply( cohorts, function(element){
  tryCatch({
    downloadTCGA( cancerTypes = element, 
                  dataSet = "rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level",
                  destDir = rsemPath, 
                  date = releaseDate )},
    error = function(cond){
      cat("Error: Maybe there weren't rnaseq data for ", element, " cancer.\n")
    }
  )
})

# Change permissions for all files
Sys.chmod(list.files( rsemPath), "777", use_umask = FALSE)

# Reading downloaded RNAseq dataset
# Shortening paths and directories


list.files( "rsem") %>% 
  file.path( "rsem", .) %>%
  file.rename( to = substr(.,start=1,stop=55))
# Removing NA files from rsem folder
# If there were not RNAseq data for some cohorts we should remove corresponding NA files.

list.files( "rsem") %>%
  file.path( "rsem", .) %>%
  sapply(function(x){
    if (x == "rsem/NA")
      file.remove(x)      
  })
# Paths to RNAseq data
# Below is the code that removes unneeded "MANIFEST.txt" file from each RNAseq cohort folder.

list.files( "rsem") %>% 
  file.path( "rsem", .) %>%
  sapply(function(x){
    file.path(x, list.files(x)) %>%
      grep(pattern = "MANIFEST.txt", x = ., value=TRUE) %>%
      file.remove()
  })

####################
### Session info ###
####################

session_info <- sessionInfo()
writeLines(capture.output(session_info), paste0(projectDir,sessionInfoPath,"A_01_download_tcga_rsem.txt"))

sessionInfo()

