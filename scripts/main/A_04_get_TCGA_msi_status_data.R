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
pacman::p_load(TCGAbiolinks,jsonlite,dplyr,plyr)
###########################
### Main analysis paths ###
###########################
projectDir <- 'D:/BIOINF/PROJECTS/2_ANALYSES/18_GIT_PAPERS/TCR_BCR_antiPD1/'
setwd(projectDir)

projectDataPath <- paste0("data/")
tcgaInputData <- paste0(projectDir, projectDataPath,"TCGA/")

# Session/dependencies
sessionInfoPath <- paste0("session_info/")

# setwd(tcgaInputData)

#################
### Functions ###
#################

getProjectSummary <- function(project, legacy = FALSE){
  baseURL <- ifelse(legacy,"https://api.gdc.cancer.gov/legacy/projects/","https://api.gdc.cancer.gov/projects/")
  url <- paste0(baseURL, project,"?expand=summary,summary.data_categories&pretty=true")
  return(fromJSON(url,simplifyDataFrame = TRUE)$data$summary)
}
checkDataCategoriesInput <- function(project,data.category, legacy = FALSE){
  for(proj in project){
    project.summary <- getProjectSummary(proj, legacy)
    # if(missing(data.category)) {
    #   #print(knitr::kable(project.summary$data_categories))
    #   print("Please set a data.category argument from the column data_category above")
    # }
    # if(!(data.category %in% project.summary$data_categories$data_category)) {
    #   #print(knitr::kable(project.summary$data_categories))
    #   print("Please set a valid data.category argument from the column data_category above. We could not validade the data.category for project ", proj)
    # }
    return(project.summary$data_categories$data_category)
  }
}
# TCGAbiolinks gives an error for some other tumor types that there are no data for msi. So I use the above TCGAbiolinks
# functions (little bit transformed) to extract msi data from the projects that have these information.
get_msi_tcga <- function(selectedcancer){
  print(paste0("Getting msi data for ", selectedcancer ))
  # Check if msi data exist for this project
  categorycheck <-checkDataCategoriesInput(selectedcancer,"Other",legacy = TRUE)
  if ("Other" %in% categorycheck){
    # Run query
    query <- GDCquery(project = selectedcancer, 
                      data.category = "Other",
                      legacy = TRUE,
                      access = "open",
                      data.type = "Auxiliary test"
    )
    if (is.null(query)){
      print(paste0("No msi data for ", selectedcancer))
    }else{
      GDCdownload(query)
      msi_results <- GDCprepare_clinic(query, "msi")
      
      msi_tcga <- rbind(msi_tcga,msi_results)
      msi_tcga
    }
  }
  
  
}
######################################### 
## Download MSI data from TCGAbiolinks ##
#########################################

tcga_projects <-TCGAbiolinks:::getGDCprojects()$project_id
mycancertypes <- tcga_projects[grepl("^TCGA-",tcga_projects)] # total 33
mycancertypes <- mycancertypes[order(mycancertypes)]

# Empty dataframe for msi data

msi_tcga <- data.frame(bcr_patient_barcode=character(),
                       bcr_aliquot_uuid = character(),
                       mononucleotide_and_dinucleotide_marker_panel_analysis_status = character(),
                       mononucleotide_marker_panel_analysis_status = character(),
                       project = character(),
                       stringsAsFactors=FALSE)

# Get msi data
msi_tcga_all <- sapply(mycancertypes,function(x) get_msi_tcga(selectedcancer =x),simplify = FALSE)
# Remove null elements- projects with no data for msi
msi_tcga_all_f <-compact(msi_tcga_all)
# Remove projects with no data for msi
msi_tcga_all_f <-msi_tcga_all_f[lengths(msi_tcga_all_f) != 1]
msi_tcga_all_f <-bind_rows(msi_tcga_all_f , .id = "column_label")
# Save data
save(msi_tcga_all_f,file=paste0(tcgaInputData,"msi_tcga_all.RData"))

####################
### Session info ###
####################

session_info <- sessionInfo()
writeLines(capture.output(session_info), paste0(projectDir,sessionInfoPath,"A_05_get_TCGA_biospecimen_RNA.txt"))

sessionInfo()
