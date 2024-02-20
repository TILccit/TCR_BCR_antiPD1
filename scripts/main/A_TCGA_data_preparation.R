library(pacman)
pacman::p_load(rmarkdown, tictoc)

# Main project Paths
projectDir <- 'D:/BIOINF/PROJECTS/2_ANALYSES/18_GIT_PAPERS/TCR_BCR_antiPD1/'
# setwd(projectDir)
mainScriptsPath <- paste0("scripts/main/")

#####################

# # Get a list of objects in the global environment
# objects_in_global <- ls()
# 
# # Identify functions among the objects
# functions_in_global <- objects_in_global[sapply(objects_in_global, function(obj) typeof(get(obj)) == "closure")]
# 
# # Remove functions from the global environment
# rm(list = functions_in_global, envir = .GlobalEnv)

remove_functions_from_global <- function() {
  # Get a list of objects in the global environment
  objects_in_global <- ls(envir = .GlobalEnv)
  
  # Identify functions among the objects and remove them
  for (obj in objects_in_global) {
    if(obj=="remove_functions_from_global"){
      #skip
    }else{
      if (is.function(get(obj, envir = .GlobalEnv))) {
        rm(list = obj, envir = .GlobalEnv)
      }
    }
   
  }
}


#####################

tic("Total time")

# # # [1] Download TCGA RSEM data - No longer working
# # tic("Download TCGA RSEM data")
# # suppressMessages(source(paste0(projectDir,mainScriptsPath,"A_01_download_tcga_rsem.R"), local = new.env())) # 345.94 sec elapsed
# # toc(log = TRUE)
#
# # # [2] Download TCGA htseq FPKM-UQ, Counds data - No longer working due to changes in TCGAbiolinks package updates, data provided in TCGA input folder
# # tic("Download TCGA htseq Counts, FPKM-UQ data")
# # source(paste0(projectDir,mainScriptsPath,"A_02_download_tcga_htseq.R"), local = T)
# # remove_functions_from_global()
# # toc(log = TRUE)
#
# [3] Get TCGA vst normalized counts
tic("Get TCGA vst normalized counts")
suppressMessages(source(paste0(projectDir,mainScriptsPath,"A_03_get_vst_normalized.R"), local = new.env())) # 422.79 sec elapsed
remove_functions_from_global()
toc(log = TRUE)
#
# # # [4] Get MSI data - No longer working due to archiving of GDC legacy data
# # tic("Get MSI status data")
# # source(paste0(projectDir,mainScriptsPath,"A_04_get_TCGA_msi_status_data.R"), local = T) # 186.72 sec elapsed
# # remove_functions_from_global()
# # toc(log = TRUE)
#
# [5] Get biospecimen data
tic("Get biospecimen data")
suppressMessages(rmarkdown::render(paste0(projectDir,mainScriptsPath,"A_05_get_TCGA_biospecimen_RNA.RMD"), quiet=T, clean = T,
                  envir = new.env(),"html_document"))# 269.51 sec elapsed
remove_functions_from_global()
toc(log = TRUE)

# rm(list=ls(all=TRUE))
# projectDir <- 'D:/BIOINF/PROJECTS/2_ANALYSES/18_GIT_PAPERS/TCR_BCR_antiPD1/'
# [6] Gather and process clinical and biomarker metadata for TCGA
tic("Gather and process clinical and biomarker metadata for TCGA",quiet = TRUE)
suppressMessages(rmarkdown::render(paste0(projectDir,mainScriptsPath,"A_06_tcga_clinical_biomarker_metadata.RMD"), quiet=T, clean = T,
                  envir = new.env(),"html_document")) #2421.66 sec elapsed
remove_functions_from_global()
toc(log = TRUE)

# End
toc(log = TRUE)
#Total time: 2986.96 sec elapsed

