library(pacman)
pacman::p_load(rmarkdown, tictoc)

# Main project Paths
projectDir <- 'D:/BIOINF/PROJECTS/2_ANALYSES/18_GIT_PAPERS/TCR_BCR_antiPD1/'
mainScriptsPath <- paste0("scripts/main/")

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
#####
tic("Total time")

# [1] Gather and merge anti-PD-L-1 datasets normalized counts, clinical information, perform immune cell deconvolution
tic("Gather and merge anti-PD-L-1 datasets normalized counts, clinical information, perform immune cell deconvolution",quiet = TRUE)
suppressMessages(rmarkdown::render(paste0(projectDir,mainScriptsPath,"B_01_aPD1_data_count_deconvolution_processing.RMD"),  quiet=T, clean = T,
                                   envir = new.env(),"html_document")) #4386.52 sec elapsed
remove_functions_from_global()
toc(log = TRUE)

# [2] Processing and analysis of TCR and BCR clones from MiXCR
tic("Processing and analysis of TCR and BCR clones from MiXCR",quiet = TRUE)
suppressMessages(rmarkdown::render(paste0(projectDir,mainScriptsPath,"B_02_aPD1_TCR_BCR_clones_analysis.RMD"),  quiet=T, clean = T,
                                   envir = new.env(),"html_document")) # 1306.16 sec elapsed
remove_functions_from_global()
toc(log = TRUE)

# [3] Merge all antiPD-1 data, clinical metadata, biomarkers
tic("Merge all antiPD-1 data, clinical metadata, biomarkers")
suppressMessages(source(paste0(projectDir,mainScriptsPath,"B_03_aPD1_merge_ready_data.R"), local = new.env())) # 
remove_functions_from_global()
toc(log = TRUE)

# End
toc(log = TRUE)# 

