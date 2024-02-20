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
tic("Differential expression analysis on TCGA data",quiet = TRUE)
suppressMessages(rmarkdown::render(paste0(projectDir,mainScriptsPath,"D_01_DEA_GOenrich_TCGA.RMD"), quiet=T, clean = T,
                                   envir = new.env(),"html_document")) #4079.98 sec elapsed

remove_functions_from_global()
toc(log = TRUE)

# End
toc(log = TRUE)# 447.76 sec elapsed

