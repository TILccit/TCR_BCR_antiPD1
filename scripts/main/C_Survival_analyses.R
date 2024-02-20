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

# [1] Uni and Multivariate survival analysis for TCGA data
tic("Uni and Multivariate survival analysis for TCGA data",quiet = TRUE)
suppressMessages(rmarkdown::render(paste0(projectDir,mainScriptsPath,"C_01_TCGA_survival_analysis.RMD"), quiet=T, clean = T,
                                   envir = new.env(),"html_document")) # 345.94 sec elapsed

remove_functions_from_global()
toc(log = TRUE)

# [2] Uni and Multivariate survival analysis for anti-PD(L)1 public and restricted data
tic("Uni and Multivariate survival analysis for anti-PD(L)1 public and restricted data",quiet = TRUE)
suppressMessages(rmarkdown::render(paste0(projectDir,mainScriptsPath,"C_02_aPD1_survival_analysis.RMD"), quiet=T, clean = T,
                                   envir = new.env(),"html_document")) # 347.59 sec elapsed
remove_functions_from_global()
toc(log = TRUE)

# [3] Univariable TCGA KM plot
tic("Univariable TCGA KM plots")
suppressMessages(source(paste0(projectDir,mainScriptsPath,"C_03_TCGA_survival_KMplots.R"), local = new.env())) 
remove_functions_from_global()
toc(log = TRUE)

# [4] Univariable anti-PD(L)1 datasets KM plots
tic("Univariable anti-PD(L)1 datasets KM plots")
suppressMessages(source(paste0(projectDir,mainScriptsPath,"C_04_aPD1_survival_KMplots.R"), local = new.env()))
remove_functions_from_global()
toc(log = TRUE)

# End
toc(log = TRUE)


