#######################
### Functions used ###
######################

###############################################
##  ##
###############################################
# Claculate FPKM-UQ normalized values from counts (htseq)
fpkm_uq <- function(counts, lengths) {
  rate <- counts / lengths
  uq <- quantile(counts)
  g_uq <-as.numeric(uq[4])
  #rate / g_uq * 1e6
  fpkm_uq <-(rate / g_uq)* 1e7
  fpkm_uq
}


##########################
## Loading RData object ##
##########################


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


get_patient_num <- function(sample_id){
  patient_num <- str_extract(sample_id, "[^_]+")
  patient_num
}

select_columns <- function(DF,columnNames){
  DF.s <- DF %>% dplyr::select(columnNames)
  colnames(DF.s) <-columnNames
  DF.s
}
############################
## Uppercase first letter ##
############################
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
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

####################
### TCR functions ##
####################

# Add chain information in dataframe
add_chain_info <- function(dataframe){
  chain_info<-dataframe[4]
  chain_info<-substr(chain_info, start = 1, stop = 3)
  if (chain_info == "TRA") {
    dataframe[9] <- "alpha"
  } else if (chain_info == "TRB") {
    dataframe[9] <- "beta"
  } else if (chain_info == "TRG") {
    dataframe[9] <- "gamma"
  } else if (chain_info == "TRD") {
    dataframe[9] <- "delta"
  }else {
    dataframe[9] <- "immunoglobulin"
  }
  dataframe
}
##    read in and do preliminary data cleaning
readAndCleanCDR3 <- function(dataFrame){
  # print(levels(as.factor(dataFrame$dataset)))
  cdr3 <- dataFrame
  print(paste("There are", nrow(cdr3), "total CDR3s"))
  nonprodCDR3.1 <- cdr3[grepl("~",cdr3$aaSeqCDR3),]
  nonprodCDR3.2 <- cdr3[grepl("[*]",cdr3$aaSeqCDR3),]
  print(paste("There are", nrow(nonprodCDR3.1), "non productive CDR3s (~)"))
  print(paste("There are", nrow(nonprodCDR3.2), "non productive CDR3s ([*])"))
  print(paste("Example:"))
  nonprodCDR3.2[1:3,]
  cdr3 <- cdr3[!grepl("~",cdr3$aaSeqCDR3),]
  cdr3 <- cdr3[!grepl("[*]",cdr3$aaSeqCDR3),]
  print(paste("There are", nrow(cdr3), "productive CDR3s"))
  
  ## There are cases where then same CDR3 aaSeq was determined to be both alpha and beta (due to ambiguous alignemnts)
  ## Need to choose the chain which has more read support.
  print("Fixing cases where single CDR3 called as both alpha and beta.")
  arg_list <- by(cdr3,cdr3[,c("sample","aaSeqCDR3")], function(x){
    df <- as.data.frame(x)
    df$cloneCount<- as.numeric(df$cloneCount)
    
    ans <- df[df$cloneCount == max(df$cloneCount),]
    if(nrow(ans) > 1){
      ## assign unknown
      ans <- ans[1,]
      ans$chain <- "Ambiguous"
      return(ans)
    }else{
      return(ans)
    }
  })
  cdr3 <- do.call(rbind,arg_list)
  print(paste("There are", nrow(cdr3), "unique CDR3s"))
  print(paste("There are", nrow(subset(cdr3,cdr3$chain == "Ambiguous")), "ambiguous CDR3s"))
  
  return(cdr3)
}

##  read in and do preliminary data cleaning on sample info
readAndCleanSample <- function(depths, cdr3){
  print(unique(depths[["dataset"]]))
  # sampleTable<-sampleTablesReads.f$Mari
  # cdr3<-tcrTable$Mari
  # depths <- as.data.table(sampleTable)
  ## get total TCR reads for each sample.
  depths$totTCR <- apply(depths, 1, function(x){
    #print(x)
    return(sum(cdr3$cloneCount[cdr3$sample==x[["run_accession"]] & cdr3$chain %in% c("alpha","beta")]))
  })
  depths$totBCR <- apply(depths, 1, function(x){
    #print(x)
    return(sum(cdr3$cloneCount[cdr3$sample==x[["run_accession"]] & cdr3$chain %in% c("immunoglobulin")]))
  })
  depths$totTCRa <- apply(depths, 1, function(x){
    return(sum(cdr3$cloneCount[cdr3$sample==x[["run_accession"]] & cdr3$chain=="alpha"]))
  })
  
  depths$totTCRb <- apply(depths, 1,function(x){
    return(sum(cdr3$cloneCount[cdr3$sample==x[["run_accession"]] & cdr3$chain=="beta"]))
  })
  
  depths$totIGH <- apply(depths, 1,function(x){
    return(sum(cdr3$cloneCount[cdr3$sample==x[["run_accession"]] & cdr3$chain=="immunoglobulin"]))
  })
  
  depths$sampleFracTCR <- depths$totTCR/depths$seqDepth
  depths$sampleFracBCR <- depths$totBCR/depths$seqDepth
  
  #depths$response <- sampleTable[match(sampleTable$run_accession, cdr3$sample),4]
  return(depths)
  
}

#apply(samples_hugo.list,1, function(x){print(x)})
##  get fraction of TCRs & BCRs that each clone is responsible for.
getCloneFraction <- function(cdr3, depths){
  # cdr3 <- tcrTable$Hugo
  # depths <- tcrTable.c$Hugo
  # # TCR
  # x <-cdr3[1,]
  cdr3$pcTotTCR <- apply(cdr3, 1, function(x){
    # return(as.numeric(x[["cloneCount"]])/depths$totTCR[depths$run_accession==x[["sample"]]])
    
    if(is.na(x[["chain"]])){
      return(NA)
    }else if(x[["chain"]] %in% c("alpha","beta")){
      return(as.numeric(x[["cloneCount"]])/depths$totTCR[depths$run_accession==x[["sample"]]])
    }else if(x[["chain"]] == "immunoglobulin"){
      return(NA)
    }else{
      ## Ambiguous chain
      return(as.numeric(x[["cloneCount"]])/depths$totTCR[depths$run_accession==x[["sample"]]])
    }
    
    
  })
  # BCR
  cdr3$pcTotBCR <- apply(cdr3, 1, function(x){
    # return(as.numeric(x[["cloneCount"]])/depths$totBCR[depths$run_accession==x[["sample"]] ])
    if(is.na(x[["chain"]])){
      return(NA)
    }else if(x[["chain"]]== "immunoglobulin" ){
      return(as.numeric(x[["cloneCount"]])/depths$totBCR[depths$run_accession==x[["sample"]] ])
    }else if(x[["chain"]] %in% c("alpha","beta")){
      return(NA)
    }else{
      ## Ambiguous chain
      return(as.numeric(x[["cloneCount"]])/depths$totBCR[depths$run_accession==x[["sample"]] ])
    }
  })
  
  # TCR
  cdr3$pcTotTCRchain <- apply(cdr3, 1, function(x){
    if(is.na(x[["chain"]])){
      return(NA)
    }else if(x[["chain"]] == "alpha"){
      return(as.numeric(x[["cloneCount"]])/depths$totTCRa[depths$run_accession==x[["sample"]]])
    }else if(x[["chain"]] == "beta"){
      return(as.numeric(x[["cloneCount"]])/depths$totTCRb[depths$run_accession==x[["sample"]]])
    }else if(x[["chain"]] == "immunoglobulin"){
      return(NA)
    }else{
      ## Ambiguous chain
      return(as.numeric(x[["cloneCount"]])/depths$totTCR[depths$run_accession==x[["sample"]]])
    }
  })
  
  # BCR
  cdr3$pcTotBCRchain <- apply(cdr3, 1, function(x){
    if(is.na(x[["chain"]])){
      return(NA)
    }else if(x[["chain"]] == "immunoglobulin"){
      return(as.numeric(x[["cloneCount"]])/depths$totIGH[depths$run_accession==x[["sample"]]])
    }else if(x[["chain"]] == "alpha"){
      return(NA)
    }else if(x[["chain"]] == "beta"){
      return(NA)
    }else{
      ## Ambiguous chain
      return(as.numeric(x[["cloneCount"]])/depths$totBCR[depths$run_accession==x[["sample"]]])
    }
  })
  
  # x <-cdr3[1,]
  
  cdr3$pcTotTCRReads <- apply(cdr3, 1, function(x){
    # return(as.numeric(x[["cloneCount"]])/depths$seqDepth[depths$run_accession==x[["sample"]]])
    
    if(x[["chain"]] %in% c("alpha","beta")){
      return(as.numeric(x[["cloneCount"]])/depths$seqDepth[depths$run_accession==x[["sample"]]])
    }else if(x[["chain"]] == "immunoglobulin"){
      return(NA)
    }else{
      ## Ambiguous chain
      return(as.numeric(x[["cloneCount"]])/depths$seqDepth[depths$run_accession==x[["sample"]]])
    }
    
  })
  
  cdr3$pcTotBCRReads <- apply(cdr3, 1, function(x){
    
    
    if(x[["chain"]]== "immunoglobulin" ){
      return(as.numeric(x[["cloneCount"]])/depths$seqDepth[depths$run_accession==x[["sample"]]])
    }else if(x[["chain"]] %in% c("alpha","beta")){
      return(NA)
    }else{
      ## Ambiguous chain
      return(as.numeric(x[["cloneCount"]])/depths$seqDepth[depths$run_accession==x[["sample"]]])
    }
    
  })
  
  # Merge two datasets
  cdr3 <- merge(depths, cdr3, by.x = "run_accession", by.y = "sample", all=TRUE)
  return(cdr3)
}


## classify clones as dominant or not
getDominantClones <- function(cdr3){
  
  # cdr3 <- tcrTable.cCloneFr$Hugo
  
  #thresh <- max(log10(datToUse$pcTotTCRchain*datToUse$pcTotReads))
  #threshAlpha <- max(log10(datToUse$pcTotTCRchain[datToUse$chain=="alpha"]*datToUse$pcTotReads[datToUse$chain=="alpha"]))
  #threshBeta <- max(log10(datToUse$pcTotTCRchain[datToUse$chain=="beta"]*datToUse$pcTotReads[datToUse$chain=="beta"]))
  
  ## metric = (pcTotReads x pcTotTCRchain)
  ##        = (abundance/seqDepth x abundance/totTCR[chain])
  ##        = (abundance)^2/(seqDepth x totTCR[chain])
  
  ## TCR ##
  #########
  cdr3$TCRdominantMetric <- (cdr3$pcTotTCRReads*cdr3$pcTotTCRchain)
  # cdr3$dominantMetricA <- (cdr3$pcTotReads*cdr3$pcTotTCRchain)
  # cdr3$dominantMetricB <- (cdr3$pcTotReads*cdr3$pcTotTCRchain)
  ## Use control data to measure background rate, and classify CDR3s as dominant or not.
  #datToUse <- subset(cdr3, grepl("Ctrl", sample))
  # ctrlAlphaValues <- datToUse$dominantMetric[datToUse$chain=="alpha"]
  # threshAlpha <- max(ctrlAlphaValues)
  # # ctrlBetaValues <- datToUse$dominantMetric[datToUse$chain=="beta"]
  # threshBeta <- max(ctrlBetaValues)
  threshAlpha <- median(cdr3$TCRdominantMetric[cdr3$chain=="alpha"], na.rm = TRUE)
  threshBeta <- median(cdr3$TCRdominantMetric[cdr3$chain=="beta"], na.rm = TRUE)
  #print(paste("Threshold being used is ",thresh,sep=""))
  print(paste("Threshold for alpha is ",threshAlpha,sep=""))
  print(paste("Threshold for beta is ",threshBeta,sep=""))
  
  cdr3$TCRdominant <- FALSE
  cdr3$TCRdominant[which(cdr3$chain =="alpha")] <- cdr3$TCRdominantMetric[cdr3$chain=="alpha"] > threshAlpha
  cdr3$TCRdominant[which(cdr3$chain =="beta")] <- cdr3$TCRdominantMetric[cdr3$chain=="beta"] > threshBeta
  
  ## BCR ##
  #########
  cdr3$BCRdominantMetric <- (cdr3$pcTotBCRReads*cdr3$pcTotBCRchain)
  # cdr3$dominantMetricA <- (cdr3$pcTotReads*cdr3$pcTotTCRchain)
  # cdr3$dominantMetricB <- (cdr3$pcTotReads*cdr3$pcTotTCRchain)
  ## Use control data to measure background rate, and classify CDR3s as dominant or not.
  #datToUse <- subset(cdr3, grepl("Ctrl", sample))
  # ctrlAlphaValues <- datToUse$dominantMetric[datToUse$chain=="alpha"]
  # threshAlpha <- max(ctrlAlphaValues)
  # # ctrlBetaValues <- datToUse$dominantMetric[datToUse$chain=="beta"]
  # threshBeta <- max(ctrlBetaValues)
  threshIGH <- median(cdr3$BCRdominantMetric[cdr3$chain=="immunoglobulin"], na.rm = TRUE)
  # threshBeta <- median(cdr3$TCRdominantMetric[cdr3$chain=="beta"], na.rm = TRUE)
  #print(paste("Threshold being used is ",thresh,sep=""))
  print(paste("Threshold for IGH is ",threshIGH,sep=""))
  # print(paste("Threshold for beta is ",threshBeta,sep=""))
  
  cdr3$BCRdominant <- FALSE
  cdr3$BCRdominant[which(cdr3$chain =="immunoglobulin")] <- cdr3$BCRdominantMetric[cdr3$chain=="immunoglobulin"] > threshIGH
  # cdr3$BCRdominant[which(cdr3$chain =="beta")] <- cdr3$BCRdominantMetric[cdr3$chain=="beta"] > threshBeta
  
  #cdr3$dominant[cdr3$chain=="Ambiguous"] <- FALSE #implicit
  return(cdr3)
}

## calculate shannon entropy for each sample
calcEntropy <- function(cdr3, depths){
  # cdr3 <- tcrTable.f$Hugo
  # depths <- tcrTable$Hugo
  colnames(cdr3)[1] <- "sample"
  ## Calculate shannon entropy for each sample
  
  # TCR
  # x <- depths[1,]
  depths$TCRshannon <- apply(depths, 1, function(x){
    # return(entropy::entropy(cdr3$cloneCount[cdr3$sample==x[["sample"]]  & cdr3$chain %in% c("alpha","beta")]))
    
    if(x[["chain"]] %in% c("alpha","beta")){
      return(entropy::entropy(cdr3$cloneCount[cdr3$sample==x[["sample"]]  & cdr3$chain %in% c("alpha","beta")]))
    }else if(x[["chain"]] == "immunoglobulin"){
      return(NA)
    }else{
      ## Ambiguous chain
      return(entropy::entropy(cdr3$cloneCount[cdr3$sample==x[["sample"]]  & cdr3$chain %in% c("alpha","beta")]))
    }
    
  })
  
  depths$TCRnumClones <- apply(depths, 1, function(x){
    
    
    if(x[["chain"]] %in% c("alpha","beta")){
      return(nrow(subset(cdr3, sample==x[["sample"]] & cdr3$chain %in% c("alpha","beta"))))
    }else if(x[["chain"]] == "immunoglobulin"){
      return(NA)
    }else{
      ## Ambiguous chain
      return(nrow(subset(cdr3, sample==x[["sample"]] & cdr3$chain %in% c("alpha","beta"))))
    }
    
  })
  depths$TCRnormShannon <- depths$TCRshannon / log2(depths$TCRnumClones)
  
  depths$TCRshannon2 <- apply(depths, 1, function(x){
    # return(diversity(cdr3$cloneCount[cdr3$sample==x[["sample"]] & cdr3$chain %in% c("alpha","beta")],index = "shannon", MARGIN = 1, base = exp(1)))
    
    if(x[["chain"]] %in% c("alpha","beta")){
      return(diversity(cdr3$cloneCount[cdr3$sample==x[["sample"]] & cdr3$chain %in% c("alpha","beta")],index = "shannon", MARGIN = 1, base = exp(1)))
    }else if(x[["chain"]] == "immunoglobulin"){
      return(NA)
    }else{
      ## Ambiguous chain
      return(diversity(cdr3$cloneCount[cdr3$sample==x[["sample"]] & cdr3$chain %in% c("alpha","beta")],index = "shannon", MARGIN = 1, base = exp(1)))
    }
    
  })
  
  # BCR
  
  depths$BCRshannon <- apply(depths, 1, function(x){
    # return(entropy::entropy(cdr3$cloneCount[cdr3$sample==x[["sample"]]  & cdr3$chain %in% c("alpha","beta")]))
    
    if(x[["chain"]] == "immunoglobulin"){
      return(entropy::entropy(cdr3$cloneCount[cdr3$sample==x[["sample"]]  & cdr3$chain %in% c("alpha","beta")]))
    }else if(x[["chain"]]%in% c("alpha","beta") ){
      return(NA)
    }else{
      ## Ambiguous chain
      return(entropy::entropy(cdr3$cloneCount[cdr3$sample==x[["sample"]]  & cdr3$chain %in% c("alpha","beta")]))
    }
    
  })
  
  depths$BCRnumClones <- apply(depths, 1, function(x){
    
    
    if(x[["chain"]]== "immunoglobulin"){
      return(nrow(subset(cdr3, sample==x[["sample"]] & cdr3$chain %in% c("alpha","beta"))))
    }else if(x[["chain"]] %in% c("alpha","beta")){
      return(NA)
    }else{
      ## Ambiguous chain
      return(nrow(subset(cdr3, sample==x[["sample"]] & cdr3$chain %in% c("alpha","beta"))))
    }
    
  })
  depths$BCRnormShannon <- depths$BCRshannon / log2(depths$BCRnumClones)
  
  depths$BCRshannon2 <- apply(depths, 1, function(x){
    # return(diversity(cdr3$cloneCount[cdr3$sample==x[["sample"]] & cdr3$chain %in% c("alpha","beta")],index = "shannon", MARGIN = 1, base = exp(1)))
    
    if(x[["chain"]] == "immunoglobulin"){
      return(diversity(cdr3$cloneCount[cdr3$sample==x[["sample"]] & cdr3$chain %in% c("alpha","beta")],index = "shannon", MARGIN = 1, base = exp(1)))
    }else if(x[["chain"]] %in% c("alpha","beta")){
      return(NA)
    }else{
      ## Ambiguous chain
      return(diversity(cdr3$cloneCount[cdr3$sample==x[["sample"]] & cdr3$chain %in% c("alpha","beta")],index = "shannon", MARGIN = 1, base = exp(1)))
    }
    
  })
  
  
  return(depths)
}


## estimate tumor purity using relative abundance of dominant clonotype.
estimatePurity <- function(cdr3, depths){
  
  depths$estPurity <- apply(depths, 1, function(x){
    
    ## get all dominant clones for this sample
    doms <- subset(cdr3, sample==x[["sample"]] & dominant)
    
    ## take abundance of max beta if present, otherwise alpha.
    
    if("beta" %in% doms$chain){
      return(max(doms$pcTotTCRchain[doms$chain=="beta"]))
    }else{
      return(max(doms$pcTotTCRchain))
    }
    
  })
  
  depths$estPurity[is.na(depths$estPurity) | depths$estPurity==-Inf] <- 0
  
  depths$hasDominant <- apply(depths, 1, function(x){
    return(TRUE %in% cdr3$dominant[cdr3$sample==x[["sample"]]])
  })
  
  return(depths)
}

####################################
## Process sample tables function ##
####################################
process_sample_table <- function(sampleDF,colName,ColElement){
  # Prepare variables
  colName.un <- as.name(colName)
  
  # Transform table
  sampleDF <- sampleDF %>% 
    mutate_at(vars(!!colName.un),as.character) %>% # if factor make it char vector to be able to change its elements
    mutate(!!colName:=ifelse(!!colName.un=="CR","CR+PR",ifelse(!!colName.un=="PR","CR+PR",ifelse(!!colName.un==ColElement,"CR+PR",!!colName.un)))) %>% # Merge CR with PR
    filter(!!colName.un %in% c("PD","CR+PR")) %>% # Select only PD, CR,PR
    mutate(!!colName:=factor(!!colName.un))%>% 
    droplevels() # Make response column a factor
  # Return table
  sampleDF
}


#########################################
## Collect MiXCR results for a dataset ##
#########################################

grabMiXcrDataset<- function(parentFolder,foldersList,datasetLabel,mixcrPath,sampleDF){
  # parentFolder<-public_aPD1Data_ParentPath
  # foldersList<-data.foldersList
  # datasetLabel<-"Hugo"
  # mixcrPath <-pathToMixcr
  # sampleDF <-survTissues.merged\samples_all.merged
  # 
  print(datasetLabel)
  # Subset to dataset defined
  sampleDF.d <- sampleDF %>% dplyr::filter(dataset==datasetLabel)
  # Extract vector of sample run accessions
  sampleRunAccession<- as.character(sampleDF.d$run_accession)
  if(datasetLabel=="Immotion150"){
    print("Fixing run accession for Immotion 150")
    sampleRunAccession<- str_after_nth(sampleRunAccession, "_", 7)
  }
  # Collect results from Mixcr in table
  datasetMixcr <-collect_mixcr_results_table(sampleRunAccession,paste0(parentFolder,foldersList[[datasetLabel]],"/",mixcrPath))
}

collect_mixcr_results_table <- function(sampleNamesVector,pathToResults){
  clones_results <- list()
  for (i in sampleNamesVector){
    clonesFile <- paste0(pathToResults,i, ".clonotypes.ALL.txt")
    
    # clonesFileTCRA <- paste0(pathToResults,i, ".clonotypes.TRA.txt")
    # clonesFileTCRB <- paste0(pathToResults,i, ".clonotypes.TRB.txt")
    # clonesFile <- rbind(clonesFileTCRA,clonesFileTCRB)
    
    #emFile <- paste0("C:/Users/Emilia/Desktop/SpecialCourse/Analysis/Bedtools/", i, ".exon_mean")
    #print("Reading clonotypes file for sample ", i ,".txt")
    # print(paste0("Reading clonotypes file for sample ", clonesFile))
    
    # Read the data if the file is non-empty.
    if(is.na(file.size(clonesFile))){
      print(paste0("File does not exist:", clonesFile))
    } else{
      if (!file.size(clonesFile) == 0){
        clones_results[[i]] <- read.delim(file= paste0(
          pathToResults,i, ".clonotypes.ALL.txt"), header=T)
        
        clones_resultsTRA <- read.delim(file= paste0(
          pathToResults,i, ".clonotypes.TRA.txt"), header=T)
        clones_resultsTRB <- read.delim(file= paste0(
          pathToResults,i, ".clonotypes.TRB.txt"), header=T)
        ## NEW
        clones_resultsIGH <- read.delim(file= paste0(
          pathToResults,i, ".clonotypes.IGH.txt"), header=T)
        ######
        clones_results[[i]] <-rbind(clones_resultsTRA,clones_resultsTRB,clones_resultsIGH)
        #colnames(clones_results_hugo[[i]]) <- c("refSeq", "length", "MappedReads", "UnmappedReads")
      }
    }
  }
  clones_results
}


grabMiXcrDatasetTCRA<- function(parentFolder,foldersList,datasetLabel,mixcrPath,sampleDF){
  # parentFolder<-public_aPD1Data_ParentPath
  # foldersList<-data.foldersList
  # datasetLabel<-"Hugo"
  # mixcrPath <-pathToMixcr
  # sampleDF <-survTissues.merged\samples_all.merged
  # 
  print(datasetLabel)
  # Subset to dataset defined
  sampleDF.d <- sampleDF %>% dplyr::filter(dataset==datasetLabel)
  # Extract vector of sample run accessions
  sampleRunAccession<- as.character(sampleDF.d$run_accession)
  if(datasetLabel=="Immotion150"){
    print("Fixing run accession for Immotion 150")
    sampleRunAccession<- str_after_nth(sampleRunAccession, "_", 7)
  }
  # Collect results from Mixcr in table
  datasetMixcr <-collect_mixcr_results_tableTCRA(sampleRunAccession,paste0(parentFolder,foldersList[[datasetLabel]],"/",mixcrPath))
}

collect_mixcr_results_tableTCRA <- function(sampleNamesVector,pathToResults){
  clones_results <- list()
  for (i in sampleNamesVector){
    clonesFile <- paste0(pathToResults,i, ".clonotypes.TRA.txt")
    
    # clonesFileTCRA <- paste0(pathToResults,i, ".clonotypes.TRA.txt")
    # clonesFileTCRB <- paste0(pathToResults,i, ".clonotypes.TRB.txt")
    # clonesFile <- rbind(clonesFileTCRA,clonesFileTCRB)
    
    #emFile <- paste0("C:/Users/Emilia/Desktop/SpecialCourse/Analysis/Bedtools/", i, ".exon_mean")
    #print("Reading clonotypes file for sample ", i ,".txt")
    # print(paste0("Reading clonotypes file for sample ", clonesFile))
    
    # Read the data if the file is non-empty.
    if(is.na(file.size(clonesFile))){
      print(paste0("File does not exist:", clonesFile))
    } else{
      if (!file.size(clonesFile) == 0){
        clones_results[[i]] <- read.delim(file= paste0(
          pathToResults,i, ".clonotypes.TRA.txt"), header=T)
        
        # clones_resultsTRA <- read.delim(file= paste0(
        #   pathToResults,i, ".clonotypes.TRA.txt"), header=T)
        # clones_resultsTRB <- read.delim(file= paste0(
        #   pathToResults,i, ".clonotypes.TRB.txt"), header=T)
        
        # clones_results[[i]] <-rbind(clones_resultsTRA,clones_resultsTRB)
        #colnames(clones_results_hugo[[i]]) <- c("refSeq", "length", "MappedReads", "UnmappedReads")
      }
    }
  }
  clones_results
}

final_clonesTable <- function(clonesDFlist,chainInfoColName,chainColName){
  # clonesDFlist <-mixcr.collect$Hugo
  # chainInfoColName <-"allVHitsWithScore"
  # chainColName<-"chain"
  
  
  # Prepare variables
  chainInfoColName.un <- as.name(chainInfoColName)
  chainColName.un <- as.name(chainColName)
  #colName.un <- as.name(colName)
  # Get columns: 1,2,3,6,24,33
  clonesDFlist <- lapply(clonesDFlist, function(x){x[,c(1,2,3,6,24,33)]})
  
  # # Add a column of sample info in each dataframe of the list
  # clonesDFlist2 <- mapply(cbind, clonesDFlist, "sample"=spls_hugo, SIMPLIFY=F)
  # Unlist list into dataframe
  clonesDF <- as.data.frame(ldply(clonesDFlist, data.frame))
  clonesDF_final <- clonesDF %>% 
    mutate_at(vars(!!chainInfoColName.un),as.character) %>% 
    dplyr::mutate(!!chainColName := ifelse(grepl("^TRA",substr(!!chainInfoColName.un, start = 1, stop = 3)),"alpha",
                                           ifelse(grepl("^TRB",substr(!!chainInfoColName.un, start = 1, stop = 3)),"beta",
                                                  ifelse(grepl("^TRG",substr(!!chainInfoColName.un, start = 1, stop = 3)),"gamma",
                                                         ifelse(grepl("^TRD",substr(!!chainInfoColName.un, start = 1, stop = 3)),"delta",
                                                                ifelse(grepl("^IGH",substr(!!chainInfoColName.un, start = 1, stop = 3)),"immunoglobulin","other")))))) %>%
    dplyr::filter(!!chainColName.un %in% c("alpha","beta","immunoglobulin")) %>% #  SELECTING ONLY TCRS alpha, beta
    dplyr::mutate(!!chainColName :=factor(!!chainColName.un)) %>% droplevels()
  clonesDF_final
}

totalReadsLoad <- function(parentFolder,foldersList,datasetLabel,dataPath,sampleDFList){
  # parentFolder<-public_aPD1Data_ParentPath
  # foldersList<-data.foldersList
  # datasetLabel<-"Hugo"
  # dataPath <-pathToSampleTables
  # sampleDFList <-samples_all.sub
  
  ##
  sampleDFdataset <- sampleDFList[[datasetLabel]]
  readsTable <-as.data.table(read.table(paste0(parentFolder,foldersList[[datasetLabel]],"/",dataPath,"counts.table.txt"),header = TRUE))
  
  # Add post reads to sample tables
  sampleDFdataset.f <- merge(sampleDFdataset,readsTable[,c(1,2,3)], by.x = "run_accession", by.y = "sample")
  # colnames(sampleDFdataset.f)
  sampleDFdataset.f
}


# 1.3) Load MiXCR data with repLoad - from all datasets
loadMiXCR <- function(apd1ParentPath, mixcrPath, datafolderL, suffix){
  # print(datafolderL)
  # apd1ParentPath <- aPD1.parentFolder
  # mixcrPath <- pathToMixcr
  # datafolderL <- data.foldersList[2]
  # suffix <-".TRA.txt"
  clonAllFiles <-list.files(paste0(apd1ParentPath, datafolderL,"/",mixcrPath), pattern=paste0("*.clonotypes",suffix), all.files=FALSE,
                            full.names=FALSE)
  
  # clonAllFiles.clean<-str_before_nth(clonAllFiles, "[.]", 1)
  # clonAllFiles.clean.df <- as.data.frame(clonAllFiles.clean)
  # colnames(clonAllFiles.clean.df)<- "Sample"
  # 
  mixcr_data <- immunarch::repLoad(paste0(apd1ParentPath, datafolderL,"/",mixcrPath,clonAllFiles),.format="mixcr")#
  names(mixcr_data$data) <-  str_before_nth(names(mixcr_data$data), "[.]", 1)
  mixcr_data$meta$Sample <-  str_before_nth(mixcr_data$meta$Sample, "[.]", 1)
  mixcr_data
}





addRestMetadata<- function(metaDFmerged, tcrDatasetList, datasetName){
  # metaDFmerged <- samples_all.merged.fix
  # tcrDatasetList <- immdata.apd1[[datasets[1]]]
  # datasetName <- datasets[1]
  print(datasetName)
  # filter metadata merged
  metaDFmerged.sub <- metaDFmerged %>% dplyr::filter(dataset == datasetName) %>% droplevels() %>% distinct()
  tcrDatasetList$meta <- merge(metaDFmerged.sub,tcrDatasetList$meta, by.x="run_accession", by.y = "Sample" )
  tcrDatasetList
}

createTcrCloneFiles <- function(apd1ParentPath, mixcrPath,mixcrPathFilter, datafolderL){
  print(datafolderL)
  
  # apd1ParentPath <- paste0(rootDir,antiPDL1publicInputData)
  # mixcrPathFilter <- pathToMixcrFilter
  # mixcrPath <- pathToMixcr
  # datafolderL <- data.foldersList[2]

  clonAllFiles <-list.files(paste0(apd1ParentPath, datafolderL,"/",mixcrPath), pattern="*.clonotypes.ALL.txt", all.files=FALSE,
                            full.names=FALSE)
  
  
  # Now load all tcr clone files for each sample of the data and write a new one.
  processCloneFiles <- sapply(clonAllFiles, function(x) readCloneProcess(paste0(apd1ParentPath, datafolderL,"/",mixcrPath), paste0(apd1ParentPath, datafolderL,"/",mixcrPathFilter), x), simplify = FALSE)
  names(processCloneFiles) <- str_replace(names(processCloneFiles), ".clonotypes.ALL.txt","")
  processCloneFiles
  
}



createBcrCloneFiles <- function(apd1ParentPath, mixcrPath,mixcrPathFilter, datafolderL, igpattern = NULL){
  print(datafolderL)
  
  # apd1ParentPath <- paste0(rootDir,antiPDL1publicInputData)
  # mixcrPathFilter <- pathToMixcrFilter
  # mixcrPath <- pathToMixcr
  # datafolderL <- data.foldersList[2]
  
  # clonAllFiles <-list.files(paste0(apd1ParentPath, datafolderL,"/",mixcrPath), pattern=igpattern, all.files=FALSE,
  #                           full.names=FALSE)
  
  clonAllFiles <-list.files(paste0(apd1ParentPath, datafolderL,"/",mixcrPath), pattern="*.clonotypes.ALL.txt", all.files=FALSE,
                            full.names=FALSE)
  
  # Now load all tcr clone files for each sample of the data and write a new one.
  processCloneFiles <- sapply(clonAllFiles, function(x) readCloneProcessBCR(paste0(apd1ParentPath, datafolderL,"/",mixcrPath), paste0(apd1ParentPath, datafolderL,"/",mixcrPathFilter), x), simplify = FALSE)
  # names(processCloneFiles) <- str_replace(names(processCloneFiles),  paste0(str_sub(igpattern,2),".txt"),"")
  names(processCloneFiles) <- str_replace(names(processCloneFiles), ".clonotypes.ALL.txt","")
  processCloneFiles
  # clonAllFiles.clean<-str_before_nth(clonAllFiles, "[.]", 1)
  # clonAllFiles.clean.df <- as.data.frame(clonAllFiles.clean)
  # colnames(clonAllFiles.clean.df)<- "Sample"
  # # 
  # mixcr_data <- repLoad(paste0(apd1ParentPath, datafolderL,"/",mixcrPath,clonAllFiles))
  # names(mixcr_data$data) <-  str_before_nth(names(mixcr_data$data), "[.]", 1)
  # mixcr_data$meta$Sample <-  str_before_nth(mixcr_data$meta$Sample, "[.]", 1)
  # mixcr_data
}
readCloneProcess <- function(pathToCloneFile, pathToCloneFileProc, cloneFile){
  ## applied to each clonotype file, of each sample ##
  # # paste0(apd1ParentPath, datafolderL,"/",mixcrPath), paste0(apd1ParentPath, datafolderL,"/",mixcrPathFilter)
  # pathToCloneFile <-paste0(apd1ParentPath, datafolderL,"/",mixcrPath)
  # pathToCloneFileProc <- paste0(apd1ParentPath, datafolderL,"/",mixcrPathFilter)
  # cloneFile <-clonAllFiles[1]

  # Check if folder for processed data exists, if not create it
  dir.create(file.path(pathToCloneFileProc), showWarnings = FALSE)
  # # Load in original file for sample
  # clones_results <- read.delim(file= paste0(
  #         pathToCloneFile, cloneFile), header=T)
  # # Process file to filter out immunoglobulin data
  # tcrCdr3s.ind <-clones_results$allVHitsWithScore %like% '^TRA.|^TRB.|^TRG.|^TRD.'
  # clones_results.tcr <- clones_results[tcrCdr3s.ind,]
  # # Write file to new mixcr folder, same name as before
  # cloneFile_new <- str_replace(cloneFile,"ALL","TCR")
  # write_delim(clones_results.tcr, 
  #             file= paste0(pathToCloneFileProc, cloneFile_new), delim = "\t")
  
  # write.table(clones_results.tcr, file = paste0(
  #         pathToCloneFileProc, cloneFile_new), sep = "\t",
  #           row.names = TRUE, col.names = TRUE)
  
  # clones_results.c <- read.delim(file= paste0(
  #         pathToCloneFileProc, cloneFile_new), header=T)
  
  # OR LOAD AND MERGE ALL TCR CLONE FILES IN ONE
  clones_tcra <- read.delim(file= paste0(
    pathToCloneFile, str_replace(cloneFile,"ALL","TRA")), header=T)
  clones_tcrb <- read.delim(file= paste0(
    pathToCloneFile, str_replace(cloneFile,"ALL","TRB")), header=T)
  clones_tcrd <- read.delim(file= paste0(
    pathToCloneFile, str_replace(cloneFile,"ALL","TRD")), header=T)
  clones_tcrg <- read.delim(file= paste0(
    pathToCloneFile, str_replace(cloneFile,"ALL","TRG")), header=T)
  
  clone_tcr <- rbind(clones_tcra,clones_tcrb, clones_tcrd, clones_tcrg)
  # clone_tcr[order(clone_tcr$cloneId),]
  clone_tcr.f <- clone_tcr %>% distinct()
  clone_tcr.f <- clone_tcr.f[order(clone_tcr.f$cloneId),]
  cloneFile_new <- str_replace(cloneFile,"ALL","TCR")
  # write_delim(clone_tcr.f, 
  #             file= paste0(pathToCloneFileProc, cloneFile_new), delim = "\t", na=" ",quote_escape=FALSE)
  clone_tcr.f
}


readCloneProcessBCR <- function(pathToCloneFile, pathToCloneFileProc, cloneFile){
  ## applied to each clonotype file, of each sample ##
  
  # pathToCloneFile <-paste0(apd1ParentPath, datafolderL,"/",mixcrPath)
  # pathToCloneFileProc <- paste0(apd1ParentPath, datafolderL,"/",mixcrPathFilter)
  # cloneFile <-clonAllFiles[11]
  
  # Check if folder for processed data exists, if not create it
  dir.create(file.path(pathToCloneFileProc), showWarnings = FALSE)
  # # Load in original file for sample
 
  # OR LOAD AND MERGE ALL TCR CLONE FILES IN ONE
  clones_igh <- read.delim(file= paste0(
    pathToCloneFile, str_replace(cloneFile,"ALL","IGH")), header=T)
  clones_igl <- read.delim(file= paste0(
    pathToCloneFile, str_replace(cloneFile,"ALL","IGL")), header=T)
  # clones_tcrb <- read.delim(file= paste0(
  #   pathToCloneFile, str_replace(cloneFile,"ALL","TRB")), header=T)
  # clones_tcrd <- read.delim(file= paste0(
  #   pathToCloneFile, str_replace(cloneFile,"ALL","TRD")), header=T)
  # clones_tcrg <- read.delim(file= paste0(
  #   pathToCloneFile, str_replace(cloneFile,"ALL","TRG")), header=T)
  
  clone_bcr <-  rbind(clones_igh, clones_igl)
  # clone_tcr[order(clone_tcr$cloneId),]
  clone_bcr.f <- clone_bcr %>% distinct()
  clone_bcr.f <- clone_bcr.f[order(clone_bcr.f$cloneId),]
  cloneFile_new <- str_replace(cloneFile,"ALL","BCR")
  # write_delim(clone_bcr.f, 
  #             file= paste0(pathToCloneFileProc, cloneFile_new), delim = "\t", na=" ",quote_escape=FALSE)
  clone_bcr.f
}


getTotalCloneReads <- function(datasetName, cloneDataList){
  # datasetName <- datasets[1:10][1]
  # cloneDataList <- process.res[[datasetName]]
  
  # Loop through all samples, and extract total clone reads
  datasetCloneReads <-sapply(cloneDataList, function(x) sum(x[["cloneCount"]]))
  datasetCloneReads <- as.data.frame(datasetCloneReads)
  datasetCloneReads <- datasetCloneReads %>% rownames_to_column(var="sampleid")
  datasetCloneReads
}


process_CDR3data <- function(datasetName, mixcrList_dataset,chainSel="all"){
  # datasetName <-datasets[1]
  # mixcrList_dataset <-immdata.apd1.f[[datasetName]]
  # print(paste0(datasetName,", TCR chain: ", chainSel))
  
  immdata_tcrCDR3 <- sapply(mixcrList_dataset$data, function(x) extract_TCR_CDR3(x,chain=chainSel), simplify = FALSE)
  immdata_tcrCDR3.final <-list()
  immdata_tcrCDR3.final$data <-immdata_tcrCDR3
  immdata_tcrCDR3.final$meta <-mixcrList_dataset$meta
  immdata_tcrCDR3.final
}

#loop through data list
extract_TCR_CDR3 <- function(sampleDF, chain="all"){
  # sampleDF <- mixcrList_dataset$data$SRR3184279
  if(chain=="all"){
    tcrCdr3s.ind <-sampleDF$V.name %like% '^TRA.|^TRB.|^TRG.|^TRD.'
  }else if(chain=="a"){
    tcrCdr3s.ind <-sampleDF$V.name %like% '^TRA.'
  }else if(chain=="b"){
    tcrCdr3s.ind <-sampleDF$V.name %like% '^TRB.'
  }else if(chain=="g"){
    tcrCdr3s.ind <-sampleDF$V.name %like% '^TRG.'
  }else if(chain=="d"){
    tcrCdr3s.ind <-sampleDF$V.name %like% '^TRD.'
  }else if(chain=="igh"){
    tcrCdr3s.ind <-sampleDF$V.name %like% '^IGH.'
  }else if(chain=="ab"){
    tcrCdr3s.ind <-sampleDF$V.name %like% '^TRA.|^TRB.'
  }else if(chain=="hl"){
    tcrCdr3s.ind <-sampleDF$V.name %like% '^IGH.|^IGL..'
  }
  
  
  
  sampleDF.tcr <- sampleDF[tcrCdr3s.ind,]
  # if(nrow(sampleDF.tcr)){
  #   print(paste0("No TCR CDR3 sequences found for sample"))
  # }
  sampleDF.tcr
}

sampleFilteringRunAcc <- function(immDF){
  # sampleDF <- samples_all.merged.fix
  # immDF <-immdata.apd1.tcr[[datasets[1:10][3]]]
  
  
  print(paste0("We have available TCR/BCR data for ",length(immDF$data), " samples, we include ", nrow(immDF$meta), " samples!"  ))
  
  # length(immDF$data)
  # nrow(immDF$meta)
  immdata <- immDF$data
  immmeta <- immDF$meta
  immdata <-immdata[immmeta$run_accession]
  immDF$data <- immdata
  
  immDF
}


plot_CountCloneReads <- function(resObj,title,yLabel,filterVar,dataseSel,yInter = 10,printLabels=FALSE,transformY=FALSE){
  
  cat('\n')  
  cat('\n') 
  cat("#### ", dataseSel, "\n") 
  cat('\n') 
  
  # resObj <- cloneReads.df
  # title <-"TCR Clone reads"
  # yLabel <-"Clone Reads"
  # filterVar <- "dataset"
  # dataseSel <- "Hugo"
  # printLabels=FALSE
  # transformY=FALSE
  title <- paste0(dataseSel," :",title)
  # Plot
  filterVar.un <- as.name(filterVar)
  # Results object to dataframe
  resObj.DF <- as.data.frame(melt(resObj))
  res.m <-resObj.DF %>% dplyr::filter(!!filterVar.un %in% c(dataseSel)) %>% droplevels()
  res.m[["variable"]]<- as.factor(res.m[["variable"]])
  # res.m <-resObj.DF
  p <- ggplot(res.m,aes(x=title, y=value,group =title, label=value, fill=variable))+ scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
  p <- p + geom_line( color="black")
  p <- p+ geom_point(shape=21, size=3)
  p <- p + geom_hline(yintercept=yInter, linetype="dashed", color = "red") + 
    geom_text(aes( 1, yInter, label = yInter, vjust = -1), size = 3)
  if(isTRUE(printLabels)){
    p <- p + geom_label_repel(fill="white",
                              colour = "black",
                              # fontface = "bold",
                              nudge_x = 0.25,
                              nudge_y = 0.25,
                              check_overlap = TRUE,
                              size=3.5)
  }
  
  if(isTRUE(transformY)){
    p <- p +  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                            labels = trans_format("log10", math_format(10^.x)))
  }
  p <- p + facet_wrap(as.formula(paste0("~ ","response")), ncol=2, strip.position = "left")
  #p <- p + scale_x_log10()
  p <- p + ylab(yLabel)
  p <- p + xlab("Sample ID")
  p <- p + ggtitle(title)
  p <- p + theme_bw(base_family = "Palatino Linotype", base_size = 12)
  p <- p + theme(legend.position = "bottom",
                 panel.border = element_rect(colour = "black", size=0.2),
                 panel.background = element_rect(fill = "white"),
                 axis.text.x = element_text(angle = 40, hjust = 1,size=10),
                 # axis.text.x = element_text(size=10),
                 axis.text.y = element_text(size = 12),
                 axis.title.x = element_text(size=14, face="bold"),
                 axis.title.y = element_text(size=14, face="bold"),
                 legend.key = element_rect(size = 2),
                 legend.key.height=unit(0.6,"line"),
                 legend.key.width=unit(1,"line"),
                 legend.text = element_text(size = 11),
                 legend.background = element_rect(colour = "black"),
                 legend.direction='horizontal',legend.box='horizontal',
                 legend.title = element_text(colour="black", size=12, face="plain"))
  p <- p + guides(fill=guide_legend(title="TCR Chain clones"))
  #p <- p + scale_color_brewer(palette="Paired")
  # p <- p + scale_color_manual(values = brewer.pal(length(nSelected), "Dark2"))
  # p <- p + scale_fill_manual(values = brewer.pal(length(nSelected), "Dark2"))
  # p <- p + labs(color="gene",fill="gene") 
  # p
  print(p)
}



getTotalCloneReads_immArch <- function(datasetName, cloneDataList){
  # datasetName <- datasets[1:10][1]
  # cloneDataList <- immdata.apd1.tcrB.f[[datasetName]]$data
  
  # Loop through all samples, and extract total clone reads
  datasetCloneReads <-sapply(cloneDataList, function(x) sum(x[["Clones"]]))
  datasetCloneReads <- as.data.frame(datasetCloneReads)
  datasetCloneReads <- datasetCloneReads %>% rownames_to_column(var="sampleid")
  datasetCloneReads
}


sampleFilteringRunAccForDowns <- function(immDF, samplesKeep, datasetName){
  # samplesKeep <- cloneReads_tcrb.filt10[[datasets[1:10][1]]]
  # immDF <-immdata.apd1.tcrB.f[[datasets[1:10][1]]]
  # datasetName <-datasets[1:10][1]

  
  # length(immDF$data)
  # nrow(immDF$meta)
  immdata <- immDF$data
  immmeta <- immDF$meta
  immdata <-immdata[samplesKeep]
  immmeta <- immmeta %>% dplyr::filter(run_accession %in% samplesKeep) %>% droplevels()
  print(paste0(datasetName, ": keeping ", length(samplesKeep), " out of ", length(immDF$data), " samples for downsampling!"))
  
  immDF$data <- immdata
  immDF$meta <- immmeta
  
  
  immDF
}


mixcrProcessing <- function(datasetName, mixcrList_dataset, filterCoding = TRUE, downSampling = FALSE, dsMethod = NULL, dsN = NULL){
  # datasetName <- datasets[1]
  # mixcrList_dataset <- immdata.apd1.tcrA.filt[[datasetName]]
  # dsMethod <- "downsample"
  # dsN <-10
  # filterCoding = TRUE
  # print(datasetName)
  if(isTRUE(filterCoding)){
    ## 1. Filter functional/non-functional/in-frame/out-of-frame clonotypes: WORKS PER SAMPLE!So needs sapply to be applied over all samples
    immdata_cod <- sapply(mixcrList_dataset$data, function(x) coding(x), simplify = FALSE)
  }else{
    immdata_cod <- mixcrList_dataset$data
  }
  
  if(isTRUE(downSampling)){
    ## 2. Downsampling: - Subsampling
    
    if(is.null(dsMethod )){
      print("Please define downsampling method!")
    }else if(dsMethod=="downsample"){
      #repSample chooses .n clones (not clonotypes!) from the input repertoires without any probabilistic simulation, but exactly computing each choosed clones
      #more consistent and biologically pleasant  
      immdata_cod.ds = repSample(immdata_cod, dsMethod, dsN)
      
    }else if(dsMethod=="sample"){
      # repSample chooses .n clonotypes (not clones!) randomly. Depending on the .prob argument, the function chooses clonotypes either according to their size (if .prob is TRUE, by default), or each clonotype has an equal chance to be choosed (if .prob is FALSE). Note that sampling is done without replacing
      immdata_cod.ds = repSample(immdata_cod, dsMethod, .n = dsN)
      
    }
  }else{
    immdata_cod.ds <- immdata_cod
  }
  
  # ## 3. Repertoire Diversity calculation
  # immdata_cod.ds.div <- repDiversity(immdata_cod.ds, "chao1")# chao1 returns 4 values: estimated number of species, standart deviation of this number and two 95
  # immdata_cod.ds.div
  immdata_cod.ds
  immDF <- list()
  immDF$data <- immdata_cod.ds
  immDF$meta <-mixcrList_dataset$meta
  immDF
}


diversityCalculations <- function(immdata,stepSelected){
  # immdata <- immdata.apd1.proc.TcrA.noDS[[datasets[1:10][2]]]
  
  # Chao1 diversity measure
  div_chao <- repDiversity(immdata$data, "chao1")
  # Hill numbers
  div_hill <- repDiversity(immdata$data, "hill")
  
  # D50
  div_d50 <- repDiversity(immdata$data, "d50")
  
  # Ecological diversity measure
  div_div <- repDiversity(immdata$data, "div",.q=3)
  
  # # # # Rarefaction
  # imm_raref <- repDiversity(immdata$data, "raref", .verbose = T,.norm=T,.step = stepSelected)
  
  # dXX
  div_dXX <- repDiversity(immdata$data, "dxx",.perc=30)
  # immdata$diversityEst <- list("chao"=div_chao, "hill"=div_hill, "d50"=div_d50,"trueDiv"=div_div,"dxx"=div_dXX,"raref" = imm_raref)
  immdata$diversityEst <- list("chao"=div_chao, "hill"=div_hill, "d50"=div_d50,"trueDiv"=div_div,"dxx"=div_dXX)
  
  immdata
}

diversityEst_vis <- function(datasetName, immdata){
  # datasetName <- datasets[1:10][2]
  # immdata<- immdata.apd1.proc.TcrB.divEst.noDS[[datasets[1:10][2]]]
  cat('\n')  
  cat('\n') 
  cat("##### ", datasetName, "\n") 
  cat('\n') 
  
  immdata$meta <- as.tibble(immdata$meta)
  immdata$meta$response <- ifelse(immdata$meta$response=="CR+PR","CRPR",immdata$meta$response)
  colnames(immdata$meta)[1] <- "Sample"
  
  # Chao
  p1 <- vis(immdata$diversityEst$chao, .by = c("response", "gender"), .meta = as.tibble(immdata$meta))
  # Hill
  p2 <- vis(immdata$diversityEst$hill, .by = c("response", "gender"), .meta = as.tibble(immdata$meta))
  # d50
  p3 <- vis(immdata$diversityEst$d50, .by = c("response"), .meta = as.tibble(immdata$meta))
  # True diversity
  p4 <- vis(immdata$diversityEst$trueDiv, .by = c("response"), .meta = as.tibble(immdata$meta))
  # dxx
  p5 <- vis(immdata$diversityEst$dxx, .by = c("response"), .meta = as.tibble(immdata$meta))
  
  print(p1+p2+p3+p4+p5)
  
}


get_modeledEntropy <- function(immdata, minAbund, colnameEntropy){
  # immdata <-immdata.apd1.proc.TcrA.divEst[[datasets[1:10][1]]]
  # minAbund <- 4
  # colnameEntropy <- "modEntr.ds"
  
  modelEntr <-sapply(immdata$data, function(x) predict_true_entropy_from_counts(x[["Clones"]],min_abundance = minAbund,should_screen_counts = TRUE), simplify = FALSE)
  modelEntr.df <-ldply(modelEntr) 
  # modelEntr.df <- modelEntr.df %>% column_to_rownames(var = ".id")
  colnames(modelEntr.df) <- c("sampleid",colnameEntropy)
  immdata$diversityEst[[colnameEntropy]] <-modelEntr.df
  immdata
}

get_Richness <- function(immdata, colnameEntropy){
  # immdata <-immdata.apd1.proc.TcrA.divEst[[datasets[1:10][1]]]
  # colnameEntropy <- "TCRArichness.ds"
  # 
  richness <-sapply(immdata$data, function(x) nrow(x), simplify = FALSE)
  richness.df <-ldply(richness) 
  # richness.df <- richness.df %>% column_to_rownames(var = ".id")
  colnames(richness.df) <- c("sampleid",colnameEntropy)
  immdata$diversityEst[[colnameEntropy]] <- richness.df
  immdata
}

#########################
## SIGNATURE FUNCTIONS ##
#########################
load_rna_sigScore.weights.allSigs <- function(pathToData,vstData,fpkmData,selectedStudy, TuTackSig, GepSig, GepNormGenes, degsOnly, samplesDF,sigTissue, addConfCoef=NA,scaleFactor=NULL){
  # selectedStudy<-datasets[1]
  # TuTackSig <- TuTack_mel
  # GepSig<- gep_weights_df
  # GepNormGenes <- normalization_genes_df
  # degsOnly <- degs.only
  # samplesDF <- samples_datasets[[1]]
  # addConfCoef=NA
  # scaleFactor=NULL
  # ######################
  print(selectedStudy)
  # print(paste0(projectRDataPath,selectedStudy,"_",vst.RData))
  # Load expression FPKM-UQ data
  # tpm <- loadRData(paste0(expressionDataDir,selectedStudy,"_",tpm.RData))
  # Add pseudo normalized count and log10
  df_all <- loadRData(paste0(pathToData,selectedStudy,"_ALLresp",vstData))
  ## PDL1 expression
  pdl1 <- as.data.frame(t(subset(df_all, rownames(df_all) %in% c("CD274"))))
  ## CD8A gene
  cd8a <- as.data.frame(t(subset(df_all, rownames(df_all) %in% c("CD8A"))))
  ## CD8B gene
  cd8b <- as.data.frame(t(subset(df_all, rownames(df_all) %in% c("CD8B"))))
  # df_all <- apply(df_all, 2, log10)
  ## Subset data to TuTACK genes
  # df_all <- as.data.frame(t(subset(df_all, rownames(df_all) %in% TuTackSig$genes)))
  ## CHECK IF ALL GENES IN TUTACK
  # setdiff(str_replace(TuTack,"_","-"),rownames(df_all))
  if(!is.null(scaleFactor)){
    # print("STEP 1")
    df_all <- as.data.frame(t(subset(df_all, rownames(df_all) %in% str_replace(scaleFactor$method$scale,"_","-"))))
    # print("STEP 2")
    colnames(df_all) <-str_replace(colnames(df_all),"-","_")
    # Correct perhaps diff order of columns
    df_all <- df_all[,scaleFactor$method$scale]
    # print("STEP 3")
    df_all.scaled <-predict(scaleFactor, df_all)
    # print("STEP 4")
    colnames(df_all.scaled) <-str_replace(colnames(df_all.scaled),"_","-")
    ## Subset data to TuTACK genes
    df_all <- df_all.scaled %>% dplyr::select(all_of(TuTackSig$genes))    # scaleFactor$method$scale
  }else{
    # df_all <- apply(df_all, 2, log10)
    ## Subset data to TuTACK genes
    # print("STEP 5")
    df_all <- as.data.frame(t(subset(df_all, rownames(df_all) %in% TuTackSig$genes)))
    # Correct perhaps diff order of columns
    df_all <- df_all[,TuTackSig$genes]
  }
  
  #####
  
  ######
  ## TUTACK
  # Multiply each normalized RNA value by corresponding scoring weight of the gene
  weighted_RNA_df <- df_all*TuTackSig$weights[match(colnames(df_all), TuTackSig$genes)][col(df_all)]
  # print(dim(weighted_RNA_df))
  # Add weighted RNA expression values and generate gene signature score for each sample
  if (!is.na(addConfCoef)){
    TuTack.scores <- rowSums(weighted_RNA_df) + addConfCoef
  }else{
    TuTack.scores <- rowSums(weighted_RNA_df)
  }
  df_all$TuTACK_sigscore <- TuTack.scores[match(rownames(df_all),names(TuTack.scores))]
  
  # ### TuTack IFNG
  
  ## GEP
  ### IMPORTANT: THIS DATA SHOULD NOT BE SCALED!!!!! ##########
  
  # Load expression FPKM-UQ data
  load(paste0(pathToData,selectedStudy,"_ALLresp",fpkmData))
  my_data_df <-as.data.frame(fpkm.uq)
  ## CHECK IF ALL GENES IN TUTACK
  setdiff(GepSig$gene,rownames(my_data_df))
  setdiff(GepNormGenes$gene,rownames(my_data_df))
  #rownames(my_data_df) <- as.character(geneid2symb$gene_name)[match(rownames(my_data_df),as.character(geneid2symb$gene_id))]
  # Ensemble IDs to gene IDs
  # subset df to TIS
  # setdiff(GepSig$gene,rownames(my_data_df))
  #setdiff(GepNormGenes$gene,rownames(my_data_df))
  tis_init <- as.data.frame(subset(my_data_df, rownames(my_data_df) %in% GepSig$gene))
  df_all_tis <- as.data.frame(t(tis_init))
  ## Fix column order
  df_all_tis <- df_all_tis[,as.character(GepSig$gene)]
  df_all_tis$Sample_ID <- rownames(df_all_tis)
  df_all_tis$dataset<-selectedStudy
  # METHOD B: as in Ayers et al, instead use log10
  # Step1: subset df to housekeeping & normalize with log2 housekeeping values
  houseKeep_init <- as.data.frame(subset(my_data_df, rownames(my_data_df) %in% GepNormGenes$gene))
  
  houseKeep_log10 <- as.data.frame(t(log10(houseKeep_init+1)))
  ### Fix column order
  houseKeep_log10 <- houseKeep_log10[,as.character(GepNormGenes$gene)]
  # Step2: Calculate the arithmetic mean of each housekeeping gene
  houseKeep_log10$Mean <- rowMeans(houseKeep_log10)
  # Create vector with samples and the housekeeping genes arithmetic means
  houseKeep_Mean <-houseKeep_log10$Mean
  names(houseKeep_Mean)<- rownames(houseKeep_log10)
  # Step3:log10 transform data to avoid extremely skewed gene expression distribtions
  tis_df_log10 <- log10(df_all_tis[,-c(dim(df_all_tis)[2]-1,dim(df_all_tis)[2])]+1)
  # Step4: re-normalize rsem data using the housekeeping genes - Follow Danaher
  # log2 count of each gene was normalized by substracting the arithmetic mean of the log2 counts of the housekeeping genes
  # function that substracts housekeeping gene arithmetic mean: log2GeneValue - log2(MeanHouseKeepGenes) - IN THE SAMPLE
  housekeeping_genes <- colnames(tis_df_log10)
  houseKeep_normalization <- function(sample,tis_df_log10, sampleHKmeans){
    sample_mean <- sampleHKmeans[sample]
    sampleIndRow <- which(rownames(tis_df_log10)==sample)
    row_update <-tis_df_log10[sampleIndRow,] - sample_mean
    tis_df_log10[sampleIndRow,] <<- row_update
  }
  housKeep_norm_res <- sapply(rownames(tis_df_log10),function(x) houseKeep_normalization(sample=x,tis_df_log10, houseKeep_Mean))
  housKeep_norm_res <-tis_df_log10
  # Now we have normalized and housekeep normalized values of the TIS genes
  # Multiply each normalized RNA value by corresponding scoring weight of the gene
  weighted_RNA_df <- housKeep_norm_res*GepSig$scoring.weight[match(colnames(housKeep_norm_res), GepSig$gene)][col(housKeep_norm_res)]
  # Add weighted RNA expression values and generate gene signature score for each sample
  Tis_scores <- rowSums(weighted_RNA_df)
  df_all$TIS_sigscore <- Tis_scores[match(rownames(df_all),names(Tis_scores))]
  
  ## ADD PD-L1
  df_all$PDL1 <- pdl1$CD274[match(rownames(df_all),rownames(pdl1))]
  ## ADD cd8A,B genes
  df_all$CD8A<- cd8a$CD8A[match(rownames(df_all),rownames(cd8a))]
  df_all$CD8B<- cd8b$CD8B[match(rownames(df_all),rownames(cd8b))]
  # #########################
  df_all <- df_all %>% rownames_to_column("run_accession")
  df.merged <-merge(df_all[,c(1,(dim(df_all)[2]-4):dim(df_all)[2])], samplesDF, by.x = c("run_accession"), by.y = "run_accession")
  df.merged <- df.merged[,c("title","run_accession","response","TuTACK_sigscore", "TIS_sigscore", "PDL1","CD8A","CD8B","dataset")]
  df.merged$signature <- sigTissue
  # df.merged
  # df.merged$dataset <- selectedStudy
  as.data.frame(df.merged)
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


compareBox <- function(DF, xVar,yVar,groupVar=NULL, facetVar =NULL,plotTitle,xLab, yLab){
  # DF <-immdata.apd1.proc.TcrA.divEst.fin
  # xVar <-"response_gender"
  # yVar <-"TCRArichness.ds"
  # # yVar <- c("")
  # groupVar <- "response_gender"
  # # groupVar <- c("response","gender")
  # 
  # facetVar <- "tissue"
  # xLab <- "response"
  # yLab <- "TCR alpha richness"
  # plotTitle <- ""
  # 
  xVar.un <- as.name(xVar)
  yVar.un <- as.name(yVar)
  
  if(length(groupVar)>1){
    DF[[paste(groupVar,collapse = "_")]]  <- interaction(DF[[groupVar[1]]],DF[[groupVar[2]]],sep = ":")
    groupVar <-paste(groupVar,collapse = "_")
  }
  
  #, label=!!labelVar.un
  if(!is.null(groupVar)){
    groupVar.un <- as.name(groupVar)
    p <-ggplot(DF, aes(x=!!xVar.un, y=!!yVar.un, fill=!!groupVar.un))
    
  }else{
    p <- ggplot(DF, aes(x=!!xVar.un, y=!!yVar.un))
    
  }
  
  p <- p + geom_boxplot()
  p <- p + geom_jitter(color="black", size=0.4, alpha=0.9)
  # p + stat_compare_means()
  # library(ggpubr)
  # library(rstatix)
  # stat.test <- DF %>% t_test(as.formula(paste0(xVar,"~ ",yVar)))
  # stat.test <- stat.test %>% add_xy_position(x = "dose")
  # 
  # p <- p + geom_point(shape=21, size=2,na.rm = TRUE,show.legend=TRUE,
  #                     alpha = 0.8)
  p <- p + theme_bw(base_family = "Palatino Linotype", base_size = 12)
  p <- p + theme(legend.position = "bottom",
                 panel.border = element_rect(colour = "black", size=0.2),
                 panel.background = element_rect(fill = "white"),
                 #axis.text.x = element_text(angle = 40, hjust = 1,size=6),
                 axis.text.x = element_text(size=12,angle = 40,hjust = 1),
                 axis.text.y = element_text(size = 12),
                 legend.key = element_rect(size = 2),
                 legend.key.height=unit(0.6,"line"),
                 legend.key.width=unit(1,"line"),
                 legend.text = element_text(size = 12),
                 legend.background = element_rect(colour = "black"),
                 legend.direction='horizontal',legend.box='horizontal'
  )
  if(!is.null(groupVar)){
    p <- p + scale_fill_brewer(palette = "Dark2")
  }
  if(!is.null(facetVar)){
    facetVar.un <- as.name(facetVar)
    p <- p + facet_wrap(as.formula(paste0("~ ",facetVar)),scales='free') 
  }
  p <- p + ylab(paste0(yLab))
  p <- p + xlab(paste0(xLab))
  p <- p + ggtitle(plotTitle)
  
  print(p)
  # p
  # 
  # ggstatsplot::ggbetweenstats(
  #   data = DF,
  #   x = response,
  #   y = modEntr.ds,
  #   grouping.var = tissue,
  #   plot.type = "boxviolin",
  #   pairwise.comparisons = TRUE,
  #   pairwise.display = "significant"
  # )
  
  
}