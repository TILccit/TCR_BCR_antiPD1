####################################################
## Create table with TCGA sample/patient info and ##
## Danaher TILs infiltration data, TMB, GEP, PD-1 ##
## expression, survival data                      ##
####################################################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##-----------##
## FUNCTIONS ##
##-----------##
##########################
## Loading RData object ##
##########################

`%notin%` <- Negate(`%in%`)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

##################################################
## Correct TCGA molecular subtypes cohort names ##
##################################################
tcga_cohortTypesCorrection <- function(DF,cohortColumn){
  DF[[cohortColumn]] <- ifelse(grepl("COAD.MSI_High",DF[[cohortColumn]]),"COAD_MSI",ifelse(grepl("COAD.MSI_Low",DF[[cohortColumn]]),"COAD_MSI.low",ifelse(grepl("BRCA.TripleNegative",DF[[cohortColumn]]),"BRCA_TN",DF[[cohortColumn]])))
  
  DF[[cohortColumn]] <- gsub("\\.","_",DF[[cohortColumn]])
  DF
}

################################################
## Get IFNG, Expanded immune signature scores ##
################################################

load_tcga_rna_add_signatureScore  <- function(selectedcancer,IFNG_rna, gs, genesDF,dataDir,dataDirVST,normal=NULL){
  # selectedcancer <- "TCGA-ACC"
  # IFNG_rna <- gene_signatures_rna
  # gs <- TuTACK_ensemble
  # genesDF <- tls_df
  # dataDir <- htseqRnaPath
  # dataDirVST <- htseqVSTRnaPath
  
  fname <- file.path(dataDir,
                     paste0(selectedcancer,"_RNASeqFPKM_UQ.RData"))
  #cat("working on ", selectedcancer, "\n")
  load(paste0(dataDir,selectedcancer,"_RNASeqFPKM_UQ.RData"))
  df_all <- df_all + 1
  df_all <- apply(df_all, 2, log10)
  
  # Add expanded ifg
  expanded_refined_imm_ges <- c("CD3D", "IDO1", "CIITA", "CD3E", 
                                "CCL5", "GZMK", "CD2", "HLA-DRA", "CXCL13",
                                "IL2RG", "NKG7", "HLA-E", "CXCR6", 
                                "LAG3", "TAGAP", "CXCL10", "STAT1", "GZMB")
  expanded_refined_imm_ges_ensembl <- c("ENSG00000167286", "ENSG00000131203", "ENSG00000179583", "ENSG00000198851", "ENSG00000271503", 
                                        "ENSG00000113088", "ENSG00000116824", "ENSG00000204287", "ENSG00000156234", "ENSG00000147168", "ENSG00000105374", "ENSG00000204592", "ENSG00000172215", 
                                        "ENSG00000089692", "ENSG00000164691", "ENSG00000169245", "ENSG00000115415", "ENSG00000100453")
  
  IFNG_ensemble <- c("ENSG00000167286", "ENSG00000131203", "ENSG00000179583", "ENSG00000198851", "ENSG00000271503", 
                     "ENSG00000113088", "ENSG00000116824", "ENSG00000204287", "ENSG00000156234", "ENSG00000147168", "ENSG00000105374", "ENSG00000204592", "ENSG00000172215", 
                     "ENSG00000089692", "ENSG00000164691", "ENSG00000169245", "ENSG00000115415", "ENSG00000100453")
  # subset df to IFNG
  df_all_exp <- as.data.frame(t(subset(df_all, rownames(df_all) %in% IFNG_ensemble)))
  df_all_exp$Sample_ID <- rownames(df_all_exp)
  df_all_exp$tcga.dataset<-selectedcancer
  colnames(df_all_exp)[1:18] <- expanded_refined_imm_ges[match(colnames(df_all_exp[1:18]),expanded_refined_imm_ges_ensembl)]
  
  df_all_exp$IFNg_expanded_sigscore <- apply(df_all_exp[, c("CD3D", "IDO1", "CIITA", "CD3E", "CCL5", "GZMK", "CD2", "HLA-DRA", "CXCL13","IL2RG", "NKG7", "HLA-E", "CXCR6", "LAG3", "TAGAP", "CXCL10", "STAT1", "GZMB")], 1, mean)
  # Select only cols with gs and sample id
  df_all_exp <- df_all_exp[,c(19,21,20)]
  df_all_exp
  
  ## Add small ifgn gs
  ifng_refined_ges <- c("IDO1", "CXCL10", "CXCL9", "HLA-DRA1", "STAT1", "IFNG")
  ifng_refined_ges_ensembl <- c("ENSG00000131203", "ENSG00000169245", "ENSG00000138755", "ENSG00000204287", "ENSG00000115415", "ENSG00000111537")
  
  IFNG_ensemble <- c("ENSG00000131203", "ENSG00000169245", "ENSG00000138755", "ENSG00000204287", "ENSG00000115415", "ENSG00000111537")
  # subset df to IFNG small
  df_all_small <- as.data.frame(t(subset(df_all, rownames(df_all) %in% IFNG_ensemble)))
  df_all_small$Sample_ID <- rownames(df_all_small)
  df_all_small$tcga.dataset<-selectedcancer
  colnames(df_all_small)[1:6] <- ifng_refined_ges[match(colnames(df_all_small[1:6]),ifng_refined_ges_ensembl)]
  # colnames(df_all) <- c("IDO1", "CXCL10", "CXCL9", "HLA-DRA1", "STAT1", "IFNg_sigscore","Sample_ID","tcga.dataset")
  df_all_small$IFNg_small_sigscore <- apply(df_all_small[, c("IDO1", "CXCL10", "CXCL9", "HLA-DRA1", "STAT1", "IFNG")], 1, mean)
  df_all_small <- df_all_small[,c(7,9,8)]
  df_all_small
  
  ## Add TuTACK gene signature sxore
  # rownames(df_all)<-genesDF$hgnc_symbol[match(rownames(df_all),genesDF$ensembl_gene_id)]
  # subset df to TuTACK
  # print(paste0(selectedcancer,"_RNASeqCountsVST.RData"))
  # fname <- file.path(dataDirVST,
  #                    paste0(selectedcancer,"_RNASeqCounts.RData"))
  # if (file.exists(fname)){
  #   load(paste0(dataDirVST,selectedcancer,"_RNASeqCountsVST.RData"))
  #   df_all<-vst.DF
  # }else{
  #   print(paste0("No htseq count&VST data for ", selectedcancer))
  # }
  # 
  if (!is.null(normal)){
    print("NORMAL TCGA DATA")
    print(paste0(selectedcancer,"_RNASeqCountsVST.RData"))
    fname <- file.path(dataDirVST,
                       paste0(selectedcancer,"_RNASeqCountsVST.RData"))
    if (file.exists(fname)){
      load(paste0(dataDirVST,selectedcancer,"_RNASeqCountsVST.RData"))
      df_all<-vst.DF
    }else{
      print(paste0("No htseq count & VST data for ", selectedcancer))
    }
    
  }else{
    print(paste0(selectedcancer,"_RNASeqCountsVST.RData"))
    fname <- file.path(dataDirVST,
                       paste0(selectedcancer,"_RNASeqCountsVST.RData"))
    if (file.exists(fname)){
      load(paste0(dataDirVST,selectedcancer,"_RNASeqCountsVST.RData"))
      df_all<-vst.DF
    }else{
      print(paste0("No htseq count & VST data for ", selectedcancer))
    }
    
  }
  #################################
  ## Subset data to TuTACK genes
  df_all_TuTACK <- as.data.frame(t(subset(df_all, rownames(df_all) %in% gs)))
  df_all_TuTACK$Sample_ID <- rownames(df_all_TuTACK)
  df_all_TuTACK$tcga.dataset<-selectedcancer
  ## Multiply each normalized RNA value by corresponding scoring weight of the gene
  weighted_RNA_df <- df_all_TuTACK[, gs]*TuTACK.DF$weights[match(colnames(df_all_TuTACK[, gs]), TuTACK.DF$ensemble)][col(df_all_TuTACK[, gs])]
  ## Add weighted RNA expression values and generate gene signature score for each sample
  print(paste0("Genes: ",dim(weighted_RNA_df)[2]))
  TuTack.scores <- rowSums(weighted_RNA_df)
  
  df_all_TuTACK$TuTACK_sigscore <- TuTack.scores[match(rownames(df_all_TuTACK),names(TuTack.scores))]
  
  df_all_TuTACK <- df_all_TuTACK[,c(dim(TuTACK.DF)[1] + 1,dim(TuTACK.DF)[1] + 2, dim(TuTACK.DF)[1] + 3)]
  ######################
  ## PDL1 expression
  pd_l1_ensemble <- "ENSG00000120217"
  # subset df to pd-L1
  df_all.pd_L1 <- as.data.frame(t(subset(df_all, rownames(df_all) %in% pd_l1_ensemble)))
  df_all.pd_L1$Sample_ID <- rownames(df_all.pd_L1)
  df_all.pd_L1$tcga.dataset<-selectedcancer
  colnames(df_all.pd_L1)[1] <- "PD_L1"
  df_all.pd_L1 <- df_all.pd_L1[,c(2,1,3)]
  ######################
  ## Calculate TLS signature score - mean gene expression across a sample
  df_all.tls <- as.data.frame(t(subset(df_all, rownames(df_all) %in% genesDF$ensemble_id)))
  df_all.tls <-rowMeans(df_all.tls, na.rm = TRUE) %>% as.data.frame()
  colnames(df_all.tls) <- "TLS_score"
  df_all.tls$Sample_ID <- rownames(df_all.tls)
  df_all.tls$tcga.dataset<-selectedcancer
  ####################
  # Merge all dataframes with sig scores
  df_sigScore <- merge(df_all_exp, df_all_small, by= c("Sample_ID","tcga.dataset"))
  df_sigScore <- merge(df_sigScore, df_all_TuTACK, by= c("Sample_ID","tcga.dataset"))
  df_sigScore <- merge(df_sigScore, df_all.pd_L1, by= c("Sample_ID","tcga.dataset"))
  df_sigScore <- merge(df_sigScore, df_all.tls, by= c("Sample_ID","tcga.dataset"))
  
  # df_sigScore <- df_sigScore[,c(1,6,2,4,7,8)]
  # df_sigScore
  IFNG_rna <- rbind(IFNG_rna,df_sigScore)
  IFNG_rna
}

######################################
## Get TIS score for TCGA datasets  ##
######################################
get_TIS_score_tcga_rsem <- function(selectedcancer, TIS_rna, gep_weights_df,normalization_genes_df){
  print(selectedcancer)
  # selectedcancer <-mycancertypes.tcgaBiolinks.ord.clean[1]
  # TIS_rna <- gene_TISsignature_rna
  # gep_weights_df <- gep_weights_df
  # normalization_genes_df <- normalization_genes_df
  # Load Each dataset separately
  readTCGA( path = pathRNA[selectedcancer], dataType = 'rnaseq' ) -> my_data
  my_data_df = do.call(rbind,my_data) 
  my_data_df <-as.data.frame(my_data_df[-c(1),])
  colnames(my_data_df)<- as.character(my_data$bcr_patient_barcode)
  my_data_df$gene <- as.character(rownames(my_data_df))
  my_data_df <-cSplit(my_data_df, "gene", "|")
  # subset df to TIS
  tis_init <- as.data.frame(subset(my_data_df, my_data_df$gene_1 %in% gep_weights_df$gene))
  rownames(tis_init) <- tis_init$gene_1
  tis_init <- tis_init[,-c(dim(tis_init)[2]-1,dim(tis_init)[2])]
  df_all_tis <- as.data.frame(t(tis_init))
  df_all_tis$Sample_ID <- rownames(df_all_tis)
  df_all_tis$tcga.dataset<-selectedcancer
  # #########################

  # #########################
  # METHOD B: as in Ayers et al, instead use log10
  # Step1: subset df to housekeeping & normalize with log2 housekeeping values
  houseKeep_init <- as.data.frame(subset(my_data_df, my_data_df$gene_1 %in% normalization_genes_df$gene))
  rownames(houseKeep_init) <- houseKeep_init$gene_1
  houseKeep_init <- houseKeep_init[,-c(dim(houseKeep_init)[2]-1,dim(houseKeep_init)[2])]
  # Issue with character, transform to numeric
  houseKeep_init <- houseKeep_init %>% mutate_if(is.character,as.numeric)
  houseKeep_log10 <- as.data.frame(t(log10(houseKeep_init+1)))
  # Step2: Calculate the arithmetic mean of each housekeeping gene
  houseKeep_log10$Mean <- rowMeans(houseKeep_log10)
  # Create vector with samples and the housekeeping genes arithmetic means
  houseKeep_Mean <-houseKeep_log10$Mean
  names(houseKeep_Mean)<- rownames(houseKeep_log10)
  # Step3:log10 transform data to avoid extremely skewed gene expression distribtions
  tis_df_log10 <- log10(df_all_tis[,-c(dim(df_all_tis)[2]-1,dim(df_all_tis)[2])]%>% mutate_if(is.character,as.numeric) +1)
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
  ## Now we have normalized and housekeep normalized values of the TIS genes
  ## Multiply each normalized RNA value by corresponding scoring weight of the gene
  weighted_RNA_df <- housKeep_norm_res*gep_weights_df$scoring.weight[match(colnames(housKeep_norm_res), gep_weights_df$gene)][col(housKeep_norm_res)]
  ## Add weighted RNA expression values and generate gene signature score for each sample
  Tis_scores <- rowSums(weighted_RNA_df)
  ########################

  # #########################
  # Add TIS score to dataset df
  df_all_tis$TIS_sigscore <- Tis_scores[match(df_all_tis$Sample_ID,names(Tis_scores))]
  # Add to TIS dataframe
  TIS_rna <- rbind(TIS_rna,df_all_tis[,c("Sample_ID","tcga.dataset","TIS_sigscore")])
  TIS_rna
}



load_obj <- function(f)
{
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

calculateEstimateScore <- function(pathToRnaData,RnaFileSuffix,pathToData,tissue,geneid2symb){
  
  # pathToRnaData <-pathtoTCGAvst
  # RnaFileSuffix <-"_RNASeqCounts.RData"
  # pathToData <-scriptDataPath
  # tissue <-paste0(mycancertypes.off[1])
  # 
  print(paste0("Calculating stromal, immune and estimate scores for: ", tissue))
  
  if (file.exists(paste0(pathToRnaData,tissue,RnaFileSuffix))){
    # Load rnaseq vst data
    # vst.df <- loadRData(paste0(pathToRnaData,tissue,RnaFileSuffix))
    vst.df <- load_obj(paste0(pathToRnaData,tissue,RnaFileSuffix))
    rownames(vst.df) <- geneid2symb$gene_name[match(rownames(vst.df),geneid2symb$gene_id)]
    # Subset to estimate genes
    vst.df<-vst.df %>% as.data.frame() %>% rownames_to_column(var="geneSymb") %>% dplyr::filter(geneSymb %in% common_genes$GeneSymbol) %>% column_to_rownames(var="geneSymb")
    
    # Write input file for estimate
    write.table(vst.df,paste0(pathToData,tissue,"_estimate",".input.txt"), sep="\t", quote=F)
    #
    inputFile <- paste0(pathToData,tissue,"_estimate",".input.txt")
    # intersect input data with 10,412 common genes
    # unifies different number of genes per platform against 10,412 common genes.
    filterCommonGenes(input.f=inputFile , output.f=paste0(pathToData,tissue,"_estimate",".gct"), id="GeneSymbol")
    #[1] "Merged dataset includes 10056 genes (356 mismatched)."
    
    # computes stromal, immune, and ESTIMATE scores per sample using gene expression data.
    estimateScore(paste0(pathToData,tissue,"_estimate",".gct"), paste0(pathToData,tissue,"_estimate",".score.gct"), platform="illumina")
    
  }else{
    print(paste0("No htseq count data for ", tissue))
  }
  
  
  
}


readEstimateFiles <- function(pathToData,tissue){
  print(paste0("Reading in stromal, immune and estimate scores for: ", tissue))
  if (file.exists(paste0(pathToData,tissue,"_estimate",".score.gct"))){
    # pathToData <-projectDataPath
    # tissue <-tissues[1]
    estFile <-read.table(paste0(pathToData,tissue,"_estimate",".score.gct"), header = TRUE, sep = "\t", dec = ".", skip =2)
    estFile <- estFile %>% dplyr::select(-Description) %>% column_to_rownames(var="NAME")
    estFile <- t(estFile)
    estFile <- as.data.frame(estFile)
    estFile <- estFile %>% rownames_to_column(var ="run_accession")
    estFile
  }else{
    print(paste0("No tumor purity data for ", tissue))
    estFile<- NULL
    
  }
  
}

# Aggregate to compact immune populations
sumColumns <- function(DF,aggregPopName,sumColsGroupList){
  # aggregPopName<-"Bcells"
  # DF <-resCiber_abs.c
  # sumColsGroupList <-compactPopGroupsList
  # print(sumColsGroupList)
  colsSelect<- as.character(sumColsGroupList)
  if(length(colsSelect)==1){
    DF.aggreg <- as.data.frame(DF[,colsSelect])
  }else{
    DF.aggreg <- as.data.frame(rowSums(DF[,colsSelect]))
  }
  
  colnames(DF.aggreg) <-aggregPopName
  DF.aggreg
}

# Function should be applied by row
newFractions <- function(row, newTotal){
  elem <- (100*row)/newTotal
  elem
}

get_clinicalTCGA <- function(project, tcgaDF, barcodeCol){
  # project <- "SKCM"
  # tcgaDF <- tcga.tp.Filt.final
  # barcodeCol <-"bcr_patient_barcode"
  print(paste0("Getting drug clinical data for: TCGA-",project))
  tcgaDF.myproj <- tcga.tp.Filt.final %>% filter(project==project) %>% droplevels()
  query <- GDCquery(project = paste0("TCGA-",project), 
                    data.category = "Clinical", 
                    file.type = "xml", 
                    barcode = tcgaDF.myproj[[barcodeCol]])
  GDCdownload(query)
  clinical.myproj <- GDCprepare_clinic(query, clinical.info = "drug")
  clinical.myproj
}