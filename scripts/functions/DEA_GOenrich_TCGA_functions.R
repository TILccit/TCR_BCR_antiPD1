#################
## FUNCITONS ##
###

test_match_order <- function(x,y) {
  
  if (all(x==y)) print('Perfect match in same order')
  
  if (!all(x==y) && all(sort(x)==sort(y))) print('Perfect match in wrong order')
  
  if (!all(x==y) && !all(sort(x)==sort(y))) print('No match')
}

# test_match_order(samples.mel$run_accession,y)
`%notin%` <- Negate(`%in%`)

# Try with TCGA htseq counts object
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

select_columns <- function(DF,columnNames){
  DF.s <- DF %>% dplyr::select(columnNames)
  colnames(DF.s) <-columnNames
  DF.s
}
get_patient_num <- function(sample_id){
  patient_num <- str_extract(sample_id, "[^_]+")
  patient_num
}

#############
findDEGs <- function(ctsDF,coldataDF,condA,condB,plotTitle,lfcThres,padjThres,fileName, degCol,
                     controlDataset=FALSE,tcga=FALSE, outpath){
  # ctsDF <-txi.mel
  # coldataDF <-samples.mel
  # lfcThres <- 0
  # padjThres <- 0.05
  # condA <-"low"
  # condB <-"high"
  # plotTitle <-"BCR high vs BCR low"
  # fileName="melDegs.dataset_lfc0"
  # controlDataset=TRUE
  # annotDF <-annot
  # annotCol.ens <-"ensembl_gene_id"
  # annotCol.mgi <- "mgi_symbol"
  
  # ctsDF <-skcm$counts
  # coldataDF = skcm$coldata
  # condA="low"
  # condB="high"
  # plotTitle="PD-L1 high vs PD-L1 low"
  # lfcThres = 0
  # padjThres=0.05
  # fileName="SKCM.Degs.dataset_lfc0"
  # degCol = "PDL1"
  # controlDataset=FALSE
  # tcga = TRUE
  
  
  degCol.un <- as.name(degCol)
  degColform <- paste0("~ ", degCol)
  controlForm <-paste0("~ ", "dataset + ",degCol)
  ####################
  if(isTRUE(tcga)){
    dds <- DESeqDataSetFromMatrix(round(ctsDF), colData = coldataDF, as.formula(degColform))
  }else{
    if(isTRUE(controlDataset)){
      dds <- DESeqDataSetFromTximport(ctsDF, colData = coldataDF, as.formula(controlForm))# + dataset
    }else{
      dds <- DESeqDataSetFromTximport(ctsDF, colData = coldataDF, as.formula(degColform))
    }
  }
  
  # dds$Condition <- factor(x = dds$Condition, levels = c(condB,condA))
  dds[[degCol]] <- factor(x = dds[[degCol]], levels = c(condA,condB))
  
  # # Normalizing with vsd
  # vsd <- vst(dds, blind = FALSE)
  # head(assay(vsd), 3)
  # colData(vsd)
  # Normalizing with rlog
  print("Normalizing data, with DESeq2 vst() function")
  vsd <- vst(dds, blind = FALSE)
  # head(assay(rld), 3)
  
  # Here I want to have gene symbol not refseq ID
  ## Running the differential expression 
  ### What DESeq function does:
  # the estimation of size factors (controlling for differences in the
  # sequencing depth of the samples), the estimation of
  # dispersion values for each gene, and fitting a generalized linear model
  print("== Running the differential expression ==")
  dds <- estimateSizeFactors(dds)
  idx <- rowSums(counts(dds)) >= 10#rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3
  dds <- dds[idx,]
  dds <- DESeq(dds)
  
  colData(dds)
  resultsNames(dds)
  # Extracting the results
  print("== Extracting results ==")
  res <- results(dds, contrast=c(degCol,condB,condA),lfcThreshold=lfcThres,alpha=padjThres)
  # since we define contrast the comparison will be condB compared to condA,
  #i.e. up/downregulated in condB compared to A
  # res
  summary(res)
  # resultsNames(dds)
  print("== LFC shrink ==")
  
  if(isTRUE(controlDataset)){
    res.shrink <- lfcShrink(dds=dds, coef=resultsNames(dds)[length(resultsNames(dds))], type="apeglm", res=res)
    
  }else{
    res.shrink <- lfcShrink(dds=dds, coef=resultsNames(dds)[2], type="apeglm", res=res)
    
  }
  
  
  # mcols(res.shrink, use.names = TRUE)
  # summary(res.shrink)
  
  res.shrink.df <- as.data.frame(res.shrink)
  ## FDR < 0.1, |LFC| > 0
  # padj.cutoff <- 0.1 ; lfc.cutoff <- 0
  res.shrink.df.degs <- res.shrink.df %>% dplyr::filter(padj < padjThres) #%>% dplyr::filter(abs(log2FoldChange) >=lfcThres)
  res.shrink.df.degs <- res.shrink.df.degs %>% rownames_to_column(var="gene")
  res.shrink.df <- res.shrink.df %>% rownames_to_column(var="gene")
 
  print("== Exporting Differential expression analysis results, ")
  # print("vst normalized expression values to use for heatmaps/PCA plots... ")
  # print(" the volcano plot ==")
  # list(res.shrink.df,vsd,res.shrink.df.degs)
  write.xlsx(res.shrink.df.degs, file = paste0(outpath,"/",fileName,".xlsx")
  )
}


# function that loads degs table and creates ordered list and table with up down
get_geneLists <- function(pathTofile,fileName, lfcThres=0){
  # xlsx filename
  # fileName <- "melDegs.dataset_lfc0"
  # pathTofile <- figuresPath
  
  degs <- read_excel(paste0(pathTofile,"/",fileName,".xlsx"))
  # Get ensembl and entrez ids
  degs$ensembl <- geneid2symb$gene_id[match(degs$gene,geneid2symb$gene_name)]
  degs$entrez <- mapIds(x = org.Hs.eg.db,
                        keys = degs$ensembl,
                        column = "ENTREZID",
                        keytype = "ENSEMBL",
                        multiVals = "first")
  # Get ordered gene list
  geneList <-  degs$log2FoldChange
  names(geneList) <- degs$gene
  geneList <- geneList[order(geneList, decreasing = TRUE)]
  # Get ordered gene list entrez
  # Get ordered gene list
  geneList.entz <-  degs$log2FoldChange
  names(geneList.entz) <- degs$entrez
  geneList.entz <- geneList.entz[order(geneList.entz, decreasing = TRUE)]
  ## Dataframe, ordered and annotated
  degs.df <- as.data.frame(degs)
  colnames(degs.df)[3] <- "FC"
  degs.df$group <- ifelse(degs.df$FC >lfcThres, "upregulated","downregulated")
  degs.df <- degs.df[order(degs.df$FC,decreasing = TRUE),]
  list(gL = geneList,
       gL.entrez = geneList.entz,
       gL.df = degs.df)
}

# Enrichment analysis function
oraGSEA <- function(resObj){
  # resObj <-mel.res
  
  dataDF <- resObj$gL.df
  gL <- resObj$gL
  gL.ent <- resObj$gL.entrez
  genes <- resObj$gL.df$gene
  
  # general enrichment GO, BP
  bp.go <- compareCluster(gene~group, data=dataDF, fun="enrichGO",keyType = "SYMBOL",ont = "BP", 
                          OrgDb = org.Hs.eg.db,pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff = 0.05)#,universe = allg.c7
  
  bp.go.dot <- dotplot(bp.go,showCategory=15) # GOOD
  # Enrich based on c2 genes
  bp.go.c2 <- compareCluster(gene~group, data=dataDF, fun="enrichGO",universe = allg.c2,keyType = "SYMBOL",ont = "BP", 
                             OrgDb = org.Hs.eg.db,pAdjustMethod = "BH",
                             pvalueCutoff  = 0.05,
                             qvalueCutoff = 0.05)#
  
  bp.go.c2.dot <- dotplot(bp.go.c2,showCategory=15)
  
  # c2.enrich <- enricher(genes, TERM2GENE=msig_c2,pvalueCutoff = 0.05,
  #                       pAdjustMethod = "BH",
  #                       universe = allg.c2,
  #                       minGSSize = 2)
  # c2.enrich.dot <-dotplot(c2.enrich ,showCategory=15)
  
  ## GSEA
  gse.go.all <- gseGO(geneList=gL, 
                      keyType = "SYMBOL",
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",
                      minGSSize = 2,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05
  )
  
  
  c2.gsea <- GSEA(gL, TERM2GENE = msig_c2)
  
  if(nrow(gse.go.all)==0 & nrow(c2.gsea)==0){
    list(go.bp = bp.go.dot,
         go.bp.c2 = bp.go.c2.dot
         # go.bp.c2.v2 = c2.enrich.dot
    )
    # gse.go.dot = gse.go.all.dot,
    # gse.go.ridge = gse.go.all.ridge,
    # c2.gsea.dot = c2.gsea.dot,
    # c2.gsea.ridge)
    
  }else if(nrow(gse.go.all)==0 | nrow(c2.gsea)==0){
    if(nrow(gse.go.all)==0){
      c2.gsea.dot <-dotplot(c2.gsea, showCategory=20)
      c2.gsea.ridge <- ridgeplot(c2.gsea,showCategory = 20)
      list(go.bp = bp.go.dot,
           go.bp.c2 = bp.go.c2.dot,
           # go.bp.c2.v2 = c2.enrich.dot,
           # gse.go.dot = gse.go.all.dot,
           # gse.go.ridge = gse.go.all.ridge,
           c2.gsea.dot = c2.gsea.dot,
           c2.gsea.ridge = c2.gsea.ridge)
    }else{
      gse.go.all.dot <- dotplot(gse.go.all, showCategory=20) # GOOD
      gse.go.all.ridge <-ridgeplot(gse.go.all,showCategory = 20)
      list(go.bp = bp.go.dot,
           go.bp.c2 = bp.go.c2.dot,
           # go.bp.c2.v2 = c2.enrich.dot,
           gse.go.dot = gse.go.all.dot,
           gse.go.ridge = gse.go.all.ridge)
      # c2.gsea.dot = c2.gsea.dot,
      # c2.gsea.ridge)
    }
  }else{
    gse.go.all.dot <- dotplot(gse.go.all, showCategory=20) # GOOD
    gse.go.all.ridge <-ridgeplot(gse.go.all,showCategory = 20)
    c2.gsea.dot <-dotplot(c2.gsea, showCategory=20)
    c2.gsea.ridge <- ridgeplot(c2.gsea,showCategory = 20)
    list(go.bp = bp.go.dot,
         go.bp.c2 = bp.go.c2.dot,
         # go.bp.c2.v2 = c2.enrich.dot,
         gse.go.dot = gse.go.all.dot,
         gse.go.ridge = gse.go.all.ridge,
         c2.gsea.dot = c2.gsea.dot,
         c2.gsea.ridge = c2.gsea.ridge)
  }
  
  
}

## load TCGA count data ##
# selectedcancer <- "SKCM"
loadTCGAcounts <- function(selectedcancer,tcgaRnaPath,dataDF, filterVar,sampleBarcode,biomarker,biomarkerName,geneid2symb){
  # selectedcancer <- "SKCM"
  # tcgaRnaPath <- paste0(rootDir,tcgaInputData,htseqCounts_Path)
  # dataDF = tcgaData.tp.sub
  # filterVar = "project"
  # sampleBarcode = "bcr_patient_barcode"
  # biomarker=tcgaOrigCol
  # biomarkerName = newCol
  
  filterVar.un <- as.name(filterVar)
  sampleBarcode.un <- as.name(sampleBarcode)
  biomarker.un <- as.name(biomarker)
  biomarkerName.un <- as.name(biomarkerName)
  print(paste0(selectedcancer,"_RNASeqCounts.RData"))
  fname <- file.path(tcgaRnaPath ,
                     paste0("TCGA-",selectedcancer,"_RNASeqCounts.RData"))
  load(fname)
  data<- df_all
  data.tr <- t(data)
  data.tr.df <- data.tr %>% as.data.frame() %>% rownames_to_column("sampleID")
  data.tr.df$patientID <- substr(data.tr.df$sampleID,1,12)
  samples.tcga <- dataDF %>% dplyr::filter(!!filterVar.un==selectedcancer) %>% droplevels()
  data.tr.df.sub <- data.tr.df %>% dplyr::filter(patientID %in% samples.tcga[[sampleBarcode.un]]) %>% droplevels()
  data.tr.df.sub.f <- data.tr.df.sub %>% column_to_rownames(var="sampleID")
  which(duplicated(data.tr.df.sub.f$patientID))
  print(paste0("Initial TCGA patients: ", nrow(data.tr.df.sub.f)))
  # dATA
  DF.filt <- dataDF %>% dplyr::filter(!!filterVar.un==selectedcancer) %>% droplevels()
  biomarker.thres <- median(DF.filt[[biomarker]], na.rm = TRUE);biomarker.thres
  DF.filt[[biomarkerName]] <-ifelse(DF.filt[[biomarker]] >=biomarker.thres, "high","low")
  # Extract patient barcodes of batients with BCR info
  biomarker.samps <- DF.filt %>% dplyr::select(!!sampleBarcode.un,!!biomarkerName.un) %>% dplyr::filter(complete.cases(.)) %>% pull(!!sampleBarcode.un)
  
  toMatch <- biomarker.samps
  matches2 <- unique (grep(paste(toMatch,collapse="|"), 
                           rownames(data.tr.df.sub.f), value=TRUE))
  
  # Subset to sampels with BCR info
  data.tr.df.sub.f.b <- data.tr.df.sub.f[matches2,]
  which(duplicated(data.tr.df.sub.f.b$patientID))
  
  if(length(which(duplicated(data.tr.df.sub.f.b$patientID)))==0){
    data.tr.df.sub.f.b1 <- data.tr.df.sub.f.b[,]
    
  }else{
    data.tr.df.sub.f.b1 <- data.tr.df.sub.f.b[-which(duplicated(data.tr.df.sub.f.b$patientID)),]
    
  }
  rownames(data.tr.df.sub.f.b1) <- NULL
  data.final <- data.tr.df.sub.f.b1 %>% column_to_rownames(var="patientID") %>% t() %>% as.data.frame()
  # rownames(data.blca) <- geneid2symb$gene_name[match(rownames(data.blca),geneid2symb$gene_id)]
  print(paste0("Final TCGA patients: ", ncol(data.final)))
  
  data.final$geneName <- geneid2symb$gene_name[match(rownames(data.final),geneid2symb$gene_id)]
  # Gene Ids to Gene names
  data.final.agg<- setDT(data.final)[, lapply(.SD, mean, na.rm=TRUE) , by = .(geneName)]
  data.final.agg <- data.final.agg %>% column_to_rownames(var="geneName")
  # data.tr.df.sub.f.b2 <- data.tr.df.sub.f.b[which(duplicated(data.tr.df.sub.f.b$patientID)),]
  # dupl <- setDT(data.tr.df.sub.f.b2)[, lapply(.SD, median, na.rm=TRUE) , by = .(patientID)]
  # 
  # data.tr.df.sub.f.agg <- data.tr.df.sub.f %>% group_by(patientID) %>% dplyr::summarize(across(everything(),median))
  # which(duplicated(data.tr.df.sub.f.agg$patientID))
  # 
  # Now make metadata table
  df.blca.final <- DF.filt %>% dplyr::select(!!sampleBarcode.un,!!biomarkerName.un) %>% dplyr::filter(complete.cases(.)) %>% droplevels()
  # Order
  df.blca.final.ord <- df.blca.final[match(colnames(data.final.agg), df.blca.final[[sampleBarcode]]),]
  
  list(counts=data.final.agg,
       coldata=df.blca.final.ord)
}
