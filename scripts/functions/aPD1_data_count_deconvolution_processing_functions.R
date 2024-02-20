###############
## FUNCTIONS ##
###############
`%notin%` <- Negate(`%in%`)

# Claculate FPKM-UQ normalized values from counts (htseq)
fpkm_uq <- function(counts, lengths) {
  rate <- counts / lengths
  uq <- quantile(counts)
  g_uq <-as.numeric(uq[4])
  #rate / g_uq * 1e6
  fpkm_uq <-(rate / g_uq)* 1e7
  fpkm_uq
}


# # Claculate FPKM-UQ normalized values from counts (htseq)
# fpkm_uq <- function(counts, lengths) {
#   rate <- counts / lengths
#   uq <- quantile(counts)
#   g_uq <-as.numeric(uq[4])
#   #rate / g_uq * 1e6
#   fpkm_uq <-(rate / g_uq)* 1e7
#   fpkm_uq
# }

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

get_mixture_from_tpm <- function(pathToRData,pathToData,fileName,RDataSuffix,signatureMatrix){
  ## Function that loads RData object of Tpm (txi abundance) and
  ## extracts the mixture matrix ready to be used in cibersort.
  print(paste0("Processing: ",fileName))
  filehandle <- paste0(pathToRData,fileName,RDataSuffix)
  # Load file
  tpm <- loadRData(filehandle)
  print(paste0("Total samples: ",dim(tpm)[2]))
  mixture <- tpm %>% as.data.frame() %>% rownames_to_column(var="gene") %>% dplyr::filter(gene %in% signatureMatrix[["Gene symbol"]])
  print(paste0("Total common genes between mixture & signature matrix: ",dim(mixture)[1], " out of ", dim(signatureMatrix)[1], " genes"))
  ## Add column
  mixture <- mixture %>% column_to_rownames(var="gene") %>% rownames_to_column(var="!Sample_title")
  # Write file
  write_tsv(mixture, file.path(paste0(pathToData, fileName,"mixture_matrix.txt")))
  mixture
}

# Change first seed in seeds
changeSeed1 <- function(SeedList,SeedVector,ind){
  SeedList[[ind]][1] <- SeedVector[ind]
  SeedList
}

run_absModeCibersort <- function(SignMatrix,MixMatrix,permutations=100,seed){
  # # set.seed(seed)
  # assign(".Random.seed", seed, envir=.GlobalEnv)
  res <-CIBERSORT(SignMatrix, MixMatrix, perm=permutations, QN=FALSE, absolute=TRUE, abs_method="sig.score")
  res
}

run_relModeCibersort <- function(SignMatrix,MixMatrix,permutations=100,seed){
  # # set.seed(seed)
  # assign(".Random.seed", seed, envir=.GlobalEnv)
  res <-CIBERSORT(SignMatrix, MixMatrix, perm=permutations, QN=FALSE, absolute=FALSE)
  res
}


# This is based on the clusterSetRNGStream function from
# the parallel package, copyrighted by The R Core Team
getseeds <- function(ntasks, iseed) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(iseed)
  seeds <- vector("list", ntasks)
  seeds[[1]] <- .Random.seed
  for (i in seq_len(ntasks - 1)) {
    seeds[[i + 1]] <- nextRNGSubStream(seeds[[i]])
  }
  seeds
}


# Function that adds rownames as column in each df
addRownames<- function(df){
  df <- df %>% as.data.frame() %>% rownames_to_column(var="sampleID")
  df <- df[,1:26]
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

niceBoxplots_pairwise <- function(DF, xVar, yVar, mytitle){
  xVar.un <- as.name(xVar)
  yVar.un <- as.name(yVar)
  
  p <-ggbetweenstats(
    data = DF,
    x = !!xVar.un,
    y = !!yVar.un,
    sort = "descending", 
    pairwise.comparisons = TRUE, # show pairwise comparison test results
    pairwise.display = "significant",
    ggtheme = hrbrthemes::theme_ipsum_tw(),
    outlier.tagging = TRUE,
    #outlier.tagging = TRUE,
    #outlier.label = education,
    ggstatsplot.layer = FALSE,
    # messages = FALSE,
    title = mytitle,
    results.subtitle = FALSE,
    bf.message = FALSE,
    caption=NULL
  )
  
  p <-p + theme_ipsum(base_family = "Palatino Linotype", base_size = 14)
  p <- p + theme(
    # panel.border = element_rect(colour = "black", size=0.2),
    # panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(size=14),
    axis.title.x= element_text(size=16),
    #axis.text.x = element_text(size=10),
    axis.text.y = element_text(size = 14),
    axis.title.y= element_text(size = 16),
    # legend.key = element_rect(size = 2),
    # # legend.key.height=unit(0.6,"line"),
    # # legend.key.width=unit(1,"line"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    # legend.background = element_rect(colour = "black"),
    # legend.direction='horizontal',legend.box='horizontal'
  )
  p
}

niceFacetBoxplots_pairwise <- function(DF, xVar, yVar, facetVar,mytitle){
  xVar.un <- as.name(xVar)
  yVar.un <- as.name(yVar)
  facetVar.un <- as.name(facetVar)
  
  p <-grouped_ggbetweenstats(
    data = DF,
    x = !!xVar.un,
    y = !!yVar.un,
    # sort = "descending", 
    grouping.var = !!facetVar.un,
    pairwise.comparisons = TRUE, # show pairwise comparison test results
    pairwise.display = "significant",
    ggtheme = hrbrthemes::theme_ipsum_tw(),
    # outlier.tagging = TRUE,
    #outlier.tagging = TRUE,
    #outlier.label = education,
    ggstatsplot.layer = FALSE,
    messages = FALSE,
    results.subtitle = FALSE,
    bf.message = FALSE,
    caption=NULL
  )
  
  p <-p + theme_ipsum(base_family = "Palatino Linotype", base_size = 14)
  p <- p + theme(
    # panel.border = element_rect(colour = "black", size=0.2),
    # panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(size=14),
    axis.title.x= element_text(size=16),
    #axis.text.x = element_text(size=10),
    axis.text.y = element_text(size = 14),
    axis.title.y= element_text(size = 16),
    # legend.key = element_rect(size = 2),
    # # legend.key.height=unit(0.6,"line"),
    # # legend.key.width=unit(1,"line"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    # legend.background = element_rect(colour = "black"),
    # legend.direction='horizontal',legend.box='horizontal'
  )
  p
}

compareScatter <- function(DF, xVar,yVar,groupVar=NULL, facetVar =NULL,plotTitle,xLab, yLab){
  # DF <-samples_meta.Abs.compact
  # xVar <-"CD8cells"
  # yVar <-"ImmuneScore"
  # groupVar <- "response"
  # facetVar <- "tissue"
  # xLab <- "CD8 T cells absolute score"
  # yLab <- "Immune score"
  # plotTitle <- ""
  # 
  xVar.un <- as.name(xVar)
  yVar.un <- as.name(yVar)
  
  
  #, label=!!labelVar.un
  if(!is.null(groupVar)){
    groupVar.un <- as.name(groupVar)
    p <-ggplot(DF, aes(x=!!xVar.un, y=!!yVar.un, fill=!!groupVar.un))
    
  }else{
    p <-ggplot(DF, aes(x=!!xVar.un, y=!!yVar.un))
    
  }
  p <- p + geom_point(shape=21, size=2,na.rm = TRUE,show.legend=TRUE,
                      alpha = 0.8)
  p <- p + theme_bw(base_family = "Palatino Linotype", base_size = 12)
  p <- p + theme(legend.position = "bottom",
                 panel.border = element_rect(colour = "black", size=0.2),
                 panel.background = element_rect(fill = "white"),
                 #axis.text.x = element_text(angle = 40, hjust = 1,size=6),
                 axis.text.x = element_text(size=12),
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
    p <- p + facet_wrap(as.formula(paste0("~ ",facetVar)))
  }
  p <- p + ylab(paste0(yLab))
  p <- p + xlab(paste0(xLab))
  p <- p + ggtitle(plotTitle)
  
  print(p)
  # p
}




# EStimate purity
calculateEstimateScore <- function(pathToRnaData,RnaFileSuffix,pathToData,tissue){
  # pathToRnaData<-antiPDL1publicIntermediateData
  # RnaFileSuffix <- vst.batch.dataset.RData
  # pathToData<-paste0(antiPDL1publicIntermediateData,"purity/")
  # tissue <- "melanoma"
  print(paste0("Calculating stromal, immune and estimate scores for: ", tissue))
  # Load rnaseq vst data
  vst.df <- loadRData(paste0(pathToRnaData,tissue,"_ALLresp",RnaFileSuffix))
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
}

readEstimateFiles <- function(pathToData,tissue){
  print(paste0("Reading in stromal, immune and estimate scores for: ", tissue))
  # pathToData <-projectDataPath
  # tissue <-tissues[1]
  estFile <-read.table(paste0(pathToData,tissue,"_estimate",".score.gct"), header = TRUE, sep = "\t", dec = ".", skip =2)
  estFile <- estFile %>% dplyr::select(-Description) %>% column_to_rownames(var="NAME")
  estFile <- t(estFile)
  estFile <- as.data.frame(estFile)
  estFile <- estFile %>% rownames_to_column(var ="run_accession")
  estFile
}