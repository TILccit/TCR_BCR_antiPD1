# library(pacman)
# pacman::p_load(ggthemes,dichromat,RColorBrewer,scales)

# display.brewer.all(colorblindFriendly = T)

#palette using grey
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#palette using black
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


myCBfriendly.dark.small <- c("#222255","#225555","#225522","#666633","#663333","#555555")
# show_col(myCBfriendly.dark.small)

myCBfriendly.pale.small <- c("#bbccee","#cceeff","#ccddaa","#eeeebb","#ffcccc","#dddddd")
# show_col(myCBfriendly.pale.small) 

myCBfriendly.light.medium <- c("#77aadd","#99ddff","#44bb99","#bbcc33","#aaaa00","#eedd88","#ee8866","#ffaabb","#dddddd")
# show_col(myCBfriendly.light.medium)

myCBfriendly.bright.small <- c("#4477aa","#66ccee","#228833","#ccbb44","#ee6677","#bbbbbb")
# show_col(myCBfriendly.bright.small)

myCBfriendly.highContr.small <- c("#ffffff","#ddaa33","#bb5566","#004488","#000000")
# show_col(myCBfriendly.highContr.small)

myCBfriendly.muted.medium<- c("#332288","#88ccee","#44aa99","#117733","#999933","#ddcc77","#cc6677","#882255","#aa4499","#dddddd")
# show_col(myCBfriendly.muted.medium)


# show_col(c(myCBfriendly.highContr.small,myCBfriendly.light.medium))

myCBfriendly.large <-c(brewer.pal(9,"Set2"),brewer.pal(8,"Dark2"),"#cc3311")
# show_col(myCBfriendly.large)


##########################
## Loading RData object ##
##########################


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


############################
## Uppercase first letter ##
############################
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


`%notin%` <- Negate(`%in%`)

### PROCESSING DATA ###

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

barplot_nominal <- function(DF,variable,angle, yLabel){
  df<-as.data.frame(table(DF[[variable]]))
  p <- ggplot(df, aes(x=Var1, y=Freq, fill=Var1))
  p <- p + geom_bar(stat="identity", color="black")
  p <- p+scale_fill_brewer(palette="Dark2")
  p <- p + theme_ipsum(base_family = "Palatino Linotype", base_size = 10)
  p <- p + theme(axis.text.x = element_text(angle = angle, hjust = 1, vjust = 1,size = 10),
                 axis.title.x = element_text(size=10, face="bold"),
                 axis.text.y = element_text( hjust = 1, vjust = 1,size = 10),
                 axis.title.y = element_text(size=11, face="bold"),
                 legend.title = element_text(size=11,color="black"),
                 legend.text = element_text(size = 10),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"))
  p <- p + labs(fill=variable,
                title = paste0(""),
                x = variable,
                y = yLabel
  )
  p
}

barplot_nominal.groupedResponse <- function(DF,variable,angle,mode,yLabel,depVar){
  variable.un <- as.name(variable)
  depVar.un <- as.name(depVar)
  df <- DF %>% dplyr::select(!!depVar.un,!!variable.un) %>% group_by(!!variable.un,!!depVar.un) %>% dplyr::summarise(counts=length(!!depVar.un))
  df <- df %>%  mutate(prop = round(counts*100/sum(counts), 1),
                       lab.ypos = cumsum(prop) - 0.5*prop)
  
  
  p <- ggplot(df, aes_string(fill=depVar, y=mode, x=variable))
  p <- p + geom_bar(position="dodge", stat="identity")
  p <- p+scale_fill_brewer(palette="Dark2")
  p <- p + theme_ipsum(base_family = "Palatino Linotype", base_size = 10)
  p <- p + theme(axis.text.x = element_text(angle = angle, hjust = 1, vjust = 1,size = 10),
                 axis.title.x = element_text(size=10, face="bold"),
                 axis.text.y = element_text( hjust = 1, vjust = 1,size = 10),
                 axis.title.y = element_text(size=11, face="bold"),
                 legend.title = element_text(size=11,color="black"),
                 legend.text = element_text(size = 10),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"))
  p <- p + labs(fill=variable,
                title = paste0(""),
                x = depVar,
                y = yLabel
  )
  p  
}


barplot_nominal.stackedResponse <- function(DF,variable,angle,mode,yLabel,depVar){
  variable.un <- as.name(variable)
  depVar.un <- as.name(depVar)
  df <- DF %>% dplyr::select(!!depVar.un,!!variable.un) %>% group_by(!!variable.un,!!depVar.un) %>% dplyr::summarise(counts=length(!!depVar.un))
  df <- df %>%  mutate(prop = round(counts*100/sum(counts), 1),
                       lab.ypos = cumsum(prop) - 0.5*prop)
  
  
  p <- ggplot(df, aes_string(fill=depVar, y=mode, x=variable))
  p <- p + geom_bar(position="stack", stat="identity")
  p <- p+scale_fill_brewer(palette="Dark2")
  p <- p + theme_ipsum(base_family = "Palatino Linotype", base_size = 10)
  p <- p + theme(axis.text.x = element_text(angle = angle, hjust = 1, vjust = 1,size = 10),
                 axis.title.x = element_text(size=10, face="bold"),
                 axis.text.y = element_text( hjust = 1, vjust = 1,size = 10),
                 axis.title.y = element_text(size=11, face="bold"),
                 legend.title = element_text(size=11,color="black"),
                 legend.text = element_text(size = 10),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"))
  p <- p + labs(fill=variable,
                title = paste0(""),
                x = depVar,
                y = yLabel
  )
  p  
  
}

hist_measurement <- function(DF, variable){
  p<-ggplot(DF, aes_string(x=variable, color=variable, fill=variable)) +
    geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
    geom_density(alpha=0.4) 
  p <- p + scale_color_brewer(palette="Dark2")
  p <- p + scale_fill_brewer(palette="Dark2")
  p <- p + theme_ipsum(base_family = "Palatino Linotype", base_size = 10)
  p <- p + theme(axis.text.x = element_text(hjust = 1, vjust = 1,size = 10),
                 axis.title.x = element_text(size=10, face="bold"),
                 axis.text.y = element_text( hjust = 1, vjust = 1,size = 10),
                 axis.title.y = element_text(size=11, face="bold"),
                 legend.title = element_text(size=11,color="black"),
                 legend.text = element_text(size = 10),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"),
                 legend.position="top")
  p
}

hist_measurement_byDependent <- function(DF, variable, depVar){
  p<-ggplot(DF, aes_string(x=variable, color=depVar, fill=depVar)) +
    geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
    geom_density(alpha=0.4) 
  p <- p + scale_color_brewer(palette="Dark2")
  p <- p + scale_fill_brewer(palette="Dark2")
  p <- p + theme_ipsum(base_family = "Palatino Linotype", base_size = 10)
  p <- p + theme(axis.text.x = element_text(hjust = 1, vjust = 1,size = 10),
                 axis.title.x = element_text(size=10, face="bold"),
                 axis.text.y = element_text( hjust = 1, vjust = 1,size = 10),
                 axis.title.y = element_text(size=11, face="bold"),
                 legend.title = element_text(size=11,color="black"),
                 legend.text = element_text(size = 10),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"),
                 legend.position="top")
  p
}


cors <- function(df) { 
  # turn all three matrices (r, n, and P into a data frame)
  M <- Hmisc::rcorr(as.matrix(df))
  # return the three data frames in a list return(Mdf)
  Mdf <- map(M, ~data.frame(.x))
}

formatted_cors <- function(df){
  cors(df) %>%
    map(~rownames_to_column(.x, var="measure1")) %>%
    map(~pivot_longer(.x, -measure1, "measure2")) %>% 
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    mutate(sig_p = ifelse(P < .05, T, F), p_if_sig = ifelse(P <.05, P, NA), r_if_sig = ifelse(P <.05, r, NA)) 
}


corHeatmap <- function(DF,groupVar,selectedGroup,selectedCols){
  # DF <- immdata.apd1.proc.TcrA.divEst.full.sub
  # groupVar <- "tissue"
  # selectedGroup <- "bladder"
  # selectedCols <- measure.vars
  cat('\n')  
  cat('\n') 
  cat("### ", selectedGroup, " \n") 
  cat('\n') 
  
  groupVar.un <- as.name(groupVar)
  
  DF.f <- DF %>% dplyr::filter(!!groupVar.un==selectedGroup) %>% droplevels()
  
  DF.f <- DF.f  %>% column_to_rownames(var="run_accession") %>% dplyr::select(all_of(selectedCols))
  
  p<-formatted_cors(DF.f) %>% 
    ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
    geom_tile() +
    labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation", title=paste0("Correlations in ",selectedGroup), subtitle="Only significant Pearson's correlation coefficients shown") + scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
    geom_text() +
    theme_bw(base_family = "Palatino Linotype", base_size = 12) +
    theme(legend.position = "bottom",
          panel.border = element_rect(colour = "black", size=0.2),
          panel.background = element_rect(fill = "white"),
          #axis.text.x = element_text(angle = 40, hjust = 1,size=6),
          axis.text.x = element_text(size=12,angle = 25, hjust = 1),
          axis.text.y = element_text(size = 12),
          legend.key = element_rect(size = 2),
          legend.key.height=unit(0.6,"line"),
          legend.key.width=unit(1,"line"),
          legend.text = element_text(size = 12),
          legend.background = element_rect(colour = "black"),
          legend.direction='horizontal',legend.box='horizontal'
    ) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0))
  print(p)
}


corHeatmap2 <- function(matDF){
  # meltedMat <- melted_mat
  # matDF <-upper_tri
  # matDF <-melted_mat
  #paste(round(value,3), c(" ","*")[(abs(value) <= .05)+1])
  p<-matDF %>% 
    # ggplot(aes(Var1, Var2, fill=value, label=paste(c(" ","***")[(abs(matDF$value) <= .05)+1]))) +
    # ggplot(aes(Var1, Var2, fill=value, label=paste(round(value,3), c(" ","*")[(abs(value) <= .05)+1]))) +
    ggplot(aes(Var1, Var2, fill=value, label=ifelse(matDF$value <= .001, "***",ifelse(matDF$value <= .01, "**", ifelse(matDF$value <= .05, "*",round(matDF$value,3)))))) +
    geom_tile() +
    labs(x = NULL, y = NULL, fill = "P-value\nPairwise\n ROC-AUC") +
    geom_text() +
    scale_fill_gradient2(mid="#0C6291",low="#FBFEF9",high="#A63446", limits=c(0,1)) + #
    theme_bw(base_family = "Palatino Linotype", base_size = 11) +
    theme(legend.position = "bottom",
          panel.border = element_rect(colour = "black", size=0.2),
          panel.background = element_rect(fill = "white"),
          #axis.text.x = element_text(angle = 40, hjust = 1,size=6),
          axis.text.x = element_text(size=12,angle = 25, hjust = 1),
          axis.text.y = element_text(size = 12),
          legend.key = element_rect(size = 8),
          # legend.key.height=unit(0.6,"line"),
          # legend.key.width=unit(1,"line"),
          legend.text = element_text(size = 9),
          legend.background = element_rect(colour = "black"),
          legend.direction='horizontal',legend.box='horizontal'
    ) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0))
  print(p)
}

pca_biPlot <- function(DF,response.var,nominal.vars,measure.vars,exclude.vars,fillVar,fillLabel,axes=NULL,title = ""){
  # DF <-data_full
  # response.var <-responseVar
  # nominal.vars <-NULL
  # measure.vars <-covariatesIn
  # exclude.vars <-NULL
  # fillVar <-as.factor(data_full$response)
  # fillLabel <-"response"
  # axes=c(1,2)
  # title = ""
  
  pca_biplot <- function(pcaObj,eig.val, fillVar, fillLabel,axes){
    mybiplot <- fviz_pca_biplot(pcaObj, 
                                geom.ind = "point",
                                fill.ind = fillVar,
                                col.ind = "black",
                                alpha.ind = 0.75,
                                pointshape = 21, 
                                pointsize = 3,
                                axes=axes,
                                #col.ind = as.factor(blca_df$ERCC2_def), 
                                palette = c("#A9A9A9", "#008000"), #"jco", 
                                # addEllipses = TRUE, 
                                label="none",
                                invisible ="var",
                                # label = "var",
                                # # #alpha.var ="contrib", 
                                # col.var = "contrib",
                                # # #col.var = "black",
                                # # #gradient.cols = "BuPu",
                                repel = TRUE,
                                # habillage="none",
                                legend.title = list(fill = fillLabel
                                                    ,
                                                    color = fillLabel
                                ),
                                title= title
                                #,
                                #col.quanti.sup="Steel Blue"
                                #,invisible = "quanti.sup"
    ) + 
      theme_minimal(base_size = 16, base_family = "Palatino Linotype") + 
      theme(panel.background = element_rect(fill = "white"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size = 16, face = "bold"),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_text(size=16, face="bold"),
            axis.text.y = element_text( size=13),#hjust = 1, vjust = 1,
            axis.title.y = element_text(size=16, face="bold"),
            legend.title = element_text(size=16,color="black",face="bold"),
            #legend.background = element_rect(size=0.1,color="black",fill=alpha("white",0.6)),
            legend.text = element_text(size=13)) + xlab(paste0("PC ",axes[1]," (",round(eig.val[axes[1],2],2)," % )")) + ylab(paste0("PC ",axes[2]," (",round(eig.val[axes[2],2],2)," % )"))
    mybiplot
  }
  
  response.index <- which(colnames(DF) %in% response.var)
  exclude.index<- which(colnames(DF) %in% exclude.vars)
  meas.index <-which(colnames(DF[,-c(response.index,exclude.index)]) %in% measure.vars)
  categ.index <- which(colnames(DF[,-c(response.index,exclude.index)]) %in% nominal.vars)
  ### biplot
  # Provide index of categorical variables as supplementary
  if(length(categ.index)==0){
    pcaObj <- PCA(DF[,-c(response.index,exclude.index)], graph = FALSE, scale.unit = TRUE, ncp=length(meas.index))
  }else{
    pcaObj <- PCA(DF[,-c(response.index,exclude.index)],quali.sup=categ.index, graph = FALSE, scale.unit = TRUE, ncp=length(meas.index))
  }
  
  
  eig.val <- get_eigenvalue(pcaObj)
  
  p<-pca_biplot(pcaObj, eig.val,fillVar, fillLabel,axes = axes)
  p
}

outlierDetect_multiple <- function(DF.X, k=4,C=1,cutoff.dens=.95, cutoff.depth=.05, cutoff.disp=.99,cutoff.m=.9,cutoff.nn=.95,cutoff.u=.95, method= "euclidean"){
  # set.seed(123)
  # dens.out <-dens(DF.X,k=k,C=C, cutoff=cutoff.dens)
  # set.seed(123)
  # depth.out <-depthout(DF.X,cutoff=cutoff.depth)
  # depth.out <-NA
  # set.seed(123)
  # disp.out <-disp(DF.X,cutoff=cutoff.disp)
  # set.seed(123)
  # maha.out <-maha(DF.X,cutoff=cutoff.m)
  # set.seed(123)
  # nn.out <- nn(DF.X,k=k, cutoff=cutoff.nn, Method=method)
  set.seed(123)
  nnk.out <-nnk(DF.X,k=k, cutoff=cutoff.nn,Method=method)
  set.seed(123)
  # inter.out <-OutlierDetection(DF.X, cutoff=cutoff.dens,Method=method, depth = FALSE,dense = FALSE, distance = FALSE, dispersion = FALSE)
  # set.seed(123)
  # pc.out <- PCOutlierDetection(DF.X, cutoff=cutoff.dens,Method=method, depth = FALSE,dense = FALSE, distance = FALSE, dispersion = FALSE)
  # uni.out <- UnivariateOutlierDetection(DF.X,cutoff=cutoff.u,Method="euclidean",rnames=FALSE)
  out.res <- list(
                  # "dens"=dens.out
                  # ,
                  # "depth" = depth.out
                  # ,
                  # "disp" = disp.out
                  # ,
                  # "maha" = maha.out,
                  # "nn" = nn.out,
                  "nnk" = nnk.out
                  # ,
                  # "inter" = inter.out,
                  # "pc" = pc.out
                  # ,
                  # "uni" = uni.out
  )
  out.res
}

pca_ObjCoords.TwoClass <- function(DF,response.var,nominal.vars,measure.vars,exclude.vars){
  # DF<- data_full
  # response.var<-"response"
  # nominal.vars<-c("dataset","tissue")
  # measure.vars<-TuTACK.genes
  # exclude.vars<- colnames(data_full.crpr)[dim(data_full.crpr)[2]]
  # class.1<- "CRPR"
  # class.2 <- "PD"
  # class.1.outlierIndex<-out.crpr$nn$`Location of Outlier`
  # class.2.outlierIndex<-out.pd$nn$`Location of Outlier`
  # axes=c("Dim.1","Dim.2")
  ###
  response.index <- which(colnames(DF) %in% response.var)
  exclude.index<- which(colnames(DF) %in% exclude.vars)
  meas.index <-which(colnames(DF[,-c(response.index,exclude.index)]) %in% measure.vars)
  categ.index <- which(colnames(DF[,-c(response.index,exclude.index)]) %in% nominal.vars)
  ### biplot
  # Provide index of categorical variables as supplementary
  pcaObj <- PCA(DF[,-c(response.index,exclude.index)],quali.sup=categ.index, graph = FALSE, scale.unit = TRUE, ncp=length(meas.index))
  
  pcDT <- as.data.frame(pcaObj$ind$coord)
  pcDT <- as.data.table(pcaObj$ind$coord)
  
  pcDT$patientID <- DF[[exclude.vars]]
  pcDT$response <- DF[[response.var]]
  pcDT
  
}


multiInsideSamplTRAINModels_LogReg.repCV_construct <- function(DFtrain,responseVar,design, ModelName,
                                                               kfolds,reps,subMethod=NULL,inoutSub="in",
                                                               metricSelected,
                                                               seedModel=NA,seedPart,seeds=NULL,
                                                               stratify=FALSE,positiveCoefs=FALSE,applyScaling=FALSE){
  
  # DFtrain = data_train2
  # responseVar = responseVar
  # design = modelDesign
  # ModelName = model.name
  # kfolds = KFOLDS 
  # reps = REPEATS
  # subMethod=subsampMethod
  # inoutSub = inoutSub
  # metricSelected =metricOpt
  # seedModel = seed.model
  # seedPart = seed.partition
  # seeds=setCVseeds 
  # stratify=stratify
  # positiveCoefs=positiveCoefs
  # applyScaling = selectedScale
  
  # DFtrain <-data_train2
  # responseVar <-responseVar
  # design <-modelDesign
  # kfolds <-KFOLDS
  # reps <-REPEATS
  # subMethod=subsampMethod
  # inoutSub=inoutSub
  # metricSelected <-metricOpt
  # seedModel <-seed.model
  # seedPart <-seed.partition
  # seeds=setCVseeds
  # stratify=stratify
  # applyScaling=selectedScale
  # positiveCoefs = FALSE
  ################################
  ## CHECK 4 : MODEL PARAMETERS ##
  ################################
  
  print("==MODEL PARAMETERS==")
  print("=====================")
  print(paste0("Folds: ",kfolds))
  print("-----------------------")
  print(paste0("Repeats: ",reps))
  print("-----------------------")
  print(paste0("Model seed: ",seedModel))
  print("-----------------------")
  print(paste0("Set CV seeds: ",seeds))
  print("-----------------------")
  print(paste0("Subsampling: ",inoutSub))
  print("-----------------------")
  print(paste0("Subsampling Method: ",subMethod))
  print("-----------------------")
  print(paste0("Model metric:",metricSelected))
  print("-----------------------")
  
  print("==========================")
  ##############################
  ## Reproducibility point A ###
  ##############################
  
  # When the models are created inside of resampling, the seeds can also be set. 
  # While setting the seed prior to calling train may guarantee that the same random numbers are used, 
  # this is unlikely to be the case when parallel processing is used (depending which technology is utilized). 
  #To set the model fitting seeds, trainControl has an additional argument called seeds that can be used. 
  # The value for this argument is a list of integer vectors that are used as seeds. 
  if(isTRUE(stratify)){
    print("~Create Stratified folds for CV~")
    set.seed(NULL)
    set.seed(seedModel)
    cvfolds= createMultiFolds(DFtrain[[responseVar]],k=kfolds, times = reps)
  }else{
    cvfolds=NULL
  }
  
  if (! is.null(seeds)){
    # if (is.null(seeds)){
    m<-reps*kfolds
    n.eval=length(lamdaSeq)*length(alpha)
    seeds <- vector(mode = "list", length = m+1)
    set.seed(seedModel)
    for(i in 1:(m)) seeds[[i]] <- sample.int(1000, n.eval)
    ## For the last model:
    seeds[[m+1]] <- sample.int(1000, 1)
  }
  set.seed(NULL)
  set.seed(seedModel)
  
  
  ## Calculate combinations number of tuning parameters.
  # gridCombos <- length(c(round(sqrt(length(feats))/2,0))) * length(round(sqrt(length(feats)),0))* length(round(sqrt(length(feats))*2,0))
  ## DO THE FOLLOWING FOR PARALLEL OPERATION-CONTROL RESAMPLES
  set.seed(NULL)
  set.seed(seedModel)
  m <- kfolds*reps
  seeds <- vector(mode = "list", length = m+1)
  # print(paste0("SEEDS LIST:", length(seeds)))
  # print(paste0("SEEDS ELEMENTS:", length(sizeVec)))
  # 
  for(i in 1:(m)) seeds[[i]] <- sample.int(1000, 1)
  ## For the last model:
  seeds[[m+1]] <- sample.int(1000, 1)
  
  
  fitControl <- trainControl(## 10-fold CV
    method="repeatedcv", 
    number=kfolds, 
    repeats=reps,
    # savePredictions = TRUE,
    savePredictions	='final',
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    seeds = seeds,
    index=cvfolds
    ,
    allowParallel = TRUE,
    # returnResamp = 'none',
    returnData=FALSE
    # ,
    # workers = 8
    # # # ## new option here:
    # # # sampling = "down"
  )
  
  
  # Design formula
  form <-as.formula(paste(responseVar,' ~ ', design))
  
  
  ##############################
  ## Reproducibility point B ###
  ##############################
  ## Logistic regression
  
  
  
  if (isTRUE(inoutSub=="in")){
    print(paste0("Selected to subsample inside CV, with ",subMethod,"!"))
    set.seed(NULL)
    set.seed(seedModel)
    fitControl$sampling <- subMethod
  }else if (isTRUE(inoutSub=="out")){
    print(paste0("Selected to subsample outside CV, with ",subMethod,"!"))
    # DFtrain <- as.data.frame(unclass(DFtrain))
    
    # set the seed to make sampling reproducible
    set.seed(NULL)
    set.seed(seedPart)
    
    if(subMethod=="smote"){
      # data_balanced <- SMOTE(form, data = DFtrain, perc.over = 1100, perc.under = 110, k = 5)
      # data_balanced <- SMOTE(form, data = DFtrain, perc.over = 900, perc.under = 100, k = 10)
      data_balanced <- SMOTE(response~., data = DFtrain, perc.over = 600, perc.under = 110, k = 5)
    }else if(subMethod=="down"){
      data_balanced <- downSample(x = as.matrix( DFtrain[, -which(colnames(DFtrain)==responseVar)]),
                                  y = DFtrain[[responseVar]])
      colnames(data_balanced)[which(colnames(data_balanced)=="Class")]<-responseVar
      
    }else if(subMethod=="up"){
      data_balanced <- upSample(x = as.matrix( DFtrain[, -which(colnames(DFtrain)==responseVar)]),
                                y = DFtrain[[responseVar]])
      colnames(data_balanced)[which(colnames(data_balanced)=="Class")]<-responseVar
    }else{
      data_balanced <-DFtrain
    }
    
    # shuffle data by rows
    set.seed(NULL)
    set.seed(seedPart)
    data_balanced <- data_balanced[sample(1:nrow(data_balanced)), ]
    data_balanced <- as.data.frame(data_balanced)
    levels(data_balanced[[responseVar]])
    # data_balanced <- as.data.table(data_balanced)
    print("Synthetic samples:")
    print(table(data_balanced$response))
    DFtrain<-data_balanced
  }
  
  # SHUFFLE DATA BY ROWS
  # shuffle data by rows
  set.seed(NULL)
  set.seed(seedModel)
  DFtrain <- DFtrain[sample(1:nrow(DFtrain)), ]
  
  set.seed(NULL)
  set.seed(seedModel)
  
  cl <- makePSOCKcluster(detectCores(all.tests = FALSE, logical = TRUE)-2)
  # cl <- makePSOCKcluster(5)
  print(paste0("Using cores:",cl))
  registerDoParallel(cl)
  
  if (isTRUE(positiveCoefs)){
    print("~Force model to get only positive coefficients!~")
    if(isTRUE(applyScaling)){
      print("~Scale data!~")
      modelFit.smote <-train(form, data = DFtrain, 
                             method="glm", 
                             family=binomial(), 
                             trControl = fitControl,
                             metric = metricSelected, # Optimize by F-measure
                             maximize = TRUE,
                             preProcess=c('knnImpute', 'center', 'scale'),
                             lower.limits=0)#
    }else{
      modelFit.smote <-train(form, data = DFtrain, 
                             method="glm", 
                             family=binomial(), 
                             trControl = fitControl,
                             metric = metricSelected, # Optimize by F-measure
                             maximize = TRUE,
                             lower.limits=0)#
    }
    
  }else{
    print("~DO NOT force model to get only positive coefficients!~")
    if(isTRUE(applyScaling)){
      print("~Scale data!~")
      modelFit.smote <-train(form, data = DFtrain, 
                             method="glm", 
                             family=binomial(), 
                             trControl = fitControl,
                             metric = metricSelected, # Optimize by F-measure
                             maximize = TRUE,
                             preProcess=c('knnImpute', 'center', 'scale'))#
    }else{
      modelFit.smote <-train(form, data = DFtrain, 
                             method="glm", 
                             family=binomial(), 
                             trControl = fitControl,
                             metric = metricSelected, # Optimize by F-measure
                             maximize = TRUE)#
    }
    
  }
  
  
  stopCluster(cl)
  registerDoSEQ()
  
  m <-list(
    # original = modelFit.imb,
    #    down = modelFit.down,
    #    up = modelFit.up,
    modelFit.smote)
  names(m)<- ModelName
  m
}



multiTRAINModels_RFE.repCV_construct <- function(DFtrain,responseVar,design,feats,
                                                 kfolds,reps,sizeVec,
                                                 subMethod=NULL,inoutSub="in",
                                                 metricSelected,
                                                 seedModel=NA,seedPart,seeds=NULL,
                                                 stratify=FALSE,positiveCoefs=FALSE,
                                                 applyScaling=FALSE,wrapperFuncs = "RF"){
  
  
  # DFtrain <-data_train2
  # responseVar <-responseVar
  # design <-modelDesign
  # kfolds <-KFOLDS
  # reps <-REPEATS
  # sizeVec <- c(1:7)
  # subMethod=subsampMethod
  # inoutSub=inoutSub
  # metricSelected <-"ROC"
  # seedModel <-seed.model
  # seedPart <-seed.partition
  # seeds=setCVseeds
  # stratify=stratify
  # applyScaling=selectedScale
  # positiveCoefs = FALSE
  # wrapperFuncs = "RF"
  # 
  
  ################################
  ## CHECK 4 : MODEL PARAMETERS ##
  ################################
  
  print("==MODEL PARAMETERS==")
  print("=====================")
  print(paste0("Folds: ",kfolds))
  print("-----------------------")
  print(paste0("Repeats: ",reps))
  print("-----------------------")
  print(paste0("Model seed: ",seedModel))
  print("-----------------------")
  print(paste0("Set CV seeds: ",seeds))
  print("-----------------------")
  print(paste0("Subsampling: ",inoutSub))
  print("-----------------------")
  print(paste0("Subsampling Method: ",subMethod))
  print("-----------------------")
  print(paste0("Model metric:",metricSelected))
  print("-----------------------")
  
  print("==========================")
  ##############################
  ## Reproducibility point A ###
  ##############################
  
  # When the models are created inside of resampling, the seeds can also be set. 
  # While setting the seed prior to calling train may guarantee that the same random numbers are used, 
  # this is unlikely to be the case when parallel processing is used (depending which technology is utilized). 
  #To set the model fitting seeds, trainControl has an additional argument called seeds that can be used. 
  # The value for this argument is a list of integer vectors that are used as seeds. 
  if(isTRUE(stratify)){
    print("~Create Stratified folds for CV~")
    set.seed(NULL)
    set.seed(seedModel)
    cvfolds= createMultiFolds(DFtrain[[responseVar]],k=kfolds, times = reps)
  }else{
    cvfolds=NULL
  }
  
  if (! is.null(seeds)){
    m<-reps*kfolds
    n.eval=length(lamdaSeq)*length(alpha)
    seeds <- vector(mode = "list", length = m+1)
    set.seed(NULL)
    set.seed(seedModel)
    for(i in 1:(m)) seeds[[i]] <- sample.int(1000, n.eval)
    ## For the last model:
    seeds[[m+1]] <- sample.int(1000, 1)
  }
  
  
  # Train control object for how the data are trained, and if scaled
  # if (isTRUE(applyScaling)){
  #   set.seed(seedModel)
  #   
  #   trainctrl <- trainControl(#number=kfolds, 
  #                             #repeats=reps,
  #                             # savePredictions = TRUE,
  #                             savePredictions	='final',
  #                             classProbs = TRUE,
  #                             summaryFunction = twoClassSummary,
  #                             seeds = seeds,
  #                             index=cvfolds,
  #                             preProc = c("center", "scale"),
  #                             allowParallel = TRUE
  #                             # returnResamp = 'none',
  #                             # returnData=FALSE
  #                             )
  #   
  # }else{
  #   set.seed(seedModel)
  #   
  #   trainctrl <- trainControl(#number=kfolds, 
  #                             #repeats=reps,
  #                             # savePredictions = TRUE,
  #                             savePredictions	='final',
  #                             classProbs = TRUE,
  #                             summaryFunction = twoClassSummary,
  #                             seeds = seeds,
  #                             index=cvfolds,
  #                             # preProc = c("center", "scale"),
  #                             allowParallel = TRUE
  #                             # returnResamp = 'none',
  #                             # returnData=FALSE
  #                             )
  #   
  # }
  ## HOW THE FINAL MODEL AFTER RFE AND FEATURE SELECTION WILL BE TRAINED
  
  
  
  trainctrl <- trainControl(
    # method = "repeatedcv",
    method = "cv",
    number=kfolds, 
    # repeats=reps,
    # savePredictions = TRUE,
    savePredictions	='final',
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    # seeds = seeds,
    # index=cvfolds,
    # preProc = c("center", "scale"),# doesnt apply when using rfe
    allowParallel = TRUE
    # returnResamp = 'none',
    # returnData=FALSE
  )
  
  # Design formula
  form <-as.formula(paste(responseVar,' ~ ', design))
  
  
  ##############################
  ## Reproducibility point B ###
  ##############################
  
  
  
  if (isTRUE(inoutSub=="in")){
    print(paste0("Selected to subsample inside CV, with ",subMethod,"!"))
    set.seed(NULL)
    set.seed(seedModel)
    trainctrl$sampling <- subMethod
  }else if (isTRUE(inoutSub=="out")){
    print(paste0("Selected to subsample outside CV, with ",subMethod,"!"))
    # DFtrain <- as.data.frame(unclass(DFtrain))
    
    # set the seed to make sampling reproducible
    set.seed(NULL)
    set.seed(seedPart)
    
    if(subMethod=="smote"){
      # data_balanced <- SMOTE(form, data = DFtrain, perc.over = 1100, perc.under = 110, k = 5)
      # data_balanced <- SMOTE(form, data = DFtrain, perc.over = 900, perc.under = 100, k = 10)
      data_balanced <- SMOTE(response~., data = DFtrain, perc.over = 600, perc.under = 110, k = 5)
    }else if(subMethod=="down"){
      data_balanced <- downSample(x = as.matrix( DFtrain[, -which(colnames(DFtrain)==responseVar)]),
                                  y = DFtrain[[responseVar]])
      colnames(data_balanced)[which(colnames(data_balanced)=="Class")]<-responseVar
      
    }else if(subMethod=="up"){
      data_balanced <- upSample(x = as.matrix( DFtrain[, -which(colnames(DFtrain)==responseVar)]),
                                y = DFtrain[[responseVar]])
      colnames(data_balanced)[which(colnames(data_balanced)=="Class")]<-responseVar
    }else{
      data_balanced <-DFtrain
    }
    
    # shuffle data by rows
    set.seed(NULL)
    set.seed(seedPart)
    data_balanced <- data_balanced[sample(1:nrow(data_balanced)), ]
    data_balanced <- as.data.frame(data_balanced)
    levels(data_balanced[[responseVar]])
    # data_balanced <- as.data.table(data_balanced)
    print("Synthetic samples:")
    print(table(data_balanced$response))
    DFtrain<-data_balanced
  }
  
  # SHUFFLE DATA BY ROWS
  # shuffle data by rows
  set.seed(NULL)
  set.seed(seedModel)
  DFtrain <- DFtrain[sample(1:nrow(DFtrain)), ]
  
  
  
  
  
  #############################
  #############################
  ## Recursive feature elimination
  
  set.seed(NULL)
  set.seed(seedModel)
  
  cl <- makePSOCKcluster(detectCores(all.tests = FALSE, logical = TRUE)-2)
  # cl <- makePSOCKcluster(5)
  print(paste0("Using cores:",cl))
  registerDoParallel(cl)
  
  if(wrapperFuncs =="RF"){
    print(">> Fitting Random Forests models")
    ## Random forest wrapper
    
    # rfFuncs$summary <- twoClassSummary
    caretFuncs$summary <- twoClassSummary
    
    
    # RUN RFE with variable selecion and final model
    rfGrid <- expand.grid(.mtry = c(round(sqrt(length(feats))/2,0),round(sqrt(length(feats)),0),round(sqrt(length(feats))*2,0)))
    ## Calculate combinations number of tuning parameters.
    gridCombos <- length(c(round(sqrt(length(feats))/2,0))) * length(round(sqrt(length(feats)),0))* length(round(sqrt(length(feats))*2,0))
    ## DO THE FOLLOWING FOR PARALLEL OPERATION-CONTROL RESAMPLES
    set.seed(NULL)
    set.seed(seedModel)
    m <- kfolds*reps
    seeds <- vector(mode = "list", length = m+1)
    # print(paste0("SEEDS LIST:", length(seeds)))
    # print(paste0("SEEDS ELEMENTS:", length(sizeVec)))
    # 
    for(i in 1:(m)) seeds[[i]] <- sample.int(1000, length(sizeVec))
    ## For the last model:
    seeds[[m+1]] <- sample.int(1000, 1)
    # print(seeds)
    # RFE CONTROL
    rctrl <- rfeControl(method = "repeatedcv", # repeated cv
                        repeats = reps, # number of repeats
                        number = kfolds, # number of folds,
                        #  method = "LGOCV",
                        # number = 30,
                        returnResamp = "all",
                        # functions = rfFuncs,# rfFuncs,caretFuncs
                        functions = caretFuncs,# rfFuncs,caretFuncs
                        saveDetails = TRUE,
                        allowParallel= TRUE,
                        seeds = seeds)
    
    if (isTRUE(applyScaling)){
      print(">>>>>>SCALE DATA<<<<<<")
      set.seed(NULL)
      set.seed(seedModel)
      model <- rfe(form, data = DFtrain,
                   sizes = sizeVec,
                   method = "rf",
                   trControl = trainctrl,
                   tuneGrid = rfGrid,
                   metric=metricSelected,
                   maximize = TRUE,
                   preProc = c("center", "scale"),
                   rfeControl = rctrl,
                   returnResamp = "all")
    }else{
      print(">>>>>>DONT SCALE DATA<<<<<<")
      set.seed(NULL)
      set.seed(seedModel)
      model <- rfe(form, data = DFtrain,
                   sizes = sizeVec,
                   method = "rf",
                   trControl = trainctrl,
                   tuneGrid = rfGrid,
                   metric=metricSelected,
                   maximize = TRUE,
                   # preProc = c("center", "scale"),
                   rfeControl = rctrl)
    }
    
    
  }else if (wrapperFuncs =="LR"){
    print(">> Fitting Linear Regression models")
    ## Linear regression wrapper
    # lrFuncs$summary <- twoClassSummary
    caretFuncs$summary <- twoClassSummary
    
    ## Calculate combinations number of tuning parameters.
    gridCombos <- 1
    ## DO THE FOLLOWING FOR PARALLEL OPERATION-CONTROL RESAMPLES
    set.seed(NULL)
    set.seed(seedModel)
    m <- kfolds*reps
    seeds <- vector(mode = "list", length = m+1)
    # print(paste0("SEEDS LIST:", length(seeds)))
    # print(paste0("SEEDS ELEMENTS:", length(sizeVec)))
    # 
    for(i in 1:(m)) seeds[[i]] <- sample.int(1000, length(sizeVec))
    ## For the last model:
    seeds[[m+1]] <- sample.int(1000, 1)
    
    # RFE CONTROL
    rctrl <- rfeControl(method = "repeatedcv", # repeated cv
                        repeats = reps, # number of repeats
                        number = kfolds, # number of folds,
                        #  method = "LGOCV",
                        # number = 30,
                        returnResamp = "all",
                        # functions = lrFuncs,# rfFuncs,caretFuncs
                        functions = caretFuncs,
                        saveDetails = TRUE,
                        allowParallel= TRUE,
                        seeds = seeds)
    
    if (isTRUE(applyScaling)){
      print(">>>>>>SCALE DATA<<<<<<")
      set.seed(NULL)
      set.seed(seedModel)
      model <- rfe(form, data = DFtrain,
                   sizes = sizeVec,
                   method = 'glm',
                   family=binomial(),
                   trControl = trainctrl,
                   # tuneGrid = data.frame(k = 1:10),
                   metric=metricSelected,
                   maximize = TRUE,
                   preProc = c("center", "scale"),
                   rfeControl = rctrl)
    }else{
      set.seed(NULL)
      set.seed(seedModel)
      model <- rfe(form, data = DFtrain,
                   sizes = sizeVec,
                   method = 'glm',
                   family=binomial(),
                   trControl = trainctrl,
                   # tuneGrid = data.frame(k = 1:10),
                   
                   metric=metricSelected,
                   maximize = TRUE,
                   rfeControl = rctrl)
    }
    
    
  }else if (wrapperFuncs =="NB"){
    print(">> Fitting Naive Bayes models")
    ## Linear regression wrapper
    # nbFuncs$summary <- twoClassSummary
    caretFuncs$summary <- twoClassSummary
    
    # nb_grid <- expand.grid(fL=c(0,0.5,1.0), usekernel = c(TRUE, FALSE), adjust=seq(0,5,0.5))# Laplace Correction (fL, numeric), Distribution Type (usekernel, logical),Bandwidth Adjustment (adjust, numeric)
    
    ## Calculate combinations number of tuning parameters.
    gridCombos <- 1
    ## DO THE FOLLOWING FOR PARALLEL OPERATION-CONTROL RESAMPLES
    set.seed(NULL)
    set.seed(seedModel)
    m <- kfolds*reps
    seeds <- vector(mode = "list", length = m+1)
    # print(paste0("SEEDS LIST:", length(seeds)))
    # print(paste0("SEEDS ELEMENTS:", length(sizeVec)))
    # 
    for(i in 1:(m)) seeds[[i]] <- sample.int(1000, length(sizeVec))
    ## For the last model:
    seeds[[m+1]] <- sample.int(1000, 1)
    
    # RFE CONTROL
    rctrl <- rfeControl(method = "repeatedcv", # repeated cv
                        repeats = reps, # number of repeats
                        number = kfolds, # number of folds,
                        #  method = "LGOCV",
                        # number = 30,
                        returnResamp = "all",
                        # functions = nbFuncs,# rfFuncs,caretFuncs
                        functions = caretFuncs,
                        saveDetails = TRUE,
                        allowParallel= TRUE,
                        seeds = seeds)
    # Grid makes it even slower
    if (isTRUE(applyScaling)){
      print(">>>>>>SCALE DATA<<<<<<")
      set.seed(NULL)
      set.seed(seedModel)
      model <- rfe(form, data = DFtrain,
                   sizes = sizeVec,
                   method = "nb",
                   trControl = trainctrl,
                   # tuneGrid = nb_grid,
                   metric=metricSelected,
                   maximize = TRUE,
                   preProc = c("center", "scale"),
                   rfeControl = rctrl)
    }else{
      set.seed(NULL)
      set.seed(seedModel)
      model <- rfe(form, data = DFtrain,
                   sizes = sizeVec,
                   method = "nb",
                   trControl = trainctrl,
                   tuneGrid = nb_grid,
                   metric=metricSelected,
                   maximize = TRUE,
                   rfeControl = rctrl)
    }
    
  }else if (wrapperFuncs =="SVMr"){
    print(">> Fitting SVM radial model")
    ## Linear regression wrapper
    # caretFuncs$summary <- twoClassSummary
    
    ## Calculate combinations number of tuning parameters.
    gridCombos <- 1
    ## DO THE FOLLOWING FOR PARALLEL OPERATION-CONTROL RESAMPLES
    set.seed(NULL)
    set.seed(seedModel)
    m <- kfolds*reps
    seeds <- vector(mode = "list", length = m+1)
    # print(paste0("SEEDS LIST:", length(seeds)))
    # print(paste0("SEEDS ELEMENTS:", length(sizeVec)))
    # 
    for(i in 1:(m)) seeds[[i]] <- sample.int(1000, length(sizeVec))
    ## For the last model:
    seeds[[m+1]] <- sample.int(1000, 1)
    
    # RFE CONTROL
    rctrl <- rfeControl(method = "repeatedcv", # repeated cv
                        repeats = reps, # number of repeats
                        number = kfolds, # number of folds,
                            returnResamp = "all",
                        functions = caretFuncs,# rfFuncs,caretFuncs
                        saveDetails = TRUE,
                        verbose=TRUE,
                        allowParallel= TRUE,
                        seeds = seeds)
    caretFuncs$summary <- twoClassSummary
    
    
    if (isTRUE(applyScaling)){
      print(">>>>>>SCALE DATA<<<<<<")
      set.seed(NULL)
      set.seed(seedModel)
      model <- rfe(form, data = DFtrain,
                   sizes = sizeVec,
                   method = "svmRadial",
                   trControl = trainctrl,
                   # tuneGrid = data.frame(k = 1:10),
                   metric=metricSelected,
                   maximize = TRUE,
                   preProc = c("center", "scale"),
                   rfeControl = rctrl)
    }else{
      set.seed(NULL)
      set.seed(seedModel)
      model <- rfe(form, data = DFtrain,
                   sizes = sizeVec,
                   method = "svmRadial",
                   trControl = trainctrl,
                   # tuneGrid = data.frame(k = 1:10),
                   metric=metricSelected,
                   maximize = TRUE,
                   rfeControl = rctrl)
    }
    
    
    
  }else if (wrapperFuncs =="SVMl"){
    print(">> Fitting SVM Linear model")
    ## Linear regression wrapper
    caretFuncs$summary <- twoClassSummary
    
    
    
    ## Calculate combinations number of tuning parameters.
    gridCombos <- 1
    ## DO THE FOLLOWING FOR PARALLEL OPERATION-CONTROL RESAMPLES
    set.seed(NULL)
    set.seed(seedModel)
    m <- kfolds*reps
    seeds <- vector(mode = "list", length = m+1)
    # print(paste0("SEEDS LIST:", length(seeds)))
    # print(paste0("SEEDS ELEMENTS:", length(sizeVec)))
    # 
    for(i in 1:(m)) seeds[[i]] <- sample.int(1000, length(sizeVec))
    ## For the last model:
    seeds[[m+1]] <- sample.int(1000, 1)
    
    # RFE CONTROL
    rctrl <- rfeControl(method = "repeatedcv", # repeated cv
                        repeats = reps, # number of repeats
                        number = kfolds, # number of folds,
                        #  method = "LGOCV",
                        # number = 30,
                        returnResamp = "all",
                        functions = caretFuncs,# rfFuncs,caretFuncs
                        saveDetails = TRUE,
                        allowParallel= TRUE,
                        seeds = seeds)
    if (isTRUE(applyScaling)){
      print(">>>>>>SCALE DATA<<<<<<")
      set.seed(NULL)
      set.seed(seedModel)
      model <- rfe(form, data = DFtrain,
                   sizes = sizeVec,
                   method = "svmLinear",
                   trControl = trainctrl,
                   # tuneGrid = data.frame(k = 1:10),
                   metric=metricSelected,
                   maximize = TRUE,
                   preProc = c("center", "scale"),
                   rfeControl = rctrl)
    }else{
      set.seed(NULL)
      set.seed(seedModel)
      model <- rfe(form, data = DFtrain,
                   sizes = sizeVec,
                   method = "svmLinear",
                   trControl = trainctrl,
                   # tuneGrid = data.frame(k = 1:10),
                   metric=metricSelected,
                   maximize = TRUE,
                   rfeControl = rctrl)
    }
    
    
  }else if (wrapperFuncs =="KNN"){
    print(">> Fitting KNN method")
    ## Linear regression wrapper
    caretFuncs$summary <- twoClassSummary
    
    ## Calculate combinations number of tuning parameters.
    gridCombos <- 1
    ## DO THE FOLLOWING FOR PARALLEL OPERATION-CONTROL RESAMPLES
    set.seed(NULL)
    set.seed(seedModel)
    m <- kfolds*reps
    seeds <- vector(mode = "list", length = m+1)
    # print(paste0("SEEDS LIST:", length(seeds)))
    # print(paste0("SEEDS ELEMENTS:", length(sizeVec)))
    # 
    for(i in 1:(m)) seeds[[i]] <- sample.int(1000, length(sizeVec))
    ## For the last model:
    seeds[[m+1]] <- sample.int(1000, 1)
    
    # RFE CONTROL
    rctrl <- rfeControl(method = "repeatedcv", # repeated cv
                        repeats = reps, # number of repeats
                        number = kfolds, # number of folds,
                        #  method = "LGOCV",
                        # number = 30,
                        returnResamp = "all",
                        functions = caretFuncs,# rfFuncs,caretFuncs
                        saveDetails = TRUE,
                        seeds = seeds)
    
    if (isTRUE(applyScaling)){
      print(">>>>>>SCALE DATA<<<<<<")
      set.seed(NULL)
      set.seed(seedModel)
      model <- rfe(form, data = DFtrain,
                   sizes = sizeVec,
                   method = "knn",
                   trControl = trainctrl,
                   metric=metricSelected,
                   maximize = TRUE,
                   tuneGrid = data.frame(k = 1:10),
                   preProc = c("center", "scale"),
                   rfeControl = rctrl)
    }else{
      print(">>>>>>NOT SCALE DATA<<<<<<")
      set.seed(NULL)
      set.seed(seedModel)
      model <- rfe(form, data = DFtrain,
                   sizes = sizeVec,
                   method = "knn",
                   trControl = trainctrl,
                   metric=metricSelected,
                   maximize = TRUE,
                   tuneGrid = data.frame(k = 1:10),
                   rfeControl = rctrl)
    }
    
  }else if (wrapperFuncs =="NNET"){
    print(">> Fitting Neural Network method")
    ## Linear regression wrapper
    caretFuncs$summary <- twoClassSummary
    
    
    nnetGrid <-  expand.grid(size = seq(from = 1, to = 10, by = 1),
                             decay = seq(from = 0.1, to = 0.5, by = 0.1))
    
    ## Calculate combinations number of tuning parameters.
    gridCombos <- length(seq(from = 1, to = 10, by = 1))* length(seq(from = 0.1, to = 0.5, by = 0.1))
    ## DO THE FOLLOWING FOR PARALLEL OPERATION-CONTROL RESAMPLES
    set.seed(NULL)
    set.seed(seedModel)
    m <- kfolds*reps
    seeds <- vector(mode = "list", length = m+1)
    # print(paste0("SEEDS LIST:", length(seeds)))
    # print(paste0("SEEDS ELEMENTS:", length(sizeVec)))
    # 
    for(i in 1:(m)) seeds[[i]] <- sample.int(1000, length(sizeVec))
    ## For the last model:
    seeds[[m+1]] <- sample.int(1000, 1)
    # RFE CONTROL
    rctrl <- rfeControl(method = "repeatedcv", # repeated cv
                        repeats = reps, # number of repeats
                        number = kfolds, # number of folds,
                        #  method = "LGOCV",
                        # number = 30,
                        returnResamp = "all",
                        functions = caretFuncs,# rfFuncs,caretFuncs
                        saveDetails = TRUE,
                        seeds = seeds)
    
    if (isTRUE(applyScaling)){
      set.seed(NULL)
      set.seed(seedModel)
      model <- rfe(form, data = DFtrain,
                   sizes = sizeVec,
                   method = 'nnet',
                   trControl = trainctrl,
                   metric=metricSelected,
                   maximize = TRUE,
                   tuneGrid = nnetGrid,
                   preProc = c("center", "scale"),
                   rfeControl = rctrl)
    }else{
      set.seed(NULL)
      set.seed(seedModel)
      model <- rfe(form, data = DFtrain,
                   sizes = sizeVec,
                   method = 'nnet',
                   trControl = trainctrl,
                   metric=metricSelected,
                   maximize = TRUE,
                   tuneGrid = nnetGrid,
                   rfeControl = rctrl)
    }
    
  }
  
  stopCluster(cl)
  registerDoSEQ()
  model
  
  
}



extract_listObj <- function(List, ObjIndex){
  List[[ObjIndex]]
}

`%notin%` <- Negate(`%in%`)

logit2prob <- function(logit){## GIVES THE PROBABILITY OF THE NEGATIVE CLASS for some reason
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}
## Assign probabilistic score to any new sample
extractProb_test <- function(DFtest,objectName){#, columnName,
  # DFtest[[objectName]] <-prob_score(DFtest[[columnName]])
  # DFtest
  d <-logit2prob(DFtest)
  k <-as.data.frame(d)
  colnames(k)<-objectName
  k
}
getProbDF <- function(DFtest,sigNames,modelNames,classes,LevelIndex){
  # DFtest<-data.testScores
  # sigNames<-signatures[1:2]
  # modelNames<- model.names[1:2]
  # classes<-ClassLevels
  # LevelIndex<-1
  pos.un <- as.name(classes[LevelIndex])
  neg.un <- as.name(classes[-LevelIndex])
  
  test_prob.DF <-lapply(1:length(sigNames), function(x) extractProb_test(DFtest[,sigNames[x]],modelNames[x]))
  
  test_prob.DF.merged <-bind_cols(test_prob.DF,.id = NULL)
  test_prob.DF.merged.all <- cbind(DFtest %>% dplyr::select(-sigNames),test_prob.DF.merged)
  
  test_prob.DF.merged.all.probs <- test_prob.DF.merged.all %>% dplyr::select(response,modelNames)
  colnames(test_prob.DF.merged.all.probs)[1] <- "obs"
  
  test_prob.DF.merged.all.probs.m <- melt(test_prob.DF.merged.all.probs)
  
  colnames(test_prob.DF.merged.all.probs.m)[3]<- classes[LevelIndex]
  
  test_prob.DF.merged.all.probs.m[[classes[-LevelIndex]]] <- 1-test_prob.DF.merged.all.probs.m[[classes[LevelIndex]]]
  
  test_prob.DF.merged.all.probs.m$model <- "glm"
  test_prob.DF.merged.all.probs.m$dataType <- "Unknown"
  colnames(test_prob.DF.merged.all.probs.m)[2]<- "object"
  
  test_prob.DF.merged.all.probs.m$pred <- ifelse(test_prob.DF.merged.all.probs.m[[classes[LevelIndex]]] >= 0.5,classes[LevelIndex],classes[-LevelIndex])
  
  test_prob.DF.merged.all.probs.m <- test_prob.DF.merged.all.probs.m %>% dplyr::select(pos.un, neg.un,obs,pred, model,dataType,object)
  test_prob.DF.merged.all.probs.m$pred<- factor(test_prob.DF.merged.all.probs.m$pred,levels = classes) 
  test_prob.DF.merged.all.probs.m$pred<- relevel(test_prob.DF.merged.all.probs.m$pred, classes[-LevelIndex])
  test_prob.DF.merged.all.probs.m$object<- as.factor(test_prob.DF.merged.all.probs.m$object)
  prob_test.f<-test_prob.DF.merged.all.probs.m
  prob_test.f<- prob_test.f[,c(classes,"obs", "pred", "model", "dataType", "object")]
  prob_test.f
}

trainRUN <- function(DF,split,ClassLevels,PositiveClass,LevelIndex,seed.partition,seed.model,setCVseeds=NULL,
                                     dummy=FALSE,dummyVars=NULL,categVars=NULL,obsColumn,responseVar,REPEATS,KFOLDS,
                                     metricOpt,subsampMethod=NULL,inoutSub="in",
                                     modelDesign,controlDum=FALSE,model.name,mainFeats,
                                     removeOutliers=TRUE,quantileClass1=.95,quantileClass2=.95,selectOutlierMethod= "meanDist",
                                     stratify=FALSE,fullResults=TRUE,
                                     positiveCoefs=FALSE, filter=FALSE,filterVar=NULL,filterSelected=NULL,CorrRemove=FALSE,corrThres=NULL,selectedScale=FALSE){
  ####################################
  # DF <-immdata.apd1.proc.TcrA.divEst.full.sub.ALL.filt 
  # split<-0.66
  # ClassLevels<-c("CRPR","PD")
  # PositiveClass<-"CRPR"
  # LevelIndex<-1
  # seed.partition=3456
  # seed.model=5627
  # setCVseeds=NULL
  # dummy=FALSE
  # dummyVars=c("tissue")
  # categVars=c("tissue","dataset")
  # obsColumn<-"run_accession"
  # responseVar<-"response"
  # REPEATS=30
  # KFOLDS=10
  # metricOpt="F"
  # subsampMethod="smote"
  # inoutSub="in"
  # modelDesign=paste("logTMB",collapse = " + ")
  # controlDum =FALSE
  # model.name<-"logTMB"
  # # subset.notIn.sig<-degs.only
  # mainFeats<-c("logTMB")
  # removeOutliers=FALSE
  # quantileClass1=.99
  # quantileClass2=.99
  # selectOutlierMethod= "none"
  # stratify=FALSE
  # fullResults=TRUE
  # positiveCoefs=FALSE
  # filter=TRUE
  # filterVar="dataset"
  # filterSelected=c("Powles","Gide","Mari","Rod","Hugo","Riaz")
  # CorrRemove= FALSE
  # corrThres=0.85
  # selectedScale<-TRUE
  #####################################
  #################################################################
  #################################################################
  #################################
  ##################################
  ################################
  ## CHECK 0 : SAMPLES, TISSUES ##
  ################################
  # DF <- DF %>% dplyr::filter(complete.cases(!!as.name(mainFeats[1]))) %>% droplevels()
  DF <-DF %>% dplyr::filter_at(vars(all_of(mainFeats)), all_vars(complete.cases(.))) %>% droplevels()
  print("==INITIAL INPUT DATA==")
  print("=======================")
  print(table(DF[[responseVar]]))
  print(table(DF[[dummyVars]]))
  print("=======================")
  
  ##########################
  ## CHECK 1 : PARAMETERS ##
  ##########################
  
  print("==MODEL/PIPELINE PARAMETERS==")
  print("==============================")
  print(paste0("Outlier filtering: ",removeOutliers))
  if(isTRUE(removeOutliers)){
    print(paste0("Outlier Method: ",selectOutlierMethod))
  }
  
  print(paste0("Dataset Filtering: ",filter))
  if(isTRUE(filter)){
    print(paste0("Datasets selected: ",filterSelected))
  }
  
  print(paste0("Higly Correlated Features Filtering: ",CorrRemove))
  if(isTRUE(CorrRemove)){
    print(paste0("Higly Correlated Features Filtering Threshold: ",corrThres))
  }
  
  if(isTRUE(selectedScale)){
    print(paste0("Scaling of train data, scaling of test data with train data scale factor: ",selectedScale))
  }
  print("===============================")
  
  
  ################################
  ## CHECK 2 : MODEL PARAMETERS ##
  ################################
  
  print("==MODEL PARAMETERS==")
  
  print("=====================")
  print(paste0("Split: ",  split))
  print("-----------------------")
  print(paste0("Folds: ",KFOLDS))
  print("-----------------------")
  print(paste0("Repeats: ",REPEATS))
  print("-----------------------")
  print(paste0("Model seed: ",seed.model))
  print("-----------------------")
  print(paste0("Set CV seeds: ",setCVseeds))
  print("-----------------------")
  print(paste0("Subsampling: ",inoutSub))
  print("-----------------------")
  print(paste0("Subsampling Method: ",subsampMethod))
  print("-----------------------")
  print(paste0("Model metric:",metricOpt))
  print("-----------------------")
  print(paste0("Positive Class: ",PositiveClass))
  print("==========================")
  
  ###########################################
  ## Detect outliers with multiple methods ##
  ###########################################
  # Make sure response var for split is not a factor#
  DF[[responseVar]] <- as.factor(DF[[responseVar]])
  # Split to responders and non-responders
  data_full.crpr <- DF %>% dplyr::filter(response==ClassLevels[LevelIndex])
  data_full.pd <- DF %>% dplyr::filter(response==ClassLevels[-LevelIndex])
  
  ## Run automatic methods from package OutlierDetection
  method.names<- c("Robust Kernal-based Outlier Factor(RKOF) algorithm","Depth based method", "Generalised dispersion method","Mahalanobis Distance","k Nearest Neighbours Distance method", "kth Nearest Neighbour Distance method", "Outlier Detection(Intersection of all the methods)","Principal Component Outlier Detection(Intersection of all the methods applied on pc’s)" )
  
  # out.crpr <- outlierDetect_multiple(as.data.frame(data_full.crpr[,mainFeats]), k=sqrt(dim(DF)[1]),C=1,cutoff.dens=.95, cutoff.depth=.05, cutoff.disp=.95,cutoff.m=.99,cutoff.u=.95)
  # 
  # out.pd <- outlierDetect_multiple(data_full.pd[,mainFeats], k=sqrt(dim(DF)[1]),C=1,cutoff.dens=.95, cutoff.depth=.05, cutoff.disp=.99,cutoff.m=.99,cutoff.u=.95)
  
  ## Mean Euclidean Distance-Miklos
  
  # pcDT <-pca_ObjCoords.TwoClass(DF,responseVar,categVars,mainFeats,colnames(DF)[dim(DF)[2]])
  # 
  # trainableData <- pcDT[,.(patientID, response, Dim.1, Dim.2)]
  # trainableData[response == ClassLevels[LevelIndex], respNum := 0]
  # trainableData[response == ClassLevels[-LevelIndex], respNum := 1]
  # 
  # dt_crpr <- trainableData[response == ClassLevels[LevelIndex]]
  # dt_pd <- trainableData[response == ClassLevels[-LevelIndex]]
  # 
  # 
  # k <- 10
  # 
  # # Caucasians:
  # for(i in 1:nrow(dt_crpr)){
  #   temp <- dt_crpr[-i,]
  #   temp[, Dim.1 := Dim.1 - dt_crpr[i,Dim.1]]
  #   temp[, Dim.2 := Dim.2 - dt_crpr[i,Dim.2]]
  #   temp[, dist := sqrt(Dim.1^2 + Dim.2^2)]
  #   temp <- temp[order(dist)]
  #   dt_crpr[i, dist:= mean(temp[1:k, dist])]
  # }
  # 
  # # AA:
  # for(i in 1:nrow(dt_pd)){
  #   temp <- dt_pd[-i,]
  #   temp[, Dim.1 := Dim.1 - dt_pd[i,Dim.1]]
  #   temp[, Dim.2 := Dim.2 - dt_pd[i,Dim.2]]
  #   temp[, dist := sqrt(Dim.1^2 + Dim.2^2)]
  #   temp <- temp[order(dist)]
  #   dt_pd[i, dist:= mean(temp[1:k, dist])]
  # }
  # 
  # 
  # trainableData[, isout := 0]
  # 
  # dt_crpr <- dt_crpr[order(-dist)]
  # dt_crpr[, isout := 0] # initially
  # 
  # 
  # dt_pd <- dt_pd[order(-dist)]
  # dt_pd[, isout := 0] # initially
  # 
  # percentileClass1 <- as.numeric(quantile(dt_crpr$dist[order(dt_crpr$dist)],c(quantileClass1))[1])
  # percentileClass2 <- as.numeric(quantile(dt_pd$dist[order(dt_pd$dist)],c(quantileClass2))[1])
  # 
  # dt_crpr[dist > percentileClass1, isout := 1]
  # trainableData[patientID %in% dt_crpr[isout==1,patientID],isout:=1]
  # 
  # dt_pd[dist > percentileClass2, isout := 1]
  # trainableData[patientID %in% dt_pd[isout==1,patientID],isout:=1]
  # 
  # ## Merge data to have the same order
  # colnames(data_full.crpr)[dim(data_full.crpr)[2]]<- "patientID"
  # data_full.crpr.m <- inner_join(data_full.crpr,dt_crpr[,.(patientID, respNum,isout)])
  # data_full.crpr.m<- as.data.table(data_full.crpr.m)
  # out.resp <- data_full.crpr.m[respNum==0 & isout==1,patientID]
  # out.resp.ind <- which(data_full.crpr.m$patientID %in% out.resp)
  # 
  # 
  # colnames(data_full.pd)[dim(data_full.pd)[2]]<- "patientID"
  # data_full.pd.m <- inner_join(data_full.pd,dt_pd[,.(patientID, respNum,isout)])
  # data_full.pd.m<- as.data.table(data_full.pd.m)
  # 
  # out.nonresp <- data_full.pd.m[respNum==1 & isout==1,patientID]
  # out.nonresp.ind <- which(data_full.pd.m$patientID %in% out.nonresp)
  # 
  # # outliers.all.ids<- c(out.resp,out.nonresp)
  # # ## Indexes in full train data
  # # outliers.all.index <- which(data_train[[subjectVar]] %in% outliers.all.ids)
  # # outlier.obs <- data_train[outliers.all.index,]
  # # Collect in list
  # out.meanDist.crpr <- list(`Outlier Observations`=data_full.crpr[out.resp.ind,],
  #                           `Location of Outlier` = out.resp.ind)
  # 
  # out.meanDist.pd <- list(`Outlier Observations`=data_full.pd[out.nonresp.ind,],
  #                         `Location of Outlier` = out.nonresp.ind)
  # 
  # out.crpr[[length(out.crpr)+1]]<- out.meanDist.crpr
  # names(out.crpr)[length(out.crpr)]<- "meanDist"
  # 
  # out.pd[[length(out.pd)+1]]<- out.meanDist.pd
  # names(out.pd)[length(out.pd)]<- "meanDist"
  # ## Extract patient IDs for full train data for selected method
  # patient.crpr <- data_full.crpr[out.crpr[[selectOutlierMethod]]$`Location of Outlier`,"patientID"]
  # patient.pd <- data_full.pd[out.pd[[selectOutlierMethod]]$`Location of Outlier`,"patientID"]
  # outliers.all.ids<- c(patient.crpr,patient.pd)
  # outliersMethod.crpr<-out.crpr[[selectOutlierMethod]]
  # outliersMethod.pd <-out.pd[[selectOutlierMethod]]
  # ## Indexes in full train data
  # outliers.all.index <- which(DF[[obsColumn]] %in% outliers.all.ids)
  # # outlier.obs <- data_train[outliers.all.index,]
  
  ####################################
  ## Clean full data from outliers ##
  ####################################
  if (isTRUE(removeOutliers)){
    print("~Removing outliers~")
    print(paste0("== Removing ",length(patient.crpr)," responders =="))
    print(paste0("== Removing ",length(patient.pd)," non-responders =="))
    print(paste0("~We will train/test the model with ", dim(DF)[1]-length(patient.crpr)-length(patient.pd),
                 " out of ",dim(DF)[1]," observations!~"))
    # data_train.full<-data_train
    
    DF <- DF[-outliers.all.index,]
  }else{
    print("~Outliers detected BUT NOT REMOVED~")
    print("~Outlier observations still reported in results!~")
    DF<-DF
  }
  
  # Make sure response var for split is not a factor#
  DF[[responseVar]] <- as.character(DF[[responseVar]])
  
  if(isTRUE(filter)){
    filterVar.un <- as.name(filterVar)
    ## Here I keep all the leftover from the filtering to be added as test data later, otherwise they are missed
    DF.left <- DF %>% dplyr::filter(!!filterVar.un %notin% filterSelected) %>% droplevels()
    DF <- DF %>% dplyr::filter(!!filterVar.un %in% filterSelected) %>% droplevels()
    print("~Filtering data~")
    print("== Including the following :==")
    print(paste(filterSelected, collapse = ", "))
    print(table(DF[[responseVar]]))
  }
  ## CHECK IF YOU NEED TO FILTER ##
  
  #######################################
  ## Remove highly correlated features ##
  #######################################
  if (isTRUE(CorrRemove)){
    print("~Identifying Correlated Features~")
    # calculate correlation matrix
    correlationMatrix <- cor(DF[,mainFeats], method = "spearman")
    #
    # find attributes that are highly corrected (ideally >0.75)
    print(paste0("~Removing Correlated Features, with thresold:",corrThres, "~"))
    highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=corrThres,names = TRUE, exact=TRUE)
    keep.cols<-setdiff(colnames(DF),highlyCorrelated)
    print("~Features removed:~")
    print(paste(highlyCorrelated,collapse = ", "))
    DF<-DF[,keep.cols]    # Re-define design of models
    if(isTRUE(filter)){
      DF.left <-DF.left[,keep.cols]
    }
    
    print("~Re-defining the model design after removing correlated features~")
    mainFeats.new <-setdiff(mainFeats,highlyCorrelated)
    modelDesign<-paste(mainFeats.new,collapse = " + ")
    mainFeats<- mainFeats.new
  }
  
  ################
  ## DATA SPLIT ##
  ################
  # Make sure response var for split is not a factor#
  DF[[responseVar]] <- as.character(DF[[responseVar]])
  
  print("== FINAL DATA AFTER FILTERING AND COMPLETE CASES FOR PREDICTORS :==")
  
  print(table(DF[[responseVar]]))
  print(table(DF[["tissue"]]))
  print(table(DF[["dataset"]]))
  print("===========================================================")
  ## Split data
  set.seed(NULL)
  set.seed(seed.partition)
  trainIndex <- createDataPartition(DF[[responseVar]], p=split, list=FALSE)
  data_train <- DF[ trainIndex,]
  data_test <- DF[-trainIndex,]
  # NOW SET FACTOR LEVELS
  data_train[[responseVar]] <- factor(data_train[[responseVar]], levels=ClassLevels)
  data_test[[responseVar]] <- factor(data_test[[responseVar]], levels=ClassLevels)
  if(isTRUE(filter)){
    DF.left[[responseVar]] <- factor(DF.left[[responseVar]], levels=ClassLevels)
  }
  
  print("~Using train data~")
  print(table(data_train[[responseVar]]))
  print("~Using test data~")
  print(table(data_test[[responseVar]]))
  
  ####################################
  ## Create dummy vars if necessary ##
  ####################################
  
  if (isTRUE(controlDum)){
    print("Creating dummy variables for categorical variable(s)")
    form.dummy <- as.formula(paste(responseVar,' ~ ', paste(dummyVars,collapse = " + ")))
    
    ## Train ##
    data_train[[responseVar]] <- as.character(data_train[[responseVar]])
    dummies.df <- dummyVars(form.dummy, data = data_train,fullRank=T,levelsOnly=TRUE)# Gide is used as reference
    
    dummies.df <- data.frame(predict(dummies.df, newdata = data_train))
    # Merge with original df
    data_train <- cbind(data_train,dummies.df)
    dummyCols <- colnames(dummies.df)
    
    ## Test ##
    if(length(levels(as.factor(data_test$tissue)))<=1){
      print("Cannot create dummy vars for test data here")
    }else{
      data_test[[responseVar]] <- as.character(data_test[[responseVar]])
      dummies.df <- dummyVars(form.dummy, data = data_test,fullRank=T,levelsOnly=TRUE)# Gide is used as reference
      
      dummies.df <- data.frame(predict(dummies.df, newdata = data_test))
      # Merge with original df
      data_test <- cbind(data_test,dummies.df)
      if(isTRUE(filter)){
        DF.left[[responseVar]] <- as.character(DF.left[[responseVar]])
        # levels(DF.left[[dummyVars]]) <- list(melanoma="melanoma",
        #                           bladder="bladder",
        #                           rcc="rcc",
        #                           gastric="gastric")
        # dummies.df <- dummyVars(form.dummy, data = DF.left,fullRank=T,levelsOnly=TRUE)# Gide is used as reference
        # 
        # dummies.df <- data.frame(predict(dummies.df, newdata = DF.left))
        # DF.left <- cbind(DF.left,dummies.df)
        # Here I mannually add this
        for (dumCol in dummyCols){
          # dumCol.un <- as.name(dummyCols[1])
          # Add a column
          DF.left <- DF.left %>% add_column(newCol = rep(NA,dim(DF.left)[1]))
          colnames(DF.left)[dim(DF.left)[2]] <-dumCol
        }
        ## Manual fix
        DF.left$tissuegastric <- ifelse(DF.left$tissue=="gastric",1,0)
        DF.left$tissuemelanoma <- ifelse(DF.left$tissue=="melanoma",1,0)
        DF.left$tissueRCC <- ifelse(DF.left$tissue=="RCC",1,0)
      }
      
    }
    
    
  }
  
  # NOW SET FACTOR LEVELS
  print("==SET FACTOR LEVELS TO RESPONSE VARIABLE==")
  data_train[[responseVar]] <- factor(data_train[[responseVar]], levels=ClassLevels)
  data_test[[responseVar]] <- factor(data_test[[responseVar]], levels=ClassLevels)
  if(isTRUE(filter)){
    print("==Filtering applied, merging leftover test data with test data==")
    DF.left[[responseVar]] <- factor(DF.left[[responseVar]], levels=ClassLevels)
    ## Now merge data_test and DF leftover
    data_test <- rbind(data_test,DF.left)
  }
  #####################################################
  ## Run  Logistic Regression model ##
  #####################################################
  
  ## Stratified Repeate CV
  if(isTRUE(controlDum)){
    print("==Adjusting for confounding categorical variable, changing design==")
    controlModel <-paste(dummyCols,collapse = " + ")
    
    modelDesign <- paste0(modelDesign, " + ",controlModel)
  }
  
  # ###########################################################################
  # ## SCALE IF SELECTED - FIRST TRAIN THEN TEST with train's scaling factor ## WRONG-SCALING SHOULD BE IN CV
  # ###########################################################################
  # if(isTRUE(scale)){
  #   ## MEthod A-from caret
  #   print("==Scaling TRAIN data==")
  #   preProcValues <- preProcess(data_train[,mainFeats], method = "scale")
  #   trainTransformed <- predict(preProcValues, data_train[,mainFeats])
  #   print("==Scaling TEST data with the scaling factor from TRAIN data==")
  #   testTransformed <- predict(preProcValues, data_test[,mainFeats])
  #   
  #   
  #   # ## MEthod B- just scale
  #   # trainTransformed2<- data_train %>% dplyr::select(all_of(mainFeats)) %>% scale()
  #   # 
  #   # testTransformed2<- data_test %>% dplyr::select(all_of(mainFeats)) %>% scale(center=attr(trainTransformed2,"scaled:center"),
  #   #                                                                        scale=attr(trainTransformed2,"scaled:scale"))
  #   
  #   data_train <- cbind(trainTransformed,data_train[,(length(mainFeats)+1): dim(data_train)[2]])
  #   data_test <- cbind(testTransformed,data_test[,(length(mainFeats)+1): dim(data_test)[2]])
  #   
  # }
  
  
  ##########################################
  ## CHECK 3 : SAMPLES, TISSUES IN TRAIN  ##
  ##########################################
  
  print("==INITIAL TRAIN DATA==")
  print("=======================")
  print(table(data_train[[responseVar]]))
  print(table(data_train[[dummyVars]]))
  print(table(data_train[["dataset"]]))
  print("=======================")
  
  # For training the model create a train object with only features and response, since other info might confound, especially
  # when out-subsampling
  if(isTRUE(controlDum)){
    data_train2 <-data_train[,c(mainFeats,dummyCols,responseVar)]
  }else{
    data_train2 <-data_train[,c(mainFeats,responseVar)]
  }
  
  
  print("==MODEL TRAIN INITIATION==")
  
  model.logreg <-multiInsideSamplTRAINModels_LogReg.repCV_construct(DFtrain = data_train2,responseVar = responseVar,
                                                                    design = modelDesign,ModelName = model.name, 
                                                                    kfolds = KFOLDS, reps = REPEATS,
                                                                    subMethod=subsampMethod,inoutSub = inoutSub,metricSelected =metricOpt,
                                                                    seedModel = seed.model,seedPart = seed.partition,seeds=setCVseeds, 
                                                                    stratify=stratify,
                                                                    positiveCoefs=positiveCoefs,
                                                                    applyScaling = selectedScale) ##
  
 
  ########################################################
  ########################################################
  ##====================================================##
  ####
  ####
  
  
  ##====================================================##
  ##########################################################
  ##########################################################
  
  
  print("====GATHERING  OUTPUT DATA===")
  
  # print("=STEP 1=")
  # Extract signature coefficients in DF
  model_weights.logreg <-coef(model.logreg[[model.name]]$finalModel)
  # print("=STEP 2=")
  model_weights.logreg.mat <- as.matrix(model_weights.logreg)
  # print("=STEP 3=")
  model_weights.logreg.df <- as.data.frame(model_weights.logreg.mat) ##
  
  # Get intercept and rest coefs sep
  # print("=STEP 4=")
  intercept.logreg <-model_weights.logreg.df[1,]
  # If dummy variable was used remove that as well
  if(isTRUE(controlDum)){
    # Remove last raw, the number is the number of dummyCols
    signature_weights<- as.data.frame(model_weights.logreg.mat[2:(dim(model_weights.logreg.mat)[1]-(length(dummyCols))),])
    
  }else{
    signature_weights<- as.data.frame(model_weights.logreg.mat[-1,])
  }
  
  colnames(signature_weights)<- "weights"
  # print("=STEP 5=")
  ## Extract probabilities on Test data with logit2prob
  # Think about how to include the second signature which is a subset of the first.
  if(isTRUE(controlDum)){
    data_test.s <-data_test[,c(mainFeats,dummyCols)]
  }else{
    # data_test.s <-data_test[,mainFeats]
    data_test.s <-data_test %>% dplyr::select(all_of(mainFeats))
  }
  
  ## HERE IS THE POINT were I would need to preprocess test data - I am not applying the model here, but instead calculating the score.
  if(isTRUE(selectedScale)){
    preProc.train <-model.logreg[[model.name]]$preProcess
    # data_test.scaled <- predict(preProcValues, data_test.s)
    data_test.s<- predict(preProc.train, data_test.s)
    ## THERE IS AN ISSUE WITH SCALING THE DATA WHEN HAVING DUMMY VARS, since preprocess takes everything in!!
    data_test.sc<-data_test.s
    if(isTRUE(controlDum)){
      # data_test.s <-data_test.s[,mainFeats]
      data_test.s <-data_test %>% dplyr::select(all_of(mainFeats))
    }
  }
  
  if(isTRUE(controlDum)){
    # data_test.sub <-data_test.s[,mainFeats]
    data_test.sub <-data_test %>% dplyr::select(all_of(mainFeats))
  }
  # weighted.df <- data_test.s*signature_weights$weights[match(colnames(data_test.s), rownames(signature_weights))][col(data_test.s)]
  # ## Add weighted RNA expression values and generate gene signature score for each sample
  # # print("=STEP 6=")
  # scores <- rowSums(weighted.df) + intercept.lasso
  # df<-as.data.frame(scores)
  # df[[responseVar]]<- data_test[[responseVar]]
  
  # print("=STEP 7=")
  ## Do below when test data are not empty
  if(isTRUE(nrow(data_test.s)<=4)){
    print("No calculations on test data, since they are empty,very few")
    
  }else{
    # inside_models.lasso
    # probsDF.logit<-getProbDF(df,"scores",model.name,ClassLevels,2)##
    # print("=STEP 8=")
    ## Confustion Matrix
    test_pred=extractProb(model.logreg,unkX = data_test.s)
    cm_log= confusionMatrix(test_pred[["pred"]], data_test[[responseVar]])
  }
  
  
  ## WHEN REMOVING OUTLIERS IN FULL DATA
  # print("=STEP 9=")
  ## Gather in a list
  if(isTRUE(nrow(data_test.s)<=4)){
    if (isTRUE(fullResults)){
      results <-list(model.logreg,model_weights.logreg.mat,data_train,data_test,seed.partition,mainFeats)#patient.crpr,patient.pd,outliersMethod.crpr,outliersMethod.pd,dt_crpr,percentileClass1,dt_pd,percentileClass2,
      names(results)<- c("model","model.weights","trainDF","testDF","seed.partition","features")#"outliers.crpr","outliers.pd","outMethodRes.crpr","outMethodRes.pd","pc.crpr","thres.crpr","pc.pd","thres.pd",
    }else{
      results <-list(model.logreg,model_weights.logreg.mat)#,patient.crpr,patient.pd,outliersMethod.crpr,outliersMethod.pdmainFeats
      names(results)<- c("model","model.weights")#,"outliers.crpr","outliers.pd","outMethodRes.crpr","outMethodRes.pd","features"
    }
  }else{
    if (isTRUE(fullResults)){
      results <-list(model.logreg,model_weights.logreg.mat,data_train,data_test,seed.partition,cm_log,mainFeats)#probsDF.logit,,patient.crpr,patient.pd,outliersMethod.crpr,outliersMethod.pd,dt_crpr,percentileClass1,dt_pd,percentileClass2
      names(results)<- c("model","model.weights","trainDF","testDF","seed.partition","confusionMatrix","features")#,"logit2prob.test","outliers.crpr","outliers.pd","outMethodRes.crpr","outMethodRes.pd","pc.crpr","thres.crpr","pc.pd","thres.pd",
    }else{
      results <-list(model.logreg,model_weights.logreg.mat,cm_log,mainFeats)#,patient.crpr,patient.pd,outliersMethod.crpr,outliersMethod.pd
      names(results)<- c("model","model.weights","confusionMatrix","features")#,"outliers.crpr","outliers.pd","outMethodRes.crpr","outMethodRes.pd"
    }
  }
  
  # if(isTRUE(scale)){
  #   results <-c(results, scalingFactor=list(preProcValues))
  # }
  # print("=STEP 10=")
  results
}


checkLogRegModelTEST <- function(DF,models.all,model.name, obsColumn,filtered = FALSE,scaled=FALSE,rfeE=FALSE){
  # DF <-immdata.apd1.proc.TcrA.divEst.full.sub.filt
  # models.all <-model_nb
  # model.name <-"NB"
  # filtered = TRUE
  # scaled=TRUE
  # obsColumn<-"run_accession"

  
  model.final.res <- models.all
  
  if (isTRUE(scaled)){
    print(">>> Extract scaling factor <<<")
    if(isTRUE(rfeE)){
      trainScaleFactor<-model.final.res$model$fit$preProcess
    }else{
      trainScaleFactor<-model.final.res$model[[model.name]]$preProcess
    }
    
  }
  
  ## FIND TEST DATA FIRST
  ## Initial train object from initial partition, AFTER OUTLIER & FILTERING , split,NO SCALING
  train.df<-models.all$trainDF
  ## Initial test object from initial partition, AFTER OUTLIER & FILTERING , split, NO SCALING
  test.df<- models.all$testDF
  # If filtering has been applied there might be extra test data left in original set. We need to include those in test data.
  if(isTRUE(filtered)){
    print("There might be left over test data in original data")
    postfiltSamples <- c(train.df[[obsColumn]],test.df[[obsColumn]]) #  test df here includes the filtered test data
    extraTestSamples<-setdiff(DF[[obsColumn]],postfiltSamples)
    finalTestSamples <-c(test.df[[obsColumn]],extraTestSamples)
    obsColumn.un <- as.name(obsColumn)
    testDF <- DF %>% dplyr::filter(!!obsColumn.un %in% finalTestSamples) %>% droplevels()
  }else{
    postfiltSamples <- c(train.df[[obsColumn]],test.df[[obsColumn]]) #  test df here includes the filtered test data
    extraTestSamples<-setdiff(DF[[obsColumn]],postfiltSamples)
    finalTestSamples <-c(test.df[[obsColumn]],extraTestSamples)
    obsColumn.un <- as.name(obsColumn)
    testDF <- DF %>% dplyr::filter(!!obsColumn.un %in% finalTestSamples) %>% droplevels()
  }
  print("==Checking if there is a leak of train in test data")
  intersect(train.df$run_accession,testDF$run_accession)
  print("===TOTAL TEST DATA==")
  print(table(testDF$tissue))
  print(table(testDF$dataset))
  testDF
}


checkLogRegModelTEST.rfe <- function(DF,models.all,model.name, obsColumn,filtered = FALSE,scaled=FALSE){
  # DF <-immdata.apd1.proc.TcrA.divEst.full.sub.filt
  # models.all <-models.all
  # model.name <-"smote"
  # filtered = TRUE
  # scaled=TRUE
  # obsColumn<-"run_accession"

  
  model.final.res <- models.all
  
  # if (isTRUE(scaled)){
  #   print(">>> Extract scaling factor <<<")
  #   trainScaleFactor<-model.final.res$model[[model.name]]$preProcess
  # }
  
  ## FIND TEST DATA FIRST
  ## Initial train object from initial partition, AFTER OUTLIER & FILTERING , split,NO SCALING
  train.df<-models.all[[model.name]]$trainDF
  ## Initial test object from initial partition, AFTER OUTLIER & FILTERING , split, NO SCALING
  test.df<- models.all[[model.name]]$testDF
  # If filtering has been applied there might be extra test data left in original set. We need to include those in test data.
  if(isTRUE(filtered)){
    print("There might be left over test data in original data")
    postfiltSamples <- c(train.df[[obsColumn]],test.df[[obsColumn]]) #  test df here includes the filtered test data
    extraTestSamples<-setdiff(DF[[obsColumn]],postfiltSamples)
    finalTestSamples <-c(test.df[[obsColumn]],extraTestSamples)
    obsColumn.un <- as.name(obsColumn)
    testDF <- DF %>% dplyr::filter(!!obsColumn.un %in% finalTestSamples) %>% droplevels()
  }else{
    postfiltSamples <- c(train.df[[obsColumn]],test.df[[obsColumn]]) #  test df here includes the filtered test data
    extraTestSamples<-setdiff(DF[[obsColumn]],postfiltSamples)
    finalTestSamples <-c(test.df[[obsColumn]],extraTestSamples)
    obsColumn.un <- as.name(obsColumn)
    testDF <- DF %>% dplyr::filter(!!obsColumn.un %in% finalTestSamples) %>% droplevels()
  }
  # print("==Checking if there is a leak of train in test data")
  # intersect(train.df$run_accession,testDF$run_accession)
  # print("===TOTAL TEST DATA==")
  # print(table(testDF$tissue))
  # print(table(testDF$dataset))
  testDF
}


extract_EvaluationMetrics<- function(finalModel, test.data, classes,LevelIndex, 
                                     modelName=NULL,var=NULL,
                                     savePlots=FALSE, 
                                     getTable=NULL,printP=TRUE, thres = NULL, rfeE = FALSE){
  
  # finalModel <-train.rfe[[1]]
  # test.data <- testDF[[1]]
  # classes <- c("CRPR","PD")
  # LevelIndex<-1
  # modelName <- names(train.rfe)[1]
  # var=sigs[1]
  # getTable = TRUE
  # savePlots=FALSE
  # rfeE = TRUE
  # thres=NULL
  
  
  # finalModel <-train.rfe[[1]]
  # test.data <- testDF
  # classes <- c("CRPR","PD")
  # LevelIndex<-1
  # modelName <- names(train.rfe)[1]
  # var=sigs[1]
  # getTable = TRUE
  # savePlots=FALSE
  # rfeE = TRUE
  # thres=NULL

  if(!isTRUE(rfeE)){
    probabilities <- predict(finalModel, newdata=test.data, type="prob")
    predicted.classes <- predict(finalModel, newdata=test.data, type="raw")
    predClassesProbs <-data.frame(obs = factor(test.data$response, levels = classes),pred = factor(predicted.classes, levels = classes))
    predClassesProbs <- cbind(predClassesProbs,probabilities)
    
  }else{
    probabilities <- predict(finalModel, newdata=test.data, type="prob")
    probabilities$obs <- factor(test.data$response, levels = classes)
    predClassesProbs <-probabilities
  }
  
  if( !is.null(thres) ){
    # Change predicted based on threshold
    thres_preds <- function(predDF,thr,classes){
      predDF[["pred"]] <- ifelse(predDF[[classes[LevelIndex]]] >=thr,classes[LevelIndex],classes[-LevelIndex])
      predDF[["pred"]] <- factor(predDF[["pred"]],levels=levels(predDF[["obs"]]))
      predDF
    }
    prob_test <- thres_preds(predClassesProbs,thres,classes)
  }
  
  ## Metrics
  # Precision-Recall
  pr <-prSummary(predClassesProbs, lev = classes)# HERE I can change levels
  pr.auc <- pr[1]
  precision <- pr[2]
  recall<- pr[3]
  f1.score <- pr[4]
  # ROC, sensitivity, specificity
  # Change levels-positive class
  
  if(isTRUE(all(classes==levels(predClassesProbs$obs)))){
    roc.sens.spec <- twoClassSummary(predClassesProbs, lev = classes)# HERE I CANNOT change levels    
    
  }else{
    predClassesProbs.c <-predClassesProbs
    predClassesProbs.c$obs <-factor(predClassesProbs.c$obs, levels = classes)
    predClassesProbs.c$pred <-factor(predClassesProbs.c$pred, levels = classes)
    roc.sens.spec <- twoClassSummary(predClassesProbs.c, lev = classes)
  }
  
  roc.auc <- roc.sens.spec[1]
  sens <- roc.sens.spec[2]
  spec <- roc.sens.spec[3]
  youden <- sens+spec-1
  # Accuracy, Balanced accuracy
  if(isTRUE(all(classes==levels(predClassesProbs$obs)))){
    acc.all <- caret::confusionMatrix(predClassesProbs$pred,predClassesProbs$obs,positive=classes[LevelIndex])
    
  }else{
    acc.all <- caret::confusionMatrix(prob_test.final.c$pred,prob_test.final.c$obs,positive=classes[LevelIndex])
  }
  
  
  ppv <- acc.all$byClass[[3]]
  npv <-acc.all$byClass[[4]]
  ################################
  
  ggplotConfusionMatrix <- function(m){
    mytitle <- paste("Accuracy", percent_format()(m$overall[1]),
                     "Kappa", percent_format()(m$overall[2]))
    conf.table <-as.data.frame(m$table)
    
    p <-ggplot(data = conf.table ,
             aes(x = Reference, y = Prediction)) +
      geom_tile(aes(fill = log(Freq)), colour = "white") +
      scale_fill_gradient(low = "white", high = "steelblue") +
      geom_text(aes(x = Reference, y = Prediction, label = Freq, size=2.5)) +
      ggtitle(modelName)
    p <- p + theme_ipsum(base_family = "Palatino Linotype", base_size=13)
    p <- p + theme(axis.text.x = element_text(hjust = 1, vjust = 1,size=13),
                   axis.title.x = element_text(size=14),#, face="bold"
                   axis.text.y = element_text( angle=90,hjust = 1, vjust = 1,size=13),
                   axis.title.y = element_text(size=14),#, face="bold"
                   plot.title = element_text(size = 13, face = "bold"),
                   legend.position = "none",
                   plot.margin = unit(c(0, 0, 0, 0), "null"))
    p
  }
  if(isTRUE(printP)){
    # print("Printing stuff!")
    print(ggplotConfusionMatrix(acc.all))
  }
  
  if(isTRUE(savePlots)){
    eval<-ggplotConfusionMatrix(acc.all)
    save(eval,file=paste0(projectRDataPath,"figures/models/",var,".eval.aPD1",saveVar,Rdata.suffix))
  }
  ################################
  conf.table <- acc.all$table
  conf.table <- as.matrix(conf.table)
  
  rownames(conf.table) <- c(paste0("Predicted Positive (",classes[LevelIndex],")" ), paste0("Predicted Negative (",classes[-LevelIndex],")" ))
  colnames(conf.table) <- c(paste0("Actual Positive(",classes[LevelIndex],")" ), paste0("Actual Negative(",classes[-LevelIndex],")" ))
  
  acc <-  acc.all$overall[1]
  acc.bal <- acc.all$byClass[11]
  # print out confusion matrix
  #pander(conf.table)
  # results all in df
  res.model <- rbind(acc,acc.bal,sens,spec,roc.auc,precision,recall, pr.auc,f1.score,youden,ppv,npv)
  colnames(res.model)<- "value"
  rownames(res.model) <- c("Accuracy", "Balanced Accuracy", "Sensitivity", "Specificity", "ROC-AUC", "Precision", "Recall", "PR-AUC", "F1-Score","Youden's J-statistic", "PPV","NPV")
  res.model
  
  if(is.null(getTable)){
    res.model
  }else{
    panderOptions("table.style", "simple")
    # print("printing table")
    invisible(pander(res.model))
  }
}

ROCplot_multiModels_sameTestdata <- function(modelList, test.data,testX,testY,classes, LevelIndex, responseVar, modelNamesVec, dsids, yEndPos,xPos, thres = NULL,yStartPos=0.1,rfeE=FALSE,sizeAnnot=5,
                                             mixedMode = FALSE,trainModelsVec = NULL, rfeModelsVec = NULL){
  # modelList <-models.gastric
  # test.data <- test.DF.gastric
  # testX = data.testX.gastric
  # testY=test.DF.gastric$response
  # classes <- c("CRPR","PD")
  # LevelIndex<-1
  # responseVar <- "response"
  # modelNamesVec <- names(models.gastric)
  # dsids <-1:length(models.gastric)
  # thres <- NULL
  # yEndPos<- 0.5
  # xPos<- 0.8
  # yStartPos=0.1
  # rfeE=FALSE
  # sizeAnnot=5
  # mixedMode = FALSE
  # trainModelsVec = NULL
  # rfeModelsVec = NULL

  # if(isTRUE(rfeE)){
  #   probabilities <- predict(modelList, newdata=test.data, type="prob")
  #   predicted.classes <- predict(modelList, newdata=test.data, type="raw")
  #   predClassesProbs <-data.frame(obs = factor(test.data[[yVar]], levels = classes))
  #   prob_test <- cbind(predClassesProbs,probabilities)
  #   prob_test$object <- rep(modelName, dim(prob_test)[1])
  #   
  # }else{
  #   prob_test<-extractProb(modelList,unkX =testX)
  #   prob_test$obs <-dataSet[[responseVar]]
  #   prob_test$obs <-factor(prob_test$obs, levels = classes)
  #   prob_test$pred <-factor(prob_test$pred, levels = classes)
  # }
  
  if(isTRUE(rfeE)){
    # b <-sapply(1: length(modelList), function(x) predict(modelList[[x]], newdata=test.data), simplify = FALSE)
    probabilities <- predict(modelList, newdata=test.data, type="prob")
    # probabilities.df <- bind_rows(probabilities, .id = "object")
    # predicted.classes <- predict(modelList, newdata=test.data, type="raw")
    
    probabilities <- lapply(probabilities, function(x) transform(x, obs = factor(test.data[[responseVar]], levels = classes)))
    probabilities <- lapply(probabilities, function(x) transform(x, pred = factor(pred, levels = classes)))
    predClassesProbs <-bind_rows(probabilities, .id = "object")
    prob_test <-predClassesProbs
    # predClassesProbs <- cbind(predClassesProbs,probabilities)
    
  }else if(isTRUE(mixedMode)){
    ## First train objects
    if( !is.null(testX) && ! is.null(testY)){
      prob_test.t<-extractProb(modelList[trainModelsVec],unkX =testX)
      prob_test.t$obs <-test.data[[responseVar]]
    }else{
      prob_test.t<-extractProb(modelList)    
    }
    # Bring those in the format of  train models
    prob_test.t <-prob_test.t %>% dplyr::select(all_of(classes),"obs", "pred","object")
    
    ## Then RFE objects
    probabilities <- predict(modelList[rfeModelsVec], newdata=test.data, type="prob")
    # probabilities.df <- bind_rows(probabilities, .id = "object")
    # predicted.classes <- predict(modelList, newdata=test.data, type="raw")
    
    probabilities <- lapply(probabilities, function(x) transform(x, obs = factor(test.data[[responseVar]], levels = classes)))
    probabilities <- lapply(probabilities, function(x) transform(x, pred = factor(pred, levels = classes)))
    predClassesProbs <-bind_rows(probabilities, .id = "object")
    prob_test.r <-predClassesProbs
    # Bring those in the format of  train models
    prob_test.r <-prob_test.r %>% dplyr::select(all_of(classes),"obs", "pred","object")
    ## Merge objects
    prob_test <- rbind(prob_test.t, prob_test.r)
    
  }else{
    if( !is.null(testX) && ! is.null(testY)){
      prob_test<-extractProb(modelList,unkX =testX)
      prob_test$obs <-test.data[[responseVar]]
    }else{
      prob_test<-extractProb(modelList)    
    }
    
  }
  
  
  if( !is.null(thres) ){
    # Change predicted based on threshold
    thres_preds <- function(predDF,thr,classes){
      predDF[["pred"]] <- ifelse(predDF[[classes[LevelIndex]]] >=thr,classes[LevelIndex],classes[-LevelIndex])
      predDF[["pred"]] <- factor(predDF[["pred"]],levels=levels(predDF[["obs"]]))
      predDF
    }
    prob_test <- thres_preds(prob_test,thres,classes)
  }
  
  # Multiple models
  prob_test.list <- split(prob_test,prob_test$object)
  modelNamesVec <-names(prob_test.list)
  # modelList <- inside_models.rcv.m1
  # testDF <- data_test
  #test.modelList.Scores <- lapply(inside_models.rcv.m1, scoresModel, data = data_test, LevelIndex=1)
  # Join scores
  scores.all <- join_scores(lapply(prob_test.list, function(x) x[[classes[LevelIndex]]]))
  names(scores.all) <- NULL
  labels <- join_labels(lapply(prob_test.list, function(x) ifelse(x$obs==classes[LevelIndex],1,0)))
  
  # Make input dataset for evalmod
  msmdat <- mmdata(scores.all, labels,modnames = modelNamesVec, dsids=dsids)
  mscurves1 <- evalmod(msmdat)
  
  ######
  auc_res <- data.table(precrec::auc(mscurves1))
  ySeqs <-seq(yEndPos,yStartPos,length=length(modelNamesVec))
  
  # Set colors
  # if(length(modelNamesVec) ==3){
  #   cols = c("#72AA7E","#FFC107","#901401")
  # }else if (length(modelNamesVec) >=12){
  #   qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  #   col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  #   cols <- sample(col_vector, length(modelNamesVec))
  # }else if (length(modelNamesVec)>8){
  #   cols = brewer.pal(10,"Set2")
  # }else{
  #   cols = brewer.pal(length(modelNamesVec),"Dark2")
  # }
  
  if(length(modelNamesVec) <= 6){
    # palette <- myCBfriendly.pale.small
    # set.seed(1562)
    # palette <-sample(myCBfriendly.muted.medium)
    palette <- pal_uchicago("default")(6)
    cols <- palette
  }else if (length(modelNamesVec) <= 9){
    set.seed(1562)
    # palette <- cbPalette
    # palette <-  sample(brewer.pal(9,"Paired"))
    # palette <- viridis_pal(option="D")(9)
    
    # palette <- sample(myCBfriendly.light.medium)
    # palette <- c(brewer.pal(9,"Set2"),brewer.pal(8,"Dark2"),"#cc3311")
    palette <- pal_uchicago("dark")(9)
    cols <- palette
  }else{
    set.seed(1562)
    palette <- distinctColorPalette(length(modelNamesVec))
    cols <- palette
    # palette <- c(brewer.pal(9,"Set2"),brewer.pal(8,"Dark2"),"#cc3311")
  }
  
  
  p1 <-autoplot(mscurves1, "ROC", color_palette(cols))+ geom_line(size=1.5)
  p1 <- p1 + scale_fill_manual("Models",values = cols)
  p1 <- p1 + scale_color_manual(values = cols) 
  
  
  for (i in 1:length(modelNamesVec)){
    p1 <- p1 + annotate("label", x = xPos, y = ySeqs[i], 
                        label = paste0(modelNamesVec[i], '\n ROC-AUC = ', round(auc_res[modnames == modelNamesVec[i] & curvetypes == "ROC", aucs], 3)),size = sizeAnnot)
  }
  p1 <- p1 + theme_ipsum(base_family = "Palatino Linotype", base_size=13)
  p1 <- p1 + theme(axis.text.x = element_text(hjust = 1, vjust = 1,size=13),
                   axis.title.x = element_text(size=16, face="bold"),
                   axis.text.y = element_text( hjust = 1, vjust = 1,size=13),
                   axis.title.y = element_text(size=16, face="bold"),
                   legend.title = element_text(size=16,color="black"),
                   legend.text = element_text(size=15),
                   
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   axis.line = element_line(colour = "black"),
                   legend.direction = 'vertical',
                   legend.box = 'vertical',
                   legend.background = element_rect(size = rel(1.5),color="white", fill=alpha("white",0.6)),
                   legend.position = "bottom",
                   legend.key.size = unit(0.4, "cm"),
                   plot.margin = unit(c(0, 0, 0, 0), "null"),
                   plot.title = element_text(size = 13))
  p1 <- p1 + labs(fill='Model', color='Model') 
  p1 <- p1 + geom_abline(colour="black", linetype = "dashed")
  p1
}

ROCcalculation_models <- function(modelList, test.data,testX,testY,classes, LevelIndex, responseVar, modelNamesVec, dsids, thres = NULL,rfeE=FALSE,mixedMode = FALSE,trainModelsVec = NULL, rfeModelsVec = NULL){
  
  # modelList = train.rfe
  # test.data = testDF[[1]]
  # testX =  data.testX
  # testY = testDF[[1]]$response
  # classes = c("CRPR","PD")
  # LevelIndex =1
  # responseVar = "response"
  # modelNamesVec = names(train.rfe)
  # dsids = 1:length(train.rfe)
  # rfeE=TRUE
  # thres = NULL
  # yEndPos =  0.5,xPos = 0.9, thres = NULL,yStartPos=0.01,rfeE=TRUE,sizeAnnot=3
  
  if(isTRUE(rfeE)){
    # b <-sapply(1: length(modelList), function(x) predict(modelList[[x]], newdata=test.data), simplify = FALSE)
    probabilities <- predict(modelList, newdata=test.data, type="prob")
    # probabilities.df <- bind_rows(probabilities, .id = "object")
    # predicted.classes <- predict(modelList, newdata=test.data, type="raw")
    
    probabilities <- lapply(probabilities, function(x) transform(x, obs = factor(test.data[[responseVar]], levels = classes)))
    probabilities <- lapply(probabilities, function(x) transform(x, pred = factor(pred, levels = classes)))
    predClassesProbs <-bind_rows(probabilities, .id = "object")
    prob_test <-predClassesProbs
    # predClassesProbs <- cbind(predClassesProbs,probabilities)
    
  }else if(isTRUE(mixedMode)){
    ## First train objects
    if( !is.null(testX) && ! is.null(testY)){
      prob_test.t<-extractProb(modelList[trainModelsVec],unkX =testX)
      prob_test.t$obs <-test.data[[responseVar]]
    }else{
      prob_test.t<-extractProb(modelList)    
    }
    # Bring those in the format of  train models
    prob_test.t <-prob_test.t %>% dplyr::select(all_of(classes),"obs", "pred","object")
    
    ## Then RFE objects
    probabilities <- predict(modelList[rfeModelsVec], newdata=test.data, type="prob")
    # probabilities.df <- bind_rows(probabilities, .id = "object")
    # predicted.classes <- predict(modelList, newdata=test.data, type="raw")
    
    probabilities <- lapply(probabilities, function(x) transform(x, obs = factor(test.data[[responseVar]], levels = classes)))
    probabilities <- lapply(probabilities, function(x) transform(x, pred = factor(pred, levels = classes)))
    predClassesProbs <-bind_rows(probabilities, .id = "object")
    prob_test.r <-predClassesProbs
    # Bring those in the format of  train models
    prob_test.r <-prob_test.r %>% dplyr::select(all_of(classes),"obs", "pred","object")
    ## Merge objects
    prob_test <- rbind(prob_test.t, prob_test.r)
    
  }else{
    if( !is.null(testX) && ! is.null(testY)){
      prob_test<-extractProb(modelList,unkX =testX)
      prob_test$obs <-test.data[[responseVar]]
    }else{
      prob_test<-extractProb(modelList)    
    }
    
  }
  
  
  if( !is.null(thres) ){
    # Change predicted based on threshold
    thres_preds <- function(predDF,thr,classes){
      predDF[["pred"]] <- ifelse(predDF[[classes[LevelIndex]]] >=thr,classes[LevelIndex],classes[-LevelIndex])
      predDF[["pred"]] <- factor(predDF[["pred"]],levels=levels(predDF[["obs"]]))
      predDF
    }
    prob_test <- thres_preds(prob_test,thres,classes)
  }
  
  # Multiple models
  prob_test.list <- split(prob_test,prob_test$object)
  modelNamesVec <-names(prob_test.list)
  # modelList <- inside_models.rcv.m1
  # testDF <- data_test
  #test.modelList.Scores <- lapply(inside_models.rcv.m1, scoresModel, data = data_test, LevelIndex=1)
  # Join scores
  scores.all <- join_scores(lapply(prob_test.list, function(x) x[[classes[LevelIndex]]]))
  names(scores.all) <- NULL
  labels <- join_labels(lapply(prob_test.list, function(x) ifelse(x$obs==classes[LevelIndex],1,0)))
  
  # Make input dataset for evalmod
  msmdat <- mmdata(scores.all, labels,modnames = modelNamesVec, dsids=dsids)
  mscurves1 <- evalmod(msmdat)
  
  ######
  auc_res <- data.table(precrec::auc(mscurves1))
  roc.only <- auc_res %>% dplyr::filter(curvetypes =="ROC")
  roc.only
}


PRCplot_multiModels_sameTestdata2 <- function(modelList, test.data,testX,testY,classes, LevelIndex, responseVar, modelNamesVec, dsids, yEndPos,xPos, thres = NULL,yStartPos=0.1,rfeE=FALSE,sizeAnnot=5,
                                              
                                              mixedMode = FALSE,trainModelsVec = NULL, rfeModelsVec = NULL){
  
  # modelList <-models.final
  # test.data <- testDF[[1]]
  # testX = data.testX
  # testY=testDF[[1]]$response
  # classes <- c("CRPR","PD")
  # LevelIndex<-1
  # responseVar <- "response"
  # modelNamesVec <- names(models.final)[1]
  # dsids <-1:length(models.final)
  # thres <- NULL
  # yEndPos<- 0.5
  # xPos<- 0.8
  # yStartPos=0.1
  # 
  if(isTRUE(rfeE)){
    # b <-sapply(1: length(modelList), function(x) predict(modelList[[x]], newdata=test.data), simplify = FALSE)
    probabilities <- predict(modelList, newdata=test.data, type="prob")
    # probabilities.df <- bind_rows(probabilities, .id = "object")
    # predicted.classes <- predict(modelList, newdata=test.data, type="raw")
    
    probabilities <- lapply(probabilities, function(x) transform(x, obs = factor(test.data[[responseVar]], levels = classes)))
    probabilities <- lapply(probabilities, function(x) transform(x, pred = factor(pred, levels = classes)))
    predClassesProbs <-bind_rows(probabilities, .id = "object")
    prob_test <-predClassesProbs
    # predClassesProbs <- cbind(predClassesProbs,probabilities)
    
  }else if(isTRUE(mixedMode)){
    ## First train objects
    if( !is.null(testX) && ! is.null(testY)){
      prob_test.t<-extractProb(modelList[trainModelsVec],unkX =testX)
      prob_test.t$obs <-test.data[[responseVar]]
    }else{
      prob_test.t<-extractProb(modelList)    
    }
    # Bring those in the format of  train models
    prob_test.t <-prob_test.t %>% dplyr::select(all_of(classes),"obs", "pred","object")
    
    ## Then RFE objects
    probabilities <- predict(modelList[rfeModelsVec], newdata=test.data, type="prob")
    # probabilities.df <- bind_rows(probabilities, .id = "object")
    # predicted.classes <- predict(modelList, newdata=test.data, type="raw")
    
    probabilities <- lapply(probabilities, function(x) transform(x, obs = factor(test.data[[responseVar]], levels = classes)))
    probabilities <- lapply(probabilities, function(x) transform(x, pred = factor(pred, levels = classes)))
    predClassesProbs <-bind_rows(probabilities, .id = "object")
    prob_test.r <-predClassesProbs
    # Bring those in the format of  train models
    prob_test.r <-prob_test.r %>% dplyr::select(all_of(classes),"obs", "pred","object")
    ## Merge objects
    prob_test <- rbind(prob_test.t, prob_test.r)
    
  }else{
    if( !is.null(testX) && ! is.null(testY)){
      prob_test<-extractProb(modelList,unkX =testX)
      prob_test$obs <-test.data[[responseVar]]
    }else{
      prob_test<-extractProb(modelList)    
    }
    
  }
  
  
  
  if( !is.null(thres) ){
    # Change predicted based on threshold
    thres_preds <- function(predDF,thr,classes){
      predDF[["pred"]] <- ifelse(predDF[[classes[LevelIndex]]] >=thr,classes[LevelIndex],classes[-LevelIndex])
      predDF[["pred"]] <- factor(predDF[["pred"]],levels=levels(predDF[["obs"]]))
      predDF
    }
    prob_test <- thres_preds(prob_test,thres,classes)
  }
  
  # Multiple models
  prob_test.list <- split(prob_test,prob_test$object)
  modelNamesVec <-names(prob_test.list)
  # modelList <- inside_models.rcv.m1
  # testDF <- data_test
  #test.modelList.Scores <- lapply(inside_models.rcv.m1, scoresModel, data = data_test, LevelIndex=1)
  # Join scores
  scores.all <- join_scores(lapply(prob_test.list, function(x) x[[classes[LevelIndex]]]))
  names(scores.all) <- NULL
  labels <- join_labels(lapply(prob_test.list, function(x) ifelse(x$obs==classes[LevelIndex],1,0)))
  
  # Make input dataset for evalmod
  msmdat <- mmdata(scores.all, labels,modnames = modelNamesVec, dsids=dsids)
  mscurves1 <- evalmod(msmdat)
  
  ######
  auc_res <- data.table(precrec::auc(mscurves1))
  ySeqs <-seq(yEndPos,yStartPos,length=length(modelNamesVec))
  
  # Set colors
  
  if (length(modelNamesVec) >=12){
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    cols <- sample(col_vector, length(modelNamesVec))
  }else if (length(modelNamesVec)>8){
    cols = brewer.pal(10,"Set2")
  }else{
    cols = brewer.pal(length(modelNamesVec),"Dark2")
  }
  
  p1 <- autoplot(mscurves1, "PRC", color_palette(cols)) + geom_line(size=1.5)
  p1 <- p1 + scale_fill_manual("Models",values = cols)
  p1 <- p1 + scale_color_manual(values = cols) 
  
  # add_annot <- function(plotObj, x, yseqObj, aucObj, ind, modelnames){
  #   plotObj <- plotObj + annotate("text", x = x, y = yseqObj[ind], 
  #                     label = paste0(modelnames[ind], '~ AUC = ', round(aucObj[modnames == modelnames[ind] & curvetypes == "ROC", aucs], 2)))
  #   plotObj
  #   #print(plotObj)
  # }
  
  for (i in 1:length(modelNamesVec)){
    p1 <- p1 + annotate("label", x = xPos, y = ySeqs[i], 
                        label = paste0(modelNamesVec[i], '\n PR-AUC = ', round(auc_res[modnames == modelNamesVec[i] & curvetypes == "PRC", aucs], 4)),size = sizeAnnot)
  }
  p1 <- p1 + theme_ipsum(base_family = "Palatino Linotype", base_size=13)
  p1 <- p1 + theme(axis.text.x = element_text(hjust = 1, vjust = 1,size=13),
                   axis.title.x = element_text(size=16, face="bold"),
                   axis.text.y = element_text( hjust = 1, vjust = 1,size=13),
                   axis.title.y = element_text(size=16, face="bold"),
                   legend.title = element_text(size=16,color="black"),
                   legend.text = element_text(size=15),
                   
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   axis.line = element_line(colour = "black"),
                   legend.direction = 'vertical',
                   legend.box = 'vertical',
                   legend.background = element_rect(size = rel(1.5),color="white", fill=alpha("white",0.6)),
                   legend.position = "bottom",
                   legend.key.size = unit(0.4, "cm"),
                   plot.margin = unit(c(0, 0, 0, 0), "null"),
                   plot.title = element_text(size = 13))
  p1 <- p1 + labs(fill='Model', color='Model') 
  #p1 <- p1 + geom_hline(colour="black", linetype = "dashed")
  # print(p1)
  p1
}


compareTwoRocCurves <- function(modelList, test.data, yVar, classes,LevelIndex,compareName, testX = NULL, testY=NULL, thres = NULL,rfeE=FALSE,
                                mixedMode = FALSE,trainModelsVec = NULL, rfeModelsVec = NULL, saveD = FALSE){
  
  # names(train.modelsSolo) <- c('M1: R ~ TCRa richness','M2: R ~ BCRigh richness', 
  #                              "M3: R ~ logTMB",'M4: R ~ PD-L1','M5: R ~ TIS GEP')
  # modelNames<- names(train.modelsSolo)
  # 
  # model.pairwiseCombos <-split(t(combn(modelNames,2)), f = 1:ncol(combn(modelNames,2)))
  # modelList <- train.modelsSolo[which(names(train.modelsSolo) %in% model.pairwiseCombos$`1`)]
  # test.data <-testDF[[1]]
  # yVar <-"response"
  # classes <-ClassLevels
  # LevelIndex <-1
  # compareName <-model.pairwiseCombos$`1`
  # testX = data.testX
  # testY=testDF[[1]]$response
  # thres = NULL
  # rfeE=FALSE
  # mixedMode = FALSE
  # trainModelsVec = NULL
  # rfeModelsVec = NULL
  # saveD=FALSE
  # 
  
  if(!isTRUE(saveD)){
    cat('\n')  
    cat('\n') 
    cat("#### ", paste(compareName, collapse = " vs. "),"\n") 
    cat('\n') 
  }
  
  
  if(isTRUE(rfeE)){
    # b <-sapply(1: length(modelList), function(x) predict(modelList[[x]], newdata=test.data), simplify = FALSE)
    probabilities <- predict(modelList, newdata=test.data, type="prob")
    # probabilities.df <- bind_rows(probabilities, .id = "object")
    # predicted.classes <- predict(modelList, newdata=test.data, type="raw")
    
    probabilities <- lapply(probabilities, function(x) transform(x, obs = factor(test.data[[yVar]], levels = classes)))
    probabilities <- lapply(probabilities, function(x) transform(x, pred = factor(pred, levels = classes)))
    predClassesProbs <-bind_rows(probabilities, .id = "object")
    prob_test <-predClassesProbs
    # predClassesProbs <- cbind(predClassesProbs,probabilities)
    
  }else if(isTRUE(mixedMode)){
    ## First train objects
    if( !is.null(testX) && ! is.null(testY)){
      prob_test.t<-extractProb(modelList[trainModelsVec],unkX =testX)
      prob_test.t$obs <-test.data[[yVar]]
    }else{
      prob_test.t<-extractProb(modelList)    
    }
    # Bring those in the format of  train models
    prob_test.t <-prob_test.t %>% dplyr::select(all_of(classes),"obs", "pred","object")
    
    ## Then RFE objects
    probabilities <- predict(modelList[rfeModelsVec], newdata=test.data, type="prob")
    # probabilities.df <- bind_rows(probabilities, .id = "object")
    # predicted.classes <- predict(modelList, newdata=test.data, type="raw")
    
    probabilities <- lapply(probabilities, function(x) transform(x, obs = factor(test.data[[yVar]], levels = classes)))
    probabilities <- lapply(probabilities, function(x) transform(x, pred = factor(pred, levels = classes)))
    predClassesProbs <-bind_rows(probabilities, .id = "object")
    prob_test.r <-predClassesProbs
    # Bring those in the format of  train models
    prob_test.r <-prob_test.r %>% dplyr::select(all_of(classes),"obs", "pred","object")
    ## Merge objects
    prob_test <- rbind(prob_test.t, prob_test.r)
    
  }else{
    if( !is.null(testX) && ! is.null(testY)){
      prob_test<-extractProb(modelList,unkX =testX)
      prob_test$obs <-test.data[[yVar]]
    }else{
      prob_test<-extractProb(modelList)    
    }
    
  }
  
  
  
  
  if( !is.null(thres) ){
    # Change predicted based on threshold
    thres_preds <- function(predDF,thr,classes){
      predDF[["pred"]] <- ifelse(predDF[[classes[LevelIndex]]] >=thr,classes[LevelIndex],classes[-LevelIndex])
      predDF[["pred"]] <- factor(predDF[["pred"]],levels=levels(predDF[["obs"]]))
      predDF
    }
    prob_test <- thres_preds(prob_test,thres,classes)
  }
  
  prob_test.final <- split(prob_test, prob_test$object)
  aucs <- lapply(prob_test.final, function(x) format(pROC::auc(ifelse(x$obs==classes[LevelIndex],1,0),x[[classes[LevelIndex]]]), digits=5, scientific=FALSE))
  rocs<- lapply(prob_test.final, function(x) pROC::roc(ifelse(x$obs==classes[LevelIndex],1,0),x[[classes[LevelIndex]]]))
  # Delongs test
  # Delongs test
  set.seed(3456)
  delong <- pROC::roc.test(rocs[[1]], rocs[[2]])
  delong <- as.data.frame(cbind(round(delong$statistic,3),round(delong$p.value,3),t(delong$estimate)))
  colnames(delong)[1:2]<- c("Test statistic", "P-value")
  invisible(pander(delong))
  # Bootstrap
  # roc.test(rocs[[1]], rocs[[2]],method="bootstrap",boot.n=10000)
  boot <- pROC::roc.test(rocs[[1]], rocs[[2]],method="bootstrap")
  boot <- as.data.frame(cbind(round(boot$statistic,3),round(boot$p.value,3),t(boot$estimate)))
  colnames(boot)[1:2]<- c("Test statistic", "P-value")
  
  if(!isTRUE(saveD)){
    invisible(pander(boot))
    
  }else{
    boot
  }
}

optimal_threshold <- function(thresSeq,modelList,dataSet,yVar, testX=NULL,testY=NULL, modelName ,classes,LevelIndex, MetrixMax="F1.Score", extractedProbDF=NULL,rfeMode = FALSE
                              # mixedMode = FALSE,
                              # trainModelsVec = NULL,
                              # rfeModelsVec = NULL
                              ){
  
  # thresSeq<-seq(0.1,0.99,length.out = 50)
  # modelList<- models.final[[1]]
  # dataSet<-testDF[[1]]
  # testX=data.testX
  # testY=testDF[[1]]$response
  # modelName<- names(models.final)[1]
  # classes<- ClassLevels
  # LevelIndex<-1
  # MetrixMax="F1.Score"
  # extractedProbDF=NULL
  # rfeMode = FALSE
  
  # Calculate Evaluation metrics for all thresholds
  metrics.thres <-lapply(thresSeq, function(x) threshold_moving(modelList,dataSet, yVar, classes, LevelIndex,testX = testX, testY=testY,x,modelName,extractedProbDF=extractedProbDF,rfeE = rfeMode
                                                                # mixedMode = mixedMode,
                                                                # trainModelsVec = trainModelsVec,
                                                                # rfeModelsVec = rfeModelsVec
                                                                ))
  
  names(metrics.thres) <- as.character(thresSeq)
  metrics.thres.df<-list.cbind(metrics.thres)
  colnames(metrics.thres.df) <-as.character(thresSeq)
  # Find classification Threshold that maximizes F-1 Score
  selectedThreshold <- max(as.numeric(names(which(metrics.thres.df[which(rownames(metrics.thres.df)==MetrixMax),]==max(metrics.thres.df[which(rownames(metrics.thres.df)==MetrixMax),],na.rm = TRUE)))))

  # print(paste0("Optimal Threshold: ",selectedThreshold))
  thresholdDT <- as.data.frame(t(metrics.thres.df))
  thresholdDT <- thresholdDT %>% rownames_to_column("threshold")
  thresholdDT.m <- melt(thresholdDT)
  thresholdDT.m <- thresholdDT.m %>% dplyr::filter(variable %in% c("Precision", "Recall"))
  thresholdDT.m$threshold <- as.numeric(thresholdDT.m$threshold)
  g <- ggplot(data = thresholdDT.m , aes(x=threshold, y=value, group=variable, color=variable))
  g <- g + geom_line(size=1) #aes(linetype=variable),
  g <- g + geom_vline(xintercept = as.numeric(selectedThreshold), color = "black",linetype="dashed")
  g <- g + geom_point()
  g <- g + ggtitle(modelName)
  g <- g + scale_color_manual(NULL, values = c("Indian red", "steelblue"))
  #g <- g + geom_vline(xintercept = as.numeric(selectedThreshold, color = "red3",linetype="dashed"))

  #g <- g + scale_x_continuous(breaks=c(0,0.25,0.5,1))
  g <- g + theme_classic(base_family = "Palatino Linotype", base_size=13)
  g <- g + theme(legend.position = c(0.99,0.02),
                 legend.justification = c(1,0))
  #                #plot.margin = margin(l = 150, r = 150))
  g <- g+ annotate("text", x=selectedThreshold - 0.15 , y=0.2, label = paste0("Optimal \nthreshold: ",round(selectedThreshold,3)))
  print(g)
}


threshold_moving <- function(modelList, dataSet, yVar, classes,LevelIndex, testX = NULL, testY=NULL,thr,modelName,optimalThres=NULL,extractedProbDF=NULL,rfeE=FALSE
                             # mixedMode = FALSE,
                             # trainModelsVec = NULL,
                             # rfeModelsVec = NULL
                             ){
  # modelList<- train.modelsSolo
  # dataSet <- testDF[[1]]
  # yVar <-"response"
  # classes<- c("CRPR","PD")# ClassLevels
  # LevelIndex <- 1
  # testX = data.testX
  # testY=testDF[[1]]$response
  # thr<-seq(0.1,0.99,length.out = 30)[1]
  # modelName <- names(train.modelsSolo)[1]
  # optimalThres=NULL
  
  
  # thresSeq<-seq(0.1,0.99,length.out = 50)
  # modelList<- models.final[[1]]
  # dataSet<-testDF[[1]]
  # testX=data.testX
  # testY=testDF[[1]]$response
  # modelName<- names(models.final)[1]
  # classes<- ClassLevels
  # LevelIndex<-1
  # MetrixMax="F1.Score"
  # extractedProbDF=NULL
  # rfeMode = FALSE
  # 
  # 
  # thresSeq<-seq(0.1,0.99,length.out = 50)
  # modelList<- modelList
  # dataSet<-testDF[[1]]
  # yVar<- "response"
  # testX=data.testX
  # testY=testDF[[1]][[yVar]]
  # modelName<- modelName
  # classes<- ClassLevels
  # LevelIndex<-1
  # MetrixMax="Youdens.J.statistic"
  # extractedProbDF=NULL
  # printP=FALSE
  # 
  # extractedProbDF=NULL
  # testX = testX
  # thr <- 0.1
  # optimalThres=NULL
  # extractedProbDF=NULL
  # rfeE=rfeMode
  #################
  # Extract probabilities, oneither train or test data
  if( !is.null(testX) && ! is.null(testY)){
    if(isTRUE(rfeE)){
      probabilities <- predict(modelList, newdata=dataSet, type="prob")
      predicted.classes <- predict(modelList, newdata=dataSet, type="raw")
      predClassesProbs <-data.frame(obs = factor(dataSet[[yVar]], levels = classes))
      prob_test <- cbind(predClassesProbs,probabilities)
      prob_test$object <- rep(modelName, dim(prob_test)[1])

    }else{
      prob_test<-extractProb(modelList,unkX =testX)
      prob_test$obs <-dataSet[[yVar]]
      prob_test$obs <-factor(prob_test$obs, levels = classes)
      prob_test$pred <-factor(prob_test$pred, levels = classes)
    }

  }else{
    prob_test<-extractProb(modelList)
  }

  
  
  # if(isTRUE(rfeE)){
  #   # b <-sapply(1: length(modelList), function(x) predict(modelList[[x]], newdata=test.data), simplify = FALSE)
  #   probabilities <- predict(modelList, newdata=dataSet, type="prob")
  #   # probabilities.df <- bind_rows(probabilities, .id = "object")
  #   # predicted.classes <- predict(modelList, newdata=test.data, type="raw")
  #   
  #   probabilities <- lapply(probabilities, function(x) transform(x, obs = factor(dataSet[[responseVar]], levels = classes)))
  #   probabilities <- lapply(probabilities, function(x) transform(x, pred = factor(pred, levels = classes)))
  #   predClassesProbs <-bind_rows(probabilities, .id = "object")
  #   prob_test <-predClassesProbs
  #   # predClassesProbs <- cbind(predClassesProbs,probabilities)
  #   
  # }else if(isTRUE(mixedMode)){
  #   ## First train objects
  #   if( !is.null(testX) && ! is.null(testY)){
  #     prob_test.t<-extractProb(modelList[trainModelsVec],unkX =testX)
  #     prob_test.t$obs <-dataSet[[responseVar]]
  #   }else{
  #     prob_test.t<-extractProb(modelList)    
  #   }
  #   # Bring those in the format of  train models
  #   prob_test.t <-prob_test.t %>% dplyr::select(all_of(classes),"obs", "pred","object")
  #   
  #   ## Then RFE objects
  #   probabilities <- predict(modelList[rfeModelsVec], newdata=test.data, type="prob")
  #   # probabilities.df <- bind_rows(probabilities, .id = "object")
  #   # predicted.classes <- predict(modelList, newdata=test.data, type="raw")
  #   
  #   probabilities <- lapply(probabilities, function(x) transform(x, obs = factor(dataSet[[responseVar]], levels = classes)))
  #   probabilities <- lapply(probabilities, function(x) transform(x, pred = factor(pred, levels = classes)))
  #   predClassesProbs <-bind_rows(probabilities, .id = "object")
  #   prob_test.r <-predClassesProbs
  #   # Bring those in the format of  train models
  #   prob_test.r <-prob_test.r %>% dplyr::select(all_of(classes),"obs", "pred","object")
  #   ## Merge objects
  #   prob_test <- rbind(prob_test.t, prob_test.r)
  #   
  # }else{
  #   if( !is.null(testX) && ! is.null(testY)){
  #     prob_test<-extractProb(modelList,unkX =testX)
  #     prob_test$obs <-dataSet[[responseVar]]
  #     prob_test$obs <-factor(prob_test$obs, levels = classes)
  #     prob_test$pred <-factor(prob_test$pred, levels = classes)
  #   }else{
  #     prob_test<-extractProb(modelList)    
  #   }
  #   
  # }
  
  
  # Select model
  prob_test.final <- prob_test %>% dplyr::filter(object %in% c(modelName))
  #thr <- 0.1
  # Change predicted based on threshold
  thres_preds <- function(predDF,thr,classes){
    predDF[["pred"]] <- ifelse(predDF[[classes[LevelIndex]]] >=thr,classes[LevelIndex],classes[-LevelIndex])
    predDF[["pred"]] <- factor(predDF[["pred"]],levels=levels(predDF[["obs"]]))
    predDF
  }
  prob_test.final.thr <- thres_preds(prob_test.final,thr,classes)
  
  
  ## Metrics
  # Precision-Recall
  pr <-prSummary(prob_test.final.thr, lev = classes)# HERE I can change levels
  pr.auc <- pr[1]
  precision <- pr[2]
  recall<- pr[3]
  f1.score <- pr[4]
  # ROC, sensitivity, specificity
  # Change levels-positive class
  
  if(isTRUE(all(classes==levels(prob_test.final.thr$obs)))){
    roc.sens.spec <- twoClassSummary(prob_test.final.thr, lev = classes)# HERE I CANNOT change levels    
    
  }else{
    prob_test.final.thr.c <-prob_test.final.thr
    prob_test.final.thr.c$obs <-factor(prob_test.final.thr.c$obs, levels = classes)
    prob_test.final.thr.c$pred <-factor(prob_test.final.thr.c$pred, levels = classes)
    roc.sens.spec <- twoClassSummary(prob_test.final.thr.c, lev = classes)
  }
  
  roc.auc <- roc.sens.spec[1]
  sens <- roc.sens.spec[2]
  spec <- roc.sens.spec[3]
  youden <- sens+spec-1
  # Accuracy, Balanced accuracy
  if(isTRUE(all(classes==levels(prob_test.final.thr$obs)))){
    acc.all <- caret::confusionMatrix(prob_test.final.thr$pred,prob_test.final.thr$obs,positive=classes[LevelIndex])
    
  }else{
    acc.all <- caret::confusionMatrix(prob_test.final.thr.c$pred,prob_test.final.thr.c$obs,positive=classes[LevelIndex])
  }
  
  
  ppv <- acc.all$byClass[[3]]
  npv <-acc.all$byClass[[4]]
  
  conf.table <- acc.all$table
  conf.table <- as.matrix(conf.table)
  
  rownames(conf.table) <- c(paste0("Predicted Positive (",classes[LevelIndex],")" ), paste0("Predicted Negative (",classes[-LevelIndex],")" ))
  colnames(conf.table) <- c(paste0("Actual Positive(",classes[LevelIndex],")" ), paste0("Actual Negative(",classes[-LevelIndex],")" ))
  
  acc <-  acc.all$overall[1]
  acc.bal <- acc.all$byClass[11]
  # print out confusion matrix
  #pander(conf.table)
  # results all in df
  res.model <- rbind(acc,acc.bal,sens,spec,roc.auc,precision,recall, pr.auc,f1.score,youden,ppv,npv)
  colnames(res.model)<- "value"
  rownames(res.model) <- c("Accuracy", "Balanced Accuracy", "Sensitivity", "Specificity", "ROC-AUC", "Precision", "Recall", "PR-AUC", "F1.Score","Youdens.J.statistic", "PPV","NPV")
  res.model
  
  
}


optimalScoreThreshold.logTMB <- function(modelList,dataSet,yVar,classes,LevelIndex,testX,testY, modelName,scorethres,var, printP=FALSE){
  # scorethres<-log(10)
  # modelList <-train.modelsSolo
  # dataSet<-testDF[[1]]
  # yVar<-"response"
  # classes<-c("CRPR","PD")
  # LevelIndex<-1
  # testX = data.testX
  # testY=testDF[[1]][["response"]]
  # modelName=names(train.modelsSolo)[1]
  # scorethres = models.optThres[1]
  # 
  # var=sigs[1]
  # extractedProbDF=NULL
  # Extract probabilities, oneither train or test data
  
  if( !is.null(testX) && ! is.null(testY)){
    prob_test<-extractProb(modelList,unkX =testX)
    prob_test$obs <-dataSet[[yVar]]
    prob_test$score <-dataSet[[var]]
  }else{
    prob_test<-extractProb(modelList)    
  }

  
  
  
  # Select model
  prob_test.final <- prob_test %>% dplyr::filter(object %in% c(modelName))
  # head(prob_test.final)
  # ## Find intersection
  AF1 = approxfun(prob_test.final$score, prob_test.final$CRPR, rule=2)
  classthress<-AF1(as.numeric(scorethres))

  # line1 <- data.frame(x = rep(scorethres,length(prob_test.final$score)), y = prob_test.final$CRPR )
  # curve2 <-data.frame(x = prob_test.final$score, y = prob_test.final$CRPR)
  #
  # intersect.point<-curve_intersect(curve2,line1)
  #

  g <- ggplot(data = prob_test.final , aes(x=score, y=CRPR))
  g <- g + geom_line(size=1) #aes(linetype=variable),
  g <- g + geom_vline(xintercept = as.numeric(scorethres), color = "red",linetype="dashed")
  g <- g + geom_hline(yintercept = classthress, color = "blue",linetype="dashed")
  # g <- g + geom_point()
  g <- g + ggtitle(modelName)

  g <- g + theme_classic(base_family = "Palatino Linotype", base_size=13)

  if(classthress >0){
    g <- g+ annotate("text", y=classthress - (classthress/10) , x=as.numeric(scorethres)-(as.numeric(scorethres)*2), label = paste0("Optimal \n Classification threshold: ",round(classthress,3)))
  }else{
    g <- g+ annotate("text", y=classthress + (classthress/10) , x=as.numeric(scorethres)-(as.numeric(scorethres)*2), label = paste0("Optimal \n Classification threshold: ",round(classthress,3)))
  }
  if(isTRUE(printP)){
    print(g)
  }


  classthress
}


optimal_threshold_value <- function(thresSeq,modelList,dataSet,yVar, testX=NULL,testY=NULL, modelName ,classes,LevelIndex, MetrixMax="F1.Score",extractedProbDF=NULL, printP=FALSE,rfeMode = FALSE
                                    # mixedMode = FALSE,
                                    # trainModelsVec = NULL,
                                    # rfeModelsVec = NULL
                                    ){
  
  # thresSeq<-seq(0.1,0.99,length.out = 50)
  # modelList<- train.modelsSolo
  # dataSet<-testDF[[1]]
  # yVar<- "response"
  # testX=data.testX
  # testY=testDF[[1]]$response
  # modelName<- names(train.modelsSolo)[1]
  # classes<- ClassLevels
  # LevelIndex<-1
  # MetrixMax="Youdens.J.statistic"
  # extractedProbDF=NULL
  # printP=FALSE
  
  # thresSeq<-seq(0.1,0.99,length.out = 50)
  # modelList<- ModelsList[[ModelIndex]]
  # dataSet<-testDF
  # yVar<- "response"
  # testX=data.testX
  # testY=testDF[[yVar]]
  # modelName<- names(train.rfe)[1]
  # classes<- ClassLevels
  # LevelIndex<-1
  # MetrixMax="Youdens.J.statistic"
  # extractedProbDF=NULL
  # printP=FALSE
  # 

  # Calculate Evaluation metrics for all thresholds
  metrics.thres <-lapply(thresSeq, function(x) threshold_moving(modelList,dataSet, yVar, classes, LevelIndex,testX = testX, testY=testY,x,modelName,extractedProbDF=NULL,rfeE=rfeMode
                                                                # mixedMode = mixedMode,
                                                                # trainModelsVec = trainModelsVec,
                                                                # rfeModelsVec = rfeModelsVec
                                                                ))
  
  names(metrics.thres) <- as.character(thresSeq)
  metrics.thres.df<-list.cbind(metrics.thres)
  colnames(metrics.thres.df) <-as.character(thresSeq)
  # Find classification Threshold that maximizes F-1 Score
  selectedThreshold <- max(as.numeric(names(which(metrics.thres.df[which(rownames(metrics.thres.df)==MetrixMax),]==max(metrics.thres.df[which(rownames(metrics.thres.df)==MetrixMax),],na.rm = TRUE)))))
  
  # print(paste0("Optimal Threshold: ",selectedThreshold))
  thresholdDT <- as.data.frame(t(metrics.thres.df))
  thresholdDT <- thresholdDT %>% rownames_to_column("threshold")
  thresholdDT$threshold <- as.numeric(thresholdDT$threshold)
  # thresholdDT.m <- melt(thresholdDT)
  # thresholdDT.m <- thresholdDT.m %>% dplyr::filter(variable %in% c("Precision", "Recall"))
  # thresholdDT.m$threshold <- as.numeric(thresholdDT.m$threshold)
  
  #selectedThreshold
  
  if (isTRUE(printP)){
    # Plot Threshold with Accuracy
    g <- ggplot(data = thresholdDT , aes(x=threshold, y=Accuracy))
    g <- g + geom_line(size=1) #aes(linetype=variable), 
    g <- g + geom_vline(xintercept = as.numeric(selectedThreshold), color = "black",linetype="dashed")
    g <- g + geom_point()
    g <- g + theme_classic(base_family = "Palatino Linotype", base_size=13)
    g <- g + scale_x_continuous(breaks=c(0,0.25,0.5,1))
    # g <- g + theme(legend.position = c(0.99,0.02),
    #                legend.justification = c(1,0))
    # #                #plot.margin = margin(l = 150, r = 150))
    g <- g+ annotate("text", x=selectedThreshold -0.2 , y=0.4, label = paste0("Optimal \nthreshold: ",round(selectedThreshold,3)))
    #print(g)
    ## Plot Threshold with F1-Score
    # Plot Threshold with Accuracy
    g.f1 <- ggplot(data = thresholdDT , aes_string(x="threshold", y="Youdens.J.statistic"))
    g.f1 <- g.f1 + geom_line(size=1) #aes(linetype=variable), 
    g.f1 <- g.f1 + geom_vline(xintercept = as.numeric(selectedThreshold), color = "black",linetype="dashed")
    g.f1 <- g.f1 + geom_point()
    g.f1 <- g.f1 + theme_classic(base_family = "Palatino Linotype", base_size=13)
    g.f1 <- g.f1 + scale_x_continuous(breaks=c(0,0.25,0.5,1))
    # g.f1 <- g.f1 + theme(leg.f1end.position = c(0.99,0.02),
    #                leg.f1end.justification = c(1,0))
    # #                #plot.marg.f1in = marg.f1in(l = 150, r = 150))
    g.f1 <- g.f1+ annotate("text", x=selectedThreshold - 0.2 , y=0.25, label = paste0("Optimal \nthreshold: ",round(selectedThreshold,3)))
    #print(g.f1)
    # print()
    fig<-ggpubr::ggarrange(g,g.f1, ncol = 2, nrow = 1)
    print(annotate_figure(fig,
                          top = text_grob(paste0(MetrixMax,": ",round(max(metrics.thres.df[which(rownames(metrics.thres.df)==MetrixMax),],na.rm = TRUE),3)), color = "black", face = "bold", size = 12)))
  }
  
  # print(paste0(MetrixMax,": ",round(max(metrics.thres.df[which(rownames(metrics.thres.df)==MetrixMax),],na.rm = TRUE),3)))
  selectedThreshold
}

optimalScoreThreshold <- function(modelList,dataSet,yVar,classes,LevelIndex,testX,testY, modelName,thres,var,extractedProbDF=NULL, printP=FALSE){
  # modelList <-train.modelsSolo
  # dataSet<-testDF[[1]]
  # yVar<-"response"
  # classes<-ClassLevels
  # LevelIndex<-1
  # testX = data.testX
  # testY=testDF[[1]][["response"]]
  # modelName=names(train.modelsSolo)[1]
  # thres = 0.5
  # var=sigs[1]
  # extractedProbDF=NULL
  # Extract probabilities, oneither train or test data
  
  if( !is.null(testX) && ! is.null(testY)){
    prob_test<-extractProb(modelList,unkX =testX)
    prob_test$obs <-dataSet[[yVar]]
    prob_test$score <-dataSet[[var]]
  }else{
    prob_test<-extractProb(modelList)    
  }
  
  ### FIX FOR INSANE PRED THRES
  if(thres==0.99){
    thres<-0.5
  }
  
  # AF1 = approxfun(prob_test.final$score, prob_test.final$CRPR)
  # AF1(as.numeric(thres))
  # Select model
  prob_test.final <- prob_test %>% dplyr::filter(object %in% c(modelName))
  
  ## Find intersection
  
  line1 <- data.frame(x = prob_test.final$score, y = rep(as.numeric(thres),length(prob_test.final$score)))
  curve2 <-data.frame(x = prob_test.final$score, y = prob_test.final$CRPR)
  
  
  # threshold_or1.fix1 <- map2_df(
  #   recall_or1_4, precision_or1_4,
  #   ~tryCatch({
  #     curve_intersect(.x, .y, empirical = TRUE, domain = NULL)
  #   }, error = function(e){
  #     return(tibble(.rows = 1))
  #   }),
  #   .id = "i"
  # )
  
  intersect.point<-curve_intersect(line1, curve2)
  
  
  g <- ggplot(data = prob_test.final , aes(x=score, y=CRPR))
  g <- g + geom_line(size=1) #aes(linetype=variable),
  g <- g + geom_hline(yintercept = as.numeric(thres), color = "red",linetype="dashed")
  g <- g + geom_vline(xintercept = intersect.point$x, color = "blue",linetype="dashed")
  # g <- g + geom_point()
  g <- g + ggtitle(modelName)
  
  g <- g + theme_classic(base_family = "Palatino Linotype", base_size=13)
  
  if(intersect.point$x >0){
    g <- g+ annotate("text", x=intersect.point$x - (intersect.point$x*2) , y=as.numeric(thres)-(as.numeric(thres)/10), label = paste0("Optimal \n Score threshold: ",round(intersect.point$x,3)))
  }else{
    g <- g+ annotate("text", x=intersect.point$x + (intersect.point$x*2) , y=as.numeric(thres)-(as.numeric(thres)/10), label = paste0("Optimal \n Score threshold: ",round(intersect.point$x,3)))
  }
  if(isTRUE(printP)){
    print(g)
  }
  
  
  intersect.point
}


scores_barplot <- function(modelList,dataSet,groupVar,xVar,yVar,labelVar,xLabel,yLabel, fillLabel, shapeLabel, title,
                           testX,testY,modelName,classes,LevelIndex,posClassName,
                           decr=FALSE,reverse=FALSE,thres=NULL,ScoreThres = NULL,scoreLabelPos="right",extractedProbDF=NULL,rfeE=FALSE){
  # modelList<-modelList.sub
  # dataSet<-testDF[[1]]
  # groupVar<-"response"
  # xVar <- "run_accession"
  # yVar <- sigs[1]
  # labelVar <- "tissue"
  # xLabel<-"Patients"
  # yLabel <- sigs_axisNames[1]
  # fillLabel<-"Response to anti-PD-1"
  # shapeLabel<-"tissue"
  # title<-sigs_titleNames[1]
  # testX <- data_testX
  # testY <- testDF[[1]]$response
  # modelName <- names(modelList.sub[1])
  # classes <- c("CRPR","PD")
  # LevelIndex <- 1
  # thres=NULL
  # posClassName <- "CRPR"
  # decr <- FALSE
  # ScoreThres = models.optScoreThres[1,1]$x
  # scoreLabelPos=scorelabelspos[1]
  # reverse = TRUE
  ####
  xVar.un <- as.name(xVar)
  yVar.un <- as.name(yVar)
  labelVar.un <- as.name(labelVar)
  posClassName.un <- as.name(posClassName)
  # Extract probabilities, oneither train or test data
  if(isTRUE(rfeE)){
    # b <-sapply(1: length(modelList), function(x) predict(modelList[[x]], newdata=test.data), simplify = FALSE)
    probabilities <- predict(modelList, newdata=dataSet, type="prob")
    # probabilities.df <- bind_rows(probabilities, .id = "object")
    # predicted.classes <- predict(modelList, newdata=test.data, type="raw")
    
    probabilities <- lapply(probabilities, function(x) transform(x, obs = factor(dataSet[[groupVar]], levels = classes)))
    probabilities <- lapply(probabilities, function(x) transform(x, pred = factor(pred, levels = classes)))
    predClassesProbs <-bind_rows(probabilities, .id = "object")
    prob_test <-predClassesProbs
    # predClassesProbs <- cbind(predClassesProbs,probabilities)
    
  }else{
    if(!is.null(extractedProbDF)){
      prob_test <- extractedProbDF
      prob_test[[yVar]] <- dataSet[[yVar]]
      prob_test[[xVar]] <- dataSet[[xVar]]
      prob_test[[labelVar]] <- dataSet[[labelVar]]
    }else{
      if( !is.null(testX) && ! is.null(testY)){
        prob_test<-extractProb(modelList,unkX =testX)
        prob_test$obs <-dataSet[[groupVar]]
        prob_test[[yVar]] <- dataSet[[yVar]]
        prob_test[[xVar]] <- dataSet[[xVar]]
        prob_test[[labelVar]] <- dataSet[[labelVar]]
      }else{
        prob_test<-extractProb(modelList)    
      }
    }
  }
  
  
  if( !is.null(thres) ){
    # Change predicted based on threshold
    thres_preds <- function(predDF,thr,classes){
      predDF[["pred"]] <- ifelse(predDF[[classes[LevelIndex]]] >=thr,classes[LevelIndex],classes[-LevelIndex])
      predDF[["pred"]] <- factor(predDF[["pred"]],levels=levels(predDF[["obs"]]))
      predDF
    }
    prob_test <- thres_preds(prob_test,thres,classes)
  }
  
  
  # Select model
  prob_test.final <- prob_test %>% dplyr::filter(object %in% c(modelName))
  
  #pdf("g.pdf", width = 12, height = 6)
  prob_test.final.ord <- prob_test.final[order(prob_test.final[[yVar]],decreasing = decr),]
  prob_test.final.ord[[xVar]] <- factor(prob_test.final.ord[[xVar]])
  
  mu_PD <- prob_test.final.ord %>% dplyr::filter(obs==classes[LevelIndex]) %>% dplyr::select(yVar.un)
  mu_PD <- mean(mu_PD[,1])
  # %>% summarize(mean(TuTACK_sigscore)) %>% pull()
  
  mu_score <- prob_test.final.ord %>% dplyr::select(yVar.un) 
  # %>% summarize(mean(TuTACK_sigscore)) %>% pull()
  mu_score <- mean(mu_score[,1])
  
  g <- ggplot(data = prob_test.final.ord)
  g <- g + geom_bar(aes(x =reorder(!!xVar.un,!!yVar.un), y = !!yVar.un, fill = obs), stat = "identity",
                    width = 1)
  #g <- g +scale_fill_brewer(palette="Dark2")
  if (isTRUE(reverse)){
    g <- g + scale_fill_manual(values=rev(c("#CD5C5C", "#004D40")))
  }else{
    g <- g + scale_fill_manual(values=c("#CD5C5C", "#004D40"))
  }
  
  g <- g + scale_x_discrete(limits = prob_test.final.ord[[xVar]])
  # g <- g + ylim(min(prob_test.final.ord[[yVar]]),max(prob_test.final.ord[[yVar]])+ 0.4)  + ylab(yLabel) + xlab(xLabel)
  g <- g + ylab(yLabel) + xlab(xLabel)
  g <- g + ggtitle(title)
  g <- g + theme_ipsum(base_family = "Palatino Linotype", base_size=13)
  g <- g +theme(plot.title = element_text(size = 16, face = "bold"),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.x = element_text(size=16, face="bold"),
                axis.text.y = element_text( hjust = 1, vjust = 1,size=13),
                axis.title.y = element_text(size=16, face="bold"),
                legend.title = element_text(size=16,color="black",face="bold"),
                #legend.background = element_rect(size=0.1,color="black",fill=alpha("white",0.6)),
                legend.text = element_text(size=13),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                legend.position="bottom",legend.box="vertical", legend.margin=margin())
  
  g <- g + geom_hline(yintercept=mu_PD, linetype="dotdash", color = "Black", size=0.8)
  if(!is.null(ScoreThres)){
    g <- g + geom_hline(yintercept=ScoreThres, linetype="dashed", color = "#4682B4", size=0.8)
  }else{
    g <- g + geom_hline(yintercept=mu_score, linetype="dashed", color = "#4682B4", size=0.8)
    
  }
  
  g <- g + geom_point(data = prob_test.final.ord, aes(x = reorder(!!xVar.un,!!yVar.un), y = ifelse(prob_test.final.ord[[yVar]] >=0 , prob_test.final.ord[[yVar]] + 0.3,prob_test.final.ord[[yVar]] - 0.3 ), color = !!labelVar.un), shape=18,size=2.2)
  g <- g +scale_color_manual(values=c("#80BD7D", "#56B4E9", "#CC79A7","#F0E442","Dark Slate Gray"))#, "#F0E442","Dark Slate Gray"
  #g <- g +scale_color_manual(values=c("#E69F00","Blue Violet","#999999", "#56B4E9"))
  #g <- g +scale_color_manual(values=RColorBrewer::brewer.pal(4,"Set2"))
  if (scoreLabelPos=="left"){
    xlen <- 11
    # xlen <- dim(dataSet)[1]-4
  }else{
    xlen <- dim(dataSet)[1]-10
  }
  
  if(!is.null(ScoreThres)){
    
    if (abs(mu_PD) > abs(ScoreThres) && mu_PD>0){
      # print("mean larger,positive")
      # print(mu_PD)
      # print(ScoreThres)
      g <- g + 
        annotate("text", y = mu_PD+0.15, x=xlen, label =bquote(mu[.(posClassName)]) , color = 'black',
                 size = 5)
      g <- g + 
        annotate("text", y = ScoreThres-0.15, x=xlen, label =expression(paste(score[optimal], "")) , color = "#4682B4",
                 size = 5)
    }else if(abs(mu_PD) > abs(ScoreThres) && mu_PD<0){
      # print("mean larger,negative")
      # print(mu_PD)
      # print(ScoreThres)
      g <- g + 
        annotate("text", y = mu_PD-0.15, x=xlen, label =bquote(mu[.(posClassName)]) , color = 'black',
                 size = 5)
      g <- g + 
        annotate("text", y = ScoreThres+0.15, x=xlen, label =expression(paste(score[optimal], "")) , color = "#4682B4",
                 size = 5)
    }else if(abs(mu_PD) < abs(ScoreThres) && mu_PD>0){
      # print("mean smaller,positive")
      # print(mu_PD)
      # print(ScoreThres)
      g <- g + 
        annotate("text", y = mu_PD-0.15, x=xlen, label =bquote(mu[.(posClassName)]) , color = 'black',
                 size = 5)
      g <- g + 
        annotate("text", y = ScoreThres+0.15, x=xlen, label =expression(paste(score[optimal], "")) , color = "#4682B4",
                 size = 5)
    }else{
      # print("mean smaller,negative")
      # print(mu_PD)
      # print(ScoreThres)
      g <- g + 
        annotate("text", y = mu_PD+0.15, x=xlen, label =bquote(mu[.(posClassName)]) , color = 'black',
                 size = 5)
      g <- g + 
        annotate("text", y = ScoreThres-0.15, x=xlen, label =expression(paste(score[optimal], "")) , color = "#4682B4",
                 size = 5)
    }
    
  }else{
    if (mu_PD > mu_score){
      g <- g + 
        annotate("text", y = mu_PD+0.15, x=xlen, label =bquote(mu[.(posClassName)]) , color = 'black',
                 size = 6)
      
      g <- g + 
        annotate("text", y = mu_score-0.15, x=xlen, label =expression(paste(mu[score], "")) , color = "#4682B4",
                 size = 6)
      
    }else{
      g <- g + 
        annotate("text", y = mu_PD-0.15, x=xlen, label =bquote(mu[.(posClassName)]) , color = 'black',
                 size = 6)
      
      g <- g + 
        annotate("text", y = mu_score+0.15, x=xlen, label =expression(paste(mu[score], "")) , color = "#4682B4",
                 size = 6)
      
    }
  }
  
  
  #g <- g + scale_color_brewer(palette="Dark2")
  # + ylab("TuTack-score") + xlab("Patients")
  g <- g + labs(x = xLabel, y=yLabel, fill= fillLabel, shape=shapeLabel)
  g <- g + guides(color = guide_legend(override.aes = list(size=11)))
  print(g)   
}

classProb_barplot <- function(modelList,dataSet,groupVar,xVar,sigVar,yVar,labelVar,xLabel,yLabel, fillLabel, shapeLabel, title,testX,testY,modelName,classes,LevelIndex,
                              decr=FALSE,reverse=FALSE,thres=NULL,extractedProbDF=NULL,rfeE=FALSE,
                              mixedMode = FALSE,
                              trainModelsVec = NULL,
                              rfeModelsVec = NULL){
  
  # modelList<-modelList.sub
  # dataSet<-testDF[[1]]
  # groupVar<-"response"
  # xVar <- "run_accession"
  # sigVar<-sigs[1]
  # yVar <- "CRPR"
  # labelVar <- "tissue"
  # xLabel<-"Patients"
  # yLabel <- sigs_axisNames[1]
  # fillLabel<-"Response to anti-PD-1"
  # shapeLabel<-"tissue"
  # title<-sigs_titleNames[1]
  # testX <- data.testX
  # testY <- testDF[[1]]$response
  # modelName <- names(modelList.sub)[1]
  # classes <- c("CRPR","PD")
  # LevelIndex <- 1
  # thres=NULL
  # posClassName <- "CRPR"
  # decr <-decr_listPlot
  # reverse=TRUE
  # rfeE=FALSE

  # modelList<-modelList.sub.rfe
  # dataSet<-testDF[[1]]
  # groupVar<-"response"
  # xVar <- "run_accession"
  # sigVar<-sigs[1]
  # yVar <- "CRPR"
  # labelVar <- "tissue"
  # xLabel<-"Patients"
  # yLabel <- sigs_axisNames[1]
  # fillLabel<-"Response to anti-PD-1"
  # shapeLabel<-"tissue"
  # title<-sigs_titleNames[1]
  # testX <- data_testX
  # testY <- testDF[[1]]$response
  # modelName <- names(modelList.sub.rfe)[1]
  # classes <- c("CRPR","PD")
  # LevelIndex <- 1
  # thres=NULL
  # posClassName <- "CRPR"
  # decr <-decr_listPlot
  # reverse=reverseParam
  # rfeE=TRUE
  
  # modelList<-modelList.final
  # dataSet<-testDF[[1]]
  # groupVar<-"gender"
  # xVar <- "run_accession"
  # sigVar<-sigs[3]
  # yVar <- "CRPR"
  # labelVar <- "gender"
  # xLabel<-"Patients"
  # yLabel <- sigs_axisNames[3]
  # fillLabel<-"Response to anti-PD-1"
  # shapeLabel<-"tissue"
  # title<-sigs_titleNames[3]
  # testX <- data.testX
  # testY <- testDF[[1]]$gender
  # modelName <- names(modelList.final)[3]
  # classes <- c("CRPR","PD")
  # LevelIndex <- 1
  # thres=models.optThres.rfe[[3]]
  # posClassName <- "CRPR"
  # decr <-decr_listPlot
  # reverse=TRUE
  # rfeE=FALSE
  # extractedProbDF=NULL
  # mixedMode = TRUE
  # trainModelsVec = c(1:7)
  # rfeModelsVec = c(8)
  ##############################
  xVar.un <- as.name(xVar)
  yVar.un <- as.name(yVar)
  sigVar.un<- as.name(sigVar)
  labelVar.un <- as.name(labelVar)
  
  # # Extract probabilities, oneither train or test data
  # if(isTRUE(rfeE)){
  #   # b <-sapply(1: length(modelList), function(x) predict(modelList[[x]], newdata=test.data), simplify = FALSE)
  #   probabilities <- predict(modelList, newdata=dataSet, type="prob")
  #   # probabilities.df <- bind_rows(probabilities, .id = "object")
  #   # predicted.classes <- predict(modelList, newdata=test.data, type="raw")
  #   
  #   probabilities <- lapply(probabilities, function(x) transform(x, obs = factor(dataSet[[groupVar]], levels = classes)))
  #   probabilities <- lapply(probabilities, function(x) transform(x, run_accession = factor(dataSet[[xVar]])))
  #   probabilities <- lapply(probabilities, function(x) transform(x, tissue = as.character(dataSet[[labelVar]])))
  #   probabilities <- lapply(probabilities, function(x) transform(x, pred = factor(pred, levels = classes)))
  #   predClassesProbs <-bind_rows(probabilities, .id = "object")
  #   prob_test <-predClassesProbs
  #   # predClassesProbs <- cbind(predClassesProbs,probabilities)
  #   
  # }else{
  #   if(!is.null(extractedProbDF)){
  #     prob_test <- extractedProbDF
  #     prob_test[[groupVar]] <- dataSet[[groupVar]]
  #     prob_test[[xVar]] <- dataSet[[xVar]]
  #     prob_test[[labelVar]] <- dataSet[[labelVar]]
  #   }else{
  #     if( !is.null(testX) && ! is.null(testY)){
  #       prob_test<-extractProb(modelList,unkX =testX)
  #       prob_test$obs <-dataSet[[groupVar]]
  #       prob_test[[groupVar]] <- dataSet[[groupVar]]
  #       prob_test[[xVar]] <- dataSet[[xVar]]
  #       prob_test[[labelVar]] <- dataSet[[labelVar]]
  #     }else{
  #       prob_test<-extractProb(modelList)    
  #     }
  #   }
  # }
  # 
  #####################
  if(isTRUE(rfeE)){
    # b <-sapply(1: length(modelList), function(x) predict(modelList[[x]], newdata=test.data), simplify = FALSE)
    probabilities <- predict(modelList, newdata=dataSet, type="prob")
    # probabilities.df <- bind_rows(probabilities, .id = "object")
    # predicted.classes <- predict(modelList, newdata=test.data, type="raw")
    
    probabilities <- lapply(probabilities, function(x) transform(x, obs = factor(dataSet[[groupVar]], levels = classes)))
    probabilities <- lapply(probabilities, function(x) transform(x, pred = factor(pred, levels = classes)))
    probabilities <- lapply(probabilities, function(x) transform(x, run_accession = factor(dataSet[[xVar]])))
    probabilities <- lapply(probabilities, function(x) transform(x, tissue = as.character(dataSet[[labelVar]])))
    predClassesProbs <-bind_rows(probabilities, .id = "object")
    prob_test <-predClassesProbs
    # predClassesProbs <- cbind(predClassesProbs,probabilities)
    
  }else if(isTRUE(mixedMode)){
    ## First train objects
    if( !is.null(testX) && ! is.null(testY)){
      prob_test.t<-extractProb(modelList[trainModelsVec],unkX =testX)
      prob_test.t$obs <-dataSet[[groupVar]]
      prob_test.t[[xVar]] <-dataSet[[xVar]]
      prob_test.t[[labelVar]] <-dataSet[[labelVar]]
    }else{
      prob_test.t<-extractProb(modelList)    
    }
    # Bring those in the format of  train models
    prob_test.t <-prob_test.t %>% dplyr::select(all_of(classes),"obs", "pred","object",xVar,labelVar)
    
    ## Then RFE objects
    probabilities <- predict(modelList[rfeModelsVec], newdata=dataSet, type="prob")
    # probabilities.df <- bind_rows(probabilities, .id = "object")
    # predicted.classes <- predict(modelList, newdata=test.data, type="raw")
    
    probabilities <- lapply(probabilities, function(x) transform(x, obs = factor(dataSet[[groupVar]], levels = classes)))
    probabilities <- lapply(probabilities, function(x) transform(x, pred = factor(pred, levels = classes)))
    probabilities <- lapply(probabilities, function(x) transform(x, run_accession = factor(dataSet[[xVar]])))
    probabilities <- lapply(probabilities, function(x) transform(x, labelVar.un= as.character(dataSet[[labelVar]])))
    
    # for(i in probabilities){
    #   i=probabilities[[1]]
    #   colnames(i)[which(colnames(i)=="labelVar.un")] <- labelVar
    # }
    # 
    changeColName<- function(DF,newName){
      colnames(DF)[which(colnames(DF)=="labelVar.un")] <- newName
      DF
    }
    keepN <- names(probabilities)
    probabilities <- sapply(1:length(probabilities), function(x) changeColName(probabilities[[x]],labelVar), simplify = FALSE)
    names(probabilities) <-keepN
    
    
    predClassesProbs <-bind_rows(probabilities, .id = "object")
    prob_test.r <-predClassesProbs
    # Bring those in the format of  train models
    prob_test.r <-prob_test.r %>% dplyr::select(all_of(classes),"obs", "pred","object",xVar,labelVar)
    ## Merge objects
    prob_test <- rbind(prob_test.t, prob_test.r)
    
  }else{
    if( !is.null(testX) && ! is.null(testY)){
      prob_test<-extractProb(modelList,unkX =testX)
      prob_test$obs <-dataSet[[groupVar]]
      prob_test[[groupVar]] <- dataSet[[groupVar]]
      prob_test[[xVar]] <- dataSet[[xVar]]
      prob_test[[labelVar]] <- dataSet[[labelVar]]
    }else{
      prob_test<-extractProb(modelList)    
    }
    
  }
  #####################
  prob_test[["obs"]] <- factor(prob_test[["obs"]],levels = classes)
  
  
  
  if( !is.null(thres) ){
    # Change predicted based on threshold
    thres_preds <- function(predDF,thr,classes){
      # predDF <-prob_test
      # thr <-thres
      # classes <-classes
      
      predDF[["pred"]] <- ifelse(predDF[[classes[LevelIndex]]] >=thr,classes[LevelIndex],classes[-LevelIndex])
      predDF[["pred"]] <- factor(predDF[["pred"]],levels=levels(predDF[["obs"]]))
      predDF
    }
    prob_test <- thres_preds(prob_test,thres,classes)
  }
  
  
  # Select model
  prob_test.final <- prob_test %>% dplyr::filter(object %in% c(modelName)) %>% droplevels()
  
  #pdf("g.pdf", width = 12, height = 6)
  prob_test.final.ord <- prob_test.final[order(prob_test.final[[yVar]],decreasing = decr),]
  prob_test.final.ord[[xVar]] <- factor(prob_test.final.ord[[xVar]])
  
  # mu_PD <- prob_test.final.ord %>% dplyr::filter(obs==classes[LevelIndex]) %>% dplyr::select(yVar.un)
  # mu_PD <- mean(mu_PD[,1])
  # # %>% summarize(mean(TuTACK_sigscore)) %>% pull()
  
  # mu_score <- prob_test.final.ord %>% dplyr::select(yVar.un) 
  # # %>% summarize(mean(TuTACK_sigscore)) %>% pull()
  # mu_score <- mean(mu_score[,1])
  # 

  
  g <- ggplot(data = prob_test.final.ord)
  if(groupVar=="gender"){
    g <- g + geom_bar(aes(x =reorder(!!xVar.un,!!yVar.un), y = !!yVar.un, fill = gender), stat = "identity",
                      width = 1)
  }else{
    g <- g + geom_bar(aes(x =reorder(!!xVar.un,!!yVar.un), y = !!yVar.un, fill = obs), stat = "identity",
                      width = 1)
  }
  
  #g <- g +scale_fill_brewer(palette="Dark2")
  if (isTRUE(reverse)){
    
    if(groupVar=="gender"){
      g <- g + scale_fill_manual(values=rev(c("#C6B6D0", "#1E88E5","#AC1833")))
      
    }else{
      g <- g + scale_fill_manual(values=rev(c("#CD5C5C", "#004D40")))
    }
    
    
  }else{
    
    if(groupVar=="gender"){
      g <- g + scale_fill_manual(values=c("#C6B6D0", "#1E88E5","#AC1833"))
      
    }else{
      g <- g + scale_fill_manual(values=c("#CD5C5C", "#004D40"))
    }
    
    
    
  }
  
  g <- g + scale_x_discrete(limits = prob_test.final.ord[[xVar]])
  g <- g + scale_y_continuous(breaks = seq(0, 1.1, by = 0.25))
  # g <- g + ylim(min(prob_test.final.ord[[yVar]]),max(prob_test.final.ord[[yVar]])+ 0.4)  + ylab(yLabel) + xlab(xLabel)
  g <- g + ylab(yLabel) + xlab(xLabel) #+ ylim(0,1.1)
  g <- g + ggtitle(title)
  g <- g + theme_classic(base_family = "Palatino Linotype", base_size=13)
  g <- g +theme(plot.title = element_text(size = 16, face = "bold"),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.x = element_text(size=16, face="bold"),
                axis.text.y = element_text( size=13),#hjust = 1, vjust = 1,
                axis.title.y = element_text(size=8, face="bold"),
                legend.title = element_text(size=16,color="black",face="bold"),
                #legend.background = element_rect(size=0.1,color="black",fill=alpha("white",0.6)),
                legend.text = element_text(size=13),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill = "white"),
                axis.line = element_line(colour = "black"),
                legend.position="bottom",
                legend.box="vertical", 
                legend.margin=margin())
  
  
  # g <- g + geom_hline(yintercept=mu_score, linetype="dashed", color = "#4682B4", size=0.8)
  
  g <- g + geom_point(data = prob_test.final.ord, 
                      aes(x = reorder(!!xVar.un,!!yVar.un), 
                          y = ifelse(prob_test.final.ord[[yVar]] >=0 , prob_test.final.ord[[yVar]] + 0.03,prob_test.final.ord[[yVar]] - 0.03 ), 
                          color = !!labelVar.un), 
                          shape=18,
                          size=2.2)
  
  # g <- g + geom_point(data = prob_test.final.ord, 
  #                     aes(x = reorder(!!xVar.un,!!yVar.un), 
  #                         y = ifelse(prob_test.final.ord[[yVar]] >=0 , prob_test.final.ord[[yVar]] + 0.03,prob_test.final.ord[[yVar]] - 0.03 ), 
  #                         color = !!labelVar.un), 
  #                     shape=18,
  #                     size=2.2)
  
  g <- g +scale_color_manual(values=c("#80BD7D", "#56B4E9", "#CC79A7", "#F0E442","Dark Slate Gray"))
  
  #g <- g + scale_color_brewer(palette="Dark2")
  # + ylab("TuTack-score") + xlab("Patients")
  g <- g + labs(x = xLabel, y=paste0("Predicted Probability \n by ",yLabel), fill= fillLabel, shape=shapeLabel)
  g <- g + guides(color = guide_legend(override.aes = list(size=7)))
  g <- g + geom_hline(yintercept=thres, linetype="dashed", color = "black")#, size=0.8
  g   
}

tilePlot_compare_signatures<- function(modelList,dataSet,yVar,classes,LevelIndex,sigNames,testX=NULL, testY=NULL, reverse=FALSE,thres=NULL,extractedProbDF=NULL,rotate=FALSE,rfeE=FALSE, mode=NULL,
                                       mixedMode = FALSE,
                                       trainModelsVec = NULL,
                                       rfeModelsVec = NULL){
  # modelList <-train.modelsSolo
  # dataSet <-testDF[[1]]
  # yVar <-"response"
  # classes <-c("CRPR","PD")
  # LevelIndex <- 1
  # sigNames <-sigs_titleNames
  # testX=data.testX
  # testY=testDF[[1]]$response
  # thres=models.optThres
  # rfeE=FALSE
  # mode='no'
  
  # modelList <-models.final
  # dataSet <-testDF[[1]]
  # yVar <-"response"
  # classes <-c("CRPR","PD")
  # LevelIndex <- 1
  # sigNames <-sigs_titleNames
  # testX=data.testX
  # testY=testDF[[1]]$response
  # # thres=models.optThres.final
  # thres=NULL
  # 
  # reverse=TRUE
  # rotate=FALSE
  # rfeE=FALSE
  # mixedMode = TRUE
  # trainModelsVec = c(1:6)
  # rfeModelsVec = c(7)
  
  names(modelList)<- sigNames

  #############################################
  if(isTRUE(rfeE)){
    # b <-sapply(1: length(modelList), function(x) predict(modelList[[x]], newdata=test.data), simplify = FALSE)
    probabilities <- predict(modelList, newdata=dataSet, type="prob")
    # probabilities.df <- bind_rows(probabilities, .id = "object")
    # predicted.classes <- predict(modelList, newdata=test.data, type="raw")
    
    probabilities <- lapply(probabilities, function(x) transform(x, obs = factor(dataSet[[yVar]], levels = classes)))
    probabilities <- lapply(probabilities, function(x) transform(x, pred = factor(pred, levels = classes)))
    probabilities <- lapply(probabilities, function(x) transform(x, patient = factor(dataSet[["run_accession"]])))
    # probabilities <- lapply(probabilities, function(x) transform(x, tissue = as.character(dataSet[[labelVar]])))
    predClassesProbs <-bind_rows(probabilities, .id = "object")
    prob_test <-predClassesProbs
    # predClassesProbs <- cbind(predClassesProbs,probabilities)
    
  }else if(isTRUE(mixedMode)){
    ## First train objects
    if( !is.null(testX) && ! is.null(testY)){
      prob_test.t<-extractProb(modelList[trainModelsVec],unkX =testX)
      prob_test.t$obs <-dataSet[[yVar]]
      prob_test.t[["patient"]] <-dataSet[["run_accession"]]
      # prob_test.t[[labelVar]] <-dataSet[[labelVar]]
    }else{
      prob_test.t<-extractProb(modelList)    
    }
    # Bring those in the format of  train models
    # prob_test.t <-prob_test.t %>% dplyr::select(all_of(classes),"obs", "pred","object","patient",labelVar)
    prob_test.t <-prob_test.t %>% dplyr::select(all_of(classes),"obs", "pred","object","patient")
    
    ## Then RFE objects
    probabilities <- predict(modelList[rfeModelsVec], newdata=dataSet, type="prob")
    # probabilities.df <- bind_rows(probabilities, .id = "object")
    # predicted.classes <- predict(modelList, newdata=test.data, type="raw")
    
    probabilities <- lapply(probabilities, function(x) transform(x, obs = factor(dataSet[[yVar]], levels = classes)))
    probabilities <- lapply(probabilities, function(x) transform(x, pred = factor(pred, levels = classes)))
    probabilities <- lapply(probabilities, function(x) transform(x, patient = factor(dataSet[["run_accession"]])))
    # probabilities <- lapply(probabilities, function(x) transform(x, tissue = as.character(dataSet[[labelVar]])))
    predClassesProbs <-bind_rows(probabilities, .id = "object")
    prob_test.r <-predClassesProbs
    # Bring those in the format of  train models
    # prob_test.r <-prob_test.r %>% dplyr::select(all_of(classes),"obs", "pred","object","patient",labelVar)
    prob_test.r <-prob_test.r %>% dplyr::select(all_of(classes),"obs", "pred","object","patient")

    ## Merge objects
    prob_test <- rbind(prob_test.t, prob_test.r)
    
  }else{
    if( !is.null(testX) && ! is.null(testY)){
      prob_test<-extractProb(modelList,unkX =testX)
      prob_test$obs <-dataSet[[yVar]]
      prob_test$patient <- dataSet[["run_accession"]]
      prob_test[["obs"]] <- factor(prob_test[["obs"]],levels = ClassLevels)
    }else{
      prob_test<-extractProb(modelList)    
    }
    
  }
  #####################
  prob_test[["obs"]] <- factor(prob_test[["obs"]],levels = classes)
  #############################################
  
  if( !is.null(thres) ){
    # Change predicted based on threshold
    thres_preds <- function(predDF,thr,classes){
      # predDF <- prob_test.l[[1]]
      # thr <-thres[[1]]
      # classes <-classes
        
      predDF[["pred"]] <- ifelse(predDF[[classes[LevelIndex]]] >=thr,classes[LevelIndex],classes[-LevelIndex])
      predDF[["pred"]] <- factor(predDF[["pred"]],levels=levels(predDF[["obs"]]))
      predDF
    }
    prob_test.l <- split(prob_test,prob_test$object)
    prob_test.l <-  prob_test.l %>% list.subset(sigNames)
    prob_test.l <- lapply(1:length(prob_test.l), function(x) thres_preds(prob_test.l[[x]],as.numeric(thres[[x]]),classes))
    
    if(isTRUE(rfeE)){
      prob_test <- ldply(prob_test.l, rbind)
    }else{
      prob_test <- ldply(prob_test.l, rbind)[,-1]
    }
    
  }
  
  # Split by signature
  prob_test.obj <- split(prob_test,prob_test$object)
  # Order by the order of signatures
  prob_test.obj <-  prob_test.obj %>% list.subset(sigNames)
  if(isTRUE(rfeE) | isTRUE(mixedMode)){
    prob_test.list <-prob_test.obj
  }else{
    # Make a list of signatures scores
    dataSet.list <-setNames(lapply(names(dataSet[,sigs]), function(x) cbind(dataSet[,sigs][x])), names(dataSet[,sigs]))
    # Name them all score
    names(dataSet.list) <- rep("score", length(dataSet.list))
    # Merge DF list with predictions and list with signatures-ORDER SHOULD BE SAME FOR SIGNATURES
    prob_test.list <- Map(cbind, prob_test.obj, dataSet.list)
    # Rename column with signature scores into score in all elements
    scoreNameChange <-sapply(1:length(prob_test.list), function(x) colnames(prob_test.list[[x]])[dim(prob_test.list[[x]])[2]] <<- "score")
    
  }
  
  # List of dataframes into one DF
  prob_test <- ldply(prob_test.list, rbind)[,-1]
  prob_test$classification <- ifelse(prob_test$obs==prob_test$pred & prob_test$obs==classes[LevelIndex], "TP",
                                     ifelse(prob_test$obs==prob_test$pred & prob_test$obs!=classes[LevelIndex], "TN",
                                            ifelse(prob_test$obs!=prob_test$pred & prob_test$obs==classes[LevelIndex], "FN",
                                                   ifelse(prob_test$obs!=prob_test$pred & prob_test$obs!=classes[LevelIndex], "FP",NA))))
  
  # Make extra DF for ACTUAL OBSERVATIONS
  prob_test.obs <- prob_test.list[[1]][,which(colnames(prob_test.list[[1]]) %in% c("pred","obs","object", "patient","classification","score"))]
  
  prob_test.obs$object <- rep("Observed", dim(prob_test.obs)[1])
  prob_test.obs$pred <- prob_test.obs$obs
  prob_test.obs$classification <- ifelse(prob_test.obs$obs==prob_test.obs$pred & prob_test.obs$obs==classes[LevelIndex], "TP",
                                         ifelse(prob_test.obs$obs==prob_test.obs$pred & prob_test.obs$obs!=classes[LevelIndex], "TN",
                                                ifelse(prob_test.obs$obs!=prob_test.obs$pred & prob_test.obs$obs==classes[LevelIndex], "FN",
                                                       ifelse(prob_test.obs$obs!=prob_test.obs$pred & prob_test.obs$obs!=classes[LevelIndex], "FP",NA))))
  prob_test.obs$classification <- as.factor(prob_test.obs$classification)
  prob_test.obs$classification <- factor(prob_test.obs$classification, levels = c("FN","FP",levels(prob_test.obs$classification)))
  
  prob_test.obs <-prob_test.obs %>% distinct()
  
  # Select columns
  prob_test <- prob_test[,which(colnames(prob_test) %in% c("pred","obs", "object", "patient","classification","score"))]
  
  prob_test$classification <- as.factor(prob_test$classification)
  # Bind them
  prob_test.all <- rbind(prob_test,prob_test.obs)
  
  
  ##
  
  ##
  # Find order of patients
  prob_test.all$object <- factor(prob_test.all$object, levels =c(sigNames,"Observed") )
  prob_test.all.list <- split(prob_test.all,prob_test.all$object)
  
  order_listDFs_byPatients<- function(x, rfeE = FALSE){
    
    y <- split(prob_test.all.list[[x]],prob_test.all.list[[x]][["classification"]])
    
    if(isTRUE(rfeE)| isTRUE(mixedMode)){
      x1 <-y[[1]]  %>% pull(patient) %>% unique()
      x2 <-y[[2]] %>% pull(patient) %>% unique()
      x3 <-y[[3]] %>% pull(patient) %>% unique()
      x4 <-y[[4]] %>% pull(patient) %>% unique()
      patients.ordered.full <-c(x1,x4,x3,x2)
    }else{
      x1 <-y[[1]] %>% dplyr::select(patient,score) %>% group_by(patient)%>% dplyr::summarise(mean=mean(score)) %>% arrange(mean) %>% pull(patient) %>% unique()
      x2 <-y[[2]] %>% dplyr::select(patient,score) %>% group_by(patient)%>% dplyr::summarise(mean=mean(score)) %>% arrange(mean)%>% pull(patient) %>% unique()
      x3 <-y[[3]] %>% dplyr::select(patient,score) %>% group_by(patient)%>% dplyr::summarise(mean=mean(score)) %>% arrange(mean)%>% pull(patient) %>% unique()
      x4 <-y[[4]] %>% dplyr::select(patient,score) %>% group_by(patient)%>% dplyr::summarise(mean=mean(score)) %>% arrange(mean)%>% pull(patient) %>% unique()
      patients.ordered.full <-c(x1,x4,x3,x2)
    }
    
    # 
    prob_test.all.list[[x]]$patient <- factor(prob_test.all.list[[x]]$patient,levels=patients.ordered.full)
    prob_test.all.list[[x]]<- prob_test.all.list[[x]][order(prob_test.all.list[[x]]$patient),]
    prob_test.all.list[[x]]
  }
  
  # order_patients <-sapply(1:length(prob_test.all.list), function(x) prob_test.all.list[[x]]$patient <<- factor(prob_test.all.list[[x]]$patient,levels=patients.ordered.full))
  
  # orderDF_by_patients <-sapply(1:length(prob_test.all.list), function(x) prob_test.all.list[[x]]<<- prob_test.all.list[[x]][order(prob_test.all.list[[x]]$patient),])
  
  names(prob_test.all.list) <-c(sigNames,"Observed")
  if(isTRUE(rfeE) | isTRUE(mixedMode)){
    prob_test.all.list.ordered <-lapply(1:length(prob_test.all.list), function(z) order_listDFs_byPatients(z, rfeE = TRUE))
  }else{
    prob_test.all.list.ordered <-lapply(1:length(prob_test.all.list), function(z) order_listDFs_byPatients(z))
  }
  names(prob_test.all.list.ordered) <-c(sigNames,"Observed")
  
  # # Bind them, ordered
  # prob_test.all.final <- ldply(prob_test.all.list, rbind)[,-1]
  ##
  
  subTilePlot <- function(DF){
    g <- ggplot(data = DF)
    
    g <- g + geom_tile(aes(x = object, y = patient, fill = classification,width=0.85, height=0.95), color = "black")
    if(isTRUE(reverse)){
      g <- g + scale_fill_manual(values=c("Gainsboro","Gainsboro",rev(c("#004D40","#CD5C5C"))),drop=FALSE)
    }else{
      g <- g + scale_fill_manual(values=c("Gainsboro","Gainsboro","#004D40","#CD5C5C"),drop=FALSE)
    }
    
    g <- g + theme_ipsum(base_family = "Palatino Linotype", base_size=11)
    if(isTRUE(rotate)){
      g <- g +theme(plot.title = element_blank(),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.title.y = element_text(size=13, face="bold"),
                    axis.text.x = element_text(size=12,angle=45,hjust = 1, vjust = 1),#angle=25,angle=25,hjust = 1, vjust = 1,
                    axis.title.x = element_text(size=13, face="bold"),
                    # axis.text.y = element_blank(),
                    #   axis.ticks.y = element_blank(),
                    legend.title = element_blank(),
                    legend.text = element_blank(),
                    legend.position = "none",
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(colour = "black")
                    ,
                    plot.margin = unit(c(1, 0, 1, 0), "cm")
      )
    }else{
      g <- g +theme(plot.title = element_blank(),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.title.y = element_text(size=13, face="bold"),
                    axis.text.x = element_text(size=12),#angle=25,angle=25,hjust = 1, vjust = 1,
                    axis.title.x = element_text(size=13, face="bold"),
                    # axis.text.y = element_blank(),
                    #   axis.ticks.y = element_blank(),
                    legend.title = element_blank(),
                    legend.text = element_blank(),
                    legend.position = "none",
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(colour = "black")
                    ,
                    plot.margin = unit(c(1, 0, 1, 0), "cm")
      )
    }
    
    g <- g + ylab("Patients") + xlab("")#xlab("Observed response and signature predictions")
    g
  }
  
  tileplots <-lapply(1:length(prob_test.all.list.ordered), function(x) subTilePlot(prob_test.all.list.ordered[[x]]))
  
  # ggpubr::ggarrange(plotlist=tileplots, nrow = 1, common.legend = TRUE)
  # egg::ggarrange(plots=tileplots, nrow = 1)
  # 
  # cowplot::plot_grid(tileplots)
  if(mode =="rfe"){
    p <-cowplot::plot_grid(
      
      tileplots[[1]]+ theme(#axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank()),
      
      tileplots[[2]]+ theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank(),
                            axis.line.y=element_blank()
                            ,
                            # axis.text.x = element_blank(),
                            # axis.ticks.x = element_blank(),
                            axis.title.x = element_blank()
      ),
      tileplots[[3]]+ theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank() ,
                            axis.line.y=element_blank()
                            ,
                            # axis.text.x = element_blank(),
                            # axis.ticks.x = element_blank(),
                            axis.title.x = element_blank()
      ),
      tileplots[[4]]+ theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank() ,
                            axis.line.y=element_blank()
                            ,
                            # axis.text.x = element_blank(),
                            # axis.ticks.x = element_blank(),
                            axis.title.x = element_blank()
      ),
      tileplots[[5]]+ theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank() ,
                            axis.line.y=element_blank()
                            ,
                            # axis.text.x = element_blank(),
                            # axis.ticks.x = element_blank(),
                            axis.title.x = element_blank()
      ),
      tileplots[[6]]+ theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank() ,
                            axis.line.y=element_blank()
                            ,
                            # axis.text.x = element_blank(),
                            # axis.ticks.x = element_blank(),
                            axis.title.x = element_blank()
      ),
      tileplots[[7]]+ theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank() ,
                            axis.line.y=element_blank()
                            ,
                            # axis.text.x = element_blank(),
                            # axis.ticks.x = element_blank(),
                            axis.title.x = element_blank()
      ),
      # tileplots[[8]]+ theme(axis.text.y = element_blank(),
      #                         axis.ticks.y = element_blank(),
      #                         axis.title.y = element_blank() ,
      #                         axis.line.y=element_blank()
      #                         ,
      #                         # axis.text.x = element_blank(),
      #                         # axis.ticks.x = element_blank(),
      #                         axis.title.x = element_blank()
      # ),tileplots[[9]]+ theme(axis.text.y = element_blank(),
      #                         axis.ticks.y = element_blank(),
      #                         axis.title.y = element_blank() ,
      #                         axis.line.y=element_blank()
      #                         ,
      #                         # axis.text.x = element_blank(),
      #                         # axis.ticks.x = element_blank(),
      #                         axis.title.x = element_blank()
      # ),tileplots[[10]]+ theme(axis.text.y = element_blank(),
      #                         axis.ticks.y = element_blank(),
      #                         axis.title.y = element_blank() ,
      #                         axis.line.y=element_blank()
      #                         ,
      #                         # axis.text.x = element_blank(),
      #                         # axis.ticks.x = element_blank(),
      #                         axis.title.x = element_blank()
      # ),tileplots[[11]]+ theme(axis.text.y = element_blank(),
      #                         axis.ticks.y = element_blank(),
      #                         axis.title.y = element_blank() ,
      #                         axis.line.y=element_blank()
      #                         ,
      #                         # axis.text.x = element_blank(),
      #                         # axis.ticks.x = element_blank(),
      #                         axis.title.x = element_blank()
      # ),tileplots[[12]]+ theme(axis.text.y = element_blank(),
      #                         axis.ticks.y = element_blank(),
      #                         axis.title.y = element_blank() ,
      #                         axis.line.y=element_blank()
      #                         ,
      #                         # axis.text.x = element_blank(),
      #                         # axis.ticks.x = element_blank(),
      #                         axis.title.x = element_blank()
      # ),tileplots[[13]]+ theme(axis.text.y = element_blank(),
      #                         axis.ticks.y = element_blank(),
      #                         axis.title.y = element_blank() ,
      #                         axis.line.y=element_blank()
      #                         ,
      #                         # axis.text.x = element_blank(),
      #                         # axis.ticks.x = element_blank(),
      #                         axis.title.x = element_blank()
      # ),tileplots[[14]]+ theme(axis.text.y = element_blank(),
      #                         axis.ticks.y = element_blank(),
      #                         axis.title.y = element_blank() ,
      #                         axis.line.y=element_blank()
      #                         ,
      #                         # axis.text.x = element_blank(),
      #                         # axis.ticks.x = element_blank(),
      #                         axis.title.x = element_blank()
      # ),tileplots[[15]]+ theme(axis.text.y = element_blank(),
      #                         axis.ticks.y = element_blank(),
      #                         axis.title.y = element_blank() ,
      #                         axis.line.y=element_blank()
      #                         ,
      #                         # axis.text.x = element_blank(),
      #                         # axis.ticks.x = element_blank(),
      #                         axis.title.x = element_blank()
      # ),
      
      nrow=1, align = "h")
  }else if(mode =="mixed"){
    p <-cowplot::plot_grid(
      
      tileplots[[1]]+ theme(#axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=8),
        axis.title.x = element_blank()),
      
      tileplots[[2]]+ theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank(),
                            axis.line.y=element_blank()
                            ,
                            # axis.text.x = element_blank(),
                            # axis.ticks.x = element_blank(),
                            axis.text.x = element_text(size=8),
                            axis.title.x = element_blank()
      ),
      tileplots[[3]]+ theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank() ,
                            axis.line.y=element_blank()
                            ,
                            # axis.text.x = element_blank(),
                            # axis.ticks.x = element_blank(),
                            axis.text.x = element_text(size=8),
                            axis.title.x = element_blank()
      ),
      tileplots[[4]]+ theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank() ,
                            axis.line.y=element_blank()
                            ,
                            # axis.text.x = element_blank(),
                            # axis.ticks.x = element_blank(),
                            axis.text.x = element_text(size=8),
                            axis.title.x = element_blank()
      ),
      tileplots[[5]]+ theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank() ,
                            axis.line.y=element_blank()
                            ,
                            # axis.text.x = element_blank(),
                            # axis.ticks.x = element_blank(),
                            axis.text.x = element_text(size=8),
                            axis.title.x = element_blank()
      ),
      tileplots[[6]]+ theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank() ,
                            axis.line.y=element_blank()
                            ,
                            # axis.text.x = element_blank(),
                            # axis.ticks.x = element_blank(),
                            axis.text.x = element_text(size=8),
                            axis.title.x = element_blank()
      ),
      tileplots[[7]]+ theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank() ,
                            axis.line.y=element_blank()
                            ,
                            # axis.text.x = element_blank(),
                            # axis.ticks.x = element_blank(),
                            axis.text.x = element_text(size=8),
                            axis.title.x = element_blank()
      ),
      tileplots[[8]]+ theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank() ,
                            axis.line.y=element_blank()
                            ,
                            # axis.text.x = element_blank(),
                            # axis.ticks.x = element_blank(),
                            axis.text.x = element_text(size=8),
                            axis.title.x = element_blank()
      ),
      tileplots[[9]]+ theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank() ,
                            axis.line.y=element_blank()
                            ,
                            # axis.text.x = element_blank(),
                            # axis.ticks.x = element_blank(),
                            axis.text.x = element_text(size=8),
                            axis.title.x = element_blank()
      ),tileplots[[10]]+ theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_blank() ,
                              axis.line.y=element_blank()
                              ,
                              # axis.text.x = element_blank(),
                              # axis.ticks.x = element_blank(),
                              axis.text.x = element_text(size=8),
                              axis.title.x = element_blank()
      ),tileplots[[11]]+ theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_blank() ,
                              axis.line.y=element_blank()
                              ,
                              # axis.text.x = element_blank(),
                              # axis.ticks.x = element_blank(),
                              axis.text.x = element_text(size=8),
                              axis.title.x = element_blank()
      ),tileplots[[12]]+ theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_blank() ,
                              axis.line.y=element_blank()
                              ,
                              # axis.text.x = element_blank(),
                              # axis.ticks.x = element_blank(),
                              axis.text.x = element_text(size=8),
                              axis.title.x = element_blank()
      ),tileplots[[13]]+ theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_blank() ,
                              axis.line.y=element_blank()
                              ,
                              # axis.text.x = element_blank(),
                              # axis.ticks.x = element_blank(),
                              axis.text.x = element_text(size=8),
                              axis.title.x = element_blank()
      ),tileplots[[14]]+ theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_blank() ,
                              axis.line.y=element_blank()
                              ,
                              # axis.text.x = element_blank(),
                              # axis.ticks.x = element_blank(),
                              axis.text.x = element_text(size=8),
                              axis.title.x = element_blank()
      ),tileplots[[15]]+ theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_blank() ,
                              axis.line.y=element_blank()
                              ,
                              # axis.text.x = element_blank(),
                              # axis.ticks.x = element_blank(),
                              axis.text.x = element_text(size=8),
                              axis.title.x = element_blank()
      ),tileplots[[16]]+ theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_blank() ,
                              axis.line.y=element_blank()
                              ,
                              # axis.text.x = element_blank(),
                              # axis.ticks.x = element_blank(),
                              axis.text.x = element_text(size=8),
                              axis.title.x = element_blank()
      ),tileplots[[17]]+ theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_blank() ,
                              axis.line.y=element_blank()
                              ,
                              # axis.text.x = element_blank(),
                              # axis.ticks.x = element_blank(),
                              axis.text.x = element_text(size=8),
                              axis.title.x = element_blank()
      ),
      nrow=1, align = "h")
  }else{
    p <-cowplot::plot_grid(
      
      tileplots[[1]]+ theme(#axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.x = element_blank()),
      
      tileplots[[2]]+ theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank(),
                            axis.line.y=element_blank()
                            ,
                            # axis.text.x = element_blank(),
                            # axis.ticks.x = element_blank(),
                            axis.title.x = element_blank()
      ),
      tileplots[[3]]+ theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank() ,
                            axis.line.y=element_blank()
                            ,
                            # axis.text.x = element_blank(),
                            # axis.ticks.x = element_blank(),
                            axis.title.x = element_blank()
      ),
      tileplots[[4]]+ theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank() ,
                            axis.line.y=element_blank()
                            ,
                            # axis.text.x = element_blank(),
                            # axis.ticks.x = element_blank(),
                            axis.title.x = element_blank()
      ),
      tileplots[[5]]+ theme(axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title.y = element_blank() ,
                            axis.line.y=element_blank()
                            ,
                            # axis.text.x = element_blank(),
                            # axis.ticks.x = element_blank(),
                            axis.title.x = element_blank()
      ),
      nrow=1, align = "h")
  }
 
  # p <- p + theme(axis.text.x = element_text(size=20,angle=45,hjust = 1, vjust = 1))
  p
}


varimpPlot <- function(model,modelName){
  
  cat('\n')  
  cat('\n') 
  cat("#### ", modelName, " \n") 
  cat('\n') 
  
  varimp_data <- data.frame(feature = row.names(varImp(model)),
                            importance = varImp(model)[["Overall"]])
  
  g <-ggplot(data = varimp_data, 
             aes(x = reorder(feature, -importance), y = importance, fill = feature))
  g <- g +geom_bar(stat="identity") + labs(x = "Features", y = "Variable Importance")
  g <- g +geom_text(aes(label = round(importance, 2)), vjust=1.6, color="white", size=4) 
  g <- g + theme_classic(base_family = "Palatino Linotype", base_size=13)
  g <- g +theme(plot.title = element_text(size = 16, face = "bold"),
                axis.text.x = element_text(size=12,angle=25,hjust = 1, vjust = 1),
                axis.ticks.x = element_blank(),
                axis.title.x = element_text(size=16, face="bold"),
                axis.text.y = element_blank(),
                axis.title.y = element_text(size=14, face="bold"),
                legend.title = element_text(size=14,color="black",face="bold"),
                #legend.background = element_rect(size=0.1,color="black",fill=alpha("white",0.6)),
                legend.text = element_text(size=13),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                legend.position="none",legend.box="vertical", legend.margin=margin())
  print(g)
}

repCVrocPlot <- function(model,modelName){
  cat('\n')  
  cat('\n') 
  cat("#### ", modelName, " \n") 
  cat('\n') 
  
  # print(plot(model, type=c("g", "o")))
  print(ggplot(data = model, metric = "ROC") + theme_bw())
}


trainRUN.RFE <- function(DF,split,ClassLevels,PositiveClass,LevelIndex,seed.partition,seed.model,setCVseeds=NULL,
                         dummy=FALSE,dummyVars=NULL,categVars=NULL,obsColumn,responseVar,REPEATS,KFOLDS,
                         metricOpt,subsampMethod=NULL,inoutSub="in",
                         modelDesign,controlDum=FALSE,model.name,mainFeats,
                         removeOutliers=TRUE,quantileClass1=.95,quantileClass2=.95,selectOutlierMethod= "meanDist",
                         stratify=FALSE,fullResults=TRUE,
                         positiveCoefs=FALSE, filter=FALSE,filterVar=NULL,filterSelected=NULL,CorrRemove=FALSE,corrThres=NULL,
                         selectedScale=FALSE,sizeVector,rfeModel){
  ####################################
  # DF <-immdata.apd1.proc.TcrA.divEst.full.sub.ALL.filt 
  # split<-0.66
  # ClassLevels<-c("CRPR","PD")
  # PositiveClass<-"CRPR"
  # LevelIndex<-1
  # seed.partition=3456
  # seed.model=5627
  # setCVseeds=NULL
  # dummy=FALSE
  # dummyVars=c("tissue")
  # categVars=c("tissue","dataset")
  # obsColumn<-"run_accession"
  # responseVar<-"response"
  # REPEATS=30
  # KFOLDS=10
  # metricOpt="ROC"
  # subsampMethod="smote"
  # inoutSub="in"
  # modelDesign=paste(covariatesIn[c(1:3,9:11,13)],collapse = " + ")
  # controlDum =FALSE
  # model.name<-"multi"
  # # subset.notIn.sig<-degs.only
  # mainFeats<-covariatesIn[c(1:3,9:11,13)]
  # removeOutliers=FALSE
  # quantileClass1=.99
  # quantileClass2=.99
  # selectOutlierMethod= "none"
  # stratify=FALSE
  # fullResults=TRUE
  # positiveCoefs=FALSE
  # filter=TRUE
  # filterVar="dataset"
  # filterSelected=c("Powles","Gide","Mari","Rod","Hugo","Riaz")
  # CorrRemove= FALSE
  # corrThres=0.85
  # selectedScale<-TRUE
  # sizeVector = c(1:7)
  # rfeModel ="RF"
  #####################################
  #################################################################
  #################################################################
  #################################
  ##################################
  ################################
  ## CHECK 0 : SAMPLES, TISSUES ##
  ################################
  # DF %>% dplyr::filter(complete.cases(all_of(mainFeats))) %>% droplevels()
  DF <-DF %>% dplyr::filter_at(vars(all_of(mainFeats)), all_vars(complete.cases(.)))
  
  # DF <- DF %>% dplyr::filter(complete.cases(!!as.name(mainFeats[1]))) %>% droplevels()
  print("==INITIAL INPUT DATA==")
  print("=======================")
  print(table(DF[[responseVar]]))
  print(table(DF[[dummyVars]]))
  print("=======================")
  
  ##########################
  ## CHECK 1 : PARAMETERS ##
  ##########################
  
  print("==MODEL/PIPELINE PARAMETERS==")
  print("==============================")
  print(paste0("Outlier filtering: ",removeOutliers))
  if(isTRUE(removeOutliers)){
    print(paste0("Outlier Method: ",selectOutlierMethod))
  }
  
  print(paste0("Dataset Filtering: ",filter))
  if(isTRUE(filter)){
    print(paste0("Datasets selected: ",filterSelected))
  }
  
  print(paste0("Higly Correlated Features Filtering: ",CorrRemove))
  if(isTRUE(CorrRemove)){
    print(paste0("Higly Correlated Features Filtering Threshold: ",corrThres))
  }
  
  if(isTRUE(selectedScale)){
    print(paste0("Scaling of train data, scaling of test data with train data scale factor: ",selectedScale))
  }
  print("===============================")
  
  
  ################################
  ## CHECK 2 : MODEL PARAMETERS ##
  ################################
  
  print("==MODEL PARAMETERS==")
  
  print("=====================")
  print(paste0("Split: ",  split))
  print("-----------------------")
  print(paste0("Folds: ",KFOLDS))
  print("-----------------------")
  print(paste0("Repeats: ",REPEATS))
  print("-----------------------")
  print(paste0("Model seed: ",seed.model))
  print("-----------------------")
  print(paste0("Set CV seeds: ",setCVseeds))
  print("-----------------------")
  print(paste0("Subsampling: ",inoutSub))
  print("-----------------------")
  print(paste0("Subsampling Method: ",subsampMethod))
  print("-----------------------")
  print(paste0("Model metric:",metricOpt))
  print("-----------------------")
  print(paste0("Positive Class: ",PositiveClass))
  print("==========================")
  
  ###########################################
  ## Detect outliers with multiple methods ##
  ###########################################
  # Make sure response var for split is not a factor#
  DF[[responseVar]] <- as.factor(DF[[responseVar]])
  # Split to responders and non-responders
  data_full.crpr <- DF %>% dplyr::filter(response==ClassLevels[LevelIndex])
  data_full.pd <- DF %>% dplyr::filter(response==ClassLevels[-LevelIndex])
  
  # ## Run automatic methods from package OutlierDetection
  # method.names<- c("Robust Kernal-based Outlier Factor(RKOF) algorithm","Depth based method", "Generalised dispersion method","Mahalanobis Distance","k Nearest Neighbours Distance method", "kth Nearest Neighbour Distance method", "Outlier Detection(Intersection of all the methods)","Principal Component Outlier Detection(Intersection of all the methods applied on pc’s)" )
  # 
  
  ####################################
  ## Clean full data from outliers ##
  ####################################
  if (isTRUE(removeOutliers)){
    print("~Removing outliers~")
    print(paste0("== Removing ",length(patient.crpr)," responders =="))
    print(paste0("== Removing ",length(patient.pd)," non-responders =="))
    print(paste0("~We will train/test the model with ", dim(DF)[1]-length(patient.crpr)-length(patient.pd),
                 " out of ",dim(DF)[1]," observations!~"))
    # data_train.full<-data_train
    
    DF <- DF[-outliers.all.index,]
  }else{
    print("~Outliers detected BUT NOT REMOVED~")
    print("~Outlier observations still reported in results!~")
    DF<-DF
  }
  
  # Make sure response var for split is not a factor#
  DF[[responseVar]] <- as.character(DF[[responseVar]])
  
  if(isTRUE(filter)){
    filterVar.un <- as.name(filterVar)
    ## Here I keep all the leftover from the filtering to be added as test data later, otherwise they are missed
    DF.left <- DF %>% dplyr::filter(!!filterVar.un %notin% filterSelected) %>% droplevels()
    DF <- DF %>% dplyr::filter(!!filterVar.un %in% filterSelected) %>% droplevels()
    print("~Filtering data~")
    print("== Including the following :==")
    print(paste(filterSelected, collapse = ", "))
    print(table(DF[[responseVar]]))
  }
  ## CHECK IF YOU NEED TO FILTER ##
  
  #######################################
  ## Remove highly correlated features ##
  #######################################
  if (isTRUE(CorrRemove)){
    print("~Identifying Correlated Features~")
    # calculate correlation matrix
    correlationMatrix <- cor(DF[,mainFeats], method = "spearman")
    #
    # find attributes that are highly corrected (ideally >0.75)
    print(paste0("~Removing Correlated Features, with thresold:",corrThres, "~"))
    highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=corrThres,names = TRUE, exact=TRUE)
    keep.cols<-setdiff(colnames(DF),highlyCorrelated)
    print("~Features removed:~")
    print(paste(highlyCorrelated,collapse = ", "))
    DF<-DF[,keep.cols]    # Re-define design of models
    if(isTRUE(filter)){
      DF.left <-DF.left[,keep.cols]
    }
    
    print("~Re-defining the model design after removing correlated features~")
    mainFeats.new <-setdiff(mainFeats,highlyCorrelated)
    modelDesign<-paste(mainFeats.new,collapse = " + ")
    mainFeats<- mainFeats.new
  }
  
  ################
  ## DATA SPLIT ##
  ################
  # Make sure response var for split is not a factor#
  DF[[responseVar]] <- as.character(DF[[responseVar]])
  
  ## Split data
  set.seed(NULL)
  set.seed(seed.partition)
  trainIndex <- createDataPartition(DF[[responseVar]], p=split, list=FALSE)
  data_train <- DF[ trainIndex,]
  data_test <- DF[-trainIndex,]
  # NOW SET FACTOR LEVELS
  data_train[[responseVar]] <- factor(data_train[[responseVar]], levels=ClassLevels)
  data_test[[responseVar]] <- factor(data_test[[responseVar]], levels=ClassLevels)
  if(isTRUE(filter)){
    DF.left[[responseVar]] <- factor(DF.left[[responseVar]], levels=ClassLevels)
  }
  
  print("~Using train data~")
  print(table(data_train[[responseVar]]))
  print("~Using test data~")
  print(table(data_test[[responseVar]]))
  
  ####################################
  ## Create dummy vars if necessary ##
  ####################################
  
  if (isTRUE(controlDum)){
    print("Creating dummy variables for categorical variable(s)")
    form.dummy <- as.formula(paste(responseVar,' ~ ', paste(dummyVars,collapse = " + ")))
    
    ## Train ##
    data_train[[responseVar]] <- as.character(data_train[[responseVar]])
    dummies.df <- dummyVars(form.dummy, data = data_train,fullRank=T,levelsOnly=TRUE)# Gide is used as reference
    
    dummies.df <- data.frame(predict(dummies.df, newdata = data_train))
    # Merge with original df
    data_train <- cbind(data_train,dummies.df)
    dummyCols <- colnames(dummies.df)
    
    ## Test ##
    if(length(levels(as.factor(data_test$tissue)))<=1){
      print("Cannot create dummy vars for test data here")
    }else{
      data_test[[responseVar]] <- as.character(data_test[[responseVar]])
      dummies.df <- dummyVars(form.dummy, data = data_test,fullRank=T,levelsOnly=TRUE)# Gide is used as reference
      
      dummies.df <- data.frame(predict(dummies.df, newdata = data_test))
      # Merge with original df
      data_test <- cbind(data_test,dummies.df)
      if(isTRUE(filter)){
        DF.left[[responseVar]] <- as.character(DF.left[[responseVar]])
        # levels(DF.left[[dummyVars]]) <- list(melanoma="melanoma",
        #                           bladder="bladder",
        #                           rcc="rcc",
        #                           gastric="gastric")
        # dummies.df <- dummyVars(form.dummy, data = DF.left,fullRank=T,levelsOnly=TRUE)# Gide is used as reference
        # 
        # dummies.df <- data.frame(predict(dummies.df, newdata = DF.left))
        # DF.left <- cbind(DF.left,dummies.df)
        # Here I mannually add this
        for (dumCol in dummyCols){
          # dumCol.un <- as.name(dummyCols[1])
          # Add a column
          DF.left <- DF.left %>% add_column(newCol = rep(NA,dim(DF.left)[1]))
          colnames(DF.left)[dim(DF.left)[2]] <-dumCol
        }
        ## Manual fix
        DF.left$tissuegastric <- ifelse(DF.left$tissue=="gastric",1,0)
        DF.left$tissuemelanoma <- ifelse(DF.left$tissue=="melanoma",1,0)
        DF.left$tissueRCC <- ifelse(DF.left$tissue=="RCC",1,0)
      }
      
    }
    
    
  }
  
  # NOW SET FACTOR LEVELS
  print("==SET FACTOR LEVELS TO RESPONSE VARIABLE==")
  data_train[[responseVar]] <- factor(data_train[[responseVar]], levels=ClassLevels)
  data_test[[responseVar]] <- factor(data_test[[responseVar]], levels=ClassLevels)
  if(isTRUE(filter)){
    print("==Filtering applied, merging leftover test data with test data==")
    DF.left[[responseVar]] <- factor(DF.left[[responseVar]], levels=ClassLevels)
    ## Now merge data_test and DF leftover
    data_test <- rbind(data_test,DF.left)
  }
  #####################################################
  ## Run  Logistic Regression model ##
  #####################################################
  
  ## Stratified Repeate CV
  if(isTRUE(controlDum)){
    print("==Adjusting for confounding categorical variable, changing design==")
    controlModel <-paste(dummyCols,collapse = " + ")
    
    modelDesign <- paste0(modelDesign, " + ",controlModel)
  }
  
  # ###########################################################################
  # ## SCALE IF SELECTED - FIRST TRAIN THEN TEST with train's scaling factor ## WRONG-SCALING SHOULD BE IN CV
  # ###########################################################################
  # if(isTRUE(scale)){
  #   ## MEthod A-from caret
  #   print("==Scaling TRAIN data==")
  #   preProcValues <- preProcess(data_train[,mainFeats], method = "scale")
  #   trainTransformed <- predict(preProcValues, data_train[,mainFeats])
  #   print("==Scaling TEST data with the scaling factor from TRAIN data==")
  #   testTransformed <- predict(preProcValues, data_test[,mainFeats])
  #   
  #   
  #   # ## MEthod B- just scale
  #   # trainTransformed2<- data_train %>% dplyr::select(all_of(mainFeats)) %>% scale()
  #   # 
  #   # testTransformed2<- data_test %>% dplyr::select(all_of(mainFeats)) %>% scale(center=attr(trainTransformed2,"scaled:center"),
  #   #                                                                        scale=attr(trainTransformed2,"scaled:scale"))
  #   
  #   data_train <- cbind(trainTransformed,data_train[,(length(mainFeats)+1): dim(data_train)[2]])
  #   data_test <- cbind(testTransformed,data_test[,(length(mainFeats)+1): dim(data_test)[2]])
  #   
  # }
  
  
  ##########################################
  ## CHECK 3 : SAMPLES, TISSUES IN TRAIN  ##
  ##########################################
  
  print("==INITIAL TRAIN DATA==")
  print("=======================")
  print(table(data_train[[responseVar]]))
  print(table(data_train[[dummyVars]]))
  print(table(data_train[["dataset"]]))
  print("=======================")
  
  # For training the model create a train object with only features and response, since other info might confound, especially
  # when out-subsampling
  if(isTRUE(controlDum)){
    data_train2 <-data_train[,c(mainFeats,dummyCols,responseVar)]
  }else{
    data_train2 <-data_train[,c(mainFeats,responseVar)]
  }
  
  
  print("==MODEL TRAIN INITIATION==")
  
  model_rfe <-multiTRAINModels_RFE.repCV_construct(data_train2,responseVar,modelDesign, feats = mainFeats,
                                                   kfolds = KFOLDS,reps = REPEATS,sizeVec = sizeVector,
                                                   subMethod=subsampMethod,inoutSub="in",
                                                   metricSelected = metricOpt,
                                                   seedModel=seed.model,seedPart = seed.partition,seeds=setCVseeds,
                                                   stratify=stratify,positiveCoefs=FALSE,
                                                   applyScaling=selectedScale,wrapperFuncs = rfeModel)
  
  
  # model_rfe
  #################################

  print("====GATHERING  OUTPUT DATA===")

  # print("=STEP 1=")
  # Extract signature coefficients in DF
  model_optVars <-model_rfe$optVariables


  ## Extract probabilities on Test data with logit2prob
  # Think about how to include the second signature which is a subset of the first.
  if(isTRUE(controlDum)){
    data_test.s <-data_test[,c(mainFeats,dummyCols)]
  }else{
    # data_test.s <-data_test[,mainFeats]
    data_test.s <-data_test %>% dplyr::select(all_of(mainFeats))
  }

  # ## HERE IS THE POINT were I would need to preprocess test data - I am not applying the model here, but instead calculating the score.
  # if(isTRUE(selectedScale)){
  #   preProc.train <-model.logreg[[model.name]]$preProcess
  #   # data_test.scaled <- predict(preProcValues, data_test.s)
  #   data_test.s<- predict(preProc.train, data_test.s)
  #   ## THERE IS AN ISSUE WITH SCALING THE DATA WHEN HAVING DUMMY VARS, since preprocess takes everything in!!
  #   data_test.sc<-data_test.s
  #   if(isTRUE(controlDum)){
  #     # data_test.s <-data_test.s[,mainFeats]
  #     data_test.s <-data_test %>% dplyr::select(all_of(mainFeats))
  #   }
  # }

  if(isTRUE(controlDum)){
    # data_test.sub <-data_test.s[,mainFeats]
    data_test.sub <-data_test %>% dplyr::select(all_of(mainFeats))
  }

  ## Do below when test data are not empty
  if(isTRUE(nrow(data_test.s)<=4)){
    print("No calculations on test data, since they are empty,very few")

  }else{

    ## Confustion Matrix
    # test_pred=extractProb(model_rfe,unkX = data_test.s)
    # cm_log= confusionMatrix(test_pred[["pred"]], data_test[[responseVar]])

    test_pred <- predict(model_rfe, newdata=data_test.s, type="prob")
    cm_log= confusionMatrix(test_pred[["pred"]], data_test[[responseVar]])
  }


  ## WHEN REMOVING OUTLIERS IN FULL DATA
  # print("=STEP 9=")
  ## Gather in a list
  if(isTRUE(nrow(data_test.s)<=4)){
    if (isTRUE(fullResults)){
      results <-list(model_rfe,model_optVars,data_train,data_test,seed.partition,mainFeats)#patient.crpr,patient.pd,outliersMethod.crpr,outliersMethod.pd,dt_crpr,percentileClass1,dt_pd,percentileClass2,
      names(results)<- c("model","model.optVars","trainDF","testDF","seed.partition","features")#"outliers.crpr","outliers.pd","outMethodRes.crpr","outMethodRes.pd","pc.crpr","thres.crpr","pc.pd","thres.pd",
    }else{
      results <-list(model_rfe,model_optVars)#,patient.crpr,patient.pd,outliersMethod.crpr,outliersMethod.pdmainFeats
      names(results)<- c("model","model.optVars")#,"outliers.crpr","outliers.pd","outMethodRes.crpr","outMethodRes.pd","features"
    }
  }else{
    if (isTRUE(fullResults)){
      results <-list(model_rfe,model_optVars,data_train,data_test,seed.partition,cm_log,mainFeats)#probsDF.logit,,patient.crpr,patient.pd,outliersMethod.crpr,outliersMethod.pd,dt_crpr,percentileClass1,dt_pd,percentileClass2
      names(results)<- c("model","model.optVars","trainDF","testDF","seed.partition","confusionMatrix","features")#,"logit2prob.test","outliers.crpr","outliers.pd","outMethodRes.crpr","outMethodRes.pd","pc.crpr","thres.crpr","pc.pd","thres.pd",
    }else{
      results <-list(model_rfe,model_optVars,cm_log,mainFeats)#,patient.crpr,patient.pd,outliersMethod.crpr,outliersMethod.pd
      names(results)<- c("model","model.optVars","confusionMatrix","features")#,"outliers.crpr","outliers.pd","outMethodRes.crpr","outMethodRes.pd"
    }
  }

  # if(isTRUE(scale)){
  #   results <-c(results, scalingFactor=list(preProcValues))
  # }
  # print("=STEP 10=")
  results

}

## PARAMETER OPTIMIZATION

fullOptimalCutoffModelResult.RFE <- function(df,df.extra,selectedSplit,ClassLevels,PositiveClass,r,f,METRIC,
                                         subSamplingParams1,subSamplingParams2,
                                         cD,rO,quant,outMethod,pCl,filt,filtVar,filtSel,cr,crt,
                                         tissue="none",doScaling=FALSE,modelName='',applyStratify = FALSE,
                                         mainPredictors,rfeSizeVector,selectedRFEModel,yVar,sampleID ){

  #####################################
  model_rfe <-trainRUN.RFE(DF = df,
                           split = selectedSplit, ClassLevels = ClassLevels, PositiveClass = PositiveClass,LevelIndex = 1,
                           seed.partition = 3456,seed.model = 5627,setCVseeds=NULL,
                           dummy=cD,dummyVars=c("tissue"),categVars=c("tissue","dataset"),
                           obsColumn = sampleID,responseVar = yVar,
                           REPEATS=r,KFOLDS=f,
                           metricOpt = METRIC,subsampMethod=subSamplingParams1,inoutSub=subSamplingParams2,
                           modelDesign=paste(mainPredictors,collapse = " + "),controlDum=cD,
                           model.name = "multi",mainFeats=mainPredictors,
                           removeOutliers=rO,quantileClass1=quant,quantileClass2=quant,selectOutlierMethod= outMethod,
                           stratify=applyStratify,fullResults=TRUE,
                           positiveCoefs=pCl,
                           filter=filt,filterVar=filtVar,filterSelected=filtSel,
                           CorrRemove=cr,corrThres=crt,
                           selectedScale=doScaling,
                           sizeVector = rfeSizeVector,
                           rfeModel = selectedRFEModel)

  # model_rfe <- loadRData(paste0(projectRDataPath,"testing_rf",Rdata.suffix))
  subMethod <-as.character(subSamplingParams1)
  inoutSubMethod <- as.character(subSamplingParams2)

  models.all <- list(model_rfe)
  names(models.all)<- c(subMethod)
  model.final.res <- models.all[[1]]


  ## Initial train object from initial partition, AFTER OUTLIER & FILTERING , split,NO SCALING
  train.df<-model.final.res$trainDF
  ## Initial test object from initial partition, AFTER OUTLIER & FILTERING , split, NO SCALING
  test.df<- model.final.res$testDF

  ####################################
  ## Order DF by original test data ##
  ####################################

  testDF <-checkLogRegModelTEST.rfe(DF = df.extra,models.all = models.all,model.name = subMethod,
                                    obsColumn = "run_accession",
                                    filtered = TRUE,scaled=TRUE)

  data.testX <- as.data.frame(testDF[,which(colnames(testDF) %in% mainPredictors)])


  train.rfe <- list(rfeModel = model_rfe$model
                    #                   ,
                    #                   lrModel = model_lr,
                    #                   nbModel = model_nb
                    #                   # svmRadModel = model_rfe3,
                    #                   # nnetModel = model_rfe4
  )

  sigs <- c(selectedRFEModel)
  sigs_axisNames<- c(selectedRFEModel)
  sigs_titleNames<- c(selectedRFEModel)
  names(train.rfe) <- c(paste0(selectedRFEModel, ': R ~ ',paste(model_rfe$model.optVars, collapse = " + ")))
  ##################################
  # # Evaluate - no optimal threshold
  # eval_thres <- function(ModelsList, ModelIndex,responseVar, classes,LevelIndex){
  #
  #   # cat('\n')
  #   # cat('\n')
  #   # cat("#### ", names(ModelsList[ModelIndex]), "\n")
  #   # cat('\n')
  #
  #   # opt.thres <- optimal_threshold_value(seq(0.1,0.99,length.out = 30),ModelsList,data_test, testX=data_testX,testY=data_test[[responseVar]], names(ModelsList[ModelIndex]))
  #
  #   extract_EvaluationMetrics(ModelsList[[ModelIndex]],testDF, classes,LevelIndex,
  #                             modelName=names(ModelsList)[ModelIndex],var=sigs[ModelIndex],
  #                             savePlots = FALSE,
  #                             getTable=NULL,printP=FALSE, thres = NULL,rfeE = TRUE)
  # }
  #
  # eval.all.thres<- sapply(1:length(train.rfe), function(x) eval_thres(train.rfe,x,yVar, ClassLevels,1) )


  eval_thres <- function(ModelsList, ModelIndex,responseVar, classes,LevelIndex){

    # ModelsList <-train.rfe
    # ModelIndex <-1
    # responseVar <-yVar
    # classes <-ClassLevels
    # LevelIndex <-1

    opt.thres <- optimal_threshold_value(seq(0.1,0.99,length.out = 50),ModelsList[[ModelIndex]],
                                         testDF,responseVar, testX=data.testX,testY=testDF[[yVar]],
                                         names(ModelsList)[ModelIndex],classes,LevelIndex,
                                         MetrixMax="Youdens.J.statistic",extractedProbDF=NULL, printP=FALSE,rfeMode = TRUE)
    # if (ModelIndex==2){
    #   opt.thres <- logTMBClassThres
    # }
    print(sigs[ModelIndex])
    extract_EvaluationMetrics(ModelsList[[ModelIndex]],testDF, classes,LevelIndex,
                              modelName=names(ModelsList)[ModelIndex],var=sigs[ModelIndex],
                              savePlots = FALSE,
                              getTable=NULL,printP=FALSE, thres = opt.thres,rfeE = TRUE)
    #
  }

  eval.all.thres<- sapply(1:length(train.rfe), function(x) eval_thres(train.rfe,x,yVar, ClassLevels,1) ,simplify = TRUE)

  rownames(eval.all.thres) <- c("Accuracy", "Balanced Accuracy", "Sensitivity", "Specificity", "ROC-AUC", "Precision", "Recall", "PR-AUC", "F1-Score","Youden's J-statistic", "PPV","NPV")
  colnames(eval.all.thres) <-c(sigs )
  eval.all.thres.t <-as.data.frame(t(eval.all.thres))
  # eval.sub <- subset(eval.all.thres, rownames(eval.all.thres) %in% c("ROC-AUC","PPV","NPV"))
  # Add optimal lamda and alpha and gene number
  # optimal.params <- model.final.res$optimal.model.Params
  # final_res <- list(t(eval.sub),optimal.params)
  # final_res
  # optimal.params
  ## Make a character vector

  ##  Kfolds, repeats, lamdaOpt, alphaOpt, ROC_Tut,PPV_Tut, NPV_TuT, ROC_TuTifng, PPV_TuTifng, NPV_TuTifng, ROC_TIS, PPV_TIS, NPV_TIS
  res <- c(rfeModel=selectedRFEModel,NOfeaturesIn=length(mainPredictors),NOfeaturesOut=length(model_rfe$model.optVars),kfolds=f, repeats=r,
           subsampling=subSamplingParams1,subsamplingMethod=subSamplingParams2,
           split=selectedSplit,dataset = paste(filtSel,collapse = ", "),
           outlierDetection=rO,outlierMethod=outMethod,
           removeCorr=cr,corThres=crt,scaling=doScaling,optPredictors=paste(model_rfe$model.optVars, collapse = ", "),
           ROC=eval.all.thres.t$`ROC-AUC`, PPV=eval.all.thres.t$PPV,NPV=eval.all.thres.t$NPV)

  # res <- list(results = res,
  #             model = model_rfe)
  res

}

myCBfriendly.dark.small <- c("#222255","#225555","#225522","#666633","#663333","#555555")
# show_col(myCBfriendly.dark.small)

myCBfriendly.pale.small <- c("#bbccee","#cceeff","#ccddaa","#eeeebb","#ffcccc","#dddddd")
# show_col(myCBfriendly.pale.small) 

myCBfriendly.light.medium <- c("#77aadd","#99ddff","#44bb99","#bbcc33","#aaaa00","#eedd88","#ee8866","#ffaabb","#dddddd")
# show_col(myCBfriendly.light.medium)

myCBfriendly.bright.small <- c("#4477aa","#66ccee","#228833","#ccbb44","#ee6677","#bbbbbb")
# show_col(myCBfriendly.bright.small)

myCBfriendly.highContr.small <- c("#ffffff","#ddaa33","#bb5566","#004488","#000000")
# show_col(myCBfriendly.highContr.small)

myCBfriendly.muted.medium<- c("#332288","#88ccee","#44aa99","#117733","#999933","#ddcc77","#cc6677","#882255","#aa4499","#dddddd")
# show_col(myCBfriendly.muted.medium)


# show_col(c(myCBfriendly.highContr.small,myCBfriendly.light.medium))

myCBfriendly.large <-c(brewer.pal(9,"Set2"),brewer.pal(8,"Dark2"),"#cc3311")
# show_col(myCBfriendly.large)
simpleBarchart <- function(meltedDF, xVar, yVar, fillVar, facetVar=NULL, xLabel, yLabel,logYtrans=FALSE,seed, percentLabel=FALSE,legendTitle=NULL,ratioLine=FALSE, type="reps"){
  xVar.un <-  as.name(xVar)
  yVar.un <- as.name(yVar)
  fillVar.un <- as.name(fillVar)
  if (!is.null(facetVar)){
    facetVar.un <- as.name(facetVar)
    facetForm <- as.formula(paste(' ~ ', facetVar))
  }
  
  
  # # colors
  # if(length(levels(meltedDF[[fillVar]])) <= 6){
  #   palette <- myCBfriendly.pale.small
  #   set.seed(seed)
  #   palette <-sample(myCBfriendly.muted.medium)
  # }else if (length(levels(meltedDF[[fillVar]])) <= 9){
  #   set.seed(seed)
  #   # palette <- cbPalette
  #   # palette <-  sample(brewer.pal(9,"Paired"))
  #   # palette <- viridis_pal(option="D")(9)
  #   
  #   palette <- sample(myCBfriendly.light.medium)
  # }else{
  #   set.seed(seed)
  #   palette <- distinctColorPalette(length(levels(meltedDF[[fillVar]])))
  #   palette <- c(brewer.pal(9,"Set2"),brewer.pal(8,"Dark2"),"#cc3311")
  # }
  # 
  ############################
  if(length(levels(meltedDF[[fillVar]])) <= 6){
    palette <- myCBfriendly.pale.small
    set.seed(seed)
    if(type=="reps"){
      palette <-sample(myCBfriendly.muted.medium)
      
    }else{
      palette <- c(brewer.pal(9,"Set2"),brewer.pal(8,"Dark2"),"#cc3311")
    }
    
  }else if (length(levels(meltedDF[[fillVar]])) <= 9){
    set.seed(seed)
    # palette <- cbPalette
    # palette <-  sample(brewer.pal(9,"Paired"))
    # palette <- viridis_pal(option="D")(9)
    if(type=="reps"){
      palette <- sample(myCBfriendly.light.medium)
      
    }else{
      palette <- c(brewer.pal(9,"Set2"),brewer.pal(8,"Dark2"),"#cc3311")
    }
    
    
  }else{
    set.seed(seed)
    palette <- distinctColorPalette(length(levels(meltedDF[[fillVar]])))
    # palette <- c(brewer.pal(9,"Set2"),brewer.pal(8,"Dark2"),"#cc3311")
  }
  ########################################################
  g <- ggplot(meltedDF, aes(fill=!!fillVar.un, y=!!yVar.un, x=!!xVar.un)) 
  g <- g + geom_bar(stat="identity", width = 0.7,colour="black",size=0.1)#position="stack", 
  # g <- g + geom_bar(position="dodge", stat="identity")
  g <- g + ylab(yLabel) + xlab(xLabel)
  # g <- g + scale_fill_viridis(discrete = T)
  # g <- g + ggtitle(title)
  g <- g + scale_fill_manual(values=palette)
  if (!is.null(facetVar)){
    # g <- g + facet_wrap(~Var1_1, ncol = 1,strip.position = "left")
    g <- g + facet_grid(facetForm,  space  = "free", scales = "free", drop = TRUE, shrink = TRUE)  
    
  }
  
  if (isTRUE(percentLabel)){
    g <- g + geom_text(aes(x = !!xVar.un, label = paste0(round(!!yVar.un,1),'%')),
                       colour = 'black', position=position_stack(vjust=0.5), size=3)
  }
  
  
  g <- g + theme_ipsum(base_family = "Palatino Linotype", base_size=13)
  
  g <- g +theme(plot.title = element_text(size = 16, face = "bold"),
                axis.text.x = element_text(size=12,angle=25,hjust = 1, vjust = 1),
                axis.ticks.x = element_blank(),
                axis.title.x = element_text(size=16, face="bold"),
                axis.text.y = element_text( hjust = 1, size=13),#vjust = 1,
                axis.title.y = element_text(size=16, face="bold"),
                legend.title = element_text(size=16,color="black",face="bold"),
                #legend.background = element_rect(size=0.1,color="black",fill=alpha("white",0.6)),
                legend.text = element_text(size=13),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                legend.position="bottom",
                legend.box = 'vertical',
                legend.background = element_rect(size = rel(1.6),color="white", fill=alpha("white",0.6)),
                legend.margin=margin())
  if(is.null(legendTitle)){
    g <- g +labs(fill="Models")
  }else{
    g <- g +labs(fill=legendTitle)
  }
  if(isTRUE(ratioLine)){
    g<-g + geom_abline(slope=0, intercept=0.7,  col = "black",lty=2)
  }
  
  
  g
  
}


get_probabilityObj <- function(modelList,dataSet,responseVar,
                               xVar,yVar,
                               testX,testY,
                               classes, predClasses,
                               LevelIndex,posClassName,
                               thres,
                               extractedProbDF = NULL, rfeE=FALSE,mixedMode = FALSE,
                               trainModelsVec = NULL,
                               rfeModelsVec = NULL){
  # modelList<-modelList.final
  # dataSet<-testDF.survival
  # responseVar<-"response"
  # xVar <- "run_accession"
  # yVar <- "CRPR"
  # testX <- data.testX.survival
  # testY <- testDF.survival$response
  # # modelName <- names(modelList.final)[3]
  # classes <- c("CRPR","PD","MR","SD","NE","UNK")
  # predClasses <- c("CRPR","PD")
  # LevelIndex <- 1
  # thres=models.optThres.final
  # posClassName <- "CRPR"
  #   rfeE=FALSE
  # extractedProbDF=NULL
  # mixedMode = TRUE
  # trainModelsVec = c(1:6)
  # rfeModelsVec = c(7)
  
  xVar.un <- as.name(xVar)
  yVar.un <- as.name(yVar)
  
  
  if(isTRUE(rfeE)){
    # b <-sapply(1: length(modelList), function(x) predict(modelList[[x]], newdata=test.data), simplify = FALSE)
    probabilities <- predict(modelList, newdata=dataSet, type="prob")
    # probabilities.df <- bind_rows(probabilities, .id = "object")
    # predicted.classes <- predict(modelList, newdata=test.data, type="raw")
    
    probabilities <- lapply(probabilities, function(x) transform(x, obs = factor(dataSet[[responseVar]], levels = predClasses)))
    probabilities <- lapply(probabilities, function(x) transform(x, pred = factor(pred, levels = predClasses)))
    probabilities <- lapply(probabilities, function(x) transform(x, run_accession = factor(dataSet[[xVar]])))
    # probabilities <- lapply(probabilities, function(x) transform(x, tissue = as.character(dataSet[[labelVar]])))
    predClassesProbs <-bind_rows(probabilities, .id = "object")
    prob_test <-predClassesProbs
    # predClassesProbs <- cbind(predClassesProbs,probabilities)
    
  }else if(isTRUE(mixedMode)){
    ## First train objects
    if( !is.null(testX) && ! is.null(testY)){
      prob_test.t<-extractProb(modelList[trainModelsVec],unkX =testX)
      prob_test.t$obs <-dataSet[[responseVar]]
      prob_test.t[[xVar]] <-dataSet[[xVar]]
      # prob_test.t[[labelVar]] <-dataSet[[labelVar]]
    }else{
      prob_test.t<-extractProb(modelList)    
    }
    # Bring those in the format of  train models
    prob_test.t <-prob_test.t %>% dplyr::select(all_of(predClasses),"obs", "pred","object",xVar)
    
    ## Then RFE objects
    probabilities <- predict(modelList[rfeModelsVec], newdata=dataSet, type="prob")
    # probabilities.df <- bind_rows(probabilities, .id = "object")
    # predicted.classes <- predict(modelList, newdata=test.data, type="raw")
    
    probabilities <- lapply(probabilities, function(x) transform(x, obs = factor(dataSet[[responseVar]], levels = classes)))
    probabilities <- lapply(probabilities, function(x) transform(x, pred = factor(pred, levels = predClasses)))
    probabilities <- lapply(probabilities, function(x) transform(x, run_accession = factor(dataSet[[xVar]])))
    # probabilities <- lapply(probabilities, function(x) transform(x, tissue = as.character(dataSet[[labelVar]])))
    predClassesProbs <-bind_rows(probabilities, .id = "object")
    prob_test.r <-predClassesProbs
    # Bring those in the format of  train models
    prob_test.r <-prob_test.r %>% dplyr::select(all_of(predClasses),"obs", "pred","object",xVar)
    ## Merge objects
    prob_test <- rbind(prob_test.t, prob_test.r)
    
  }else{
    if( !is.null(testX) && ! is.null(testY)){
      prob_test<-extractProb(modelList,unkX =testX)
      prob_test$obs <-dataSet[[responseVar]]
      prob_test[[groupVar]] <- dataSet[[groupVar]]
      prob_test[[xVar]] <- dataSet[[xVar]]
      # prob_test[[labelVar]] <- dataSet[[labelVar]]
    }else{
      prob_test<-extractProb(modelList)    
    }
    
  }
  #####################
  prob_test[["obs"]] <- factor(prob_test[["obs"]],levels = classes)
  
  
  
  if( !is.null(thres) ){
    # Change predicted based on threshold
    thres_preds <- function(predDF,thr,classes){
      # predDF <-prob_test
      # thr <-thres
      # classes <-classes
      
      predDF[["pred"]] <- ifelse(predDF[[classes[LevelIndex]]] >=thr,classes[LevelIndex],classes[-LevelIndex])
      predDF[["pred"]] <- factor(predDF[["pred"]],levels=levels(predDF[["obs"]]))
      predDF
    }
    prob_test.2 <- thres_preds(prob_test,thres,predClasses)
    prob_test.2 <- prob_test.2 %>% droplevels()
  }
  
  prob_test.2 
}



# runMultiSurvival <- function(DF,var,varSelection,varSelected,survType,sampleID,varPlotName){
#   # DF <- immdata.apd1.proc.TcrA.divEst.fin
#   # var <-survVars[[4]]
#   # varSelection <-datasets.OS
#   # typeVar <- "cont"
#   # varSelected <-"dataset"
#   # survType <- "OS"
#   # sampleID <- "run_accession"
#   # varPlotName<-"Modeled Entropy"
#   
#   
#   # DF <- testDF.survival.prob.thres.split.l3[[2]]
#   # var <-survVars[[1]]
#   # varSelection <-datasets.OS
#   # # typeVar <- "cat"
#   # varSelected <-"dataset"
#   # survType <- "OS"
#   # sampleID <- "run_accession"
#   # varPlotName<-names(survVars)[2]
#   
#   cat('\n')  
#   cat('\n') 
#   cat("### ", varPlotName, "{.tabset .tabset-fade .tabset-pills} \n") 
#   cat('\n') 
#   
#   # The var will be the main covariate.
#   # res.cont<-performSurvivalCalcs(varSelection, DF,var, varSelected, survType, sampleVar=sampleID,covName=varPlotName,covType="cont",stratifyMethod="median",sLevelA = "high",sLevelB = "low",onlyForest=TRUE)
#   res.cat<-performSurvivalCalcs(varSelection, DF,var, varSelected, survType, sampleVar=sampleID,covName=varPlotName,covType="cat",stratifyMethod=NULL,sLevelA = "highProb",sLevelB = "lowProb")
#   
#   # cat('\n')  
#   # cat('\n') 
#   # cat("#### ", "Forest plot - continuous", "\n") 
#   # cat('\n') 
#   # 
#   # print(res.cont$forestplot)
#   
#   cat('\n')  
#   cat('\n') 
#   cat("#### ", "Forest plot - categorical", "\n") 
#   cat('\n') 
#   
#   print(res.cat$forestplot)
#   
#   cat('\n')  
#   cat('\n') 
#   cat("#### ", "KM plots", "\n") 
#   cat('\n') 
#   
#   print(res.cat$kmPlot)
#   res.cat
# }
# 
# ## WRAPPER WITH ALL SURV ANALYSIS ##

performSurvivalCalcs <- function(datasetsEvent, dataDF,mainVar, filterVar, survivalVar, sampleVar,covName,covType="cont",
                                 stratifyMethod=NULL,sLevelA = "high",sLevelB = "low",
                                 printPlots= FALSE,onlyForest=FALSE, onlyKM = FALSE){
  
  # DF <- testDF.survival.prob.thres.split.l3[[1]]
  # var <-survVars[[1]]
  # varSelection <-tissues.OS
  # # typeVar <- "cont"
  # varSelected <-"tissue"
  # survType <- "OS"
  # sampleID <- "run_accession"
  # varPlotName<-names(testDF.survival.prob.thres.split.l3)[1]
  # 
  # datasetsEvent <- varSelection
  # dataDF <- DF
  # mainVar <-var
  # filterVar <-varSelected
  # survivalVar <-survType
  # sampleVar=sampleID
  # covName=varPlotName
  # covType="cat"
  # stratifyMethod=NULL
  # sLevelA = "highProb"
  # sLevelB = "lowProb"
  # printPlots= FALSE

  

  # varSelection, DF,var, varSelected, survType, sampleVar=sampleID,covName=varPlotName,covType="cat",stratifyMethod="median",sLevelA = "highProb",sLevelB = "lowProb"
  
  # Transform some vars
  filterVar.un <- as.name(filterVar)
  
  ## Get number of points depending on filterVar
  # Empty df
  survivalPointsNumb <- data.frame(dataset=factor(),
                                   no_rows = numeric(),
                                   comboVar = character(),
                                   stringsAsFactors=FALSE)
  
  survivalPointsNumb <-get_tissueDataNumb(dataDF, survivalPointsNumb, sampleVar,paste0(survivalVar,".time"),mainVar,filterVar)
  
  ## Cox hazard ratio, HR calculations
  
  ######################
  ## As continuous var ##
  #########################
  survival_HR_logrankP <- data.frame(filterVar=character(),
                                     covariate=character(),
                                     HR= character(),
                                     HR_lower = character(),
                                     HR_upper = character(),
                                     logrankP = character(),
                                     stringsAsFactors=FALSE)
  
  colnames(survival_HR_logrankP)[which(colnames(survival_HR_logrankP)=="filterVar")] <- filterVar
  
  
  if(covType=="cont"){
    survival_HR_logrankP<- sapply(datasetsEvent,function(x) cox_hazard(dataDF, x,sampleVar, paste0(survivalVar,".time"),
                                                                       survivalVar,filterVar,c(mainVar),survival_HR_logrankP,
                                                                       cutPointsTable=NULL,signature=NULL,
                                                                       dichot = FALSE, stratMethod = stratifyMethod,stratVariable=NULL,StratLevel1=sLevelA,StratLevel2=sLevelB,
                                                                       globalLocal="local"),simplify = FALSE)#, "NK_cells"
  }else if(covType == "cat"){
    ##########################
    # Dichotomize by median ##
    ##########################
    survival_HR_logrankP<- sapply(datasetsEvent,function(x) cox_hazard(dataDF, x,sampleVar, paste0(survivalVar,".time"),
                                                                       survivalVar,filterVar,c(mainVar),survival_HR_logrankP,
                                                                       cutPointsTable=NULL,signature=NULL,
                                                                       dichot = TRUE, stratMethod = stratifyMethod,stratVariable=mainVar,StratLevel1=sLevelA,StratLevel2=sLevelB,
                                                                       globalLocal="local"),simplify = FALSE)#, "NK_cells"
  }
    
  # Merge data
  survival_HR_logrankP <- bind_rows(survival_HR_logrankP)
    ####################################
  # Process cox prop model results ##
  ###################################
  survival_HR_logrankP <- survival_HR_logrankP %>% dplyr::mutate_at(vars(HR,HR_lower,HR_upper),as.numeric)%>% dplyr::mutate_at(vars(logrankP),as.character)%>% dplyr::mutate_at(vars(logrankP),as.numeric) #%>% dplyr::mutate_at(vars(HR,HR_lower,HR_upper),levels)
  ## Cancers in which TCR richness is statistically significant(p<0.05) associated with good prognosis(HR <1, log10(HR)<0) are highlighted in green and significant associations with poor prognosis (HR>1,log10(HR)>0)are in indian red.
  survival_HR_logrankP$significant <- ifelse(survival_HR_logrankP$logrankP <= 0.05,TRUE,FALSE)
  survival_HR_logrankP$prognosis <- ifelse(survival_HR_logrankP$HR <1, "good",ifelse(survival_HR_logrankP$HR >1,"poor","noeffect"))
  survival_HR_logrankP$HR_eval <- ifelse(survival_HR_logrankP$significant==TRUE & survival_HR_logrankP$prognosis=="good","signGood",ifelse(survival_HR_logrankP$significant==TRUE & survival_HR_logrankP$prognosis=="poor","signPoor","notSign"))
  survival_HR_logrankP$HR_eval <- as.factor(survival_HR_logrankP$HR_eval)
  survival_HR_logrankP$logHR <- log10(survival_HR_logrankP$HR)
  survival_HR_logrankP$logHR_lower <- log10(survival_HR_logrankP$HR_lower)
  survival_HR_logrankP$logHR_upper <- log10(survival_HR_logrankP$HR_upper)
  
  # Merge with data points size
  survival_HR_logrankP2 <- merge(survival_HR_logrankP,survivalPointsNumb[,1:2],by.x = filterVar, by.y =filterVar, all.x=TRUE)
  
  
  survival_HR_logrankP2$HR_eval <-factor(survival_HR_logrankP2$HR_eval,
                                         levels = c("notSign","signGood","signPoor"))

  if(isTRUE(onlyForest)){
    # NOT logarithmic HR
    p <-survival_forest_plot(survival_HR_logrankP2,filterVar, predictorName=covName,facetVariable="covariate",facetsCols=2,logHR= FALSE)
    # print(p)
    res = list(hrDF=survival_HR_logrankP2,forestplot=p)
  }else if(isTRUE(onlyKM)){
    global_legend <-fit_all(as.character(datasetsEvent[1]), dataDF, survival_HR_logrankP2,
                            cutPointsTable = NULL,
                            covName = mainVar, stratMethod=stratifyMethod,signature=NULL,
                            stratVariable = mainVar,SampleCol = sampleVar,
                            paste0(survivalVar,".time"), survivalVar,facetVariable = filterVar,predictor = mainVar,
                            StratLevel1 = sLevelA,StratLevel2 = sLevelB, globalLocal="global",getLegend=TRUE)
    
    kmPlots = sapply(datasetsEvent,function(x) fit_all(as.character(x), dataDF, survival_HR_logrankP2,
                                                       cutPointsTable = NULL,
                                                       covName = mainVar, stratMethod=stratifyMethod,signature=NULL,
                                                       stratVariable = mainVar,SampleCol = sampleVar,
                                                       paste0(survivalVar,".time"), survivalVar,facetVariable = filterVar,predictor = mainVar,
                                                       StratLevel1 = sLevelA,StratLevel2 = sLevelB, globalLocal="local",covType=covType,getLegend=FALSE,HRAnnot = TRUE),simplify = FALSE)
    
    res = list(hrDF=survival_HR_logrankP2,kmPlot = kmPlots,legend=global_legend)
    
  }else{
    #######################################################
    # KM PLOTS: Kaplan-Meier estimator and log-rank test ##
    #######################################################
    
    n<-dim(survival_HR_logrankP2)[1]/2;
    graphCols <-ceiling(n/floor(sqrt(n)));
    graphRows<-floor(n)
    #
    
    global_legend <-fit_all(as.character(datasetsEvent[1]), dataDF, survival_HR_logrankP2,
                            cutPointsTable = NULL,
                            covName = mainVar, stratMethod=stratifyMethod,signature=NULL,
                            stratVariable = mainVar,SampleCol = "run_accession",
                            paste0(survivalVar,".time"), survivalVar,facetVariable = filterVar,predictor = mainVar,
                            StratLevel1 = sLevelA,StratLevel2 = sLevelB, globalLocal="global",getLegend=TRUE)
    
    
    kmPlots = sapply(datasetsEvent,function(x) fit_all(as.character(x), dataDF, survival_HR_logrankP2,
                                                       cutPointsTable = NULL,
                                                       covName = mainVar, stratMethod=stratifyMethod,signature=NULL,
                                                       stratVariable = mainVar,SampleCol = "run_accession",
                                                       paste0(survivalVar,".time"), survivalVar,facetVariable = filterVar,predictor = mainVar,
                                                       StratLevel1 = sLevelA,StratLevel2 = sLevelB, globalLocal="local",covType=covType,getLegend=FALSE,HRAnnot = TRUE),simplify = FALSE)
    
    p.km<- ggarrange(plotlist = kmPlots,ncol=graphCols, nrow =  graphRows,common.legend = TRUE, legend="bottom")
    p <-survival_forest_plot(survival_HR_logrankP2,filterVar, predictorName=covName,facetVariable="covariate",facetsCols=2,logHR= FALSE)
    res = list(hrDF=survival_HR_logrankP2,forestplot=p, kmPlot = p.km,legend=global_legend)
  }
    res
  
}




########################################################
## Get Number of patients of tissue for variable ##
########################################################
get_tissueDataNumb <- function(tcgaDF, numb_DataPoints_df, SampleCol,corVarx,corVary,facetVariable){
  # Prepare variables
  SampleCol <- as.name(SampleCol)
  corVarx.un <- as.name(corVarx)
  corVary.un <- as.name(corVary)
  facetVariable.un <- as.name(facetVariable)
  tcgaDF[[corVarx]] <- gsub(",","",tcgaDF[[corVarx]],fixed=TRUE)
  tcgaDF[[corVarx]] <- gsub("#N/A",NA,tcgaDF[[corVarx]],fixed=TRUE)
  tcgaDF[tcgaDF=="NA"] = NA
  tcgaDF[tcgaDF=="NA"] = NA
  tcgaDF.clean <- tcgaDF %>% dplyr::mutate(!!facetVariable:=factor(!!facetVariable.un)) %>% dplyr::filter(complete.cases(!!corVary.un, !!corVarx.un)) 
  # Find number of datapoints for each cancer type
  tissue_dataPoints <- tcgaDF.clean  %>% group_by(!!facetVariable.un) %>% 
    distinct(!!SampleCol,.keep_all = TRUE) %>%
    dplyr::filter(!any(is.na(!!corVarx.un))) %>%
    dplyr::summarize(no_rows=length(!!facetVariable.un)) %>% dplyr::mutate(comboVar=paste0(!!facetVariable.un, "-", corVarx))
  # tissue_dataPoints
  # tissues_fewData <- tissue_dataPoints %>% dplyr::filter(no_rows<10) %>% dplyr::select(!!facetVariable.un) %>% pull(!!facetVariable.un) %>% droplevels()
  tissue_dataPoints
  
  numb_DataPoints_df <- rbind(numb_DataPoints_df,tissue_dataPoints)
  numb_DataPoints_df
}


cox_hazard <- function(DF, selectedcancer, SampleCol, endopoint, event,facetVariable,predictor,hrDF,cutPointsTable=NULL,signature=NULL, 
                       dichot = FALSE, stratMethod = "median",stratVariable=NULL,StratLevel1=NULL,StratLevel2=NULL,globalLocal="global"){#dichotLabels=NULL,
  
  # hrDF<- survival_HR_logrankP.cat
  # DF<- dataDF
  # selectedcancer<- "Hugo"
  # SampleCol<- sampleVar
  # 
  # endopoint<- paste0(survivalVar,".time")
  # event<- survivalVar
  # facetVariable<- filterVar
  # predictor<-mainVar
  # cutPointsTable <-NULL
  # signature <- NULL
  # dichot = TRUE
  # globalLocal <- "global"
  # StratLevel1 <- "high"# The one defined as above median/cutoff set
  # StratLevel2 <-"low"# The one defined as below the median/cutoff set
  # stratVariable<- mainVar
  # 
  # stratMethod = "median"
  ###################################
  # DF =dataDF
  # selectedcancer <- datasetsEvent[1]
  # SampleCol <-sampleVar
  # endopoint <-paste0(survivalVar,".time")
  # event <- survivalVar
  # facetVariable<-filterVar
  # predictor<-c(mainVar)
  # hrDF<-survival_HR_logrankP
  # cutPointsTable=NULL
  # signature=NULL
  # dichot = TRUE
  # stratMethod = stratifyMethod
  # stratVariable=mainVar
  # StratLevel1=sLevelA
  # StratLevel2=sLevelB
  # globalLocal="local"


  # print(selectedcancer)
  # # 
  # Prepare variables
  SampleCol <- as.name(SampleCol)
  endopoint.un <- as.name(endopoint)
  event.un <- as.name(event)
  facetVariable.un <- as.name(facetVariable)
  predictor.un <- as.name(predictor)
  if(!is.null(stratVariable)){
    stratVariable.un <- as.name(stratVariable)
    
  }
  
  # Clean DF, make columns factors, remove NA rows for x-z variables and keep one row per sample ID
  ## Tumor
  DF.sub <- DF %>% dplyr::mutate(!!facetVariable:=factor(!!facetVariable.un)) %>% 
    dplyr::filter(!!facetVariable.un == selectedcancer) %>% 
    distinct(!!SampleCol,.keep_all = TRUE) %>% 
    dplyr::mutate_at(vars(!!endopoint.un,!!event.un),as.numeric)  %>% dplyr::filter(complete.cases(!!event.un,!!endopoint.un))
  #filter(complete.cases(!!endopoint.un,!!event.un,!!predictor.un))
  
  # 
  # if(isTRUE(dichot)){
  #   # DICHOTOMIZE-- HERE I HAVE TO CHECK IF TISSUE OR DATASET
  #   if(stratMethod=="median"){
  #     # Is it a bi-modal distribution?How to select the cutoff?
  #     cutglobal <- median(DF[[predictor]][!is.na(DF[[predictor]])], na.rm = TRUE)
  #     # Calculate median for cancer
  #     cutlocal <- median(DF.sub[[predictor]][!is.na(DF.sub[[predictor]])], na.rm = TRUE)
  #   }else if(stratMethod=="mean"){
  #     cutglobal <- mean(DF[[predictor]][!is.na(DF[[predictor]])], na.rm = TRUE)
  #     cutlocal <- mean(DF.sub[[predictor]][!is.na(DF.sub[[predictor]])], na.rm = TRUE)
  #   }else if(stratMethod=="cut"){
  #     cutglobal <-as.numeric(optCutpointContinuous(DF,endopoint,event,predictor,filterVar=NULL, filterVarSelected=NULL,gatherTable=TRUE)[2])
  #     cutlocal <- as.numeric(cutPointsTable[[selectedcancer]][2])
  #   }else if(stratMethod=="opt"){
  #     # check if dataset/tissue
  #     if (facetVariable=="tissue"){
  #       selectedTissue=selectedcancer
  #     }else if(facetVariable=="dataset"){
  #       # Figure out the type of tissue it is
  #       selectedTissue=as.character(unique(DF.sub[["tissue"]]))
  #     }
  #     
  #     if(signature=="TuTack"){
  #       cutglobal <-as.numeric(cutPointsTable[[selectedTissue]][1])
  #       cutlocal <-as.numeric(cutPointsTable[[selectedTissue]][1])
  #     }else if(signature=="TuTack.ifng"){
  #       cutglobal <-as.numeric(cutPointsTable[[selectedTissue]][2])
  #       cutlocal <-as.numeric(cutPointsTable[[selectedTissue]][2])
  #     }else if(signature=="GEP"){
  #       cutglobal <-as.numeric(cutPointsTable[[selectedTissue]][3])
  #       cutlocal <-as.numeric(cutPointsTable[[selectedTissue]][3])
  #     }
  #   }
  #   
  #   
  #   
  #   if(globalLocal=="global"){
  #     cutVar <- cutglobal
  #     #subtitlePlot<-""
  #   }else{
  #     cutVar <- cutlocal
  #     #subtitlePlot <-paste0("Richness Threshold (median): ", medianVar, "\n", "# Samples: ", nrow(tcgaDF_sign))
  #   }
  #   DF.sub2 <- DF.sub %>% dplyr::mutate(!!stratVariable.un := ifelse(!!predictor.un >= cutVar,StratLevel1,StratLevel2)) %>%mutate(!!stratVariable:=factor(!!stratVariable.un)) %>% dplyr::filter(complete.cases(!!event.un,!!endopoint.un))
  #   # DF.sub2[[stratVariable]] <- factor(DF.sub2[[stratVariable]],levels = c())
  # }else{
  #   stratVariable <- predictor
  #   DF.sub2 <-DF.sub
  # }
  # 
  stratVariable <- predictor
  DF.sub2 <-DF.sub
  
  DF.sub2[[predictor]] <- factor(DF.sub2[[predictor]], levels = c(StratLevel1,StratLevel2))
  # HERE I HAVE TO DEFINE THE LEVEL-looking at low or high
  if(isTRUE(dichot)){
    DF.sub2[[predictor]] = relevel(DF.sub2[[predictor]], ref = StratLevel2)
  }
  # tcgaDF_tp.clean
  # paste0("coxph(Surv(",endopoint,",",event,") ~ ",paste(predictor, sep="", collapse=" + "),", data=tcgaDF_tp.clean)")
  res.cox <-eval(parse(text=paste0("coxph(Surv(",endopoint,", ",event,") ~ ",paste(stratVariable, sep="", collapse=" + ") ,", data=DF.sub2)")))
  res.cox
  if(isTRUE(dichot)){
    # print("CATEGORICAL")
    if (all(is.na(res.cox$coefficients))){
      # HERE WE DEFINE WHICH LEVEL TO TAKE!!!!!
      cov.res <-sapply(stratVariable, function(x) getNACoxPHres(x,facetVariable,selectedcancer,varType="cat",varLabel=levels(DF.sub2[[stratVariable]])[2]),simplify = FALSE)
      # Merge data
      cov.res.merged  <-bind_rows(cov.res)
    }else{
      cov.res <-sapply(stratVariable, function(x) getCoxPHres(res.cox,x,facetVariable,selectedcancer,varType="cat",varLabel=levels(DF.sub2[[stratVariable]])[2]),simplify = FALSE)
      # Merge data
      # print("merge 1")
      cov.res.merged  <-bind_rows(cov.res)
      # cov.res
    } 
  }else{
    # print("CONTINUOUS")
    if (all(is.na(res.cox$coefficients))){
      cov.res <-sapply(stratVariable, function(x) getNACoxPHres(x,facetVariable,selectedcancer,varType="quant",varLabel=NULL),simplify = FALSE)
      
      # getNACoxPHres <- function(covar,variable,selectedcancer,varType="quant",varLabel=NULL)
      
      # Merge data
      cov.res.merged  <-bind_rows(cov.res)
    }else{
      cov.res <-sapply(stratVariable, function(x) getCoxPHres(res.cox,x,facetVariable,selectedcancer,varType="quant",varLabel=NULL),simplify = FALSE)
      # Merge data
      # print("merge 2")
      cov.res.merged  <-bind_rows(cov.res)
    }
  }
  
  hrDF <- rbind(hrDF,cov.res.merged)
  hrDF
}

survival_forest_plot <-function(survivalDF,selectedCancer,predictorName=NULL,facetVariable=NULL,facetsCols=NULL,logHR= FALSE){
  # survivalDF<-survival_HR_logrankP.s8.sub
  # selectedCancer<-"tisData"
  # predictorName="TuTack score"
  # facetVariable="varLabel"
  # facetsCols=2
  # logHR= FALSE
  
  
  if (logHR==TRUE){
    
    if(selectedCancer=="tissue"){
      p <- ggplot(survivalDF, aes(x=tissue, y=logHR, ymin=logHR_lower, ymax=logHR_upper))
    }else if(selectedCancer=="dataset"){
      p <- ggplot(survivalDF, aes(x=dataset, y=logHR, ymin=logHR_lower, ymax=logHR_upper))
    }else if(selectedCancer=="tisData"){
      p <- ggplot(survivalDF, aes(x=tisData, y=logHR, ymin=logHR_lower, ymax=logHR_upper)) 
    }
    
    # geom_linerange(size=8, colour="#a6d8f0") +
    p <- p + geom_hline(aes(yintercept=0), lty=2)
    p <- p + geom_pointrange(data=survivalDF,aes(color=HR_eval,size = no_rows),shape = 15) #,fatten = 1
    p <- p + scale_size_continuous(range = c(0, 1))
    #geom_point(aes_string(size = "no_rows",color="HR_eval"),shape = 15, show.legend = TRUE, alpha = .4) +
    #geom_errorbar(aes(ymin=HR_lower, ymax=HR_upper,color=HR_eval, width=0.5,cex=1)) + 
    # scale_color_manual("Significant HR", breaks=c("signGood", "notSign", "notSign"),values=c("black", "Medium sea green", "Indian red"))+
    p <- p + scale_color_manual(values = c('black', "Medium sea green", "Indian red")) 
    #scale_y_continuous(limits = c(0.5, 2)) +
    #scale_y_log10(breaks=c(0.5,1,2),position="top",limits=c(0.5,2)) +
    p <- p + guides(col = guide_legend(reverse = TRUE))
    p <- p + coord_flip()
    
    p <- p + ggtitle("")
    p <- p + xlab('')
    # p <- p + xlim(round(min(survivalDF[["HR_lower"]]),digits=2), round(min(survivalDF[["HR_upper"]]),digits=2))
    
    if(!is.null(facetVariable)){
      p <- p + facet_wrap(facetVariable,ncol=facetsCols,scales='free')
      
    }
    if(!is.null(predictorName)){
      p <- p + ylab(paste0("log Hazard Ratio (95% Confidence Interval) \nassociated with unit increase in ",predictorName))
      # p <- p + ylab(paste0("log Hazard Ratio (95% Confidence Interval) \nassociated with TuTack positive score"))
    }else{
      p <- p + ylab(paste0("log Hazard Ratio (95% Confidence Interval) \nassociated with unit increase in covariate"))
      
    }
    p <- p + theme_ipsum(base_family = "Palatino Linotype", base_size = 13)
    p <- p + theme(axis.text.x = element_text(size = 13),#angle = 90,#hjust = 1, vjust = 1,
                   axis.title.x = element_text(size=13),
                   axis.text.y = element_text( size = 13),#hjust = 1, vjust = 1,
                   axis.title.y = element_text(size=13),
                   legend.text = element_text(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   axis.line = element_line(colour = "black"),
                   legend.title = element_blank(),
                   legend.position = "none"
    )
  }else{
    if(selectedCancer=="tissue"){
      p <- ggplot(survivalDF, aes(x=tissue, y=HR, ymin=HR_lower, ymax=HR_upper)) 
    }else if(selectedCancer=="dataset"){
      p <- ggplot(survivalDF, aes(x=dataset, y=HR, ymin=HR_lower, ymax=HR_upper)) 
    }else if(selectedCancer=="tisData"){
      p <- ggplot(survivalDF, aes(x=tisData, y=HR, ymin=HR_lower, ymax=HR_upper)) 
    }
    p <- p + ylim(floor(min(survivalDF[["HR_lower"]])), ceiling(max(survivalDF[["HR_upper"]])))
    # geom_linerange(size=8, colour="#a6d8f0") +
    p <- p + geom_hline(aes(yintercept=1), lty=2)
    p <- p + geom_pointrange(data=survivalDF,aes(color=HR_eval,size = no_rows),shape = 15) #,fatten = 1
    p <- p + scale_size_continuous(range = c(0, 1))
    myColors <- c('black', "Medium sea green", "Indian red")
    names(myColors) <- levels(survivalDF$HR_eval)
    
    #geom_point(aes_string(size = "no_rows",color="HR_eval"),shape = 15, show.legend = TRUE, alpha = .4) +
    #geom_errorbar(aes(ymin=HR_lower, ymax=HR_upper,color=HR_eval, width=0.5,cex=1)) + 
    # scale_color_manual("Significant HR", breaks=c("signGood", "notSign", "notSign"),values=c("black", "Medium sea green", "Indian red"))+
    p <- p + scale_color_manual(name = "HR_eval",values = myColors)
    # p + scale_color_manual(values = c('black', "Medium sea green", "Indian red")) 
    #scale_y_continuous(limits = c(0.5, 2)) +
    # p <- p + ylim(round(min(survivalDF[["HR_lower"]]),digits=2), round(min(survivalDF[["HR_upper"]]),digits=2))
    #scale_y_continuous(limits = c(0.5, 2)) +
    #scale_y_log10(breaks=c(0.5,1,2),position="top",limits=c(0.5,2)) +
    p <- p + guides(col = guide_legend(reverse = TRUE))
    p <- p + coord_flip()
    
    p <- p + ggtitle("")
    p <- p + xlab('')
    
    if(!is.null(facetVariable)){
      p <- p + facet_wrap(facetVariable,ncol=facetsCols,scales='free')
      
    }
    # if(!is.null(predictorName)){
    #   p <- p + ylab(paste0("Hazard Ratio (95% Confidence Interval) \nassociated with unit increase in ",predictorName))
    #   # p <- p + ylab(paste0("Hazard Ratio (95% Confidence Interval) \nassociated with TuTack/GEP positive score"))
    # }else{
    #   p <- p + ylab(paste0("Hazard Ratio (95% Confidence Interval) \nassociated with unit increase in covariate"))
    #   
    # }
    
    p <- p + ylab(paste0("Hazard Ratio (95% Confidence Interval)"))
    
    p <- p + theme_ipsum(base_family = "Palatino Linotype", base_size = 13)
    p <- p + theme(axis.text.x = element_text(size = 13),#angle = 90,#hjust = 1, vjust = 1,
                   axis.title.x = element_text(size=13),
                   axis.text.y = element_text( size = 13),#hjust = 1, vjust = 1,
                   axis.title.y = element_text(size=13),
                   legend.text = element_text(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   axis.line = element_line(colour = "black"),
                   legend.title = element_blank(),
                   legend.position = "none"
    )
  }
  
  #p <- p + scale_y_discrete()
  p
}

fit_all <- function(selectedcancer,tcgaDF,survHRdf,cutPointsTable,covName,
                    stratMethod,signature,stratVariable,SampleCol,endopoint, event,
                    facetVariable,predictor,StratLevel1,StratLevel2,
                    globalLocal,covType="cont",getLegend=FALSE,HRAnnot=FALSE){

  # selectedcancer <-datasets.OS[1]
  # tcgaDF <-immdata.apd1.proc.TcrB.divEst.fin
  # survHRdf <-r.datasets$hrDF
  # cutPointsTable <- NULL
  # covName <-"modEntr.ds"
  # stratMethod <-"median"
  # signature <-NULL
  # stratVariable <-"modEntr.ds"
  # SampleCol <-"run_accession"
  # endopoint <-"OS.time"
  # event <-"OS"
  # facetVariable <- "dataset"
  # predictor <-"modEntr.ds"
  # StratLevel1 <-"high"
  # StratLevel2 <-"low"
  # globalLocal <-"local"
  # covType="cont"
  # getLegend=FALSE
  # HRAnnot=FALSE


  ##########################################
  # selectedcancer <-datasets.OS[1]
  # tcgaDF<-dataDF
  # survHRdf<-survival_HR_logrankP2
  # cutPointsTable = NULL
  # covName = mainVar
  # stratMethod=stratifyMethod
  # signature=NULL
  # stratVariable = mainVar
  # SampleCol = "run_accession"
  # paste0(survivalVar,".time")
  # event<-survivalVar
  # facetVariable = filterVar
  # predictor = mainVar
  # StratLevel1 = sLevelA
  # StratLevel2 = sLevelB
  # globalLocal="local"
  # covType=covType
  # getLegend=FALSE

  # selectedcancer <-datasets.OS[1]
  # tcgaDF <- dataDF
  # survHRdf <-survival_HR_logrankP2
  # cutPointsTable <-NULL
  # covName <-mainVar
  # stratMethod <-stratifyMethod
  # signature <- NULL
  # stratVariable <-mainVar
  # SampleCol <-"run_accession"
  # endopoint <-paste0(survivalVar,".time")
  # event <-survivalVar
  # facetVariable <- filterVar
  # predictor <-mainVar
  # StratLevel1 = sLevelA
  # StratLevel2 = sLevelB
  # globalLocal="local"
  # covType=covType
  # getLegend=FALSE
  # HRAnnot = TRUE
  #

  # selectedcancer <-as.character(datasetsEvent[1])
  # tcgaDF <- dataDF
  # survHRdf <-survival_HR_logrankP2
  # cutPointsTable <-NULL
  # covName <-mainVar
  # stratMethod <-stratifyMethod
  # signature <- NULL
  # stratVariable <-mainVar
  # SampleCol <-"run_accession"
  # endopoint <-paste0(survivalVar,".time")
  # event <-survivalVar
  # facetVariable <- filterVar
  # predictor <-mainVar
  # StratLevel1 = sLevelA
  # StratLevel2 = sLevelB
  # globalLocal="local"
  # covType=covType
  # getLegend=FALSE
  # HRAnnot = TRUE
  #


  ##########################################

  # print(selectedcancer)

  # Prepare variables
  SampleCol <- as.name(SampleCol)
  endopoint.un <- as.name(endopoint)
  event.un <- as.name(event)
  facetVariable.un <- as.name(facetVariable)
  predictor.un <- as.name(predictor)
  stratVariable.un <- as.name(stratVariable)
  # Survival plots
  # First subset data to the cancers with significant associations
  tcgaDF_sign <- filter(tcgaDF, !!facetVariable.un %in% dplyr::selectedcancer) %>% droplevels()
  tcgaDF_sign <- tcgaDF_sign %>% distinct(!!SampleCol,.keep_all = TRUE)
  if (selectedcancer %in% c("bladder","melanoma","RCC","gastric")){
    tissueSelected <- ''
  }else{
    tissueSelected <- tcgaDF_sign %>% pull(tissue) %>% unique()
    tissueSelected <- paste0("-",tissueSelected)
  }
  # tissueSelected <- tcgaDF_sign %>% pull(tissue) %>% unique()
  # # DICHOTOMIZE-- HERE I HAVE TO CHECK IF TISSUE OR DATASET
  # if(stratMethod=="median"){
  #   # Is it a bi-modal distribution?How to select the cutoff?
  #   cutglobal <- median(tcgaDF[[predictor]][!is.na(tcgaDF[[predictor]])],na.rm = TRUE)
  #   # Calculate median for cancer
  #   cutlocal <- median(tcgaDF_sign[[predictor]][!is.na(tcgaDF_sign[[predictor]])],na.rm = TRUE)
  # }else if(stratMethod=="mean"){
  #   cutglobal <- mean(tcgaDF[[predictor]][!is.na(tcgaDF[[predictor]])],na.rm = TRUE)
  #   cutlocal <- mean(tcgaDF_sign[[predictor]][!is.na(tcgaDF_sign[[predictor]])],na.rm = TRUE)
  # }else if(stratMethod=="cut"){
  #   cutglobal <-as.numeric(optCutpointContinuous(tcgaDF,endopoint,event,predictor,filterVar=NULL, filterVarSelected=NULL,gatherTable=TRUE)[2])
  #   cutlocal <- as.numeric(cutPointsTable[[selectedcancer]][2])
  # }else if(stratMethod=="opt"){
  #   # check if dataset/tissue
  #   if (facetVariable=="tissue"){
  #     selectedTissue=selectedcancer
  #   }else if(facetVariable=="dataset"){
  #     # Figure out the type of tissue it is
  #     selectedTissue=unique(tcgaDF_sign[["tissue"]])
  #   }
  #
  #   if(signature=="TuTack"){
  #     cutglobal <-as.numeric(cutPointsTable[[selectedTissue]][1])
  #     cutlocal <-as.numeric(cutPointsTable[[selectedTissue]][1])
  #   }else if(signature=="TuTack.ifng"){
  #     cutglobal <-as.numeric(cutPointsTable[[selectedTissue]][2])
  #     cutlocal <-as.numeric(cutPointsTable[[selectedTissue]][2])
  #   }else if(signature=="GEP"){
  #     cutglobal <-as.numeric(cutPointsTable[[selectedTissue]][3])
  #     cutlocal <-as.numeric(cutPointsTable[[selectedTissue]][3])
  #   }
  # }
  #


  # if(globalLocal=="global"){
  #   cutVar <- cutglobal
  #   #subtitlePlot<-""
  # }else{
  #   cutVar <- cutlocal
  #   #subtitlePlot <-paste0("Richness Threshold (median): ", medianVar, "\n", "# Samples: ", nrow(tcgaDF_sign))
  # }
  # tcgaDF_sign2 <- tcgaDF_sign %>% dplyr::mutate(!!stratVariable.un := ifelse(!!predictor.un >= cutVar,StratLevel1,StratLevel2)) %>%mutate(!!stratVariable:=factor(!!stratVariable.un))

  tcgaDF_sign2 <- tcgaDF_sign
  # Cleanup df
  tcgaDF_sign2[[endopoint]] <- gsub(",","",tcgaDF_sign2[[endopoint]],fixed=TRUE)
  tcgaDF_sign2[[endopoint]] <- gsub("#N/A",NA,tcgaDF_sign2[[endopoint]],fixed=TRUE)
  tcgaDF_sign2[tcgaDF_sign2=="NA"] = NA
  tcgaDF_sign2[tcgaDF_sign2=="NA"] = NA

  tcgaDF_signFinal <- tcgaDF_sign2 %>%
    mutate(!!facetVariable:=factor(!!facetVariable.un)) %>%
    filter(complete.cases(!!endopoint.un,!!event.un,!!predictor.un,!!stratVariable.un)) %>%
    filter(!!facetVariable.un == selectedcancer) %>%
    distinct(!!SampleCol,.keep_all = TRUE) %>%
    mutate_at(vars(!!endopoint.un,!!event.un),as.numeric)

  if(globalLocal=="global"){
    #medianVar <- median_global
    subtitlePlot<-""
  }else{
    #medianVar <- median_local
    #subtitlePlot <-paste0("Richness Threshold (median): ", medianVar, "\n", "# Samples: ", nrow(tcgaDF_signFinal))
    subtitlePlot<-""
  }
  #tcgaDF_signFinal
  attach(tcgaDF_signFinal)
  if (event=="OS"){
    ylabel = "Survival probability"
  }else if(event=="PFS"){
    ylabel = "Progression-free Survival probability"
  }
  # surv_object <- Surv(endopoint.un/30.5, event.un)
  # surv_object
  #fits<-eval(parse(text=paste0("survfit(Surv(OS.time/30.5, OS) ~ ",stratVariable,", data = tcgaDF_signFinal)")))
  # Here I turn weeks to months!
  fits<-eval(parse(text=paste0("survfit(Surv(",endopoint,", ",event,") ~ ",stratVariable,", data = tcgaDF_signFinal)")))
  survpval <-eval(parse(text=paste0("surv_pvalue(fits",", data = tcgaDF_signFinal)")))
  survPlot <- ggsurvplot(
    fits,
    data = tcgaDF_signFinal,
    # title = paste0(selectedcancer," (n = ",nrow(tcgaDF_signFinal) ,")"),
    title = paste0(selectedcancer,"\n(n.",StratLevel1," = ",table(tcgaDF_signFinal[[predictor]])[2] ,", n.",StratLevel2," = ", table(tcgaDF_signFinal[[predictor]])[1],")"),
    subtitle = subtitlePlot,
    legend = "none",
    pval = FALSE,
    pval.method = FALSE,
    #pval.size = 3,
    #pval.coord = c(200, 0.10),
    conf.int = FALSE,
    xlab = "Time in months",
    ylab = ylabel,
    ggtheme = theme_minimal(base_family = "Palatino Linotype", base_size = 13),
    surv.median.line = "hv",
    size = 0.2,
    #legend.labs = c("Female","Male"),
    #legend.title = "",
    palette = c("#8C3F4D","#3E606F"),
    font.main=13,
    font.legend = 10
  )
  fits
  if (getLegend==TRUE){
    survPlot <- ggsurvplot(
      fits,
      data = tcgaDF_signFinal,
      title = paste0(selectedcancer),
      subtitle = subtitlePlot,
      legend = "none",
      pval = FALSE,
      #pval.size = 3,
      #pval.method = FALSE,
      #pval.coord = c(200, 0.10),
      conf.int = FALSE,
      xlab = "Time in months",
      ggtheme = theme_minimal(base_family = "Palatino Linotype", base_size = 13),
      surv.median.line = "hv",
      #legend.labs = c("Female","Male"),
      #legend.title = "",
      palette = c("#8C3F4D","#3E606F"),
      font.main=13,
      font.legend = 10
    )
    survPlot <- survPlot + labs(color = "Biomarker score")
    legend <- get_legend(survPlot$plot + theme(legend.position = "top"))
    legend
  }else{


    # res.cox <-eval(parse(text=paste0("coxph(Surv(",endopoint,",",event,") ~ ",predictor,", data=tcgaDF_signFinal)")))
    # cox.sum <- summary(res.cox)
    # hr <- cox.sum$conf.int[1]
    # ci <- cox.sum$conf.int[3:4]


    # hr <-survHRdf %>% dplyr::filter(!!facetVariable.un==selectedcancer) %>% dplyr::filter(covariate==paste0(covName,StratLevel2)) %>% pull(HR)
    # ci.l <-survHRdf %>% dplyr::filter(!!facetVariable.un==selectedcancer) %>% dplyr::filter(covariate==paste0(covName,StratLevel2)) %>% pull(HR_lower)
    # ci.u <-survHRdf %>% dplyr::filter(!!facetVariable.un==selectedcancer) %>% dplyr::filter(covariate==paste0(covName,StratLevel2)) %>% pull(HR_upper)
    # ci <- c(ci.l,ci.u)
    # ## logrankp <- round(survpval$pval,5)
    # # logrankp <- formatC(survpval$pval,format = "e", digits = 3)
    # logrankp <- format(survHRdf %>% dplyr::filter(!!facetVariable.un==selectedcancer) %>% dplyr::filter(covariate==paste0(covName,StratLevel2)) %>% pull(logrankP),format = "e", digits = 3)
    #
    if(covType=="cont"){
      hr <-survHRdf %>% dplyr::filter(!!facetVariable.un==selectedcancer) %>% dplyr::filter(covariate==covName) %>% pull(HR)
      ci.l <-survHRdf %>% dplyr::filter(!!facetVariable.un==selectedcancer) %>% dplyr::filter(covariate==covName) %>% pull(HR_lower)
      ci.u <-survHRdf %>% dplyr::filter(!!facetVariable.un==selectedcancer) %>% dplyr::filter(covariate==covName) %>% pull(HR_upper)
      ci <- c(ci.l,ci.u)
      ## logrankp <- round(survpval$pval,5)
      # logrankp <- formatC(survpval$pval,format = "e", digits = 3)
      logrankp <- format(survHRdf %>% dplyr::filter(!!facetVariable.un==selectedcancer) %>% dplyr::filter(covariate==covName) %>% pull(logrankP),format = "e", digits = 3)
    }else if(covType=="cat"){
      # HERE WE DEFINE BASED ON WHICH LEVEL OF COVARIATE IS SET AS REFERENCE_NEED TO CHANGE IF NECESSARY, , if ref level is strat level 2(low) we show results for strat level 1 (high)
      hr <-survHRdf %>% dplyr::filter(!!facetVariable.un==selectedcancer) %>% dplyr::filter(covariate==paste0(covName,StratLevel1)) %>% pull(HR)
      ci.l <-survHRdf %>% dplyr::filter(!!facetVariable.un==selectedcancer) %>% dplyr::filter(covariate==paste0(covName,StratLevel1)) %>% pull(HR_lower)
      ci.u <-survHRdf %>% dplyr::filter(!!facetVariable.un==selectedcancer) %>% dplyr::filter(covariate==paste0(covName,StratLevel1)) %>% pull(HR_upper)
      ci <- c(ci.l,ci.u)
      ## logrankp <- round(survpval$pval,5)
      # logrankp <- formatC(survpval$pval,format = "e", digits = 3)
      logrankp <- format(survHRdf %>% dplyr::filter(!!facetVariable.un==selectedcancer) %>% dplyr::filter(covariate==paste0(covName,StratLevel1)) %>% pull(logrankP),format = "e", digits = 3)
    }

    annotations <- data.frame(
      # xpos = c(-Inf,-Inf,Inf,Inf),
      # ypos =  c(-Inf, Inf,-Inf,Inf),
      # annotateText = c("Text","tExt","teXt",paste0("HR = ",round(cox$hazard_ratio,2),"(",round(cox$confidence_intervals[1],2),"-",round(cox$confidence_intervals[2],2),")\n","logrankP = ", round(cox$logrankP[1],2))),
      # hjustvar = c(0,0,1,1) ,
      # vjustvar = c(0,1.0,0,1))
      xpos = c(-Inf),
      ypos =  c(-Inf),
      annotateText = c(paste0("HR = ",round(hr,2),"(",round(ci[1],4),"-",round(ci[2],4),")\n","logrankP = ", logrankp)),
      hjustvar = c(0) ,
      vjustvar = c(0))
    #cox$hazard_ratio
    if(isTRUE(HRAnnot)){
      survPlot$plot <- survPlot$plot+
        ggplot2::geom_label(data = annotations, aes(x=xpos,y=ypos,hjust=hjustvar,
                                                    vjust=vjustvar,label=annotateText),family="Palatino Linotype",size=3.5)
    }else{
      annotations <- data.frame(
        # xpos = c(-Inf,-Inf,Inf,Inf),
        # ypos =  c(-Inf, Inf,-Inf,Inf),
        # annotateText = c("Text","tExt","teXt",paste0("HR = ",round(cox$hazard_ratio,2),"(",round(cox$confidence_intervals[1],2),"-",round(cox$confidence_intervals[2],2),")\n","logrankP = ", round(cox$logrankP[1],2))),
        # hjustvar = c(0,0,1,1) ,
        # vjustvar = c(0,1.0,0,1))
        xpos = c(-Inf),
        ypos =  c(-Inf),
        annotateText = c(paste0("logrankP = ", logrankp)),
        hjustvar = c(0) ,
        vjustvar = c(0))
      survPlot$plot <- survPlot$plot+
        ggplot2::geom_label(data = annotations, aes(x=xpos,y=ypos,hjust=hjustvar,
                                                    vjust=vjustvar,label=annotateText),family="Palatino Linotype",size=3.5)
    }

    survPlot$plot
  }

}


getCoxPHres <- function(coxObj,covar,variable,selectedcancer,varType="quant",varLabel=NULL){
  # coxObj <-res.cox
  # covar<-stratVariable
  # variable<-facetVariable
  # selectedcancer<-selectedcancer
  # varType <-"cat"  #cat,quant
  # varLabel <-StratLevel2
  # 
  coxObjDf <-as.data.frame(t(as.data.frame(coxObj$coefficients)))
  if (varType=="cat"){
    # print("CATEGORICAL")
    covar <- paste0(covar,varLabel)
    if(is.na(coxObjDf[[covar]])){
      cox_res <- getNACoxPHres(covar,variable,selectedcancer)
      cox_res
    }else{
      coxSumObj<-summary(coxObj)
      coxCIdf <- as.data.frame(coxSumObj$conf.int)
      coxCIdf <-coxCIdf %>% rownames_to_column(var="covariate")
      coxCIdf.cov <- coxCIdf %>% dplyr::filter(covariate==covar)
      hr <- coxCIdf.cov[2]
      ci <- coxCIdf.cov[4:5]
      logrankp <- coxSumObj$logtest[3]
      # logrankp <- round(coxSumObj$sctest[3], digits = 9)
      logrankp <- formatC(coxSumObj$sctest[3], format = "e", digits = 5)
      cox_res <- list("hazard_ratio"=hr,"confidence_intervals"=ci, "logrankP"=logrankp)
      cox_res <- as.data.frame(t(c(selectedcancer,covar, unlist(cox_res))))
      colnames(cox_res) <- c(variable,"covariate","HR","HR_lower","HR_upper","logrankP")
      cox_res
      
    }
  }else{
    # print("CONTINUOUS")
    if(is.na(coxObjDf[[covar]])){
      cox_res <- getNACoxPHres(covar,selectedcancer)
      cox_res
    }else{
      coxSumObj<-summary(coxObj)
      coxCIdf <- as.data.frame(coxSumObj$conf.int)
      coxCIdf <-coxCIdf %>% rownames_to_column(var="covariate")
      coxCIdf.cov <- coxCIdf %>% dplyr::filter(covariate==covar)
      hr <- coxCIdf.cov[2]
      ci <- coxCIdf.cov[4:5]
      logrankp <- coxSumObj$logtest[3]
      # logrankp <- round(coxSumObj$sctest[3], digits = 9)
      logrankp <- formatC(coxSumObj$sctest[3], format = "e", digits = 5)
      cox_res <- list("hazard_ratio"=hr,"confidence_intervals"=ci, "logrankP"=logrankp)
      cox_res <- as.data.frame(t(c(selectedcancer,covar, unlist(cox_res))))
      colnames(cox_res) <- c(variable,"covariate","HR","HR_lower","HR_upper","logrankP")
      cox_res
      
    }
  }
  
}


getNACoxPHres <- function(covar,variable,selectedcancer,varType="quant",varLabel=NULL){
  if (varType=="cat"){
    # print("CATEGORICAL")
    covar<- paste0(covar,varLabel)
  }
  cox_res <- as.data.frame(t(c(selectedcancer,covar,NA,NA,NA,NA)))
  colnames(cox_res) <- c(variable,"covariate","HR","HR_lower","HR_upper","logrankP")
  cox_res
}


boxplotsSection <- function(modelList,dataSet,groupVar,xVar,sigVar,yVar,labelVar,xLabel,yLabel, fillLabel, shapeLabel, title,testX,testY,
                            modelName,classes,LevelIndex,
                            thres=NULL,extractedProbDF=NULL,rfeE=FALSE,
                            mixedMode = FALSE,
                            trainModelsVec = NULL,
                            rfeModelsVec = NULL, Plotmode="pancancer"){
  #####################################
  # modelList<-modelList.final
  # dataSet<-testDF[[1]]
  # groupVar<-"response"
  # xVar <- "run_accession"
  # sigVar<-sigs[3]
  # yVar <- "CRPR"
  # labelVar <- "tissue"
  # xLabel<-"Patients"
  # yLabel <- sigs_axisNames[3]
  # fillLabel<-"Response to anti-PD-1"
  # shapeLabel<-"tissue"
  # title<-sigs_titleNames[3]
  # testX <- data.testX
  # testY <- testDF[[1]]$response
  # modelName <- names(modelList.final)[3]
  # classes <- c("CRPR","PD")
  # LevelIndex <- 1
  # thres=models.optThres.rfe[[3]]
  # posClassName <- "CRPR"
  # rfeE=FALSE
  # extractedProbDF=NULL
  # mixedMode = TRUE
  # trainModelsVec = c(1:c(length(models.final)-1))#c(1:7)
  # rfeModelsVec = c(length(models.final))#c(8)
  ##############################
  
  ######################################
  cat('\n')  
  cat('\n') 
  cat("###### ", modelName, " \n") 
  cat('\n') 
  
  xVar.un <- as.name(xVar)
  yVar.un <- as.name(yVar)
  sigVar.un<- as.name(sigVar)
  labelVar.un <- as.name(labelVar)
  
  if(isTRUE(rfeE)){
    # b <-sapply(1: length(modelList), function(x) predict(modelList[[x]], newdata=test.data), simplify = FALSE)
    probabilities <- predict(modelList, newdata=dataSet, type="prob")
    # probabilities.df <- bind_rows(probabilities, .id = "object")
    # predicted.classes <- predict(modelList, newdata=test.data, type="raw")
    
    probabilities <- lapply(probabilities, function(x) transform(x, obs = factor(dataSet[[groupVar]], levels = classes)))
    probabilities <- lapply(probabilities, function(x) transform(x, pred = factor(pred, levels = classes)))
    probabilities <- lapply(probabilities, function(x) transform(x, run_accession = factor(dataSet[[xVar]])))
    probabilities <- lapply(probabilities, function(x) transform(x, tissue = as.character(dataSet[[labelVar]])))
    predClassesProbs <-bind_rows(probabilities, .id = "object")
    prob_test <-predClassesProbs
    # predClassesProbs <- cbind(predClassesProbs,probabilities)
    
  }else if(isTRUE(mixedMode)){
    ## First train objects
    if( !is.null(testX) && ! is.null(testY)){
      prob_test.t<-extractProb(modelList[trainModelsVec],unkX =testX)
      prob_test.t$obs <-dataSet[[groupVar]]
      prob_test.t[[xVar]] <-dataSet[[xVar]]
      prob_test.t[[labelVar]] <-dataSet[[labelVar]]
    }else{
      prob_test.t<-extractProb(modelList)    
    }
    # Bring those in the format of  train models
    prob_test.t <-prob_test.t %>% dplyr::select(all_of(classes),"obs", "pred","object",xVar,labelVar)
    
    ## Then RFE objects
    probabilities <- predict(modelList[rfeModelsVec], newdata=dataSet, type="prob")
    # probabilities.df <- bind_rows(probabilities, .id = "object")
    # predicted.classes <- predict(modelList, newdata=test.data, type="raw")
    
    probabilities <- lapply(probabilities, function(x) transform(x, obs = factor(dataSet[[groupVar]], levels = classes)))
    probabilities <- lapply(probabilities, function(x) transform(x, pred = factor(pred, levels = classes)))
    probabilities <- lapply(probabilities, function(x) transform(x, run_accession = factor(dataSet[[xVar]])))
    probabilities <- lapply(probabilities, function(x) transform(x, tissue = as.character(dataSet[[labelVar]])))
    predClassesProbs <-bind_rows(probabilities, .id = "object")
    prob_test.r <-predClassesProbs
    # Bring those in the format of  train models
    prob_test.r <-prob_test.r %>% dplyr::select(all_of(classes),"obs", "pred","object",xVar,labelVar)
    ## Merge objects
    prob_test <- rbind(prob_test.t, prob_test.r)
    
  }else{
    if( !is.null(testX) && ! is.null(testY)){
      prob_test<-extractProb(modelList,unkX =testX)
      prob_test$obs <-dataSet[[groupVar]]
      prob_test[[groupVar]] <- dataSet[[groupVar]]
      prob_test[[xVar]] <- dataSet[[xVar]]
      prob_test[[labelVar]] <- dataSet[[labelVar]]
    }else{
      prob_test<-extractProb(modelList)    
    }
    
  }
  #####################
  prob_test[["obs"]] <- factor(prob_test[["obs"]],levels = classes)
  
  
  
  if( !is.null(thres) ){
    # Change predicted based on threshold
    thres_preds <- function(predDF,thr,classes){
      # predDF <-prob_test
      # thr <-thres
      # classes <-classes
      
      predDF[["pred"]] <- ifelse(predDF[[classes[LevelIndex]]] >=thr,classes[LevelIndex],classes[-LevelIndex])
      predDF[["pred"]] <- factor(predDF[["pred"]],levels=levels(predDF[["obs"]]))
      predDF
    }
    prob_test <- thres_preds(prob_test,thres,classes)
  }
  
  
  # Select model
  prob_test.final <- prob_test %>% dplyr::filter(object %in% c(modelName)) %>% droplevels()
  # prob_test.final.m <- melt(prob_test.final)
  if(Plotmode=="pancancer"){
    p<-ggstatsplot::ggbetweenstats(
      data = prob_test.final,
      x = obs,
      y = CRPR,
      sort = "descending", 
      pairwise.comparisons = TRUE, # show pairwise comparison test results
      pairwise.annotation ="asterisk",
      ggtheme = hrbrthemes::theme_ipsum_tw(),
      #outlier.tagging = TRUE,
      #outlier.label = education,
      ggstatsplot.layer = FALSE,
      messages = FALSE,
      xlab = "Observed Response",
      ylab = "Predicted Probability by model",
      mean.label.args = list(size = 5),
      bf.message = FALSE#"#CD5C5C", "#004D40"
    ) + ggplot2::scale_color_manual(values = c("#004D40","#CD5C5C" ))
    p <- p + coord_cartesian(ylim = c(0,1))
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
    print(p)
  }else{
    p <-ggboxplot(prob_test.final, x = "obs", y = "CRPR",notch=TRUE,add=c("jitter"),
                  color = "obs",palette = "uchicago",
                  order = classes,
                  ylab = "Predicted Probability by model", xlab = "Observed Response",facet.by = "tissue")
    
    #  Add p-value
    #p <- p + stat_compare_means(method = "t.test", label.x = 1.3, label.y = 0)
    p <- p + ylim(0,1)
    p <-p + stat_compare_means(label = "p.format",label.x.npc = 'center',label.y.npc = 'top', method="t.test")#, label.x = 1.3, label.y = 1.8
    p <- p + labs(fill = "Observed Response", color="Observed Response") + ggplot2::scale_color_manual(values = c("#CD5C5C","#004D40" ))
    
    p <-p + theme_bw(base_family = "Palatino Linotype", base_size = 12)
    p <- p + theme(legend.position = "top",
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.text.x = element_text(size=12),
                   axis.title.x= element_text(size=12),
                   #axis.text.x = element_text(size=10),
                   axis.text.y = element_text(size = 10),
                   axis.title.y= element_text(size = 12),
                   # legend.key = element_rect(size = 2),
                   # # legend.key.height=unit(0.6,"line"),
                   # # legend.key.width=unit(1,"line"),
                   legend.text = element_text(size = 12),
                   legend.title = element_text(size = 12),
                   # legend.background = element_rect(colour = "black"),
                   # legend.direction='horizontal',legend.box='horizontal'
                   strip.text.x = element_text(size = 12))
    # p <- p + coord_trans(y = "log10")       #,clip="on",expand = FALSE
    # p <- p + scale_y_log10(formatter=numClones)
    # p <-p + scale_y_continuous(trans = "log10")
    # p <- p  +
    #      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
    #               labels = trans_format("log10", math_format(10^.x)))
    print(p)
  }
  
}






runMultiSurvival <- function(DF,var,varSelection,varSelected,survType,sampleID,varPlotName, saveP = FALSE){
  # DF <- testDF.survival.prob.thres.split.l3[[1]]
  # var <-survVars[[1]]
  # varSelection <-tissues.OS
  # # typeVar <- "cont"
  # varSelected <-"tissue"
  # survType <- "OS"
  # sampleID <- "run_accession"
  # varPlotName<-names(testDF.survival.prob.thres.split.l3)[1]
  # 
  if(isTRUE(saveP)){
    # res.cont<-performSurvivalCalcs(varSelection, DF,var, varSelected, survType, sampleVar=sampleID,covName=varPlotName,covType="cont",stratifyMethod="median",sLevelA = "high",sLevelB = "low",onlyForest=TRUE)
    res.cat.km<-performSurvivalCalcs(varSelection, DF,var, varSelected, survType, sampleVar=sampleID,covName=varPlotName,covType="cat",
                                  stratifyMethod=NULL,sLevelA = "highProb",sLevelB = "lowProb", onlyKM = TRUE)
    res.cat.f<-performSurvivalCalcs(varSelection, DF,var, varSelected, survType, sampleVar=sampleID,covName=varPlotName,covType="cat",
                                  stratifyMethod=NULL,sLevelA = "highProb",sLevelB = "lowProb", onlyForest = TRUE)
    
    
    res = list(
               forestCat = res.cat.f$forestplot,
               km =  res.cat.km$kmPlot,
               hrDF = res.cat.f$hrDF
    )
    res
  }else{
    cat('\n')  
    cat('\n') 
    cat("### ", varPlotName, "{.tabset .tabset-fade .tabset-pills} \n") 
    cat('\n') 
    
    # The var will be the main covariate.
    # res.cont<-performSurvivalCalcs(varSelection, DF,var, varSelected, survType, sampleVar=sampleID,covName=varPlotName,covType="cont",stratifyMethod="median",sLevelA = "high",sLevelB = "low",onlyForest=TRUE)
    res.cat<-performSurvivalCalcs(varSelection, DF,var, varSelected, survType, sampleVar=sampleID,covName=varPlotName,covType="cat",stratifyMethod=NULL,sLevelA = "highProb",sLevelB = "lowProb")
    
    # cat('\n')  
    # cat('\n') 
    # cat("#### ", "Forest plot - continuous", "\n") 
    # cat('\n') 
    # 
    # print(res.cont$forestplot)
    # 
    cat('\n')  
    cat('\n') 
    cat("#### ", "Forest plot - categorical", "\n") 
    cat('\n') 
    
    print(res.cat$forestplot)
    
    cat('\n')  
    cat('\n') 
    cat("#### ", "KM plots", "\n") 
    cat('\n') 
    
    print(res.cat$kmPlot)
  }
  
}



findBestReduced <- function(rocDF,compDF){
  # rocDF <- rocDF.f
  # compDF <- melted_mat
  # Find model with highest ROC AUC
  highestROC <- rocDF$modnames[which.max(rocDF$aucs)]
  # Go to ROC AUC comparisons
  # Now grab those rows that include the best model
  rowsKeep <- c(which(compDF$Var1==highestROC),which(compDF$Var2==highestROC))[order(c(which(compDF$Var1==highestROC),which(compDF$Var2==highestROC)))]
  df.bestROC <-compDF[rowsKeep ,]
  # Select rows that have a value above 0.05
  df.bestROC <- df.bestROC %>% dplyr::filter(value > 0.05) %>% droplevels()
  # Extract other models names
  # REmove row of current best model
  rocDF.f <- rocDF[-which(rocDF$modnames==highestROC),]
  # get models that have less features than current best ROC model
  noFeats.best <- rocDF %>% dplyr::filter(modnames %in% highestROC) %>% pull(NOfeaturesOut)
  # How many models do we have left - those with no sign diff and that have less feats than current model
  redModels <-rocDF.f  %>% dplyr::filter(NOfeaturesOut<noFeats.best) %>% pull(modnames) %>% length()
  if(redModels==1){
    bestReducedModel <-rocDF.f  %>% dplyr::filter(NOfeaturesOut<noFeats.best) %>% pull(modnames)
    bestReducedModelFeats <-rocDF.f %>% dplyr::filter(NOfeaturesOut<noFeats.best) %>% pull(NOfeaturesOut)
    bestReducedModelFeats2 <-rocDF.f %>% dplyr::filter(NOfeaturesOut<noFeats.best) %>% pull(optPredictors)
    res <-list(bestReducedModel,bestReducedModelFeats,bestReducedModelFeats2)
    res
  }else if(redModels==0){
    bestReducedModel <-rocDF  %>% dplyr::filter(modnames==highestROC) %>% pull(modnames)
    bestReducedModelFeats <-rocDF %>% dplyr::filter(modnames==highestROC) %>% pull(NOfeaturesOut)
    bestReducedModelFeats2 <-rocDF %>% dplyr::filter(modnames==highestROC) %>% pull(optPredictors)
    res <-list(bestReducedModel,bestReducedModelFeats,bestReducedModelFeats2)
    res
  }else{
    findBestReduced(rocDF.f,melted_mat)
  }
  
}
