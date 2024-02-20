##########################
## Loading RData object ##
##########################
`%notin%` <- Negate(`%in%`)

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


##########################
## Histogram with facet ##
##########################

hist_measurement <- function(DF, variable,facetVariable=NULL,plotCols=NULL){
  # DF <- tcgaData.tp.compl
  # variable <- "TCR_Richness"
  # facetVariable <- "project"
  # plotCols=graphCols
  # facetVariable.un <- as.name(facetVariable)
  p<-ggplot(DF, aes_string(x=variable, color=variable, fill=variable)) +
    geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
    geom_density(alpha=0.4) 
  p <- p + scale_color_brewer(palette="Dark2")
  p <- p + scale_fill_brewer(palette="Dark2")
  
  if(!is.null(facetVariable)){
    p <- p + facet_wrap(facetVariable,ncol=plotCols,scales='free')
  }
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

####################################################
## Find optimal cutpoint for continuous variables ##
####################################################

optCutpointContinuous <- function(DF,timeVar,eventVar,varVector,filterVar=NULL, filterVarSelected=NULL, gatherTable=FALSE){
  # cat('\n')  
  # cat('\n') 
  # cat("### ", paste(filterVarSelected),"\n") 
  # cat('\n')
  # DF <-tcgaData.tp.compl
  # timeVar<-"OS.time"
  # eventVar<-"OS"
  # varVector <-c("TCR_Richness", "TCR_Evenness")
  # filterVar="project"
  # filterVarSelected="ACC"
  # Filter DF if needed
  if(!isTRUE(gatherTable)){
    if (!is.null(filterVar)){
      cat('\n')  
      cat('\n') 
      cat("#### ", paste(filterVarSelected),"\n") 
      cat('\n') 
      filterVar.un <- as.name(filterVar)
      DF <- DF %>%dplyr::filter(!!filterVar.un %in% filterVarSelected) %>% droplevels()
      DF
    }
  }else{
    if (!is.null(filterVar)){
      filterVar.un <- as.name(filterVar)
      DF <- DF %>%dplyr::filter(!!filterVar.un %in% filterVarSelected) %>% droplevels()
      DF
    }
  }
  # 1.Determine the optimal cutpoint of variables
  df.cut <- surv_cutpoint(DF, time = timeVar, event = eventVar,
                          variables = varVector, minprop=0.3)
  # 2. Categorize variables
  df.cat <- surv_categorize(df.cut)
  
  if(isTRUE(gatherTable)){
    #df.cat
    df.cutTable<-df.cut$cutpoint %>% rownames_to_column(var="variable")
    df.cutTable
  }else{
    
    if(!is.null(filterVar)){
      # get first var
      firstVar <-varVector[1]
      # print summary for each cancer
      # print(filterVarSelected)
      print(summary(df.cut))
      
      print(plot(df.cut, firstVar, palette = "npg"))
      table(df.cat[[firstVar]])
      
    }else{
      # get first var
      firstVar <-varVector[1]
      # print summary for each cancer
      # print(filterVarSelected)
      print(summary(df.cut))
      print(plot(df.cut, firstVar, palette = "npg"))
      table(df.cat[[firstVar]])
      
    }
  }
  
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
  # tissues_fewData <- tissue_dataPoints %>%dplyr::filter(no_rows<10) %>% select(!!facetVariable.un) %>% pull(!!facetVariable.un) %>% droplevels()
  tissue_dataPoints
  
  numb_DataPoints_df <- rbind(numb_DataPoints_df,tissue_dataPoints)
  numb_DataPoints_df
}




######################
## COX HAZARD RATIO ##
######################



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
      coxCIdf.cov <- coxCIdf %>%dplyr::filter(covariate==covar)
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
      coxCIdf.cov <- coxCIdf %>%dplyr::filter(covariate==covar)
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



cox_hazard <- function(DF, selectedcancer, SampleCol, endopoint, event,facetVariable,predictor,hrDF,cutPointsTable=NULL,signature=NULL, 
                       dichot = FALSE, stratMethod = "median",stratVariable=NULL,StratLevel1=NULL,StratLevel2=NULL,globalLocal="global"){#dichotLabels=NULL,
  # print(selectedcancer)
  
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
  # selectedcancer <- varSelection[14]
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
    dplyr::mutate_at(vars(!!endopoint.un,!!event.un),as.numeric)  %>%dplyr::filter(complete.cases(!!event.un,!!endopoint.un))
  #filter(complete.cases(!!endopoint.un,!!event.un,!!predictor.un))
  
  
  if(isTRUE(dichot)){
    # DICHOTOMIZE-- HERE I HAVE TO CHECK IF TISSUE OR DATASET
    if(stratMethod=="median"){
      # Is it a bi-modal distribution?How to select the cutoff?
      cutglobal <- median(DF[[predictor]][!is.na(DF[[predictor]])], na.rm = TRUE)
      # Calculate median for cancer
      cutlocal <- median(DF.sub[[predictor]][!is.na(DF.sub[[predictor]])], na.rm = TRUE)
    }else if(stratMethod=="mean"){
      cutglobal <- mean(DF[[predictor]][!is.na(DF[[predictor]])], na.rm = TRUE)
      cutlocal <- mean(DF.sub[[predictor]][!is.na(DF.sub[[predictor]])], na.rm = TRUE)
    }else if(stratMethod=="cut"){
      cutglobal <-as.numeric(optCutpointContinuous(DF,endopoint,event,predictor,filterVar=NULL, filterVarSelected=NULL,gatherTable=TRUE)[2])
      cutlocal <- as.numeric(cutPointsTable[[selectedcancer]][2])
    }else if(stratMethod=="opt"){
      # check if dataset/tissue
      if (facetVariable=="tissue"){
        selectedTissue=selectedcancer
      }else if(facetVariable=="dataset"){
        # Figure out the type of tissue it is
        selectedTissue=as.character(unique(DF.sub[["tissue"]]))
      }
      
      if(signature=="TuTack"){
        cutglobal <-as.numeric(cutPointsTable[[selectedTissue]][1])
        cutlocal <-as.numeric(cutPointsTable[[selectedTissue]][1])
      }else if(signature=="TuTack.ifng"){
        cutglobal <-as.numeric(cutPointsTable[[selectedTissue]][2])
        cutlocal <-as.numeric(cutPointsTable[[selectedTissue]][2])
      }else if(signature=="GEP"){
        cutglobal <-as.numeric(cutPointsTable[[selectedTissue]][3])
        cutlocal <-as.numeric(cutPointsTable[[selectedTissue]][3])
      }
    }
    
    
    
    if(globalLocal=="global"){
      cutVar <- cutglobal
      #subtitlePlot<-""
    }else{
      cutVar <- cutlocal
      #subtitlePlot <-paste0("Richness Threshold (median): ", medianVar, "\n", "# Samples: ", nrow(tcgaDF_sign))
    }
    if(cutVar==0){
      DF.sub2 <- DF.sub %>% dplyr::mutate(!!stratVariable.un := ifelse(!!predictor.un > cutVar,StratLevel1,StratLevel2)) %>%mutate(!!stratVariable:=factor(!!stratVariable.un)) %>%dplyr::filter(complete.cases(!!event.un,!!endopoint.un))
      
    }else{
      DF.sub2 <- DF.sub %>% dplyr::mutate(!!stratVariable.un := ifelse(!!predictor.un >= cutVar,StratLevel1,StratLevel2)) %>%mutate(!!stratVariable:=factor(!!stratVariable.un)) %>%dplyr::filter(complete.cases(!!event.un,!!endopoint.un))
      
    }
    # DF.sub2[[stratVariable]] <- factor(DF.sub2[[stratVariable]],levels = c())
  }else{
    stratVariable <- predictor
    DF.sub2 <-DF.sub
  }
  
  # HERE I HAVE TO DEFINE THE LEVEL-looking at low or high
  if(isTRUE(dichot)){
    DF.sub2[[predictor]] = stats::relevel(DF.sub2[[predictor]], ref = StratLevel2)
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

#############################
## HazardRatio forest plot ##
#############################


survival_forest_plot <-function(survivalDF,selectedCancer,predictorName=NULL,facetVariable=NULL,facetsCols=NULL,logHR= FALSE){
  # survivalDF<-survival_HR_logrankP2
  # selectedCancer<-filterVar
  # predictorName=covName
  # facetVariable="covariate"
  # facetsCols=2
  # logHR= FALSE
  # 
  
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
    }else if(selectedCancer=="project"){
      p <- ggplot(survivalDF, aes(x=project, y=HR, ymin=HR_lower, ymax=HR_upper)) 
    }else if(selectedCancer=="dataset"){
      p <- ggplot(survivalDF, aes(x=dataset, y=HR, ymin=HR_lower, ymax=HR_upper)) 
    }else if(selectedCancer=="tisData"){
      p <- ggplot(survivalDF, aes(x=tisData, y=HR, ymin=HR_lower, ymax=HR_upper)) 
    }
    p <- p + ylim(floor(min(survivalDF[["HR_lower"]][!is.infinite(survivalDF[["HR_lower"]])])), 
                  ceiling(max(survivalDF[["HR_upper"]][!is.infinite(survivalDF[["HR_upper"]])])))
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
    if(!is.null(predictorName)){
      p <- p + ylab(paste0("Hazard Ratio (95% Confidence Interval) \nassociated with unit increase in ",predictorName))
      # p <- p + ylab(paste0("Hazard Ratio (95% Confidence Interval) \nassociated with TuTack/GEP positive score"))
    }else{
      p <- p + ylab(paste0("Hazard Ratio (95% Confidence Interval) \nassociated with unit increase in covariate"))
      
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
  }
  
  #p <- p + scale_y_discrete()
  p
}

#########################
## Kaplan Meier curves ##
#########################


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
  # selectedcancer,tcgaDF,survHRdf,cutPointsTable,covName,
  # stratMethod,signature,stratVariable,SampleCol,endopoint, event,
  # facetVariable,predictor,StratLevel1,StratLevel2
  
  # selectedcancer <- datasetsEvent[14]
  # tcgaDF <- dataDF
  # survHRdf <-survival_HR_logrankP2
  # cutPointsTable = NULL 
  # covName = mainVar 
  # stratMethod=stratifyMethod
  # signature=NULL 
  # stratVariable = mainVar
  # SampleCol = sampleVar
  # endopoint <-paste0(survivalVar,".time")
  # event <-survivalVar
  # facetVariable = filterVar
  # predictor = mainVar
  # StratLevel1 = sLevelA
  # StratLevel2 = sLevelB 
  # globalLocal="local"
  # covType=covType
  # getLegend=FALSE
  # HRAnnot = TRUE
  ##########################################
  
  # print(selectedcancer)
  
  # Prepare variables
  SampleCol <- as.name(SampleCol)
  endopoint.un <- as.name(endopoint)
  event.un <- as.name(event)
  facetVariable.un <- as.name(facetVariable)
  predictor.un <- as.name(predictor)
  stratVariable.un <- as.name(stratVariable)
  # facetVariable.un <- as.name(facetVariable)
  # Survival plots
  # First subset data to the cancers with significant associations
  tcgaDF_sign <-dplyr::filter(tcgaDF, !!facetVariable.un %in% selectedcancer) %>% dplyr::filter(complete.cases(!!event.un,!!endopoint.un)) %>% droplevels()
  tcgaDF_sign <- tcgaDF_sign %>% distinct(!!SampleCol,.keep_all = TRUE)
  if (selectedcancer %in% c("bladder","melanoma","RCC","gastric")){
    tissueSelected <- ''
  }else{
    tissueSelected <- tcgaDF_sign %>% pull(!!facetVariable.un) %>% unique()
    tissueSelected <- paste0("-",tissueSelected)
  }
  # tissueSelected <- tcgaDF_sign %>% pull(tissue) %>% unique()
  # DICHOTOMIZE-- HERE I HAVE TO CHECK IF TISSUE OR DATASET
  if(stratMethod=="median"){
    # Is it a bi-modal distribution?How to select the cutoff?
    cutglobal <- median(tcgaDF[[predictor]][!is.na(tcgaDF[[predictor]])],na.rm = TRUE)
    # Calculate median for cancer
    cutlocal <- median(tcgaDF_sign[[predictor]][!is.na(tcgaDF_sign[[predictor]])],na.rm = TRUE)
  }else if(stratMethod=="mean"){
    cutglobal <- mean(tcgaDF[[predictor]][!is.na(tcgaDF[[predictor]])],na.rm = TRUE)
    cutlocal <- mean(tcgaDF_sign[[predictor]][!is.na(tcgaDF_sign[[predictor]])],na.rm = TRUE)
  }else if(stratMethod=="cut"){
    cutglobal <-as.numeric(optCutpointContinuous(tcgaDF,endopoint,event,predictor,filterVar=NULL, filterVarSelected=NULL,gatherTable=TRUE)[2])
    cutlocal <- as.numeric(cutPointsTable[[selectedcancer]][2])
  }else if(stratMethod=="opt"){
    # check if dataset/tissue
    if (facetVariable=="tissue"){
      selectedTissue=selectedcancer
    }else if(facetVariable=="dataset"){
      # Figure out the type of tissue it is
      selectedTissue=unique(tcgaDF_sign[["tissue"]])
    }
    
    if(signature=="TuTack"){
      cutglobal <-as.numeric(cutPointsTable[[selectedTissue]][1])
      cutlocal <-as.numeric(cutPointsTable[[selectedTissue]][1])
    }else if(signature=="TuTack.ifng"){
      cutglobal <-as.numeric(cutPointsTable[[selectedTissue]][2])
      cutlocal <-as.numeric(cutPointsTable[[selectedTissue]][2])
    }else if(signature=="GEP"){
      cutglobal <-as.numeric(cutPointsTable[[selectedTissue]][3])
      cutlocal <-as.numeric(cutPointsTable[[selectedTissue]][3])
    }
  }
  
  
  
  if(globalLocal=="global"){
    cutVar <- cutglobal
    #subtitlePlot<-""
  }else{
    cutVar <- cutlocal
    #subtitlePlot <-paste0("Richness Threshold (median): ", medianVar, "\n", "# Samples: ", nrow(tcgaDF_sign))
  }
  
  if(cutVar==0){
    tcgaDF_sign2 <- tcgaDF_sign %>% dplyr::mutate(!!stratVariable.un := ifelse(!!predictor.un > cutVar,StratLevel1,StratLevel2)) %>%
      dplyr::mutate(!!stratVariable:=factor(!!stratVariable.un)) %>% 
     dplyr::filter(complete.cases(!!event.un,!!endopoint.un))
    
  }else{
    tcgaDF_sign2 <- tcgaDF_sign %>% dplyr::mutate(!!stratVariable.un := ifelse(!!predictor.un >= cutVar,StratLevel1,StratLevel2)) %>%
      dplyr::mutate(!!stratVariable:=factor(!!stratVariable.un)) %>% 
     dplyr::filter(complete.cases(!!event.un,!!endopoint.un))
    
  }
  
  # tcgaDF_sign2 <- tcgaDF_sign %>% dplyr::mutate(!!stratVariable.un := ifelse(!!predictor.un >= cutVar,StratLevel1,StratLevel2)) %>%mutate(!!stratVariable:=factor(!!stratVariable.un))
  # Cleanup df
  tcgaDF_sign2[[endopoint]] <- gsub(",","",tcgaDF_sign2[[endopoint]],fixed=TRUE)
  tcgaDF_sign2[[endopoint]] <- gsub("#N/A",NA,tcgaDF_sign2[[endopoint]],fixed=TRUE)
  tcgaDF_sign2[tcgaDF_sign2=="NA"] = NA
  tcgaDF_sign2[tcgaDF_sign2=="NA"] = NA
  
  tcgaDF_signFinal <- tcgaDF_sign2 %>% 
    dplyr::mutate(!!facetVariable:=factor(!!facetVariable.un)) %>%
   dplyr::filter(complete.cases(!!endopoint.un,!!event.un,!!predictor.un,!!stratVariable.un)) %>% 
   dplyr::filter(!!facetVariable.un == selectedcancer) %>%
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
    # palette = c("#8C3F4D","#3E606F"),
    palette = c("black","grey"),
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
      # palette = c("#8C3F4D","#3E606F"),
      palette = c("black","grey"),
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
    
    
    # hr <-survHRdf %>%dplyr::filter(!!facetVariable.un==selectedcancer) %>%dplyr::filter(covariate==paste0(covName,StratLevel2)) %>% pull(HR)
    # ci.l <-survHRdf %>%dplyr::filter(!!facetVariable.un==selectedcancer) %>%dplyr::filter(covariate==paste0(covName,StratLevel2)) %>% pull(HR_lower)
    # ci.u <-survHRdf %>%dplyr::filter(!!facetVariable.un==selectedcancer) %>%dplyr::filter(covariate==paste0(covName,StratLevel2)) %>% pull(HR_upper)
    # ci <- c(ci.l,ci.u)
    # ## logrankp <- round(survpval$pval,5)
    # # logrankp <- formatC(survpval$pval,format = "e", digits = 3)
    # logrankp <- format(survHRdf %>%dplyr::filter(!!facetVariable.un==selectedcancer) %>%dplyr::filter(covariate==paste0(covName,StratLevel2)) %>% pull(logrankP),format = "e", digits = 3)
    # 
    if(covType=="cont"){
      hr <-survHRdf %>%dplyr::filter(!!facetVariable.un==selectedcancer) %>%dplyr::filter(covariate==covName) %>% pull(HR)
      ci.l <-survHRdf %>%dplyr::filter(!!facetVariable.un==selectedcancer) %>%dplyr::filter(covariate==covName) %>% pull(HR_lower)
      ci.u <-survHRdf %>%dplyr::filter(!!facetVariable.un==selectedcancer) %>%dplyr::filter(covariate==covName) %>% pull(HR_upper)
      ci <- c(ci.l,ci.u)
      ## logrankp <- round(survpval$pval,5)
      # logrankp <- formatC(survpval$pval,format = "e", digits = 3)
      logrankp <- format(survHRdf %>%dplyr::filter(!!facetVariable.un==selectedcancer) %>%dplyr::filter(covariate==covName) %>% pull(logrankP),format = "e", digits = 3)
    }else if(covType=="cat"){
      # HERE WE DEFINE BASED ON WHICH LEVEL OF COVARIATE IS SET AS REFERENCE_NEED TO CHANGE IF NECESSARY, , if ref level is strat level 2(low) we show results for strat level 1 (high)
      hr <-survHRdf %>%dplyr::filter(!!facetVariable.un==selectedcancer) %>%dplyr::filter(covariate==paste0(covName,StratLevel1)) %>% pull(HR)
      ci.l <-survHRdf %>%dplyr::filter(!!facetVariable.un==selectedcancer) %>%dplyr::filter(covariate==paste0(covName,StratLevel1)) %>% pull(HR_lower)
      ci.u <-survHRdf %>%dplyr::filter(!!facetVariable.un==selectedcancer) %>%dplyr::filter(covariate==paste0(covName,StratLevel1)) %>% pull(HR_upper)
      ci <- c(ci.l,ci.u)
      ## logrankp <- round(survpval$pval,5)
      # logrankp <- formatC(survpval$pval,format = "e", digits = 3)
      logrankp <- format(survHRdf %>%dplyr::filter(!!facetVariable.un==selectedcancer) %>%dplyr::filter(covariate==paste0(covName,StratLevel1)) %>% pull(logrankP),format = "e", digits = 3)
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




# ##############
# DF <- tcgaData.tp.compl.filt
# var <-"TCR_Richness"
# varSelection <-unique(tcgaData.tp.compl.filt$project)
# typeVar <- "cont"
# varSelected <-"project"
# survType <- "OS"
# sampleID <- "bcr_patient_barcode"
# varPlotName<-"TCR Richness"
# 
# res.cont<-performSurvivalCalcs(datasetsEvent=unique(tcgaData.tp.compl.filt$project), dataDF = tcgaData.tp.compl.filt,
#                                mainVar = "TCR_Richness", filterVar = "project",
#                                survivalVar = "OS", sampleVar="bcr_patient_barcode",covName="TCR Richness",covType="cont",
#                                stratifyMethod="median",sLevelA = "high",sLevelB = "low",onlyForest=TRUE)
# 
# res.cat<-performSurvivalCalcs(datasetsEvent=unique(tcgaData.tp.compl.filt$project), dataDF = tcgaData.tp.compl.filt,
#                               mainVar = "TCR_Richness", filterVar = "project",
#                               survivalVar = "OS", sampleVar="bcr_patient_barcode",covName="TCR Richness",covType="cat",
#                               stratifyMethod="median",sLevelA = "high",sLevelB = "low")
# 
# ## WRAPPER WITH ALL SURV ANALYSIS ##
performSurvivalCalcs <- function(datasetsEvent, dataDF,mainVar, filterVar, survivalVar, sampleVar,covName,covType="cont",
                                 stratifyMethod=NULL,sLevelA = "high",sLevelB = "low",printPlots= FALSE,onlyForest=FALSE, onlyKM = FALSE){
  
  
  # dataDF <-tcgaData.tp.compl.filt
  # datasetsEvent<- unique(tcgaData.tp.compl.filt$project)
  # # datasetTissueName <- "Hugo"
  # # dataPointsDF <-
  # mainVar <- "TCR_Richnesso"
  # filterVar <- "project"
  # survivalVar <- "OS"
  # sampleVar <- "bcr_patient_barcode"
  # covName <- "TCR Richness"
  # covType <- "cont" # cont or cat
  # stratifyMethod <- "median"
  # sLevelA = "high"
  # sLevelB = "low"
  # printPlots= FALSE
  # ##########################
  # datasetsEvent <- varSelection
  # dataDF <- DF
  # mainVar <-var[[1]]
  # filterVar <-varSelected
  # survivalVar <-survType
  # sampleVar=sampleID
  # covName=varPlotName
  # covType="cont"
  # stratifyMethod="median"
  # sLevelA = "high"
  # sLevelB = "low"
  # printPlots= FALSE
  # 
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
  survival_HR_logrankP  <-bind_rows(survival_HR_logrankP)
  
  
  # survival_HR_logrankP.cat <- data.frame(filterVar=character(),
  #                                        covariate=character(),
  #                                        HR= character(),
  #                                        HR_lower = character(),
  #                                        HR_upper = character(),
  #                                        logrankP = character(),
  #                                        stringsAsFactors=FALSE)
  # 
  # colnames(survival_HR_logrankP.cat)[which(colnames(survival_HR_logrankP.cat)=="filterVar")] <- filterVar
  # 
  # survival_HR_logrankP.cat<- sapply(datasetsEvent,function(x) cox_hazard(dataDF, x,sampleVar, paste0(survivalVar,".time"), 
  #                                                                        survivalVar,filterVar,c(mainVar),survival_HR_logrankP.cat,
  #                                                                        cutPointsTable=NULL,signature=NULL,
  #                                                                        dichot = TRUE, stratMethod = "median",stratVariable=mainVar,StratLevel1="high",StratLevel2="low",globalLocal="global"),simplify = FALSE)#, "NK_cells"
  # # Merge data
  # survival_HR_logrankP.cat  <-bind_rows(survival_HR_logrankP.cat)
  
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
  
  
  # # CHOOSE
  # if(covType=="cont"){
  #   survival_HR_logrankP <- survival_HR_logrankP
  # }else if(covType == "cat"){
  #   survival_HR_logrankP <-survival_HR_logrankP.cat
  # }
  
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
    # NOT logarithmic HR
    p <-survival_forest_plot(survival_HR_logrankP2,filterVar, predictorName=covName,facetVariable="covariate",facetsCols=2,logHR= FALSE)
    # print(p)
    n<-dim(survival_HR_logrankP2)[1]/2;
    graphCols <-ceiling(n/floor(sqrt(n)));
    graphRows<-floor(n)
    # 
    
    global_legend <-fit_all(as.character(datasetsEvent[2]), dataDF, survival_HR_logrankP2,
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
    
    # p.km<- as_ggplot(arrangeGrob(grobs = kmPlots, ncol = graphCols,newpage = FALSE))
    
    p.km<- ggarrange(plotlist = kmPlots,ncol=graphCols, nrow =  graphRows,common.legend = TRUE, legend="bottom")
    
    res = list(hrDF=survival_HR_logrankP2,forestplot=p, kmPlot = p.km,legend=global_legend)
  }
  
  res
}


runMultiSurvival_TCGA <- function(DF,var,varSelection,varSelected,survType,sampleID,varPlotName,printP=FALSE,saveP = FALSE,savePath,onlySign=TRUE,selectFacets=NULL){
  # DF <- tcgaData.tp.sub.filt
  # var <-"TCR_Richness"
  # varSelection <-mycancertypes.final
  # varSelected <-"project"
  # survType <- "OS"
  # sampleID <- "bcr_patient_barcode"
  # varPlotName<-"TCR Richness"
  # onlySign=FALSE
  # selectFacets=mycancertypes.final
  # typeVar <- "cont"
  # 
  
  
  
  cat('\n')
  cat('\n')
  cat("## ", varPlotName, "{.tabset .tabset-fade .tabset-pills} \n")
  cat('\n')
  
  # The var will be the main covariate.
  res.cont<-performSurvivalCalcs(varSelection, DF,var, varSelected, survType, sampleVar=sampleID,covName=varPlotName,covType="cont",stratifyMethod="median",sLevelA = "high",sLevelB = "low",onlyForest=TRUE)
  res.cat<-performSurvivalCalcs(varSelection, DF,var, varSelected, survType, sampleVar=sampleID,covName=varPlotName,covType="cat",stratifyMethod="median",sLevelA = "high",sLevelB = "low",onlyForest=TRUE)
  
  res.cat2<-performSurvivalCalcs(varSelection, DF,var, varSelected, survType, sampleVar=sampleID,covName=varPlotName,
                                 covType="cat",stratifyMethod="median",
                                 sLevelA = "high",sLevelB = "low",
                                 onlyForest=FALSE, onlyKM = TRUE)
  
  
  # Extract Only cancers with significant HR
  if(isTRUE(onlySign)){
    hrDF.sign <- res.cat2$hrDF %>% as.data.frame() %>% dplyr::filter(HR_eval %in% c("signGood", "signPoor")) %>% droplevels()
  }else{
    if(!is.null(selectFacets)){
      varSelected.un <- as.name(varSelected)
      hrDF.sign <- res.cat2$hrDF %>% as.data.frame() %>% dplyr::filter(!!varSelected.un %in% selectFacets) %>% droplevels()
      hrDF.sign <- hrDF.sign[match(selectFacets, hrDF.sign$project),]
    }else{
      hrDF.sign <- res.cat2$hrDF %>% as.data.frame
    }
  }
  
  if(isTRUE(saveP)){
    res = list(forestCont = res.cont$forestplot,
               forestCat = res.cat$forestplot,
               km =  res.cat2$kmPlot,
               hrDF = res.cat2$hrDF)
    # print(paste0(projectRDataPath, paste0("tcga.SURV","_",survType,"_",var),Rdata.suffix))
    save(res, file = paste0(savePath, paste0("tcga.SURV","_",survType,"_",var),Rdata.suffix))
  }
  
  if(isTRUE(printP)){
    cat('\n')  
    cat('\n') 
    cat("### ", "Forest plot - continuous", "\n") 
    cat('\n') 
    
    print(res.cont$forestplot)
    
    cat('\n')  
    cat('\n') 
    cat("### ", "Forest plot - categorical", "\n") 
    cat('\n') 
    
    print(res.cat$forestplot)
    
    cat('\n')  
    cat('\n') 
    cat("### ", "KM plots", "\n") 
    cat('\n') 
    
    
    n<-dim(hrDF.sign)[1]/2;
    graphCols <-ceiling(n/floor(sqrt(n)));
    graphRows<-floor(n)
    
    p.km<- ggarrange(plotlist = res.cat2$kmPlot[hrDF.sign$project], nrow =  graphRows,ncol=graphCols,common.legend = TRUE, legend="bottom")#
    
    print(p.km)
  }
  # res = list(forestCont = res.cont$forestplot,
  #            forestCat = res.cat$forestplot,
  #            km =  res.cat2$kmPlot,
  #            hrDF = res.cat2$hrDF)
  res.cat2$hrDF
}
