#######################
### Functions used ###
######################

###############################################
##  ##
###############################################
`%notin%` <- Negate(`%in%`)


##########################
## Loading RData object ##
##########################


loadRData <- function(fileName){
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
########################################################
## Multi Correlation plot-facet variable for all tcga ##
########################################################
ks <- function (x) { number_format(accuracy = 1,
                                   scale = 1/1000,
                                   suffix = "k",
                                   big.mark = ",")(x) }

multi_corrPlots_tcga <- function(tcgaDF,selectedcancers, xVariable, yVariable, xLabel, yLabel, facetVariable,facetVariable2,SampleCol, plotCols, logxLabel=FALSE, xUnitScale=FALSE){
  # Prepare variables
  SampleCol <- as.name(SampleCol)
  xVariable.un <- as.name(xVariable)
  yVariable.un <- as.name(yVariable)
  facetVariable.un <- as.name(facetVariable)
  facetVariable2.un <- as.name(facetVariable2)
  # Clean DF, make columns factors, remove NA rows for x-z variables and keep one row per sample ID
  tcgaDF.clean <- tcgaDF %>% mutate(!!facetVariable:=factor(!!facetVariable.un), !!facetVariable2:=factor(!!facetVariable2.un)) %>% dplyr::filter(complete.cases(!!yVariable.un, !!xVariable.un)) %>% dplyr::filter(!!facetVariable.un %in% selectedcancers) %>% 
    droplevels()
  # Melt tibble
  tcgaDFmelted <-melt(tcgaDF.clean[,c(xVariable, yVariable,facetVariable)],c(xVariable, yVariable,facetVariable))
  # Find number of datapoints for each cancer type
  tissue_dataPoints <- tcgaDF.clean  %>% group_by(!!facetVariable.un) %>% 
    distinct(!!SampleCol,.keep_all = TRUE) %>%
    filter(!any(is.na(!!xVariable.un))) %>%
    summarize(no_rows=length(!!facetVariable.un))
  
  ggtheme = ggplot2::theme_minimal(base_family = "Palatino Linotype", base_size = 10)
  p1 <- ggplot(tcgaDFmelted,aes_string(xVariable,yVariable))
  p1 <- p1 + geom_point(size=1,shape = 16, show.legend = TRUE, alpha = .4)#,color="#69b3a2"
  # p1 <- p1 + geom_text_repel(size = 3,
  #                            segment.size=0.2,
  #                            nudge_x = 0.015,
  #                            point.padding = unit(8.1, "points"))
  p1 <- p1 + xlab(xLabel)
  p1 <- p1 + ylab(yLabel)
  p1 <- p1 + geom_smooth(method=lm,
                         se=FALSE,
                         color="Gray",
                         linetype = "dashed",size=0.3)
  p1 <- p1 + ggtheme
  p1 <- p1 + theme(axis.text.x = element_text(hjust = 1, vjust = 1,size = 6),#angle = 90,
                   axis.title.x = element_text(size=12),
                   axis.text.y = element_text( hjust = 1, vjust = 1,size = 6),
                   axis.title.y = element_text(size=12),
                   legend.title = element_text(size=12,color="black"),
                   legend.text = element_text(size = 10),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   axis.line = element_line(colour = "black"),
                   legend.box.background = element_rect(colour = "black"),
                   legend.position="right", legend.box = "vertical")
  if (xUnitScale==TRUE){
    p1 <- p1 + scale_x_continuous(labels = unit_format(unit = "K", scale=1e-3))
  }
  
  if (logxLabel==TRUE){
    p1 <- p1 +scale_x_continuous(trans='log10')
  }
  #
  p1 <- p1 + facet_wrap(facetVariable,ncol=plotCols,scales='free')
  print(p1)
}


corPlot <- function(DF, Xvar, Yvar, xLabel,yLabel, xPlus, yPlus, yMinus,Logx= FALSE, filterVar=NULL, filterVarSelected="", shapeVar=NULL, shapeVarLabel="", facetVar=NULL, plotCols=NULL, addCorLabel=NULL){
  # DF <-tcgaData.tp.sub.filt
  # Xvar <-"BCR_Richness"
  # Yvar <-"B_cells"
  # xLabel <-"BCR Richness"
  # yLabel <-"B cells"
  # xPlus <-2
  # yPlus <-0.35
  # yMinus <-0.17
  # Logx= FALSE
  # filterVar="project"
  # filterVarSelected=c("SKCM","KIRC","BLCA","STAD")
  # shapeVar=NULL
  # shapeVarLabel=""
  # facetVar="project"
  # plotCols=2
  # addCorLabel=TRUE
  
  if (!is.null(filterVar)){
    filterVar.un <- as.name(filterVar)
    DF <- DF %>% dplyr::filter(!!filterVar.un %in% filterVarSelected) %>% droplevels()
    DF[[filterVar]] <- factor(DF[[filterVar]], levels=filterVarSelected)
    DF
  }
  
  # if(!is.null(facetVar)){
  #   DF[[facetVar]] <- factor(DF[[facetVar]], levels=filterVarSelected)
  # }
  ######################################
  Xvar.un <- as.name(Xvar)
  Yvar.un <- as.name(Yvar)
  if(!is.null(facetVar)){
    facetVar.un <- as.name(facetVar)
  }
  ## PLOT
  
  if(is.null(shapeVar)){
    if(!is.null(facetVar)){
      DF <- DF %>% dplyr::select(!!Xvar.un,!!Yvar.un,!!facetVar.un) %>% dplyr::filter(complete.cases(.))
      
    }else{
      DF <- DF %>% dplyr::select(!!Xvar.un,!!Yvar.un) %>% dplyr::filter(complete.cases(.))
      
    }
    
    p1 <- ggplot(DF,aes_string(Xvar,Yvar))
  }else{
    shapeVar.un <- as.name(shapeVar)
    if(!is.null(facetVar)){
      DF <- DF %>% dplyr::select(!!Xvar.un,!!Yvar.un,!!facetVar.un) %>% dplyr::filter(complete.cases(.))
      
      
    }else{
      DF <- DF %>% dplyr::select(!!Xvar.un,!!Yvar.un) %>% dplyr::filter(complete.cases(.))
      
      
    }
    
    
    p1 <- ggplot(DF,aes_string(Xvar,Yvar, shape=shapeVar))
  }
  
  if (Logx==TRUE){
    p1 <- p1 + scale_x_continuous(trans='log10')
  }
  
  p1 <- p1 + geom_point(size = 1, show.legend = TRUE, alpha = .4)#,color="#69b3a2"
  p1 <- p1 + xlab(paste0(xLabel))+ylab(yLabel)
  p1 <- p1 + geom_smooth(method=lm, se=FALSE, color="black",linetype = "dashed",size=0.5)
  
  p1 <- p1 + theme_ipsum(base_family = "Palatino Linotype", base_size = 10)
  
  if(is.null(shapeVar)){
    p1 <- p1
  }else{
    p1 + scale_color_discrete(name=shapeVarLabel)
  }
  
  if(!is.null(facetVar)){
    p1 <- p1 + theme(axis.text.x = element_text(hjust = 1, vjust = 1,size = 10),
                     axis.title.x = element_text(size=16,vjust = 0.5, hjust=0.5),
                     axis.text.y = element_text( hjust = 1, vjust = 1,size = 10), 
                     axis.title.y = element_text(size=16,angle=90, vjust = 0.5,hjust=0.5),
                     legend.title = element_text(size=13,color="black"),
                     legend.text = element_text(size = 10),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"))
    p1 <- p1 + facet_wrap(facetVar,ncol=plotCols,scales='free') 
    # p1 <- p1 + facet_grid(facetVar,cols=plotCols,scales='free') 
    
    
  }
  
  if (!is.null(addCorLabel)){
    p1 <- p1 + theme(axis.text.x = element_text(hjust = 1, vjust = 1,size = 10),
                     axis.title.x = element_text(size=12,),
                     axis.text.y = element_text( hjust = 1, vjust = 1,size = 10), 
                     axis.title.y = element_text(size=12,angle=90),
                     legend.title = element_text(size=12,color="black"),
                     legend.text = element_text(size = 10),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"))
    # # Pearson
    # corTest.pearson <- cor.test(DF[[Xvar]],DF[[Yvar]],method = "pearson")
    # 
    # l.pearson <- corTest.pearson$p.value
    # l.pearson <- format(round(l.pearson,4), scientific = TRUE)
    # l.pearson <-strsplit(l.pearson,"e")
    # # quote the part before the exponent to keep all the digits
    # l.pearson1 <- l.pearson[[1]][1]
    # l.pearson2 <- gsub("-0","-",l.pearson[[1]][2])
    # 
    # annot.pearson <- annotate("text", x = min(DF[[Xvar]]) + xPlus, y = max(DF[[Yvar]])+yPlus, label = bquote(paste("Pearson R" == .(round( corTest.pearson$estimate,2))," (P" < .(l.pearson1) %*% 10^ .(l.pearson2),")")),fontface ="plain",size = 4,family= "Palatino Linotype")
    # 
    # # Spearman
    # corTest.spearman <- cor.test(DF[[Xvar]],DF[[Yvar]],method = "spearman")
    # 
    # l.spearman <- corTest.spearman$p.value
    # l.spearman <- format(round(l.spearman,4), scientific = TRUE)
    # l.spearman <-strsplit(l.spearman,"e")
    # # quote the part before the exponent to keep all the digits
    # l.spearman1 <- l.spearman[[1]][1]
    # l.spearman2 <- gsub("-0","-",l.spearman[[1]][2])
    # 
    # annot.spearman <- annotate("text", x = min(DF[[Xvar]]) + xPlus, y = max(DF[[Yvar]])+yPlus-yMinus, label = bquote(paste("Spearman R" == .(round( corTest.spearman$estimate,2))," (P" < .(l.spearman1) %*% 10^ .(l.spearman2),")")),fontface ="plain",size = 4,family= "Palatino Linotype")
    # 
    # p1 <- p1 + annot.pearson
    p1 <- p1 +
      stat_cor(aes(label = ..r.label..),cor.coef.name="R",label.y.npc="bottom", label.x.npc = "right",hjust = 1, method = "spearman")##,vjust=0
  }
  p1 <- p1 + scale_y_continuous(breaks = scales::pretty_breaks(n = 3))
  p1 <- p1 + scale_x_continuous(breaks = scales::pretty_breaks(n = 3))
  print(p1)
}


corr_calculation <- function(DF, Xvar,Yvar, corMethod,filterVar, filterVarSelected){
  if (!is.null(filterVar)){
    filterVar.un <- as.name(filterVar)
    DF <- DF %>% dplyr::filter(!!filterVar.un %in% filterVarSelected) %>% droplevels() %>% distinct()
    DF
  }
  
  ## Calculate correlations
  if(nrow(DF)>10){
    corTest <- cor.test(DF[[Xvar]],DF[[Yvar]],method = corMethod)
    corTest_se <- round(unname(sqrt((1 - corTest$estimate^2)/corTest$parameter)),3)
    data.label <- data.frame(label=paste("r = ", round(corTest$estimate,3),"\u00B1",corTest_se), x=0.1*max(DF[[Xvar]],na.rm=T), y=max(DF[[Yvar]], na.rm=T)*.8)
    
    cor <- round(corTest$estimate,3)
    cor
  }else{
    cor <- NULL
  }
  
  if(is.null(cor)==TRUE){
    cor <- NA
    cor
  }else{
    cor <- round(corTest$estimate,3)
  }
  cor
}

cor_heatmap <- function(corDF,DF,xVar,yVar,fillVar,facetVar,SampleCol,xLabel){
  # corDF <-cor_all.tis.tutack
  # DF <-orr.tutack.tmb.tcga.all
  # xVar <-"correlation"
  # yVar <-"project"
  # fillVar <-"correlation"
  # facetVar <-"project"
  # SampleCol <-"bcr_patient_barcode"
  # xLabel <-"Spearman Correlation"
  
  facetVar.un <- as.name(facetVar)
  SampleCol.un <- as.name(SampleCol)
  xVar.un <- as.name(xVar)
  # Find number of datapoints for each cancer type
  tissue_dataPoints <- DF %>% dplyr::select(!!facetVar.un,!!SampleCol.un) %>% group_by(!!facetVar.un) %>% distinct() %>% dplyr::summarize(n())
  colnames(tissue_dataPoints)[2] <- "no_rows" 
  tissues_fewData <- tissue_dataPoints %>% dplyr::filter(no_rows<15) %>% dplyr::select(!!facetVar.un) %>% dplyr::pull(!!facetVar.un)
  
  
  summary_cor <- as.data.frame(corDF)
  colnames(summary_cor) <- "correlation"
  
  summary_cor_T <- summary_cor %>% rownames_to_column("project") %>% as_tibble()  
  colnames(summary_cor_T) <- c("project","correlation")
  
  corr_data <-summary_cor_T %>% as.data.frame() %>%  as_tibble() %>% mutate(project=factor(project))
  
  ggtheme = ggplot2::theme_minimal(base_family = "Palatino Linotype", base_size = 10)
  corr_data.m <- melt(corr_data)
  corr_data.m[[yVar]] <-factor(corr_data.m[[yVar]],
                               levels = rev(levels(factor(corr_data.m[[yVar]]))))
  
  # Remove rows where tcr richness is NA
  corr_data.m.c <- corr_data.m[complete.cases(corr_data.m),]
  if (nrow(corr_data.m.c)==0){
    corr_data.m.c <- corr_data.m
    maxCor_ind <-which(abs(corr_data.m.c$value) == max(abs(corr_data.m.c$value)))
    maxCor_tissue <- corr_data.m.c[maxCor_ind ,1]
    maxCor_tissue <- as.character(corr_data.m.c)
    # Color for tissue with max correlation
    a <- rep("black",length(levels(corr_data.m.c[[yVar]])))
    b <- rep("plain",length(levels(corr_data.m.c[[yVar]])))
  }else{
    maxCor_ind <-which(abs(corr_data.m.c$value) == max(abs(corr_data.m.c$value)))
    maxCor_tissue <- corr_data.m.c[maxCor_ind ,1]
    maxCor_tissue <- as.character(maxCor_tissue)
    # Color for tissue with max correlation
    a <- ifelse(levels(corr_data.m[[yVar]]) == maxCor_tissue,"Light Slate Blue",ifelse(levels(corr_data.m[[yVar]]) %in% tissues_fewData, "Chocolate 3","black"))
    b <- ifelse(levels(corr_data.m[[yVar]]) == maxCor_tissue,"bold","plain")
  }
  g <- ggplot(corr_data.m,aes_string("variable",yVar))
  g <- g +  geom_tile(aes_string(fill = "value"), color = "black", size = 0.1)
  g <- g + geom_text(data=corr_data.m,
                     aes_string(x = "variable",
                                y = yVar,
                                label = round(corr_data.m[["value"]], digits = 2)),
                     family = "Palatino Linotype",
                     size = 3.5)
  g <- g + scale_fill_gradient2(low = 'steel blue', mid = 'white', high = 'Indian red', midpoint = 0,limit = c(-1,1))
  
  # +
  #   scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
  #                        midpoint = 0, limit = c(-1,1), space = "Lab", 
  #                        name="Pearson\nCorrelation") 
  # 
  g <- g + xlab(xLabel)
  g <- g + ylab("")
  #g <-g + scale_x_discrete(limits = rev(levels(the_factor)))
  g <- g + theme_minimal(base_family = "Palatino Linotype", base_size = 10)
  g <- g + theme(axis.text.x = element_blank(),#angle = 45, hjust = 1, vjust = 1,
                 axis.text.y = element_text(size = 10,color=a, face = b),
                 axis.title = element_text(size = 12),
                 #legend.position = c(0.33,-0.4),
                 legend.position = "bottom",
                 legend.direction = "horizontal",
                 #legend.justification = c(0,0),
                 plot.margin = margin(t=10,b=20,l=10,r=10, unit="pt"),
                 legend.title = element_text(size=8,color="black")
  )
  
  
  g <- g + ggtitle("")
  g <- g + labs(fill = "")
  print(g)
}



circleNetGrad <- function(corGraph, colorVar, sizeVar){
  
  colorVar.un <- as.name(colorVar)
  sizeVar.un <- as.name(sizeVar)
  p <-ggplot(data = ggnetwork(asNetwork(corGraph), layout = "circle"), aes(x = x, y = y, xend = xend, yend = yend)) + #, cell.jitter = 0.75
    geom_edges(aes(colour=!!colorVar.un,size = abs(!!colorVar.un)),curvature = 0.1) + 
    scale_colour_gradient2(low="navyblue",mid = "white",high = "indianred", midpoint=0,limits = c(-1,1)) +
    scale_size_continuous(range  = c(0.1, 5), 
                          limits = c(0, 1), 
                          breaks = c(0, .25, .5, .75, 1)) + 
    # scale_size_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c("0","0.25","0.5","0.75","1")) + 
    geom_nodes(aes(x, y)) +#, size = 8.5, color = "black", fill="white"
    # geom_nodetext(aes(label = vertex.names), size = 2.5,fontface="bold",nudge_x = +0,nudge_y = -0) + 
    geom_nodelabel(aes(label =vertex.names ),
                   fontface = "bold",nudge_x = +0.02) + 
    guides(size = guide_legend(order = 2),col = guide_colorbar(order = 1), col=guide_legend(title.position = "top")) + #,title.position = "right"
    labs(colour="Spearman's Correlation", size="") +# guides(size = "none") +
    # guides(size = guide_legend(order = 2),col = guide_legend(order = 1)) + 
    theme_blank()
  p
}


circleNetNograd <- function(corGraph, colorVar, sizeVar){
  
  colorVar.un <- as.name(colorVar)
  sizeVar.un <- as.name(sizeVar)
  p <-ggplot(data = ggnetwork(asNetwork(corGraph), layout = "circle"), aes(x = x, y = y, xend = xend, yend = yend)) + #, cell.jitter = 0.75
    geom_edges(aes(size = abs(!!sizeVar.un), color=ifelse(!!colorVar.un>0, 'Indian Red', 'navyblue')),curvature = 0) +
    # scale_color_identity() +
    scale_color_manual(name="Spearman Correlation", labels = c("> 0","< 0"),values=c('Indian Red','navyblue'))+ #, breaks = c("> 0", "< 0")
    scale_size_continuous(range  = c(0.1, 5), 
                          limits = c(0, 1), 
                          breaks = c(0, .25, .5, .75, 1)) + 
    geom_nodes(aes(x, y)) +#, size = 8.5, color = "black", fill="white"
    # geom_nodetext(aes(label = vertex.names), size = 2.5,fontface="bold",nudge_x = +0,nudge_y = -0) + 
    geom_nodelabel(aes(label =vertex.names ),
                   fontface = "bold",nudge_x = +0,nudge_y = -0) +
    guides(size = guide_legend(order = 2),col = guide_legend(order = 1,override.aes = list(size = 5)), col=guide_legend(title.position = "top")) + #,title.position = "right"
    labs(colour="Spearman's Correlation", size="") +# guides(size = "none") +
    geom_edgetext(aes(label = round(r_if_sig,2)), color = "black", fill = "white", size=2.5) + #,box.padding = unit(1, "lines")#,label.size = 0.2
    geom_edgelabel(aes(label = round(r_if_sig,2)),position = "identity",size=2.2) + #,label.size = 0.1
    theme_blank()
  p
}

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
  tissue_dataPoints <- tcgaDF.clean  %>% dplyr::group_by(!!facetVariable.un) %>% 
    distinct(!!SampleCol,.keep_all = TRUE) %>%
    dplyr::filter(!any(is.na(!!corVarx.un))) %>%
    dplyr::summarize(no_rows=length(!!facetVariable.un)) %>% dplyr::mutate(comboVar=paste0(!!facetVariable.un, "-", corVarx))
  # tissue_dataPoints
  # tissues_fewData <- tissue_dataPoints %>% dplyr::filter(no_rows<10) %>% dplyr::select(!!facetVariable.un) %>% pull(!!facetVariable.un) %>% droplevels()
  tissue_dataPoints
  
  numb_DataPoints_df <- rbind(numb_DataPoints_df,tissue_dataPoints)
  numb_DataPoints_df
}


summaryCor_bubblePlot <- function(corrDFmelted,xVar, yVar, bubbleSizeVar,groupVar,groupVarSelection, xLabel,legendSizeLabel){ #xvar:value, yvar=project
  # Prepare variables
  xVar.un <- as.name(xVar)
  yVar.un <- as.name(yVar)
  bubbleSizeVar.un <- as.name(bubbleSizeVar)
  groupVar.un <- as.name(groupVar)
  ## Max Correlation Tissues
  corrDFmelted <-corrDFmelted%>% mutate(!!yVar:=factor(!!yVar.un))
  maxCorInd <-corrDFmelted %>% dplyr::filter(complete.cases(!!xVar.un)) %>% group_by(!!groupVar.un) %>% mutate(project:=as.character(!!yVar.un)) %>% summarize(maxCor = max(abs(!!xVar.un)), TissueRowInd=which.max(abs(!!xVar.un)))
  
  maxCorInSpecific <- maxCorInd %>% dplyr::filter(!!groupVar.un==groupVarSelection) %>% pull() 
  
  maxCor_tissue <-corrDFmelted %>% dplyr::filter(complete.cases(!!xVar.un)) %>% group_by(!!groupVar.un) %>% dplyr::filter(!!groupVar.un==groupVarSelection, row_number()==maxCorInSpecific) %>% pull(!!yVar.un)
  ## Tissues with few values
  fewData_Df <-corrDFmelted %>% group_by(!!groupVar.un) %>% mutate(project:=as.character(!!yVar.un)) %>% mutate(fewData = ifelse(!!bubbleSizeVar.un<10,1,0)) %>% dplyr::filter(fewData==1)
  fewDataTissue <-fewData_Df %>% dplyr::filter(!!groupVar.un==groupVarSelection) %>% dplyr::select(!!yVar.un) %>% pull()
  
  # Color for tissue with max correlation
  a <- ifelse(levels(corrDFmelted[[yVar]]) == maxCor_tissue,"Light Slate Blue",ifelse(levels(corrDFmelted[[yVar]]) %in% fewDataTissue, "Chocolate 3","black"))
  b <- ifelse(levels(corrDFmelted[[yVar]]) == maxCor_tissue,"bold","plain")
  
  ggtheme = ggplot2::theme_minimal(base_family = "Palatino Linotype", base_size = 11)
  
  g <- ggplot(filter(corrDFmelted,!!groupVar.un==groupVarSelection), aes(!!xVar.un,!!yVar.un))
  g <- g + geom_point(aes(size=!!bubbleSizeVar.un),shape=21,fill="white")
  g <- g + geom_text(aes(label=!!xVar.un),size=3)
  g <- g + scale_size(range = c(6.5, 20.5),breaks = c(0,25,50,100,200,400)) #range = c(4, 18)#+ scale_size_area(breaks = c(0,25,50,100,200,400,800))
  g <- g + ggtheme
  g <- g + theme(axis.text.x = element_text(hjust = 1, vjust = 1,size = 11),#angle = 90,
                 axis.title.x = element_text(size=11),
                 #axis.text.y = element_text( hjust = 1, vjust = 1,size = 6),
                 axis.text.y = element_text(size = 11,color=a, face = b),
                 axis.title.y = element_text(size=11),
                 legend.title = element_text(size=11,color="black"),
                 legend.text = element_text(size = 11),
                 #panel.grid.major.y=element_line(linetype=2,color="black"),
                 #panel.grid.major = element_blank(),
                 #panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"),
                 legend.box.background = element_rect(colour = "black"),
                 legend.position="right", legend.box = "vertical")
  g <- g + coord_cartesian(xlim = c(-1,1))
  g <- g + xlab(xLabel)
  g <- g + ylab("")
  g <- g + ggtitle("")
  g <- g + labs(size = legendSizeLabel)
  #g <- g +facet_wrap(~variable,scales = "free", ncol = 5,strip.position = "bottom")
  #g <- g +facet_wrap(~variable, ncol = 5,strip.position = "bottom")
  #g <- g +facet_grid(~variable,scales = "free")
  g
  
}


corPlot2 <- function(DF, Xvar, Yvar, xLabel,yLabel, xPlus, yPlus, yMinus,Logx= FALSE, filterVar=NULL, filterVarSelected="", shapeVar=NULL, shapeVarLabel="", facetVar=NULL, plotCols=NULL, addCorLabel=NULL){
  # x <- 8
  # DF <-tcgaData.nt.compl.filt
  # Xvar <-as.character(topCorr_nt[x,2])
  # Yvar <-"TCR_Richness"
  # xLabel <-as.character(topCorr_nt[x,4])
  # yLabel <-"Richness"
  # xPlus <-0
  # yPlus <-0
  # yMinus <-0
  # Logx= FALSE
  # filterVar="project"
  # filterVarSelected=as.character(topCorr_nt[x,1])
  # shapeVar=NULL
  # shapeVarLabel=""
  # facetVar=NULL
  # plotCols=NULL
  # addCorLabel=NULL
  
  if (!is.null(filterVar)){
    filterVar.un <- as.name(filterVar)
    DF <- DF %>% dplyr::filter(!!filterVar.un %in% filterVarSelected) %>% droplevels()
    DF
  }
  
  
  xTick <- mean(DF[[Xvar]][!is.na(DF[[Xvar]])])
  yPlus <-max(DF[[Yvar]])/10
  yTick<- max(DF[[Yvar]])+yPlus
  yMinus <-max(DF[[Yvar]])/20
  
  # print(xLabel)
  # print(min(DF[[Xvar]]))
  # print(mean(DF[[Xvar]]))
  # print(xTick)
  ## PLOT
  # plot.new()
  if(is.null(shapeVar)){
    p1 <- ggplot(DF,aes_string(Xvar,Yvar))
  }else{
    p1 <- ggplot(DF,aes_string(Xvar,Yvar, shape=shapeVar))
  }
  
  if (Logx==TRUE){
    p1 <- p1 + scale_x_continuous(trans='log10')
  }
  
  p1 <- p1 + geom_point(size = 2, show.legend = TRUE, alpha = .4)#,color="#69b3a2"
  p1 <- p1 + xlab(paste0(xLabel))+ylab(yLabel)
  p1 <- p1 + geom_smooth(method=lm, se=FALSE, color="black",linetype = "dashed",size=0.5)
  p1 <- p1 + ggtitle(filterVarSelected)
  p1 <- p1 + theme_ipsum(base_family = "Palatino Linotype", base_size = 10)
  
  if(is.null(shapeVar)){
    p1 <- p1
  }else{
    p1 + scale_color_discrete(name=shapeVarLabel)
  }
  
  if(!is.null(facetVar)){
    p1 <- p1 + theme(axis.text.x = element_text(hjust = 1, vjust = 1,size = 10),
                     axis.title.x = element_text(size=16,vjust = 0.5, hjust=0.5),
                     axis.text.y = element_text( hjust = 1, vjust = 1,size = 10), 
                     axis.title.y = element_text(size=16,angle=90, vjust = 0.5,hjust=0.5),
                     legend.title = element_text(size=13,color="black"),
                     legend.text = element_text(size = 13),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"))
    p1 <- p1 + facet_wrap(facetVar,ncol=plotCols,scales='free')
  }
  
  if (is.null(addCorLabel)){
    p1 <- p1 + theme(axis.text.x = element_text(hjust = 1, vjust = 1,size = 13),
                     axis.title.x = element_text(size=12,),
                     axis.text.y = element_text( hjust = 1, vjust = 1,size = 13), 
                     axis.title.y = element_text(size=13,angle=90),
                     legend.title = element_text(size=13,color="black"),
                     legend.text = element_text(size = 13),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"))
    # Pearson
    corTest.pearson <- cor.test(DF[[Xvar]],DF[[Yvar]],method = "pearson")
    
    l.pearson <- corTest.pearson$p.value
    l.pearson <- format(round(l.pearson,4), scientific = TRUE)
    l.pearson <-strsplit(l.pearson,"e")
    # quote the part before the exponent to keep all the digits
    l.pearson1 <- l.pearson[[1]][1]
    l.pearson2 <- gsub("-0","-",l.pearson[[1]][2])
    
    # annot.pearson <- annotate("text", x = min(DF[[Xvar]]) + xPlus, y = max(DF[[Yvar]])+yPlus, label = bquote(paste("Pearson R" == .(round( corTest.pearson$estimate,2))," (P" < .(l.pearson1) %*% 10^ .(l.pearson2),")")),fontface ="plain",size = 4,family= "Palatino Linotype")
    annot.pearson <- ggplot2::annotate("text", x = xTick + xPlus, y = yTick, label = bquote(paste("Pearson R" == .(round( corTest.pearson$estimate,2))," (P" < .(l.pearson1) %*% 10^ .(l.pearson2),")")),fontface ="plain",size = 4,family= "Palatino Linotype")
    
    # Spearman
    corTest.spearman <- cor.test(DF[[Xvar]],DF[[Yvar]],method = "spearman")
    
    l.spearman <- corTest.spearman$p.value
    l.spearman <- format(round(l.spearman,4), scientific = TRUE)
    l.spearman <-strsplit(l.spearman,"e")
    # quote the part before the exponent to keep all the digits
    l.spearman1 <- l.spearman[[1]][1]
    l.spearman2 <- gsub("-0","-",l.spearman[[1]][2])
    
    annot.spearman <- ggplot2::annotate("text", x = xTick + xPlus, y = yTick-yMinus, label = bquote(paste("Spearman R" == .(round( corTest.spearman$estimate,2))," (P" < .(l.spearman1) %*% 10^ .(l.spearman2),")")),fontface ="plain",size = 4,family= "Palatino Linotype")
    
    p1 <- p1 + annot.pearson
    p1 <- p1 + annot.spearman
  }
  p1 <- p1 + scale_y_continuous(breaks = scales::pretty_breaks(n = 3))
  p1 <- p1 + scale_x_continuous(breaks = scales::pretty_breaks(n = 3))
  p1
  # #p1 <- recordPlot()
  # recordedPlot = recordPlot()
}

###########################################################
## Correlation plot,1 TCGA dataset ##
###########################################################
corPlot_tcga_var <- function(tcgaDF,SampleCol,tcgaVar, tcgaProject,xVar,yVar, xLabel,yLabel,xPlus,yPlus,yMinus,logxLabel=FALSE){
  
  # Prepare variables
  SampleCol <- as.name(SampleCol)
  tcgaVar.un <- as.name(tcgaVar)
  xVar.un <- as.name(xVar)
  yVar.un <- as.name(yVar)
  #tcgaDF <- subset(tcgaDF, tcgaDF[[tcgaVar]]== tcgaProject)
  
  tcgaDF <- tcgaDF %>% dplyr::select(!!SampleCol,!!tcgaVar.un,!!xVar.un,!!yVar.un)%>%filter(!!tcgaVar.un == tcgaProject) %>% dplyr::filter(complete.cases(!!xVar.un,!!yVar.un))%>% distinct(!!SampleCol,.keep_all = TRUE)
  
  # Plot
  p1 <- ggplot(tcgaDF,aes_string(xVar,yVar))
  p1 <- p1 + geom_point(size = 2, show.legend = TRUE, alpha = .4)#,color="#69b3a2"
  p1 <- p1 + xlab(xLabel)+ylab(yLabel)
  p1 <- p1 + ggtitle(tcgaProject)
  # p1 <- p1 + geom_smooth(method=lm, se=FALSE, color="black",linetype = "dashed",size=0.5)
  # p1 <- p1 + geom_smooth(se=FALSE, color="black",linetype = "dashed",size=0.5)
  p1 <- p1 + theme_ipsum(base_family = "Palatino Linotype", base_size = 10)
  p1 <- p1 + theme(axis.text.x = element_text(hjust = 1, vjust = 1,size = 10),
                   axis.title.x = element_text(size=13),
                   axis.text.y = element_text( hjust = 1, vjust = 1,size = 10), 
                   axis.title.y = element_text(size=13),
                   legend.title = element_text(size=10,color="black"),
                   legend.text = element_text(size = 10),
                   plot.title = element_text(hjust = 0.5,size = 10),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"))
  p1 <- p1 + scale_color_discrete(name="TCGA project")
  # Add correlation Info
  ## Spearman
  # Calculate correlation
  cor2 <- cor.test(tcgaDF[[xVar]],tcgaDF[[yVar]],method = "spearman",exact=FALSE)
  # cor2_se <- round(unname(sqrt((1 - cor2$estimate^2)/cor2$parameter)),4)
  # data.label <- data.frame(label=paste("r = ", round(cor2$estimate,4),"\u00B1",cor2_se), x=0.1*max(tcgaDF[[xVar]],na.rm=T), y=max(tcgaDF[[yVar]], na.rm=T)*.8)
  
  l <- as.character(cor2$p.value)
  #l <- format(round(l), scientific = TRUE)
  l <-strsplit(l,"e")
  # quote the part before the exponent to keep all the digits
  l1 <- as.character(round(as.numeric(l[[1]][1]),2))
  l2 <- gsub("-0","-",l[[1]][2])
  if(is.na(l2)){
    l2 <- 0
  }
  p1<- p1 + ggplot2::annotate("text", x = min(tcgaDF[[xVar]]) + xPlus, y = max(tcgaDF[[yVar]])+yPlus, label = bquote(paste("Spearman R" == .(round(cor2$estimate,2))," (P" == .(l1) %*% 10^ .(l2),")")),fontface ="plain",size = 4,family= "Palatino Linotype")
  
  ## Pearson
  # Calculate correlation
  cor2 <- cor.test(tcgaDF[[xVar]],tcgaDF[[yVar]],method = "pearson")
  # cor2_se <- round(unname(sqrt((1 - cor2$estimate^2)/cor2$parameter)),4)
  # data.label <- data.frame(label=paste("r = ", round(cor2$estimate,4),"\u00B1",cor2_se), x=0.1*max(tcgaDF[[xVar]],na.rm=T), y=max(tcgaDF[[yVar]], na.rm=T)*.8)
  l <- as.character(cor2$p.value)
  #l <- format(round(l), scientific = TRUE)
  l <-strsplit(l,"e")
  # quote the part before the exponent to keep all the digits
  l1 <- as.character(round(as.numeric(l[[1]][1]),2))
  l2 <- gsub("-0","-",l[[1]][2])
  if(is.na(l2)){
    l2 <- 0
  }
  p1<- p1 + ggplot2::annotate("text", x = min(tcgaDF[[xVar]]) + xPlus, y = max(tcgaDF[[yVar]])+yPlus-yMinus, label = bquote(paste("Pearson R" == .(round(cor2$estimate,2))," (P" == .(l1) %*% 10^ .(l2),")")),fontface ="plain",size = 4,family= "Palatino Linotype")
  
  if (logxLabel==TRUE){
    p1 <- p1 +scale_x_continuous(trans='log10')
  }
  p1
  
}

##############################################################
## Boxplot in all tcga datasets, Anything Kvs Any signature ##
##############################################################

TwoVars_RankTCGA_boxplot<-function(tcgaDF,selectedCancers,SampleCol,tcgaVar,xVar,yVar,aggregateMethod,aggregateMethodChar, xLabel, yLabel,ymin1,ymin2,xmin1,xmin2,logyLabel=FALSE){
  # tcgaDF <-tcgaData.tp.compl.filt
  # selectedCancers <-mycancertypes.reduced
  # SampleCol <-"bcr_patient_barcode"
  # tcgaVar <-"project"
  # xVar <-"TCR_Richness"
  # yVar <-"TCR_Richness"
  # aggregateMethod <-median
  # aggregateMethodChar <-"median"
  # xLabel <-"Richness"
  # yLabel <-"Richness"
  # ymin1 <- -1200
  # ymin2 <- -1400
  # xmin1 <- 17.5
  # xmin2 <- 35
  # logyLabel=FALSE
  
  # Prepare variables
  SampleCol <- as.name(SampleCol)
  xVar.un <- as.name(xVar)
  yVar.un <- as.name(yVar)
  tcgaVar.un <- as.name(tcgaVar)
  # Clean DF, make columns factors, remove NA rows for x-z variables and keep one row per sample ID
  tcgaDF.clean <- tcgaDF %>% mutate(!!tcgaVar:=factor(!!tcgaVar.un)) %>% distinct(!!SampleCol,.keep_all = TRUE)%>% dplyr::filter(!!tcgaVar.un %in% selectedCancers)%>% 
    droplevels()
  # tcgaDF_sub <- subset(tcgaDF.clean,tcgaDF.clean[[tcgaVar]] %in% selectedCancers)
  
  # order factors -project by the mean/median value of variable
  # Add mean/median of variable for each cancer
  # Have richness on x axis always
  tcga_aggregDF <- aggregate(tcgaDF.clean[[xVar]]~tcgaDF.clean[[tcgaVar]], FUN=aggregateMethod)
  colnames(tcga_aggregDF) <- c("type", "xvariable")
  tcga_aggregDF.ord <- tcga_aggregDF[order(tcga_aggregDF$xvariable),]
  #Order the factor levels in the big df - order the project
  tcgaDF.clean[[tcgaVar]] <- factor(tcgaDF.clean[[tcgaVar]],tcga_aggregDF.ord[["type"]])
  # Plot
  p <- ggplot(tcgaDF.clean, aes(x = !!tcgaVar.un, y = !!yVar.un))
  p <- p + geom_boxplot(alpha=0.3,outlier.shape = NA, size = 0.1, fill = "black", varwidth = FALSE)
  # p <- p + geom_beeswarm(data = tcgaDF,aes(x = !!tcgaVar.un, y = !!yVar.un), dodge.width=.2,cex=.1,alpha=0.5,size=.4)
  p <- p + geom_point(size = 0.4)
  p <- p + labs(
    title = paste0(""),
    x = "",
    y = yLabel
  )
  p <- p + theme_ipsum(base_family = "Palatino Linotype", base_size = 10)
  p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size = 10),
                 axis.title.x = element_text(size=13, face="bold"),
                 axis.text.y = element_text( hjust = 1, vjust = 1,size = 10),
                 axis.title.y = element_text(size=13, face="bold"),
                 legend.title = element_text(size=13,color="black"),
                 legend.text = element_text(size = 10),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"))
  p <- p + scale_color_discrete(name="TCGA project")
  if (logyLabel==TRUE){
    p <- p +scale_y_continuous(trans='log10')
  }
  # Add annotation labels on x-axix
  Text1 <- textGrob("low",gp = gpar(fontsize = 11,fontface="bold"))
  Text2 <- textGrob("medium",gp = gpar(fontsize = 11,fontface="bold"))
  Text3 <- textGrob("high",gp = gpar(fontsize = 11,fontface="bold"))
  Text4 <- textGrob(paste0(firstup(aggregateMethodChar)," ",xLabel),gp = gpar(fontsize = 11,fontface="bold"))
  (g1 <- p + annotation_custom(grob = Text1,  xmin = 1, xmax = 1, ymin = ymin2, ymax = ymin2) +
      annotation_custom(grob = Text2,  xmin = xmin1, xmax = xmin1, ymin = ymin2, ymax = ymin2) +
      annotation_custom(grob = Text3,  xmin = xmin2, xmax = xmin2, ymin = ymin2, ymax = ymin2)+
      annotation_custom(grob = Text4,  xmin = xmin1, xmax = xmin1, ymin = ymin1, ymax = ymin1))
  
  gg_table <- ggplot_gtable(ggplot_build(g1))
  gg_table$layout$clip[gg_table$layout$name=="panel"] <- "off"
  #grid.draw(gg_table)
  gg_table
}


TwosigScores_boxplot_tcga_old <- function(tcgaDF,selectedCancers,sigName,aggregateMethod,aggregateMethodChar, yVar,yLabel,ymin1,ymin2){
  # Subset to allowed datasets/subsets
  tcgaDF <- subset(tcgaDF,tcgaDF$type %in% selectedCancers)
  
  tcgaDF$type <- as.factor(tcgaDF$type)
  # order factors -project by the mean/median value of signature
  # Add mean/median of signature for each cancer
  # Have expanded immune signature on x-axis
  tcga_sig_aggregDF <- aggregate(tcgaDF[[sigName]]~tcgaDF$type, FUN=aggregateMethod)
  colnames(tcga_sig_aggregDF) <- c("type", "sigScore")
  tcga_sig_aggregDF.ord <- tcga_sig_aggregDF[order(tcga_sig_aggregDF$sigScore),]
  #Order the factor levels in the big df - order the project
  tcgaDF$type <- factor(tcgaDF$type,tcga_sig_aggregDF.ord$type)
  # Plot
  p <- ggplot(tcgaDF, aes_string(x = "type", y = yVar))
  p <- p + geom_boxplot(alpha=0.3,outlier.shape = NA, size = 0.1, fill = "black", varwidth = FALSE)
  p <- p + geom_beeswarm(data = tcgaDF,aes_string(x = "type", y = yVar), dodge.width=.2,cex=.1,alpha=0.5,size=.4)
  p <- p + geom_point(size = 0.4)
  p <- p + labs(
    title = paste0(""),
    x = "",
    y = yLabel
  )
  p <- p + theme_ipsum(base_family = "Palatino Linotype", base_size = 10)
  p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size = 10),
                 axis.title.x = element_text(size=13, face="bold"),
                 axis.text.y = element_text( hjust = 1, vjust = 1,size = 10),
                 axis.title.y = element_text(size=13, face="bold"),
                 legend.title = element_text(size=13,color="black"),
                 legend.text = element_text(size = 10),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"))
  p <- p + scale_color_discrete(name="TCGA project")
  
  # Add annotation labels on x-axix
  Text1 <- textGrob("low",gp = gpar(fontsize = 11,fontface="bold", fontfamily="Palatino Linotype"))
  Text2 <- textGrob("medium",gp = gpar(fontsize = 11,fontface="bold", fontfamily="Palatino Linotype"))
  Text3 <- textGrob("high",gp = gpar(fontsize = 11,fontface="bold", fontfamily="Palatino Linotype"))
  Text4 <- textGrob(paste0(firstup(aggregateMethodChar)," ",yLabel),gp = gpar(fontsize = 11,fontface="bold"))
  (g1 <- p + annotation_custom(grob = Text1, xmin = 1, xmax = 1, ymin = ymin2, ymax = ymin2) +
      annotation_custom(grob = Text2, xmin = 16, xmax = 16, ymin = ymin2, ymax = ymin2) +
      annotation_custom(grob = Text3, xmin = 31, xmax = 31, ymin = ymin2, ymax = ymin2)+
      annotation_custom(grob = Text4, xmin = 16, xmax = 16, ymin = ymin1, ymax = ymin1))
  
  gg_table <- ggplot_gtable(ggplot_build(g1))
  gg_table$layout$clip[gg_table$layout$name=="panel"] <- "off"
  #grid.draw(gg_table)
  gg_table
  
}

###############################################################
## Correlation plot, ORR vsANything, ALL avail ORR TCGA data ##
###############################################################
orr_corrPlot_tcga <- function(tcgaDF,xVariable, aggregateMethod,aggregateMethodChar, 
                              selectedORRCancers, orrDF, xLabel, 
                              xPlus,yPlus, yMinus,CorrMethod,Logx=FALSE){
  # tcgaDF <- tcgaData.tp.compl.filt
  # xVariable<-"CD8_T_cells"
  # aggregateMethod<-median
  # aggregateMethodChar<-"median"
  # selectedORRCancers<-orr_selected_cancers
  # orrDF<-tcgaData.tp.compl.filt.orr
  # xLabel<-"CD8+ T cells"
  # xPlus<-0.3
  # yPlus<-0.08
  # yMinus<-0.015
  # CorrMethod<- "pearson"
  # Logx=FALSE
  # Add mean of TuTACK signature for each cancer
  tcga_orr_aggregDF <-aggregate(tcgaDF[[xVariable]]~tcgaDF$project, FUN=aggregateMethod)
  colnames(tcga_orr_aggregDF) <- c("type","variable")
  tcga_orr_aggregDF.ord <- tcga_orr_aggregDF[order(tcga_orr_aggregDF$variable),]
  # subset to selected TCGA datasets
  tcga_orr_aggregDF_final <- subset(tcga_orr_aggregDF.ord,tcga_orr_aggregDF.ord$type %in% selectedORRCancers)
  colnames(tcga_orr_aggregDF_final) <- c("project", xVariable)
  # Merge variable with ORR
  orr_var_df <- merge(orrDF, tcga_orr_aggregDF_final, by="project")
  ## PLOT
  p1 <- ggplot(orr_var_df,aes_string(xVariable,"ORR", label="project"))
  if (Logx==TRUE){
    p1 <- p1 + scale_x_continuous(trans='log10')
  }
  p1 <- p1 + geom_point(aes_string(size = "PD1L1_totalTreated"),shape = 16, show.legend = TRUE, alpha = .4)#,color="#69b3a2"
  p1 <- p1 + geom_text_repel(size = 3, 
                             segment.size=0.2, 
                             nudge_x = 0.015,
                             point.padding = unit(8.1, "points"))
  p1 <- p1 + xlab(paste0(firstup(aggregateMethodChar)," ",xLabel))+ ylab("ORR")
  # p1 <- p1 + geom_smooth(method=lm, 
  #                        se=FALSE, 
  #                        color="Gray",
  #                        linetype = "dashed",size=0.3)
  p1 <- p1 + theme_ipsum(base_family = "Palatino Linotype", base_size = 10)
  p1 <- p1 + theme(axis.text.x = element_text(hjust = 1, vjust = 1,size = 10),
                   axis.title.x = element_text(size=13),
                   axis.text.y = element_text( hjust = 1, vjust = 1,size = 10), 
                   axis.title.y = element_text(size=13),
                   legend.title = element_text(size=13,color="black"),
                   legend.text = element_text(size = 10),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   legend.box.background = element_rect(colour = "black"),
                   legend.position="right", legend.box = "vertical")
  p1 <- p1 + scale_color_discrete(name="TCGA project")
  p1 <- p1 + scale_size_area(breaks = c(0,25,50,100,200,400,800))
  p1 <- p1 + labs(size = "No. of \npatients \nevaluated \n(ORR)")
  
  ## Spearman
  # Calculate correlation
  cor2 <- cor.test(orr_var_df[,xVariable],orr_var_df[,"ORR"],method = "spearman")
  cor2_se <- round(unname(sqrt((1 - cor2$estimate^2)/cor2$parameter)),2)
  data.label <- data.frame(label=paste("r = ", round(cor2$estimate,2),"\u00B1",cor2_se), x=0.1*max(orr_var_df[[xVariable]],na.rm=T), y=max(orr_var_df$ORR, na.rm=T)*.8)
  
  l <- cor2$p.value
  l <- format(round(l,4), scientific = TRUE)
  l <-strsplit(l,"e")
  # quote the part before the exponent to keep all the digits
  l1 <- l[[1]][1]
  l2 <- gsub("-0","-",l[[1]][2])
  
  p1<- p1 + ggplot2::annotate("text", x = min(orr_var_df[[xVariable]]) + xPlus, y = max(orr_var_df$ORR)+ yPlus, label = bquote(paste("Spearman R" == .(round(cor2$estimate,2))," (P" < .(l1) %*% 10^ .(l2),")")),fontface ="plain",size = 4,family= "Palatino Linotype")
  
  ## Pearson
  # Calculate correlation
  cor2 <- cor.test(orr_var_df[,xVariable],orr_var_df[,"ORR"],method = "pearson")
  cor2_se <- round(unname(sqrt((1 - cor2$estimate^2)/cor2$parameter)),2)
  data.label <- data.frame(label=paste("r = ", round(cor2$estimate,2),"\u00B1",cor2_se), x=0.1*max(orr_var_df[[xVariable]],na.rm=T), y=max(orr_var_df$ORR, na.rm=T)*.8)
  
  l <- cor2$p.value
  l <- format(round(l,4), scientific = TRUE)
  l <-strsplit(l,"e")
  # quote the part before the exponent to keep all the digits
  l1 <- l[[1]][1]
  l2 <- gsub("-0","-",l[[1]][2])
  
  p1<- p1 + ggplot2::annotate("text", x = min(orr_var_df[[xVariable]]) + xPlus, y = max(orr_var_df$ORR)+ yPlus-yMinus, label = bquote(paste("Pearson R" == .(round(cor2$estimate,2))," (P" < .(l1) %*% 10^ .(l2),")")),fontface ="plain",size = 4,family= "Palatino Linotype")
  p1
  # print(p1) 
}


###############################################################
## Correlation plot, ANYTHING vsANything, TCGA data, mean/median          ##
###############################################################
corrPlot_tcga <- function(tcgaDF,selectedcancers,xVariable,yVariable,SampleCol,facetVariable, xVarAggregateMethod,xVarAggregateMethodChar, yVarAggregateMethod,yVarAggregateMethodChar, xLabel,yLabel, xPlus,yPlus, yMinus,Logx=FALSE){
  # Prepare variables
  SampleCol <- as.name(SampleCol)
  xVariable.un <- as.name(xVariable)
  yVariable.un <- as.name(yVariable)
  facetVariable.un <- as.name(facetVariable)
  # Clean DF, make columns factors, remove NA rows for x-z variables and keep one row per sample ID
  tcgaDF.clean <- tcgaDF %>% mutate(!!facetVariable:=factor(!!facetVariable.un)) %>% dplyr::filter(complete.cases(!!yVariable.un, !!xVariable.un)) %>% dplyr::filter(!!facetVariable.un %in% selectedcancers) %>% 
    droplevels()%>% distinct(!!SampleCol,.keep_all = TRUE)
  
  # Find number of datapoints for each cancer type
  tissue_dataPoints <- tcgaDF.clean  %>% group_by(!!facetVariable.un) %>% 
    distinct(!!SampleCol,.keep_all = TRUE) %>%
    filter(!any(is.na(!!xVariable.un))) %>%
    summarize(no_rows=length(!!facetVariable.un))
  
  tcga_aggDF.xVar <-aggregate(tcgaDF.clean[[xVariable]]~tcgaDF.clean[[facetVariable]], FUN=xVarAggregateMethod)
  colnames(tcga_aggDF.xVar) <- c("type", "xVar")
  tcga_aggDF.yVar <-aggregate(tcgaDF.clean[[yVariable]]~tcgaDF.clean[[facetVariable]], FUN=yVarAggregateMethod)
  colnames(tcga_aggDF.yVar) <- c("type", "yVar")
  tcga_aggDF <- merge(tcga_aggDF.xVar,tcga_aggDF.yVar, by="type")
  ## PLOT
  p1 <- ggplot(tcga_aggDF,aes_string("xVar","yVar", label="type"))
  p1 <- p1 + geom_point(size=3,shape = 16, show.legend = TRUE, alpha = .4)#,color="#69b3a2"
  p1 <- p1 + geom_text_repel(size = 3,
                             segment.size=0.2,
                             nudge_x = 0.015,
                             point.padding = unit(8.1, "points"))
  p1 <- p1 + xlab(paste0(firstup(xVarAggregateMethodChar)," ",xLabel))
  p1 <- p1 + ylab(paste0(firstup(yVarAggregateMethodChar)," ",yLabel))
  p1 <- p1 + geom_smooth(method=lm,
                         se=FALSE,
                         color="Gray",
                         linetype = "dashed",size=0.3)
  p1 <- p1 + theme_ipsum(base_family = "Palatino Linotype", base_size = 10)
  p1 <- p1 + theme(axis.text.x = element_text(hjust = 1, vjust = 1,size = 10),
                   axis.title.x = element_text(size=12),
                   axis.text.y = element_text( hjust = 1, vjust = 1,size = 10),
                   axis.title.y = element_text(size=12),
                   legend.title = element_text(size=12,color="black"),
                   legend.text = element_text(size = 10),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   axis.line = element_line(colour = "black"),
                   legend.box.background = element_rect(colour = "black"),
                   legend.position="right", legend.box = "vertical")
  p1 <- p1 + scale_color_discrete(name="TCGA project")
  
  if (Logx==TRUE){
    p1 <- p1 + scale_x_continuous(trans='log10')
  }
  ## Spearman
  # Calculate correlation
  cor2 <- cor.test(tcga_aggDF[["xVar"]],tcga_aggDF[["yVar"]],method = "spearman")
  cor2_se <- round(unname(sqrt((1 - cor2$estimate^2)/cor2$parameter)),2)
  data.label <- data.frame(label=paste("r = ", round(cor2$estimate,2),"\u00B1",cor2_se), x=0.1*max(tcga_aggDF[["xVar"]],na.rm=T), y=max(tcga_aggDF[["yVar"]], na.rm=T)*.8)
  
  l <- cor2$p.value
  l <- format(round(l,4), scientific = TRUE)
  l <-strsplit(l,"e")
  # quote the part before the exponent to keep all the digits
  l1 <- l[[1]][1]
  l2 <- gsub("-0","-",l[[1]][2])
  
  p1<- p1 + ggplot2::annotate("text", x = min(tcga_aggDF[["xVar"]]) + xPlus, y = max(tcga_aggDF[["yVar"]])+yPlus, label = bquote(paste("Spearman R" == .(round(cor2$estimate,2))," (P" < .(l1) %*% 10^ .(l2),")")),fontface ="plain",size = 4,family= "Palatino Linotype")
  
  ## Pearson
  # Calculate correlation
  cor2 <- cor.test(tcga_aggDF[["xVar"]],tcga_aggDF[["yVar"]],method = "pearson")
  cor2_se <- round(unname(sqrt((1 - cor2$estimate^2)/cor2$parameter)),2)
  data.label <- data.frame(label=paste("r = ", round(cor2$estimate,2),"\u00B1",cor2_se), x=0.1*max(tcga_aggDF[["xVar"]],na.rm=T), y=max(tcga_aggDF[["yVar"]], na.rm=T)*.8)
  
  l <- cor2$p.value
  l <- format(round(l,4), scientific = TRUE)
  l <-strsplit(l,"e")
  # quote the part before the exponent to keep all the digits
  l1 <- l[[1]][1]
  l2 <- gsub("-0","-",l[[1]][2])
  
  p1<- p1 + ggplot2::annotate("text", x = min(tcga_aggDF[["xVar"]]) + xPlus, y = max(tcga_aggDF[["yVar"]])+yPlus-yMinus, label = bquote(paste("Pearson R" == .(round(cor2$estimate,2))," (P" < .(l1) %*% 10^ .(l2),")")),fontface ="plain",size = 4,family= "Palatino Linotype")
  
  print(p1)
}


correlationsSectionALL <- function(DF, paramA, paramB, labelA, labelB, facets,facetParam, filterParam, filterSelect, shapeParam,
                                   plotNumCols, addCorr, logXaxis,SampleCol,saveP=FALSE){
  
  # paramA <- "TCR_Richness" # fixed
  # paramB <- parameters[8]
  # labelA <-"TCR Richness"
  # labelB <-labels[8]
  # DF <-tcgaData.tp.sub.filt
  # facetParam<- "project"
  # filterParam <- "project"
  # filterSelect <- mycancertypes.final
  # shapeParam <- NULL
  # plotNumCols <-2
  # addCorr <-TRUE
  # logXaxis <- FALSE
  # facets <- mycancertypes.final
  # SampleCol <- "bcr_patient_barcode"
  # 
  cat('\n')  
  cat('\n') 
  cat("### ", paste(c(labelA,labelB), collapse = " vs. "),"{.tabset .tabset-fade .tabset-pills} \n") 
  cat('\n') 
  
  cat('\n')  
  cat('\n') 
  cat("#### ", "Correlation Plots","\n") 
  cat('\n') 
  
  pcor<-corPlot(DF, paramA, paramB, labelA, labelB, 2, 0.35, 0.17,Logx= logXaxis, filterVar=filterParam, filterVarSelected=filterSelect, shapeVar=shapeParam, shapeVarLabel="", facetVar=facetParam, plotCols=plotNumCols,addCorLabel=addCorr)
  # print(pcor)
  cat('\n')  
  cat('\n') 
  cat("#### ", "Correlation Heatmap","\n") 
  cat('\n') 
  
  cor.res <- sapply(facets, function(x) corr_calculation(DF, paramA, paramB, "spearman",facetParam, x))
  names(cor.res) <- gsub(".rho","",names(cor.res))
  
  
  pheat.cor <-cor_heatmap(cor.res,DF,"correlation",facetParam,"correlation",facetParam,SampleCol,"Spearman Correlation")
  
  # print(pheat.cor)
  res <- list(multiPlots = pcor, corDF = cor.res, corHeat =pheat.cor)
  #save plot
  if(isTRUE(saveP)){
    save(pcor, file=paste0(projectRDataPath,"tcga.Corr.",paramA,".",paramB,Rdata.suffix))
  }
  
  res
}


corFeatsMultiPlot <- function(DF, idVars,Yvar, xLabel,yLabel ,Logx= FALSE, filterVar=NULL, filterVarSelected="", shapeVar=NULL, shapeVarLabel="", facetVar=NULL, plotCols=NULL, addCorLabel=NULL){
  
  # DF <-tcgaData.tp.sub.filt[,2:11]
  # idVars = c("BCR_Richness","project")
  # yvar <-"BCR_Richness"
  # yLabel <-"BCR Richness"
  # xLabel <- "Feature"
  # Logx= FALSE
  # filterVar="project"
  # filterVarSelected=mycancertypes.final
  # shapeVar=NULL
  # shapeVarLabel=""
  # facetVar="variable"
  # plotCols=2
  # addCorLabel=TRUE
  
  # DF.melted <-split_columns(DF)$continuous %>% data.table::melt(id.vars=idVars)
  
  if (!is.null(filterVar)){
    filterVar.un <- as.name(filterVar)
    DF <- DF %>% dplyr::filter(!!filterVar.un %in% filterVarSelected) %>% droplevels()
    DF[[filterVar]] <- factor(DF[[filterVar]], levels=filterVarSelected)
    DF
  }
  DF.melted <- DF %>% data.table::melt(id.vars=idVars)
  
  
  if(is.null(shapeVar)){
    p1 <- ggplot(DF.melted,aes_string("value",Yvar))
  }else{
    p1 <- ggplot(DF,aes_string("value",Yvar, shape=shapeVar))
  }
  
  if (Logx==TRUE){
    p1 <- p1 + scale_x_continuous(trans='log10')
  }
  
  p1 <- p1 + geom_point(size = 1, show.legend = TRUE, alpha = .4)#,color="#69b3a2"
  p1 <- p1 + xlab(paste0(xLabel))+ylab(yLabel)
  p1 <- p1 + geom_smooth(method=lm, se=FALSE, color="black",linetype = "dashed",size=0.5)
  
  p1 <- p1 + theme_ipsum(base_family = "Palatino Linotype", base_size = 10)
  
  if(is.null(shapeVar)){
    p1 <- p1
  }else{
    p1 + scale_color_discrete(name=shapeVarLabel)
  }
  
  if(!is.null(facetVar)){
    p1 <- p1 + theme(axis.text.x = element_text(hjust = 1, vjust = 1,size = 10),
                     axis.title.x = element_text(size=16,vjust = 0.5, hjust=0.5),
                     axis.text.y = element_text( hjust = 1, vjust = 1,size = 10), 
                     axis.title.y = element_text(size=16,angle=90, vjust = 0.5,hjust=0.5),
                     legend.title = element_text(size=13,color="black"),
                     legend.text = element_text(size = 10),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"))
    p1 <- p1 + facet_wrap(facetVar,ncol=plotCols,scales='free') 
    # p1 <- p1 + facet_grid(facetVar,cols=plotCols,scales='free') 
    
    
  }
  
  if (!is.null(addCorLabel)){
    p1 <- p1 + theme(axis.text.x = element_text(hjust = 1, vjust = 1,size = 10),
                     axis.title.x = element_text(size=12,),
                     axis.text.y = element_text( hjust = 1, vjust = 1,size = 10), 
                     axis.title.y = element_text(size=12,angle=90),
                     legend.title = element_text(size=12,color="black"),
                     legend.text = element_text(size = 10),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"))
    
    p1 <- p1 +
      stat_cor(aes(label = ..r.label..),cor.coef.name="R",label.y.npc="top", label.x.npc = "right",hjust = 1, method = "spearman",size=5)##,vjust=0
  }
  p1 <- p1 + scale_y_continuous(breaks = scales::pretty_breaks(n = 3))
  p1 <- p1 + scale_x_continuous(breaks = scales::pretty_breaks(n = 3))
  #save plot
  save(p1, file=paste0(projectRDataPath,"tcga.CorrMulti.",Yvar,Rdata.suffix))
  print(p1)
}

corFeatsMultiPlotSection <- function(DF, idVars,Yvar, xLabel,yLabel ,Logx= FALSE, filterVar=NULL, filterVarSelected="", 
                                     shapeVar=NULL, shapeVarLabel="", facetVar=NULL, plotCols=NULL, addCorLabel=NULL,sectionVar,section,saveP=FALSE){
  
  
  # 
  cat('\n')  
  cat('\n') 
  cat("#### ", paste(section)," \n") 
  cat('\n') 
  
  sectionVar.un <- as.name(sectionVar)
  DF <- DF %>% dplyr::filter(!!sectionVar.un %in% section) %>% droplevels()
  
  p <-corFeatsMultiPlot(DF, idVars,Yvar, xLabel,yLabel,Logx, filterVar=NULL, filterVarSelected="", 
                        shapeVar, shapeVarLabel, facetVar, plotCols, addCorLabel)
  
  #save plot
  if(isTRUE(saveP)){
    save(p, file=paste0(projectRDataPath,"tcga.CorrMulti.",Yvar,".",section,Rdata.suffix))
    
  }
  # print(p)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

cors <- function(df, methodSel=NULL) {
  if(is.null(methodSel)){
    methodSel <- "spearman"
  }
  # turn all three matrices (r, n, and P into a data frame)
  df.c <-na.omit(as.matrix(df))
  M <- Hmisc::rcorr(df.c,type = methodSel)
  # M <- cor(as.matrix(df))
  # return the three data frames in a list return(Mdf)
  Mdf <- map(M, ~data.frame(.x))
}

formatted_cors <- function(df, methodSel=NULL){
  # df <- DF.f
  # cors(df)  %>% map(~get_upper_tri(.x)) %>% #
  #   map(~rownames_to_column(.x, var="measure1")) %>%
  #   map(~pivot_longer(.x, -measure1, "measure2")) %>% map(~drop_na(.x)) %>% #
  #   bind_rows(.id = "id") %>% 
  #   pivot_wider(names_from = id, values_from = value) %>%
  #   mutate(sig_p = ifelse(P < .05, T, F), p_if_sig = ifelse(P <.05, P, NA), r_if_sig = ifelse(P <.05, r, NA)) 
  
  # res <- cors(df)  %>% map(~get_upper_tri(.x)) %>% #
  #   map(~rownames_to_column(.x, var="measure1")) %>%
  #   map(~pivot_longer(.x, -measure1, "measure2"))
  # 
  if(is.null(methodSel)){
    methodSel <- "spearman"
  }
  res <- cors(df, method=methodSel)  %>% map(~get_upper_tri(.x)) %>% #
    map(~rownames_to_column(.x, var="measure1")) %>%
    map(~melt(.x,na.rm = TRUE))
  
  colnames(res$r)[2]<- "measure2"
  colnames(res$n)[2]<- "measure2"
  colnames(res$P)[2]<- "measure2"
  
  # nas <-which(is.na(res$r$value))
  # res$r <- res$r[-nas,]
  # res$n <- res$n[-nas,]
  # res$P <- res$P[-nas,]
  res %>% bind_rows(.id = "id") %>% 
    pivot_wider(names_from = id, values_from = value) %>%
    dplyr::mutate(sig_p = ifelse(P < .05, T, F), p_if_sig = ifelse(P <.05, P, NA), r_if_sig = ifelse(P <.05, r, NA)) 
}

corHeatmap <- function(DF,sampleCol,selectedCols,labels,methodSel=NULL,justPlot=FALSE, formCors=NULL){
  
  if(!isTRUE(justPlot)){
    # DF <- tcgaData.tp.sub.filt  %>% dplyr::filter(project %in% mycancertypes.reduced)
    # sampleCol <-"bcr_patient_barcode"
    # selectedCols <- featuresIn
    # labels <- featuresInLabels
    # 
    if(is.null(methodSel)){
      methodSel <- "spearman"
    }
    #else{
    #   methodSel <- "pearson"
    # }
    
    # DF.f <- DF  %>% column_to_rownames(var=sampleCol) %>% dplyr::select(all_of(selectedCols))
    # colnames(DF.f) <- labels
    # 
    # 
    # cormat <- round(cor(na.omit(as.matrix(DF.f)),method = methodSel),2)
    
    DF.f <- DF %>% 
      dplyr::select(all_of(selectedCols)) %>% 
      setNames(labels) %>%
      correlate(use = "pairwise.complete.obs", # This is very important
                method = "spearman", diagonal=1) #%>% 
    #stretch()
    
    # Check significance of correlations
    cormat <- DF.f %>% column_to_rownames(var="term") %>% as.matrix()
    
    # head(cormat)
    upper_tri <- get_upper_tri(cormat)
    melted_cormat <- melt(upper_tri, na.rm = TRUE)
    # head(melted_cormat)
    
    
    formCors <- formatted_cors(DF %>% 
                                 dplyr::select(all_of(selectedCols)) %>% 
                                 setNames(labels), methodSel = methodSel)
    formCors$measure2 <- gsub('\\.\\.',"+ ",formCors$measure2 )
    formCors$measure2 <- gsub('\\.'," ",formCors$measure2 )
    formCors$measure2 <- ifelse(formCors$measure2=="PD L1 expression","PD-L1 expression",
                                   ifelse(formCors$measure2=="T cell inflamed GEP","T-cell inflamed GEP",
                                          ifelse(formCors$measure2=="CD8+ T cells infiltration","CD8+ T-cells infiltration",
                                                 ifelse(formCors$measure2=="CD4+ T cells infiltration","CD4+ T-cells infiltration",formCors$measure2))))
    formCors$measure1 <- factor(formCors$measure1, levels = levels(melted_cormat$Var1))
    formCors$measure2 <- factor(formCors$measure2, levels = levels(melted_cormat$Var2))
  }
  
  # head(formCors)
  # 
  # 
  # 
  p<-formCors %>%
    ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
    geom_tile() +
    labs(x = NULL, y = NULL, fill = "Spearman's\nCorrelation",  caption="Only significant Spearman's correlation coefficients shown") + #title=paste0("Correlations in ",selectedGroup),
    scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
    geom_label(size = 6,label.size = NA)+
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0))
  # p1 <- p1 + theme_ipsum(base_family = "Palatino Linotype", base_size = 10)
  p <-p + theme_bw(base_family = "Palatino Linotype", base_size = 12) +
    theme(legend.position = "bottom",
          panel.border = element_rect(colour = "black", size=0.2),
          panel.background = element_rect(fill = "white"),
          #axis.text.x = element_text(angle = 40, hjust = 1,size=6),
          axis.text.x = element_text(size=12,angle = 25, hjust = 1),
          axis.text.y = element_text(size = 12),
          legend.key = element_rect(color = NA),#element_rect(size = 10,margin(2, 4, 2, 2, "cm")),
          legend.key.height=unit(0.8,"line"),
          legend.key.width=unit(3.2,"line"),
          legend.text = element_text(size = 12),
          legend.background = element_rect(colour = NA),legend.spacing.x = unit(1.0, 'cm'),
          legend.direction='horizontal',legend.box='horizontal'
    )
  print(p)
}




corHeatmapSection <- function(DF,groupVar,selectedGroup,selectedCols,sampleColumn, labels){
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
  
  corHeatmap(DF.f,sampleColumn,selectedCols,labels)
}



corPlot <- function(DF, Xvar, Yvar, xLabel,yLabel, xPlus, yPlus, yMinus,Logx= FALSE, filterVar=NULL, filterVarSelected="", 
                    shapeVar=NULL, shapeVarLabel="", facetVar=NULL, plotCols=NULL, addCorLabel=NULL,colorVar=NULL){
  if (!is.null(filterVar)){
    filterVar.un <- as.name(filterVar)
    DF <- DF %>% dplyr::filter(!!filterVar.un %in% filterVarSelected) %>% droplevels()
    DF[[filterVar]] <- factor(DF[[filterVar]],levels=filterVarSelected) 
    DF
  }
  
  ######################################
  Xvar.un <- as.name(Xvar)
  Yvar.un <- as.name(Yvar)
  if(!is.null(facetVar)){
    facetVar.un <- as.name(facetVar)
  }
  ## PLOT
  
  if(is.null(shapeVar)){
    if(!is.null(facetVar)){
      DF <- DF %>% dplyr::select(!!Xvar.un,!!Yvar.un,!!facetVar.un) %>% dplyr::filter(complete.cases(.))
      
    }else{
      DF <- DF %>% dplyr::select(!!Xvar.un,!!Yvar.un) %>% dplyr::filter(complete.cases(.))
      
    }
    
    # if(is.null(colorVar)){
    #   p1 <- ggplot(DF,aes_string(Xvar,Yvar))
    # }else{
    #   p1 <- ggplot(DF,aes_string(Xvar,Yvar, color=colorVar))
    # }
    
    p1 <- ggplot(DF,aes_string(Xvar,Yvar))
  }else{
    shapeVar.un <- as.name(shapeVar)
    if(!is.null(facetVar)){
      DF <- DF %>% dplyr::select(!!Xvar.un,!!Yvar.un,!!facetVar.un) %>% dplyr::filter(complete.cases(.))
      
      
    }else{
      DF <- DF %>% dplyr::select(!!Xvar.un,!!Yvar.un) %>% dplyr::filter(complete.cases(.))
      
      
    }
    
    # if(is.null(colorVar)){
    #   p1 <- ggplot(DF,aes_string(Xvar,Yvar, shape=shapeVar))
    # }else{
    #   p1 <- ggplot(DF,aes_string(Xvar,Yvar, shape=shapeVar, color=colorVar))
    # }
    p1 <- ggplot(DF,aes_string(Xvar,Yvar, shape=shapeVar))
  }
  ## PLOT
  
  # if(is.null(shapeVar)){
  #   if(is.null(colorVar)){
  #     p1 <- ggplot(DF,aes_string(Xvar,Yvar))
  #   }else{
  #     p1 <- ggplot(DF,aes_string(Xvar,Yvar, color=colorVar))
  #   }
  #   
  #   
  # }else{
  #   if(is.null(colorVar)){
  #     p1 <- ggplot(DF,aes_string(Xvar,Yvar, shape=shapeVar))
  #   }else{
  #     p1 <- ggplot(DF,aes_string(Xvar,Yvar, shape=shapeVar, color=colorVar))
  #   }
  #   
  # }
  # 
  if (Logx==TRUE){
    p1 <- p1 + scale_x_continuous(trans='log10')
  }
  
  p1 <- p1 + geom_point(size = 2, show.legend = TRUE, alpha = .4)#,color="#69b3a2"
  p1 <- p1 + xlab(paste0(xLabel))+ylab(yLabel)
  p1 <- p1 + geom_smooth(method=lm, se=FALSE, color="black",linetype = "dashed",size=0.5)
  
  p1 <- p1 + theme_ipsum(base_family = "Palatino Linotype", base_size = 10)
  
  if(is.null(shapeVar)){
    p1 <- p1
  }else{
    p1 + scale_color_discrete(name=shapeVarLabel)
  }
  
  if(!is.null(facetVar)){
    p1 <- p1 + theme(axis.text.x = element_text(hjust = 1, vjust = 1,size = 10),
                     axis.title.x = element_text(size=16,vjust = 0.5, hjust=0.5),
                     axis.text.y = element_text( hjust = 1, vjust = 1,size = 10), 
                     axis.title.y = element_text(size=16,angle=90, vjust = 0.5,hjust=0.5),
                     legend.title = element_text(size=13,color="black"),
                     legend.text = element_text(size = 10),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"))
    p1 <- p1 + facet_wrap(facetVar,ncol=plotCols,scales='free')
  }
  
  if (is.null(addCorLabel)){
    p1 <- p1 + theme(axis.text.x = element_text(hjust = 1, vjust = 1,size = 10),
                     axis.title.x = element_text(size=12,),
                     axis.text.y = element_text( hjust = 1, vjust = 1,size = 10), 
                     axis.title.y = element_text(size=12,angle=90),
                     legend.title = element_text(size=12,color="black"),
                     legend.text = element_text(size = 10),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"))
    # Pearson
    corTest.pearson <- cor.test(DF[[Xvar]],DF[[Yvar]],method = "pearson")
    
    l.pearson <- corTest.pearson$p.value
    l.pearson <- format(round(l.pearson,4), scientific = TRUE)
    l.pearson <-strsplit(l.pearson,"e")
    # quote the part before the exponent to keep all the digits
    l.pearson1 <- l.pearson[[1]][1]
    l.pearson2 <- gsub("-0","-",l.pearson[[1]][2])
    
    annot.pearson <- ggplot2::annotate("text", x = min(DF[[Xvar]]) + xPlus, y = max(DF[[Yvar]])+yPlus, label = bquote(paste("Pearson R" == .(round( corTest.pearson$estimate,2))," (P" < .(l.pearson1) %*% 10^ .(l.pearson2),")")),fontface ="plain",size = 4,family= "Palatino Linotype")
    
    # Spearman
    corTest.spearman <- cor.test(DF[[Xvar]],DF[[Yvar]],method = "spearman")
    
    l.spearman <- corTest.spearman$p.value
    l.spearman <- format(round(l.spearman,4), scientific = TRUE)
    l.spearman <-strsplit(l.spearman,"e")
    # quote the part before the exponent to keep all the digits
    l.spearman1 <- l.spearman[[1]][1]
    l.spearman2 <- gsub("-0","-",l.spearman[[1]][2])
    
    annot.spearman <- ggplot2::annotate("text", x = min(DF[[Xvar]]) + xPlus, y = max(DF[[Yvar]])+yPlus-yMinus, label = bquote(paste("Spearman R" == .(round( corTest.spearman$estimate,2))," (P" < .(l.spearman1) %*% 10^ .(l.spearman2),")")),fontface ="plain",size = 4,family= "Palatino Linotype")
    
    p1 <- p1 + annot.pearson
    p1 <- p1 + annot.spearman
  }
  p1 <- p1 + scale_y_continuous(breaks = scales::pretty_breaks(n = 3))
  p1 <- p1 + scale_x_continuous(breaks = scales::pretty_breaks(n = 3))
  print(p1)
}

corr_calculation <- function(DF, Xvar,Yvar, corMethod,filterVar, filterVarSelected){
  if (!is.null(filterVar)){
    filterVar.un <- as.name(filterVar)
    DF <- DF %>% dplyr::filter(!!filterVar.un %in% filterVarSelected) %>% droplevels()
    DF
  }
  
  ## Calculate correlations
  if(nrow(DF)>10){
    corTest <- cor.test(DF[[Xvar]],DF[[Yvar]],method = corMethod)
    corTest_se <- round(unname(sqrt((1 - corTest$estimate^2)/corTest$parameter)),3)
    data.label <- data.frame(label=paste("r = ", round(corTest$estimate,3),"\u00B1",corTest_se), x=0.1*max(DF[[Xvar]],na.rm=T), y=max(DF[[Yvar]], na.rm=T)*.8)
    
    cor <- round(corTest$estimate,3)
    cor
  }else{
    cor <- NULL
  }
  
  if(is.null(cor)==TRUE){
    cor <- NA
    cor
  }else{
    cor <- round(corTest$estimate,3)
  }
  cor
}

cor_heatmap <- function(corDF,DF,xVar,yVar,fillVar,facetVar,SampleCol,xLabel,absoluteCor=FALSE,factorLevels=NULL){
  # corDF <-cor_all.tis.tutack
  # DF <-orr.tutack.tmb.tcga.all
  # xVar <-"correlation"
  # yVar <-"project"
  # fillVar <-"correlation"
  # facetVar <-"project"
  # SampleCol <-"bcr_patient_barcode"
  # xLabel <-"Spearman Correlation"
  
  facetVar.un <- as.name(facetVar)
  SampleCol.un <- as.name(SampleCol)
  xVar.un <- as.name(xVar)
  # Find number of datapoints for each cancer type
  tissue_dataPoints <- DF %>% dplyr::select(!!facetVar.un,!!SampleCol.un) %>% group_by(!!facetVar.un) %>% distinct() %>% dplyr::summarize(n())
  colnames(tissue_dataPoints)[2] <- "no_rows" 
  tissues_fewData <- tissue_dataPoints %>% dplyr::filter(no_rows<15) %>% dplyr::select(!!facetVar.un) %>% dplyr::pull(!!facetVar.un)
  
  
  summary_cor <- as.data.frame(corDF)
  colnames(summary_cor) <- "correlation"
  
  summary_cor_T <- summary_cor %>% rownames_to_column("project") %>% as_tibble()  
  colnames(summary_cor_T) <- c("project","correlation")
  
  if(!is.null(factorLevels)){
    corr_data <-summary_cor_T %>% as.data.frame() %>%  as_tibble() %>% mutate(project=factor(project,levels = factorLevels))
    
  }else{
    corr_data <-summary_cor_T %>% as.data.frame() %>%  as_tibble() %>% mutate(project=factor(project))
    
  }
  
  if (nrow(corr_data)==0){
    
    maxCor_ind <-which(abs(corr_data[[fillVar]]) == max(abs(corr_data[[fillVar]])))
    maxCor_tissue <- corr_data[maxCor_ind ,1]
    maxCor_tissue <- as.character(maxCor_tissue[[yVar]])
    # # Color for tissue with max correlation
    # a <- rep("black",length(levels(corr_data[[yVar]])))
    # b <- rep("plain",length(levels(corr_data[[yVar]])))
  }else{
    maxCor_ind <-which(abs(corr_data[[fillVar]]) == max(abs(corr_data[[fillVar]])))
    maxCor_tissue <- corr_data[maxCor_ind ,1]
    maxCor_tissue <- as.character(maxCor_tissue[[yVar]])
    # # Color for tissue with max correlation
    # a <- ifelse(levels(corr_data[[yVar]]) == maxCor_tissue,"Light Slate Blue",ifelse(levels(corr_data[[yVar]]) %in% tissues_fewData, "Chocolate 3","black"))
    # b <- ifelse(levels(corr_data[[yVar]]) == maxCor_tissue,"bold","plain")
  }
  
  ggtheme = ggplot2::theme_minimal(base_family = "Palatino Linotype", base_size = 10)
  corr_data.m <- melt(corr_data)
  corr_data.m[[yVar]] <-factor(corr_data.m[[yVar]], 
                               levels = rev(levels(factor(corr_data.m[[yVar]]))))
  
  if(isTRUE(absoluteCor)){
    corr_data.m[["value"]]<- abs(corr_data.m[["value"]])
  }
  g <- ggplot(corr_data.m,aes_string("variable",yVar))
  g <- g +  geom_tile(aes_string(fill = "value"), color = "black", size = 0.1)
  g <- g + geom_text(data=corr_data.m,
                     aes_string(x = "variable",
                                y = yVar,
                                label = round(corr_data.m[["value"]], digits = 2)),
                     family = "Palatino Linotype",
                     size = 3.5)
  if(isTRUE(absoluteCor)){
    g <- g + scale_fill_gradient2(low = 'steel blue', mid = 'white', high = 'Indian red', midpoint = 0.2,limits=c(0,1))
  }else{
    g <- g + scale_fill_gradient2(low = 'steel blue', mid = 'white', high = 'Indian red', midpoint = 0,limits=c(-1,1))
  }
  
  g <- g + xlab(xLabel)
  g <- g + ylab("")
  #g <-g + scale_x_discrete(limits = rev(levels(the_factor)))
  g <- g + theme_minimal(base_family = "Palatino Linotype", base_size = 10)
  g <- g + theme(axis.text.x = element_blank(),#angle = 45, hjust = 1, vjust = 1,
                 axis.text.y = element_text(size = 10),#,color=a, face = b
                 axis.title = element_text(size = 12),
                 #legend.position = c(0.33,-0.4),
                 legend.position = "bottom",
                 legend.direction = "horizontal",
                 #legend.justification = c(0,0),
                 plot.margin = margin(t=10,b=20,l=10,r=10, unit="pt"),
                 legend.title = element_text(size=8,color="black"))
  
  
  g <- g + ggtitle("")
  g <- g + labs(fill = "")
  print(g)    
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
    geom_text(size=6.5, family="Palatino Linotype") +
    # geom_label(size=3.5) +
    scale_fill_gradient2(mid="#0C6291",low="#FBFEF9",high="#A63446", limits=c(0,1)) + #
    theme_bw(base_family = "Palatino Linotype", base_size = 14) +
    theme(legend.position = "bottom",
          panel.border = element_rect(colour = "black", size=0.2),
          panel.background = element_rect(fill = "white"),
          #axis.text.x = element_text(angle = 40, hjust = 1,size=6),
          axis.text.x = element_text(size=15,angle = 25, hjust = 1),
          axis.text.y = element_text(size = 15),
          legend.key = element_rect(size = 8),
          # legend.key.height=unit(0.6,"line"),
          # legend.key.width=unit(1,"line"),
          legend.text = element_text(size = 12),
          legend.background = element_rect(colour = NA),legend.spacing.x = unit(1.0, 'cm'),
          legend.direction='horizontal',legend.box='horizontal'
    ) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) + guides(fill = guide_colorbar(order = 1,barwidth = 15,title.vjust = 1,title.hjust = 1))
  print(p)
}
###############################################################
## Correlation plot, ORR vsANything, ALL avail ORR TCGA data ##
###############################################################
orr_corrPlot_tcga <- function(tcgaDF,xVariable, aggregateBy,aggregateMethod,aggregateMethodChar,mergeByCol, selectedORRCancers, orrDF, xLabel, xPlus,yPlus, yMinus,CorrMethod,Logx=FALSE){

  # tcgaDF <- tcgaData.tp.sub.filt.orr
  # xVariable <- "TCR_Richness"
  # aggregateBy <- "project"
  # aggregateMethod <- median
  # aggregateMethodChar <-"median"
  # mergeByCol <- "project"
  # selectedORRCancers <-orr_selected_cancers
  # orrDF <-orrDf
  # xLabel <- "TCR Richness"
  # xPlus <-0.4
  # yPlus <-0.08
  # yMinus <-0.015
  # CorrMethod <-"spearman"
  # Logx=FALSE
  # Add mean of TuTACK signature for each cancer
  tcga_orr_aggregDF <-aggregate(tcgaDF[[xVariable]]~tcgaDF[[aggregateBy]], FUN=aggregateMethod)
  colnames(tcga_orr_aggregDF) <- c(mergeByCol,"variable")
  tcga_orr_aggregDF.ord <- tcga_orr_aggregDF[order(tcga_orr_aggregDF$variable),]
  # subset to selected TCGA datasets
  tcga_orr_aggregDF_final <- subset(tcga_orr_aggregDF.ord,tcga_orr_aggregDF.ord[[mergeByCol]] %in% selectedORRCancers)
  colnames(tcga_orr_aggregDF_final) <- c(mergeByCol, xVariable)
  # Merge variable with ORR
  orr_var_df <- merge(orrDF, tcga_orr_aggregDF_final, by=mergeByCol)
  ## PLOT
  p1 <- ggplot(orr_var_df,aes_string(xVariable,"ORR", label=mergeByCol))
  if (Logx==TRUE){
    p1 <- p1 + scale_x_continuous(trans='log10')
  }
  p1 <- p1 + geom_point(aes_string(size = "PD1L1_totalTreated"),shape = 16, show.legend = TRUE, alpha = .4)#,color="#69b3a2"
  p1 <- p1 + geom_text_repel(size = 3, 
                             segment.size=0.2, 
                             nudge_x = 0.015,
                             point.padding = unit(8.1, "points"))
  p1 <- p1 + xlab(paste0(firstup(aggregateMethodChar)," ",xLabel))+ ylab("ORR")
  # p1 <- p1 + geom_smooth(method=lm, 
  #                        se=FALSE, 
  #                        color="Gray",
  #                        linetype = "dashed",size=0.3)
  p1 <- p1 + theme_ipsum(base_family = "Palatino Linotype", base_size = 10)
  p1 <- p1 + theme(axis.text.x = element_text(hjust = 1, vjust = 1,size = 10),
                   axis.title.x = element_text(size=13),
                   axis.text.y = element_text( hjust = 1, vjust = 1,size = 10), 
                   axis.title.y = element_text(size=13),
                   legend.title = element_text(size=13,color="black"),
                   legend.text = element_text(size = 10),
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(), 
                   axis.line = element_line(colour = "black"),
                   legend.box.background = element_rect(colour = "black"),
                   legend.position="right", legend.box = "vertical")
  p1 <- p1 + scale_color_discrete(name="TCGA project")
  p1 <- p1 + scale_size_area(breaks = c(0,25,50,100,200,400,800))
  p1 <- p1 + labs(size = "No. of \npatients \nevaluated \n(ORR)")
  
  ## Spearman
  # Calculate correlation
  cor2 <- cor.test(orr_var_df[,xVariable],orr_var_df[,"ORR"],method = "spearman")
  cor2_se <- round(unname(sqrt((1 - cor2$estimate^2)/cor2$parameter)),2)
  data.label <- data.frame(label=paste("r = ", round(cor2$estimate,2),"\u00B1",cor2_se), x=0.1*max(orr_var_df[[xVariable]],na.rm=T), y=max(orr_var_df$ORR, na.rm=T)*.8)
  
  l <- cor2$p.value
  l <- format(round(l,4), scientific = TRUE)
  l <-strsplit(l,"e")
  # quote the part before the exponent to keep all the digits
  l1 <- l[[1]][1]
  l2 <- gsub("-0","-",l[[1]][2])
  
  p1<- p1 + ggplot2::annotate("text", x = min(orr_var_df[[xVariable]]) + xPlus, y = max(orr_var_df$ORR)+ yPlus, label = bquote(paste("Spearman R" == .(round(cor2$estimate,2))," (P" < .(l1) %*% 10^ .(l2),")")),fontface ="plain",size = 4,family= "Palatino Linotype")
  
  ## Pearson
  # Calculate correlation
  cor2 <- cor.test(orr_var_df[,xVariable],orr_var_df[,"ORR"],method = "pearson")
  cor2_se <- round(unname(sqrt((1 - cor2$estimate^2)/cor2$parameter)),2)
  data.label <- data.frame(label=paste("r = ", round(cor2$estimate,2),"\u00B1",cor2_se), x=0.1*max(orr_var_df[[xVariable]],na.rm=T), y=max(orr_var_df$ORR, na.rm=T)*.8)
  
  l <- cor2$p.value
  l <- format(round(l,4), scientific = TRUE)
  l <-strsplit(l,"e")
  # quote the part before the exponent to keep all the digits
  l1 <- l[[1]][1]
  l2 <- gsub("-0","-",l[[1]][2])
  
  p1<- p1 + ggplot2::annotate("text", x = min(orr_var_df[[xVariable]]) + xPlus, y = max(orr_var_df$ORR)+ yPlus-yMinus, label = bquote(paste("Pearson R" == .(round(cor2$estimate,2))," (P" < .(l1) %*% 10^ .(l2),")")),fontface ="plain",size = 4,family= "Palatino Linotype")
  
  # print(p1)
  p1
}


violinPlotsGroupsFacets <- function(DF, filterVar=NULL, filterSelectCol, facetVar=NULL, xVar, yVar, xLabel, yLabel,Logy=FALSE,showStat=FALSE){
  # DF <- tcgaData.tp.sub.filt.m
  # filterVar <- "project"
  # filterSelectCol <- mycancertypes.final
  # facetVar <- "feats"
  # xVar<-"project"
  # yVar<-"value"
  # xLabel<-"Cancer type"
  # yLabel <-"Value"
  # Logy=FALSE
  
  if(!is.null(filterVar)){
    filterVar.un <- as.name(filterVar)
    DF <- DF %>% dplyr::filter(!!filterVar.un %in% filterSelectCol) %>% distinct() %>% droplevels()
    DF[[filterVar]] <- factor(DF[[filterVar]], levels = filterSelectCol)
  }
  
  ## Statistics calc
  stat.test <- DF  %>%
    group_by(feats) %>%
    t_test(value ~ project) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance()#tukey_hsd
  stat.test <- stat.test %>%
    add_xy_position(fun = "max", scales = "free")#x = "Condition",group="cells",
  
  ###
  
  p <-ggviolin(DF, x = xVar, y = yVar,#notch=TRUE,add=c("jitter"),
               color = xVar,palette = "uchicago",
               # order = c("PD","CR+PR"),
               ylab = yLabel, xlab = xLabel,add = c("boxplot", "mean_sd"))#,facet.by = #,"jitter"
  
  if(!is.null(facetVar)){
    facetVar.un <- as.name(facetVar)
    p <- p + facet_wrap(facetVar, scales = "free")
    # p <-p + facet_nested(.~feats + project, scales = "free_y", independent = "y")
  }
  
  p + facet_nested(.~feats + project, scales="free")
  
  # p <- p + stat_pvalue_manual(stat.test, label = "p = {p}",step.group.by="feats")#,step.group.by="cells"
  if(isTRUE(showStat)){
    p <-p + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE)
    
  }
  if (Logy==TRUE){
    p <- p + scale_y_continuous(trans='log10')
  }
  
  p <- p + labs(fill = "Cancer type", color= "Cancer type") 
  p <-p + theme_bw(base_family = "Palatino Linotype", base_size = 14)
  p <- p + theme(legend.position = "bottom",
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 axis.text.x = element_text(size=14),
                 axis.title.x= element_text(size=16),
                 #axis.text.x = element_text(size=10),
                 axis.text.y = element_text(size = 16),
                 axis.title.y= element_text(size = 16),
                 # legend.key = element_rect(size = 2),
                 # # legend.key.height=unit(0.6,"line"),
                 # # legend.key.width=unit(1,"line"),
                 legend.text = element_text(size = 14),
                 legend.title = element_text(size = 16),
                 
                 # legend.background = element_rect(colour = "black"),
                 # legend.direction='horizontal',legend.box='horizontal'
                 strip.text.x = element_text(size = 16)
  )
  list(p=p,s=stat.test)  
}


TwosigScores_boxplot_tcga <- function(tcgaDF,selectedCancers,filterCol,xVar,yVar,aggregateMethod,aggregateMethodChar,xLabel,yLabel,ymin1,ymin2,logY=FALSE, 
                                      segm = NULL, ylims=NULL,ticks=NULL, relHeights=NULL,rollSum=FALSE,filterBY=NULL,rollingStep=3){#
  # tcgaDF <-tcgaData.tp.sub.filt
  # selectedCancers <-mycancertypes.reduced
  # filterCol <- "project"
  # xVar <-"B_cells"# The one that defines ORDER
  # yVar<-"TLS_score"
  # aggregateMethod <-median
  # aggregateMethodChar<-"median"
  # xLabel <- "B cell infiltration"
  # yLabel<-"TLS_score"
  # ymin1<- -10
  # ymin2<- -8
  # logY=FALSE
  # segm = NULL
  # ylims=NULL
  # ticks = NULL
  # relHeights = NULL
  # rollSum=TRUE
  # filterBY = NULL#c("BCR_Richness","TCR_Richness","CD8cells","TIS_sigscore","B_cells")
  # Subset to allowed datasets/subsets
  filterCol.un <- as.name(filterCol)
  xVar.un <- as.name(xVar)
  yVar.un <- as.name(yVar)
  if(!is.null(filterBY)){
    tcgaDF <- tcgaDF %>% dplyr::filter(!!filterCol.un %in% selectedCancers) %>% 
      dplyr::select(all_of(c(filterBY,filterCol))) %>% dplyr::filter(complete.cases(.)) %>% droplevels()
  }else{
    tcgaDF <- tcgaDF %>% dplyr::filter(!!filterCol.un %in% selectedCancers) %>% 
      dplyr::select(!!filterCol.un,!!xVar.un,!!yVar.un) %>% dplyr::filter(complete.cases(.)) %>% droplevels()
  }
  
  
  # tcgaDF%>% 
  #   dplyr::group_by(!!filterCol.un) %>% 
  #   dplyr::summarise(n = n()) %>% as.data.frame() %>% print()
  
  tcgaDF[[filterCol]] <- as.factor(tcgaDF[[filterCol]])
  # order factors -project by the mean/median value of signature
  # Add mean/median of signature for each cancer
  # Have expanded immune signature on x-axis
  tcga_sig_aggregDF <- aggregate(tcgaDF[[xVar]]~tcgaDF[[filterCol]], FUN=aggregateMethod)
  colnames(tcga_sig_aggregDF) <- c(filterCol, "sigScore")
  tcga_sig_aggregDF.ord <- tcga_sig_aggregDF[order(tcga_sig_aggregDF$sigScore),]
  #Order the factor levels in the big df - order the project
  tcgaDF[[filterCol]] <- factor(tcgaDF[[filterCol]],tcga_sig_aggregDF.ord[[filterCol]])
  tcgaDF[[yVar]] <- ifelse(tcgaDF[[yVar]]==0, 0.00001,tcgaDF[[yVar]])
  
  # Create separate DF with rolling averages for selected yVar
  tcgaDF.rollavg <- tcgaDF %>% dplyr::group_by(!!filterCol.un)  %>% dplyr::mutate(meanVar = mean(!!yVar.un, na.rm=TRUE))
  tcgaDF.rollavg <- tcgaDF.rollavg %>% dplyr::select(!!filterCol.un, meanVar) %>%dplyr:: group_by(!!filterCol.un) %>% distinct()
  # Order the projects as they were
  tcgaDF.rollavg <- tcgaDF.rollavg[match(levels(tcgaDF[[filterCol]]), tcgaDF.rollavg[[filterCol]]),]
  tcgaDF.rollavg$rollavg <- rollsum(tcgaDF.rollavg$meanVar,rollingStep, fill = TRUE, align = "right")#rollsum(tcgaDF.rollavg$meanVar,3, fill = NA, align = "center")#rollsum(tcgaDF.rollavg$meanVar,3, fill = NA)
  # Plot
  p <- ggplot(tcgaDF, aes_string(x = filterCol, y = yVar))
  p <- p + geom_boxplot(alpha=0.3,outlier.shape = NA, size = 0.1, fill = "black", varwidth = FALSE)
  p <- p + geom_beeswarm(data = tcgaDF,aes_string(x = filterCol, y = yVar), dodge.width=.01,cex=.01,alpha=0.5,size=.4,na.rm = TRUE)#,cex=.1,
  p <- p + geom_point(size = 0.4)
  # rolling average
  # p <- ggplot(tcgaDF.rollavg, aes_string(x = filterCol, y = "rollavg"))
  # p <- p  + geom_line(aes(group=1),linetype="dotted", color="red", size=2)
  # p <- p  + geom_line(data=tcgaDF.rollavg, aes_string(x = filterCol, y = "rollavg"),linetype="dotted", color="red", size=2)
  if(isTRUE(rollSum)){
    p <- p + geom_path(data = tcgaDF.rollavg, aes_string(x = filterCol, y = "rollavg", group = 1),linetype="dotted", color="red",size=1.2)
    
  }
  # p + geom_ma(ma_fun = SMA, n = 30, color = "red")
  p <- p + labs(
    title = paste0(""),
    x = "",
    y = yLabel
  )
  p <- p + theme_ipsum(base_family = "Palatino Linotype", base_size = 10)
  p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size = 12),
                 axis.title.x = element_text(size=16, face="bold"),
                 axis.text.y = element_text( hjust = 1, vjust = 1,size = 12),
                 axis.title.y = element_text(size=12.5, face="bold"),
                 legend.title = element_text(size=13,color="black"),
                 legend.text = element_text(size = 12),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"))
  p <- p + scale_color_discrete(name="TCGA project") +
    coord_cartesian( # This focuses the x-axis on the range of interest
                    clip = 'off')
  
  # if(!is.null(segm)){
  #   p <-mygg.gap(plot=p,
  #                segments=c(250,500), ylim=c(0,1800)) 
  # }
  if(isTRUE(logY)){
    p <- p + scale_y_log10()
  }
  # Add annotation labels on x-axix
  Text1 <- textGrob("low",gp = gpar(fontsize = 12.5,fontface="bold",fontfamily = "Palatino Linotype"))
  Text2 <- textGrob("medium",gp = gpar(fontsize = 12.5,fontface="bold",fontfamily = "Palatino Linotype"))
  Text3 <- textGrob("high",gp = gpar(fontsize = 12.5,fontface="bold",fontfamily = "Palatino Linotype"))
  Text4 <- textGrob(paste0(firstup(aggregateMethodChar)," ",xLabel),gp = gpar(fontsize = 12.5,fontface="bold",fontfamily = "Palatino Linotype"))#family = "serif"
  
  
  # (g1 <- p + theme(plot.margin = unit(c(1,1,2,1), "lines")) + annotation_custom(grob = Text1, xmin = 1, xmax = 1, ymin = ymin2, ymax = ymin2) +
  #     annotation_custom(grob = Text2, xmin = 16, xmax = 16, ymin = ymin2, ymax = ymin2) +
  #     annotation_custom(grob = Text3, xmin = 30, xmax = 30, ymin = ymin2, ymax = ymin2)+
  #     annotation_custom(grob = Text4, xmin = 16, xmax = 16, ymin = ymin1, ymax = ymin1))
  g1 <- p
  if(!is.null(segm)){
    if(is.null(ticks)){
      g2 <-mygg.gap(plot=g1,
                    segments=segm, ylim=ylims)
    }else{
      g2 <-mygg.gap(plot=g1,
                    segments=segm, ylim=ylims,tick_width = ticks,rel_heights = relHeights,margin = c(top = 1, right = 1, bottom = 2, left = 1))
    }

  }else{
    g2 <- g1
  }
  
  # g2 + annotate("text", x = 8, y = 25, label = "Some text")
  # g2 + annotate("text", x = 8, y = 25, label = "Some text")
  # grid.text("MADE IT", x = unit(16, "npc"), y = unit(ymin, "npc"))
  if(!is.null(segm)){
    (g3 <- g2 + theme_ipsum(base_family = "Palatino Linotype", base_size = 12) +  theme(plot.margin = unit(c(1,1,2,1), "lines")) + 
       annotation_custom(grob = Text1, xmin = 0.06, xmax = 0.06, ymin = ymin2, ymax = ymin2) +
       annotation_custom(grob = Text2, xmin = .5, xmax = .5, ymin = ymin2, ymax = ymin2) +
       annotation_custom(grob = Text3, xmin = .95, xmax = .95, ymin = ymin2, ymax = ymin2)+
       annotation_custom(grob = Text4, xmin = .5, xmax = .5, ymin = ymin1, ymax = ymin1))
  }else{
    
    (g3 <- g2 + theme(plot.margin = unit(c(1,1,2,1), "lines")) + annotation_custom(grob = Text1, xmin = .5, xmax = .5, ymin = ymin2, ymax = ymin2) +
        annotation_custom(grob = Text2, xmin = 16, xmax = 16, ymin = ymin2, ymax = ymin2) +
        annotation_custom(grob = Text3, xmin = 30, xmax = 30, ymin = ymin2, ymax = ymin2)+
        annotation_custom(grob = Text4, xmin = 16, xmax = 16, ymin = ymin1, ymax = ymin1))
  }
  
  

  # ymin1<- -0.07
  # ymin2<-4.0
  # p + theme(plot.margin = unit(c(1,1,2,1), "lines")) + annotation_custom(grob = Text1, xmin = 1, xmax = 1, ymin = ymin2, ymax = ymin2) +
  #   annotation_custom(grob = Text2, xmin = 16, xmax = 16, ymin = ymin2, ymax = ymin2) +
  #   annotation_custom(grob = Text3, xmin = 30, xmax = 30, ymin = ymin2, ymax = ymin2)+
  #   annotation_custom(grob = Text4, xmin = 16, xmax = 16, ymin = ymin1, ymax = ymin1)
  
  # gg_table <- ggplot_gtable(ggplot_build(g2))
  # gg_table$layout$clip[gg_table$layout$name=="panel"] <- "off"
  # #grid.draw(gg_table)
  # gg_table
  g3
  
  
  
}

TwosigScores_boxplot_tcga_reversed <- function(tcgaDF,selectedCancers,filterCol,xVar,yVar,aggregateMethod,aggregateMethodChar,xLabel,yLabel,ymin1,ymin2,logY=FALSE, 
                                      segm = NULL, ylims=NULL,ticks=NULL, relHeights=NULL,rollSum=FALSE,filterBY=NULL,rollingStep=3){#
  # tcgaDF <-tcgaData.tp.sub.filt
  # selectedCancers <-mycancertypes.reduced
  # filterCol <- "project"
  # xVar <-"CD8cells"# The one that defines ORDER
  # yVar<-"TCR_Richness"
  # aggregateMethod <-median
  # aggregateMethodChar<-"median"
  # xLabel <- "CD8+ T-cells infiltration"
  # yLabel<-"TCR Richness"
  # ymin1<- -10
  # ymin2<- -8
  # logY=FALSE
  # segm = c(250,500)
  # ylims=c(0,1800)
  # ticks = c(50,300)
  # relHeights = c(1,0,.4)
  # rollSum=TRUE
  # filterBY = NULL
  # rollingStep=3#c("BCR_Richness","TCR_Richness","CD8cells","TIS_sigscore","B_cells")
  # Subset to allowed datasets/subsets
  filterCol.un <- as.name(filterCol)
  xVar.un <- as.name(xVar)
  yVar.un <- as.name(yVar)
  if(!is.null(filterBY)){
    tcgaDF <- tcgaDF %>% dplyr::filter(!!filterCol.un %in% selectedCancers) %>% 
      dplyr::select(all_of(c(filterBY,filterCol))) %>% dplyr::filter(complete.cases(.)) %>% droplevels()
  }else{
    tcgaDF <- tcgaDF %>% dplyr::filter(!!filterCol.un %in% selectedCancers) %>% 
      dplyr::select(!!filterCol.un,!!xVar.un,!!yVar.un) %>% dplyr::filter(complete.cases(.)) %>% droplevels()
  }
  
  
  tcgaDF%>% 
    dplyr::group_by(!!filterCol.un) %>% 
    dplyr::summarise(n = n()) %>% as.data.frame() %>% print()
  
  tcgaDF[[filterCol]] <- as.factor(tcgaDF[[filterCol]])
  # order factors -project by the mean/median value of signature
  # Add mean/median of signature for each cancer
  # Have expanded immune signature on x-axis
  tcga_sig_aggregDF <- aggregate(tcgaDF[[xVar]]~tcgaDF[[filterCol]], FUN=aggregateMethod)
  colnames(tcga_sig_aggregDF) <- c(filterCol, "sigScore")
  tcga_sig_aggregDF.ord <- tcga_sig_aggregDF[order(tcga_sig_aggregDF$sigScore),]
  #Order the factor levels in the big df - order the project
  tcgaDF[[filterCol]] <- factor(tcgaDF[[filterCol]],tcga_sig_aggregDF.ord[[filterCol]])
  tcgaDF[[yVar]] <- ifelse(tcgaDF[[yVar]]==0, 0.00001,tcgaDF[[yVar]])
  
  # Create separate DF with rolling averages for selected yVar
  tcgaDF.rollavg <- tcgaDF %>% dplyr::group_by(!!filterCol.un)  %>% dplyr::mutate(meanVar = mean(!!yVar.un, na.rm=TRUE))
  tcgaDF.rollavg <- tcgaDF.rollavg %>% dplyr::select(!!filterCol.un, meanVar) %>%dplyr:: group_by(!!filterCol.un) %>% distinct()
  # Order the projects as they were
  tcgaDF.rollavg <- tcgaDF.rollavg[match(levels(tcgaDF[[filterCol]]), tcgaDF.rollavg[[filterCol]]),]
  tcgaDF.rollavg$rollavg <- rollsum(tcgaDF.rollavg$meanVar,rollingStep, fill = TRUE, align = "right")#rollsum(tcgaDF.rollavg$meanVar,3, fill = NA, align = "center")#rollsum(tcgaDF.rollavg$meanVar,3, fill = NA)
  # Plot
  p <- ggplot(tcgaDF, aes_string(x = filterCol, y = yVar))
  p <- p + geom_boxplot(alpha=0.3,outlier.shape = NA, size = 0.1, fill = "black", varwidth = FALSE)
  p <- p + geom_beeswarm(data = tcgaDF,aes_string(x = filterCol, y = yVar), dodge.width=.01,cex=.01,alpha=0.5,size=.4,na.rm = TRUE)#,cex=.1,
  p <- p + geom_point(size = 0.4)
  # rolling average
  # p <- ggplot(tcgaDF.rollavg, aes_string(x = filterCol, y = "rollavg"))
  # p <- p  + geom_line(aes(group=1),linetype="dotted", color="red", size=2)
  # p <- p  + geom_line(data=tcgaDF.rollavg, aes_string(x = filterCol, y = "rollavg"),linetype="dotted", color="red", size=2)
  if(isTRUE(rollSum)){
    p <- p + geom_path(data = tcgaDF.rollavg, aes_string(x = filterCol, y = "rollavg", group = 1),linetype="dotted", color="red",size=1.2)
    
  }
  # p + geom_ma(ma_fun = SMA, n = 30, color = "red")
  p <- p + labs(
    title = paste0(""),
    x = "",
    y = yLabel
  )
  p <- p + theme_ipsum(base_family = "Palatino Linotype", base_size = 10)
  p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size = 10),
                 axis.title.x = element_text(size=13, face="bold"),
                 axis.text.y = element_text( hjust = 1, vjust = 1,size = 10),
                 axis.title.y = element_text(size=13, face="bold"),
                 legend.title = element_text(size=13,color="black"),
                 legend.text = element_text(size = 10),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"))
  p <- p + scale_color_discrete(name="TCGA project") +
    coord_cartesian( # This focuses the x-axis on the range of interest
      clip = 'off')
  
  # if(!is.null(segm)){
  #   p <-mygg.gap(plot=p,
  #                segments=c(250,500), ylim=c(0,1800)) 
  # }
  if(isTRUE(logY)){
    p <- p + scale_y_log10()
  }
  # Add annotation labels on x-axix
  Text1 <- textGrob("low",gp = gpar(fontsize = 11,fontface="bold",fontfamily = "Palatino Linotype"))
  Text2 <- textGrob("medium",gp = gpar(fontsize = 11,fontface="bold",fontfamily = "Palatino Linotype"))
  Text3 <- textGrob("high",gp = gpar(fontsize = 11,fontface="bold",fontfamily = "Palatino Linotype"))
  Text4 <- textGrob(paste0(firstup(aggregateMethodChar)," ",xLabel),gp = gpar(fontsize = 11,fontface="bold",fontfamily = "Palatino Linotype"))#family = "serif"
  
  
  # (g1 <- p + theme(plot.margin = unit(c(1,1,2,1), "lines")) + annotation_custom(grob = Text1, xmin = 1, xmax = 1, ymin = ymin2, ymax = ymin2) +
  #     annotation_custom(grob = Text2, xmin = 16, xmax = 16, ymin = ymin2, ymax = ymin2) +
  #     annotation_custom(grob = Text3, xmin = 30, xmax = 30, ymin = ymin2, ymax = ymin2)+
  #     annotation_custom(grob = Text4, xmin = 16, xmax = 16, ymin = ymin1, ymax = ymin1))
  g1 <- p + coord_flip()
  
  
 
  if(!is.null(segm)){
    if(is.null(ticks)){
      g2 <-mygg.gap_X(plot=g1,
                    segments=segm, xlim=ylims)
    }else{
      g2 <-mygg.gap_X(plot=g1,
                    segments=segm, xlim=ylims,tick_width = ticks,rel_heights = relHeights,margin = c(top = 1, right = 1, bottom = 2, left = 1))
    }
    
  }else{
    g2 <- g1
  }
  
  # g2 + annotate("text", x = 8, y = 25, label = "Some text")
  # g2 + annotate("text", x = 8, y = 25, label = "Some text")
  # grid.text("MADE IT", x = unit(16, "npc"), y = unit(ymin, "npc"))
  if(!is.null(segm)){
    (g3 <- g2 + theme_ipsum(base_family = "Palatino Linotype", base_size = 10) +  theme(plot.margin = unit(c(1,1,2,1), "lines")) + annotation_custom(grob = Text1, xmin = 0.06, xmax = 0.06, ymin = ymin2, ymax = ymin2) +
       annotation_custom(grob = Text2, xmin = .5, xmax = .5, ymin = ymin2, ymax = ymin2) +
       annotation_custom(grob = Text3, xmin = .95, xmax = .95, ymin = ymin2, ymax = ymin2)+
       annotation_custom(grob = Text4, xmin = .5, xmax = .5, ymin = ymin1, ymax = ymin1))
  }else{
    
    (g3 <- g2 + theme(plot.margin = unit(c(1,1,2,1), "lines")) + annotation_custom(grob = Text1, xmin = .5, xmax = .5, ymin = ymin2, ymax = ymin2) +
       annotation_custom(grob = Text2, xmin = 16, xmax = 16, ymin = ymin2, ymax = ymin2) +
       annotation_custom(grob = Text3, xmin = 30, xmax = 30, ymin = ymin2, ymax = ymin2)+
       annotation_custom(grob = Text4, xmin = 16, xmax = 16, ymin = ymin1, ymax = ymin1))
  }
  
  
  
  # ymin1<- -0.07
  # ymin2<-4.0
  # p + theme(plot.margin = unit(c(1,1,2,1), "lines")) + annotation_custom(grob = Text1, xmin = 1, xmax = 1, ymin = ymin2, ymax = ymin2) +
  #   annotation_custom(grob = Text2, xmin = 16, xmax = 16, ymin = ymin2, ymax = ymin2) +
  #   annotation_custom(grob = Text3, xmin = 30, xmax = 30, ymin = ymin2, ymax = ymin2)+
  #   annotation_custom(grob = Text4, xmin = 16, xmax = 16, ymin = ymin1, ymax = ymin1)
  
  # gg_table <- ggplot_gtable(ggplot_build(g2))
  # gg_table$layout$clip[gg_table$layout$name=="panel"] <- "off"
  # #grid.draw(gg_table)
  # gg_table
  g3
  
  
  
}



mygg.gap <- function (plot, ylim, segments, tick_width, rel_heights, vjust = 0, 
                      margin = c(top = 1, right = 2, bottom = 1, left = 1), ...) {
  if (!is.list(segments)) {
    segments = list(segments)
  }
  if (all(missing(ylim), is.null(plot$coordinates$limits$y))) {
    stop("ylim is undefined")
  }
  else if (ylim[1] == ylim[2]) {
    stop("ylim should not be the same number")
  }
  else if (missing(ylim)) {
    ylim = plot$coordinates$limits$y
  }
  for (j in 1:length(segments)) {
    seg1 = segments[[j]][1]
    seg2 = segments[[j]][2]
    if (seg1 > seg2) {
      if (ylim[1] < ylim[2]) {
        msg = paste0("No.", j, " segment: c(", 
                     seg1, ",", seg2, ") is wrong. It should be ", 
                     "c(", seg2, ",", seg1, ")")
        stop(msg)
      }
    }
    else if (seg1 < seg2) {
      if (ylim[1] > ylim[2]) {
        msg = paste0("No.", j, " segment: c(", 
                     seg1, ",", seg2, ") is wrong. It should be ", 
                     "c(", seg2, ",", seg1, ")")
        stop(msg)
      }
    }
    else if (seg1 == seg2) {
      msg = paste0("No.", j, " segment: c(", 
                   seg1, ",", seg2, ") is wrong. tick_width should not be equal")
      stop(msg)
    }
  }
  if (length(segments) >= 2) {
    if (ylim[1] < ylim[2]) {
      for (k in 2:length(segments)) {
        pre.2 = segments[[k - 1]][2]
        suf.1 = segments[[k]][1]
        if (pre.2 > suf.1) {
          pre = paste0("c(", segments[[k - 1]][1], 
                       ",", segments[[k - 1]][2], ")")
          suf = paste0("c(", segments[[k]][1], 
                       ",", segments[[k]][2], ")")
          msg = paste0("Segments ", k - 1, " and ", 
                       k, ": ", pre, ",", suf, " are wrong. They should be ", 
                       suf, ",", pre)
          stop(msg)
        }
      }
    }
    else if (ylim[1] > ylim[2]) {
      for (k in 2:length(segments)) {
        pre.2 = segments[[k - 1]][2]
        suf.1 = segments[[k]][1]
        if (pre.2 < suf.1) {
          pre = paste0("c(", segments[[k - 1]][1], 
                       ",", segments[[k - 1]][2], ")")
          suf = paste0("c(", segments[[k]][1], 
                       ",", segments[[k]][2], ")")
          msg = paste0("Segments ", k - 1, " and ", 
                       k, ": ", pre, ",", suf, " are wrong. They should be ", 
                       suf, ",", pre)
          stop(msg)
        }
      }
    }
  }
  if (ylim[1] < ylim[2]) {
    if (min(unlist(segments)) <= ylim[1]) 
      stop("the minimum of segments must be more than the minium of ylim")
    if (max(unlist(segments)) > ylim[2]) 
      stop("the maximum of segments must be lower than maximum of ylim")
  }
  else if (ylim[1] > ylim[2]) {
    if (min(unlist(segments)) <= ylim[2]) 
      stop("the minimum of segments must be more than the minium of ylim")
    if (max(unlist(segments)) > ylim[1]) 
      stop("the maximum of segments must be lower than maximum of ylim")
  }
  if (missing(tick_width)) {
    tick_width = rep(abs(ylim[2] - ylim[1])/10, (length(segments) + 
                                                   1))
  }
  if ((length(tick_width) - length(segments)) < 1) {
    int_len = length(tick_width)
    for (m in (int_len + 1):(length(segments) + 1)) {
      tick_width[m] = tick_width[int_len]
    }
  }
  seg_heights = 0
  y_heights = 1
  if (length(seg_heights) < length(segments)) {
    seg_heights_len = length(seg_heights)
    for (m in (seg_heights_len + 1):length(segments)) {
      seg_heights[m] = seg_heights[seg_heights_len]
    }
  }
  if (length(y_heights) < (length(segments) + 1)) {
    y_heights_len = length(y_heights)
    for (m in (y_heights_len + 1):(length(segments) + 1)) {
      y_heights[m] = y_heights[y_heights_len]
    }
  }
  if (length(plot$scales$scales) == 0) {
    trans = "identity"
  }
  else if ("trans" %in% names(plot$scales$scales[[1]])) {
    trans = plot$scales$scales[[1]]$trans
  }
  else {
    trans = "identity"
  }
  if ("reverse" %in% trans) {
    if (ylim[1] < ylim[2]) {
      msg = paste0("ylim: ", "c(", ylim[1], 
                   ",", ylim[2], ")", " is wrong. It should be ", 
                   "c(", ylim[2], ",", ylim[1], ")")
      stop(msg)
    }
  }
  if ("identity" %in% trans) {
    if (ylim[1] > ylim[2]) {
      msg = paste0("ylim: ", "c(", ylim[1], 
                   ",", ylim[2], ")", " is wrong. It should be ", 
                   "c(", ylim[2], ",", ylim[1], ")")
      stop(msg)
    }
  }
  for (i in 1:length(segments)) {
    gap = unlist(segments[i])
    if (i == 1) {
      if (ylim[1] < ylim[2]) {
        breaks = seq(ylim[1], gap[1], by = tick_width[i])
      }
      else if (ylim[1] > ylim[2]) {
        breaks = seq(gap[1], ylim[1], by = tick_width[i])
      }
      
      p_segment.i <- plot + coord_cartesian(ylim = c(ylim[1], 
                                                     gap[1])) + theme(panel.border = element_blank()) + 
        theme(axis.line.y = element_line(), axis.line.x.bottom = element_line(), 
              plot.title = element_blank(), legend.position = "none") + 
        scale_y_continuous(expand = c(0, 0), trans = trans, 
                           breaks = breaks) + ylab(label = NULL)
      
      
      p_segment = list(p_segment.i)
      names(p_segment)[length(p_segment)] = i
      rel_heigh = c(y_heights[i], seg_heights[i])
    }
    else {
      if (ylim[1] < ylim[2]) {
        breaks = seq(ylim[1], gap[1], by = tick_width[i])
      }
      else if (ylim[1] > ylim[2]) {
        breaks = seq(gap[1], ylim[1], by = tick_width[i])
      }
      p_segment.i <- plot + coord_cartesian(ylim = c(unlist(segments[i - 
                                                                       1])[2], gap[1])) + theme(panel.border = element_blank()) + 
        theme(axis.line.y = element_line(), legend.position = "none", 
              axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
              title = element_blank(), axis.title.x = element_blank()) + 
        scale_y_continuous(expand = c(0, 0), breaks = breaks, 
                           trans = trans) + ylab(label = NULL)  
      p_segment = c(p_segment, list(NULL), list(p_segment.i))
      names(p_segment)[length(p_segment)] = i
      rel_heigh = c(rel_heigh, y_heights[i], seg_heights[i])
    }
    if (i == length(segments)) {
      if (ylim[1] < ylim[2]) {
        breaks = seq(gap[2], ylim[2], by = tick_width[i + 
                                                        1])
      }
      else if (ylim[1] > ylim[2]) {
        breaks = seq(ylim[2], gap[2], by = tick_width[i + 
                                                        1])
      }
      p_segment.i <- plot + coord_cartesian(ylim = c(gap[2], 
                                                     ylim[2])) + theme(panel.border = element_blank()) + 
        theme(axis.line.y = element_line(), axis.line.x.top = element_line(), 
              legend.position = "none", axis.text.x = element_blank(), 
              axis.ticks.x = element_blank(), axis.title.x = element_blank()) + 
        scale_y_continuous(expand = c(0, 0), breaks = breaks, 
                           trans = trans) + ylab(label = NULL) 
      p_segment = c(p_segment, list(NULL), list(p_segment.i))
      names(p_segment)[length(p_segment)] = i + 1
      rel_heigh = c(rel_heigh, y_heights[i])
    }
  }
  
  
  num_parts <- length(p_segment)
  sbt <- p_segment[[1]]$labels$subtitle
  
  p_segment <- purrr::map(p_segment,
                          ~ if(is.ggplot(.)) {.+labs(subtitle=NULL)} else {NULL})
  
  p_segment = rev(p_segment)
  
  p_segment[[1]]$labels$subtitle <- sbt
  
  if (missing(rel_heights)) {
    rel_heights = rev(rel_heigh)
  }
  else {
    rel_heights = rev(rel_heights)
  }
  if (is.null(plot$theme$axis.title.y$angle)) {
    angle = 90
  }
  else {
    angle = plot$theme$axis.title.y$angle
  }
  plot_grid(plotlist = p_segment, ncol = 1, align = "v", 
            rel_heights = rel_heights) + theme(plot.margin = unit(margin, 
                                                                  "cm")) + draw_label(label = plot$labels$y, x = 0, 
                                                                                      hjust = plot$theme$axis.title.y$hjust, vjust = vjust, 
                                                                                      fontfamily = plot$theme$axis.title.y$family, fontface = plot$theme$axis.title.y$face, 
                                                                                      size = plot$theme$axis.title.y$size, angle = angle, lineheight = plot$theme$axis.title.y$lineheight, 
                                                                                      colour = plot$theme$axis.title.y$colour)
}



mygg.gap_X <- function (plot, xlim, segments, tick_width, rel_heights, vjust = 0, 
                        margin = c(top = 1, right = 2, bottom = 1, left = 1), ...) {
  
  # plot=g1
  # segments=segm
  # xlim=ylims
  # tick_width = ticks
  # rel_heights = relHeights
  # margin = c(top = 1, right = 1, bottom = 2, left = 1)
  if (!is.list(segments)) {
    segments = list(segments)
  }
  if (all(missing(xlim), is.null(plot$coordinates$limits$x))) {
    stop("xlim is undefined")
  }
  else if (xlim[1] == xlim[2]) {
    stop("xlim should not be the same number")
  }
  else if (missing(xlim)) {
    xlim = plot$coordinates$limits$x
  }
  for (j in 1:length(segments)) {
    seg1 = segments[[j]][1]
    seg2 = segments[[j]][2]
    if (seg1 > seg2) {
      if (xlim[1] < xlim[2]) {
        msg = paste0("No.", j, " segment: c(", 
                     seg1, ",", seg2, ") is wrong. It should be ", 
                     "c(", seg2, ",", seg1, ")")
        stop(msg)
      }
    }
    else if (seg1 < seg2) {
      if (xlim[1] > xlim[2]) {
        msg = paste0("No.", j, " segment: c(", 
                     seg1, ",", seg2, ") is wrong. It should be ", 
                     "c(", seg2, ",", seg1, ")")
        stop(msg)
      }
    }
    else if (seg1 == seg2) {
      msg = paste0("No.", j, " segment: c(", 
                   seg1, ",", seg2, ") is wrong. tick_width should not be equal")
      stop(msg)
    }
  }
  if (length(segments) >= 2) {
    if (xlim[1] < xlim[2]) {
      for (k in 2:length(segments)) {
        pre.2 = segments[[k - 1]][2]
        suf.1 = segments[[k]][1]
        if (pre.2 > suf.1) {
          pre = paste0("c(", segments[[k - 1]][1], 
                       ",", segments[[k - 1]][2], ")")
          suf = paste0("c(", segments[[k]][1], 
                       ",", segments[[k]][2], ")")
          msg = paste0("Segments ", k - 1, " and ", 
                       k, ": ", pre, ",", suf, " are wrong. They should be ", 
                       suf, ",", pre)
          stop(msg)
        }
      }
    }
    else if (xlim[1] > xlim[2]) {
      for (k in 2:length(segments)) {
        pre.2 = segments[[k - 1]][2]
        suf.1 = segments[[k]][1]
        if (pre.2 < suf.1) {
          pre = paste0("c(", segments[[k - 1]][1], 
                       ",", segments[[k - 1]][2], ")")
          suf = paste0("c(", segments[[k]][1], 
                       ",", segments[[k]][2], ")")
          msg = paste0("Segments ", k - 1, " and ", 
                       k, ": ", pre, ",", suf, " are wrong. They should be ", 
                       suf, ",", pre)
          stop(msg)
        }
      }
    }
  }
  if (xlim[1] < xlim[2]) {
    if (min(unlist(segments)) <= xlim[1]) 
      stop("the minimum of segments must be more than the minium of xlim")
    if (max(unlist(segments)) > xlim[2]) 
      stop("the maximum of segments must be lower than maximum of xlim")
  }
  else if (xlim[1] > xlim[2]) {
    if (min(unlist(segments)) <= xlim[2]) 
      stop("the minimum of segments must be more than the minium of xlim")
    if (max(unlist(segments)) > xlim[1]) 
      stop("the maximum of segments must be lower than maximum of xlim")
  }
  if (missing(tick_width)) {
    tick_width = rep(abs(xlim[2] - xlim[1])/10, (length(segments) + 
                                                   1))
  }
  if ((length(tick_width) - length(segments)) < 1) {
    int_len = length(tick_width)
    for (m in (int_len + 1):(length(segments) + 1)) {
      tick_width[m] = tick_width[int_len]
    }
  }
  seg_heights = 0
  x_heights = 1
  if (length(seg_heights) < length(segments)) {
    seg_heights_len = length(seg_heights)
    for (m in (seg_heights_len + 1):length(segments)) {
      seg_heights[m] = seg_heights[seg_heights_len]
    }
  }
  if (length(x_heights) < (length(segments) + 1)) {
    x_heights_len = length(x_heights)
    for (m in (x_heights_len + 1):(length(segments) + 1)) {
      x_heights[m] = x_heights[x_heights_len]
    }
  }
  if (length(plot$scales$scales) == 0) {
    trans = "identity"
  }
  else if ("trans" %in% names(plot$scales$scales[[1]])) {
    trans = plot$scales$scales[[1]]$trans
  }
  else {
    trans = "identity"
  }
  if ("reverse" %in% trans) {
    if (xlim[1] < xlim[2]) {
      msg = paste0("xlim: ", "c(", xlim[1], 
                   ",", xlim[2], ")", " is wrong. It should be ", 
                   "c(", xlim[2], ",", xlim[1], ")")
      stop(msg)
    }
  }
  if ("identity" %in% trans) {
    if (xlim[1] > xlim[2]) {
      msg = paste0("xlim: ", "c(", xlim[1], 
                   ",", xlim[2], ")", " is wrong. It should be ", 
                   "c(", xlim[2], ",", xlim[1], ")")
      stop(msg)
    }
  }
  for (i in 1:length(segments)) {
    gap = unlist(segments[i])
    if (i == 1) {
      if (xlim[1] < xlim[2]) {
        breaks = seq(xlim[1], gap[1], by = tick_width[i])
      }
      else if (xlim[1] > xlim[2]) {
        breaks = seq(gap[1], xlim[1], by = tick_width[i])
      }
      
      p_segment.i <- plot + coord_cartesian(xlim = c(xlim[1],
                                                     gap[1])) + theme(panel.border = element_blank()) +
        theme(axis.line.x = element_line(), axis.line.y.bottom = element_line(),
              plot.title = element_blank(), legend.position = "none") +
        scale_x_continuous(expand = c(0, 0), trans = trans,
                           breaks = breaks) + xlab(label = NULL)
      
      
      p_segment = list(p_segment.i)
      names(p_segment)[length(p_segment)] = i
      rel_heigh = c(x_heights[i], seg_heights[i])
    }
    else {
      if (xlim[1] < xlim[2]) {
        breaks = seq(xlim[1], gap[1], by = tick_width[i])
      }
      else if (xlim[1] > xlim[2]) {
        breaks = seq(gap[1], xlim[1], by = tick_width[i])
      }
      p_segment.i <- plot + coord_cartesian(xlim = c(unlist(segments[i -
                                                                       1])[2], gap[1])) + theme(panel.border = element_blank()) +
        theme(axis.line.x = element_line(), legend.position = "none",
              axis.text.y = element_blank(), axis.ticks.y = element_blank(),
              title = element_blank(), axis.title.y = element_blank()) +
        scale_x_continuous(expand = c(0, 0), breaks = breaks,
                           trans = trans) + xlab(label = NULL)
      p_segment = c(p_segment, list(NULL), list(p_segment.i))
      names(p_segment)[length(p_segment)] = i
      rel_heigh = c(rel_heigh, y_heights[i], seg_heights[i])
    }
    if (i == length(segments)) {
      if (xlim[1] < xlim[2]) {
        breaks = seq(gap[2], xlim[2], by = tick_width[i +
                                                        1])
      }
      else if (xlim[1] > xlim[2]) {
        breaks = seq(xlim[2], gap[2], by = tick_width[i +
                                                        1])
      }
      p_segment.i <- plot + coord_cartesian(xlim = c(gap[2],
                                                     xlim[2])) + theme(panel.border = element_blank()) +
        theme(axis.line.x = element_line(), axis.line.y.top = element_line(),
              legend.position = "none", axis.text.y = element_blank(),
              axis.ticks.y = element_blank(), axis.title.y = element_blank()) +
        scale_x_continuous(expand = c(0, 0), breaks = breaks,
                           trans = trans) + ylab(label = NULL)
      p_segment = c(p_segment, list(NULL), list(p_segment.i))
      names(p_segment)[length(p_segment)] = i + 1
      rel_heigh = c(rel_heigh, x_heights[i])
    }
  }
  
  
  num_parts <- length(p_segment)
  sbt <- p_segment[[1]]$labels$subtitle
  
  p_segment <- purrr::map(p_segment,
                          ~ if(is.ggplot(.)) {.+labs(subtitle=NULL)} else {NULL})
  
  p_segment = rev(p_segment)
  
  p_segment[[1]]$labels$subtitle <- sbt
  # p_segment
  if (missing(rel_heights)) {
    rel_heights = rev(rel_heigh)
  }
  else {
    rel_heights = rev(rel_heights)
  }
  if (is.null(plot$theme$axis.title.x$angle)) {
    angle = 90
  }
  else {
    angle = plot$theme$axis.title.x$angle
  }
  plot_grid(plotlist = p_segment, ncol = 1, align = "h",
            rel_heights = rel_heights) + theme(plot.margin = unit(margin,
                                                                  "cm")) + draw_label(label = plot$labels$x, y = 0,
                                                                                      hjust = plot$theme$axis.title.x$hjust, vjust = vjust,
                                                                                      fontfamily = plot$theme$axis.title.x$family, fontface = plot$theme$axis.title.x$face,
                                                                                      size = plot$theme$axis.title.x$size, angle = angle, lineheight = plot$theme$axis.title.x$lineheight,
                                                                                      colour = plot$theme$axis.title.x$colour)
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
  resObj <-skcm.res
  
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
  
  
  
  # # Enrich & GSEA KEGG
  # 
  # kegg.enrich <- compareCluster(entrez~group, data=dataDF, fun="enrichKEGG",organism= 'hsa', 
  #                               pvalueCutoff = 0.05)
  # 
  # dotplot(kegg.enrich,showCategory=20)
  # kegg.gsea <- gseKEGG(geneList=gL.ent, 
  #                      organism = "hsa",
  #                      keyType = "kegg",
  #                      nPerm = 10000, 
  #                      minGSSize = 2, 
  #                      maxGSSize = 800, 
  #                      pvalueCutoff = 0.05,
  #                      pAdjustMethod = "none")
  # 
  # dotplot(kegg.gsea,showCategory=20)
  
}


change_legend_breaks <- function(gg, aesthetic, breaks) {
  
  ## Find the scales associated with the specifed aesthetic
  sc <- as.list(gg$scales)$scales
  all_aesthetics <- sapply(sc, function(x) x[["aesthetics"]][1]) 
  idx <- which(aesthetic == all_aesthetics) 
  
  ## Overwrite the breaks of the specifed aesthetic
  gg$scales$scales[[idx]][["breaks"]] <- breaks
  
  return(gg)
} 

# df <- data.frame(x = 1:10, y = 1:10, z1 = 1:10, z2 = 1:10)
# gg <- ggplot(df, aes(x, y, colour = z1, size = z2)) + 
#   geom_point() + 
#   scale_size(name = "Original size title") + 
#   scale_colour_distiller(palette = "Spectral", name = "Original colour title") 
# 
# change_legend_breaks(gg, "colour", breaks = c(2.5, 7.5))
# change_legend_breaks(gg, "size", breaks = c(1, 9))



mybarplot <- function(df,xVar,yVar,xLabel,yLabel,facetVar,bootFreq=0.5){
  # Basic barplot
  p<-ggplot(data=df, aes_string(x=xVar, y=yVar)) +
    geom_bar(stat="identity")
  p <- p + facet_wrap(facetVar,ncol=2,scales='free')
  p <-p + geom_hline(yintercept = bootFreq, colour="indian red", linetype="dashed")
  p <- p + ylab(yLabel)
  p <- p + xlab(xLabel)
  p <- p + ylim(0,1)
  # Horizontal bar plot
  p <-p + coord_flip()
  #p
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


plotHRratio_old <- function(DF, xVar,yVar,xLabel,yLabel,shapeVar, facetVar=NULL,fillVar,filterVar=NULL,selectedFilter="",yLim=NULL,
                        scaleYBreak=NULL,logy=FALSE,allIn=FALSE, sign3=FALSE,errorBarY=FALSE,Yminerr = NULL,Ymaxerr = NULL){
  # breaks_fun <- function(x) {
  #   # x <- DF$HR
  #   if (max(x, na.rm = T) > 10){
  #     c(6, max(x, na.rm = T)-1)
  #    } else {
  #     c(6, 10)
  #   }
  # }
  # 
  # breaks_fun <- function(x) {
  #   if (max(x) > 24) {
  #     seq(0, 120, 24)
  #   } else {
  #     seq(0, 8, 1)
  #   }
  # }
  
  # library(reshape2)
  # library(scales)
  # library(ggbreak) 
  # library(patchwork)
  # DF <- categorical.4.tcga.aPD1
  # yVar <- "HR"
  # xVar <- "factor.name"
  # xLabel<-"Hazard Ratio"
  # yLabel<-""
  # shapeVar <- "analysis"
  # facetVar <- "project"
  # fillVar <- "sign"
  # filterVar <-NULL
  # selectedFilter<-""
  # scaleYBreak <-NULL
  # logy=TRUE
  # allIn=TRUE
  # errorBarY=TRUE
  # Yminerr = "Lower_CI"  
  # Ymaxerr = "Upper_CI"
  # 
  if(!is.null(filterVar)){
    filterVar.un <- as.name(filterVar)
    DF <- DF %>% dplyr::filter(!!filterVar.un ==selectedFilter) #%>% droplevels() %>% distinct()
    DF$title <- selectedFilter
  }else{
    DF$title <- selectedFilter
  }
  
  if(isTRUE(errorBarY)){
    p<-ggplot(DF, aes_string(y=yVar,x=xVar,shape=shapeVar, ymin = Yminerr, ymax = Ymaxerr))#,color=fillVar,group=shapeVar, fill=fillVar
    p<-p + geom_errorbar(width = 0.2)
    
  }else{
    p<-ggplot(DF, aes_string(y=yVar,x=xVar,shape=shapeVar))#,color=fillVar,group=shapeVar, fill=fillVar
    
  }
  p <- p + geom_point(size=2.5,alpha = 0.65,color="black",aes_string(shape=shapeVar, fill=fillVar))#,aes_string(shape=shapeVar)
  # p <- p + scale_y_break(c(6, max(categorical.4.aPD1$HR, na.rm = T)-1))
  if(isTRUE(allIn)){
    p <- p + scale_shape_manual(values = c(21,22,24), labels = analysisLabels)
    
  }else{
    p <- p + scale_shape_manual(values = c(21,24), labels = analysisLabels)
    
  }
  # p <- p + scale_fill_manual(values=c("#FFa319FF","#350E20FF"), labels=fillLabels, drop=FALSE)#, labels=fillLabels
  if(isTRUE(sign3)){
    p <- p + scale_fill_manual(values=c("#C76024","#FFC107","#350E20FF"), labels=fillLabels, drop=FALSE)#, labels=fillLabels
    
  }else{
    p <- p + scale_fill_manual(values=c("#C76024","#350E20FF"), labels=fillLabels, drop=FALSE)#, labels=fillLabels
    
  }
  
  # p <- p + scale_fill_manual(values=c("#C76024","#FFC107","#350E20FF"), labels=fillLabels, drop=FALSE)#, labels=fillLabels
  
  p <- p + guides(fill=guide_legend(override.aes=list(shape=21)))
  if (logy==TRUE){
    p <- p + scale_y_continuous(trans='log10')
  }
  
  if(!is.null(facetVar)){
    p <- p + facet_wrap(as.formula(paste0("~ ",facetVar)), nrow = 1)#,scales = 'free_y'
  }else{
    p <- p + facet_grid(. ~ title) # as.formula(paste0(".~ ",selectedFilter))
  }
  if(!is.null(yLim)){
    p <- p + coord_cartesian(ylim = c(0, yLim))
  }
  # p <- p + facet_wrap(as.formula(paste0("~ ",facetVar)),scales = 'free_y')
  if(!is.null(scaleYBreak)){
    p <- p + scale_y_break(breaks=c(scaleYBreak, max(DF[[yVar]], na.rm = T)-1), scales = 'free',expand=TRUE,ticklabels = NULL)#round(max(DF[[yVar]], na.rm = T),0)#breaks_fun
  }
  # p<- p+ scale_y_continuous(breaks = breaks_fun)
  # p <- p + scale_y_continuous(breaks = round(seq(min(DF[[yVar]], max(DF[[yVar]], by = 0.5),1))))
  p <- p + scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
  p <- p + geom_hline(yintercept=1)
  # if (logy==TRUE){
  #   p <- p + scale_y_continuous(trans='log10')
  # }
  # p + scale_x_continuous(breaks = breaks_fun, limits = c(0, NA))
  # p <- p + scale_y_break(breaks=c(6, max(DF$HR, na.rm = T)-1), scales = 'free')#breaks_fun
  p <- p + ylab(yLabel)
  p <- p + xlab(xLabel)
  p <- p + ggtitle("")
  
  p <- p + theme_bw(base_family = "Palatino Linotype", base_size = 12)
  p <- p + theme(legend.position = "bottom",
                 panel.border = element_rect(colour = "black", size=0.2),
                 # panel.background = element_rect(fill = "white"),
                 # remove the horizontal grid lines
                 panel.grid.major.y = element_line(color = "grey",
                                                   size = 0.1,
                                                   linetype = "dashed") ,
                 panel.grid.minor.y = element_line(color = "grey",
                                                   size = 0.1,
                                                   linetype = "dashed"),
                 # panel.grid.minor.y = element_line( size=.1, color="grey", linetype = "dashed" ) ,
                 
                 # explicitly set the horizontal lines (or they will disappear too)
                 panel.grid.major.x = element_line( size=.1, color="grey", linetype = "dashed" ), 
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
  p
}

plotHRratio <- function(DF, xVar,yVar,xLabel,yLabel,shapeVar, facetVar=NULL,fillVar,filterVar=NULL,selectedFilter="",yLim=NULL,
                        scaleYBreak=NULL,logy=FALSE,allIn=FALSE, sign3=FALSE,
                        errorBarY=FALSE,Yminerr = NULL,Ymaxerr = NULL, pvalCol="p.adj"){
  DF[[xVar]] <- factor(DF[[xVar]], levels=DF[[xVar]] %>% unique())
  
  if(isTRUE(sign3)){
    DF[[fillVar]]<- as.factor(DF[[fillVar]])
    
  }else{
    DF[[fillVar]] <- ifelse(DF[[pvalCol]] <= 0.05, "p <= 0.05","ns")
    DF[[fillVar]]<- as.factor(DF[[fillVar]])
    
    
  }
  if(!is.null(filterVar)){
    filterVar.un <- as.name(filterVar)
    DF <- DF %>% dplyr::filter(!!filterVar.un ==selectedFilter) #%>% droplevels() %>% distinct()
    DF$title <- selectedFilter
  }else{
    DF$title <- selectedFilter
  }
  
  if(isTRUE(errorBarY)){
    p<-ggplot(DF, aes_string(y=yVar,x=xVar,shape=shapeVar, ymin = Yminerr, ymax = Ymaxerr))#,color=fillVar,group=shapeVar, fill=fillVar
    p<-p + geom_errorbar(width = 0.2)
    
  }else{
    p<-ggplot(DF, aes_string(y=yVar,x=xVar,shape=shapeVar))#,color=fillVar,group=shapeVar, fill=fillVar
    
  }
  
  
  
  
  p <- p + geom_point(size=2.5,alpha = 0.65,color="black",aes_string(shape=shapeVar, fill=fillVar))#,aes_string(shape=shapeVar)
  # p <- p + scale_y_break(c(6, max(categorical.4.aPD1$HR, na.rm = T)-1))
  if(isTRUE(allIn)){
    p <- p + scale_shape_manual(values = c(21,22,24,23), labels = analysisLabels)
    
  }else{
    p <- p + scale_shape_manual(values = c(21,24), labels = analysisLabels)
    
  }
  # p <- p + scale_fill_manual(values=c("#FFa319FF","#350E20FF"), labels=fillLabels, drop=FALSE)#, labels=fillLabels
  
  if(isTRUE(sign3)){
    fillLabels <-levels(DF[[fillVar]])# ns, p<0.05, p <0.1
    p <- p + scale_fill_manual(values=c("#350E20FF","#FFC107","#C76024"), labels=fillLabels, drop=FALSE)#, labels=fillLabels
    
  }else{
    DF[[fillVar]] <- ifelse(DF[[pvalCol]] <= 0.05, "p <= 0.05","ns")
    DF[[fillVar]]<- as.factor(DF[[fillVar]])
    fillLabels <- levels(DF[[fillVar]])# ns first level
    p <- p + scale_fill_manual(values=c("#350E20FF","#FFC107"), labels=fillLabels, drop=FALSE)#, labels=fillLabels
    
  }
  
  # p <- p + scale_fill_manual(values=c("#C76024","#FFC107","#350E20FF"), labels=fillLabels, drop=FALSE)#, labels=fillLabels
  
  p <- p + guides(fill=guide_legend(override.aes=list(shape=21)))
  if (logy==TRUE){
    p <- p + scale_y_continuous(trans='log10')
  }
  
  if(!is.null(facetVar)){
    p <- p + facet_wrap(as.formula(paste0("~ ",facetVar)), nrow = 1)#,scales = 'free_y'
  }else{
    p <- p + facet_grid(. ~ title) # as.formula(paste0(".~ ",selectedFilter))
  }
  if(!is.null(yLim)){
    p <- p + coord_cartesian(ylim = c(0, yLim))
  }
  # p <- p + facet_wrap(as.formula(paste0("~ ",facetVar)),scales = 'free_y')
  if(!is.null(scaleYBreak)){
    p <- p + scale_y_break(breaks=c(scaleYBreak, max(DF[[yVar]], na.rm = T)-1), scales = 'free',expand=TRUE,ticklabels = round(max(DF[[yVar]], na.rm = T),0))#breaks_fun
  }
  # p<- p+ scale_y_continuous(breaks = breaks_fun)
  p <- p + scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
  p <- p + geom_hline(yintercept=1)
  # if (logy==TRUE){
  #   p <- p + scale_y_continuous(trans='log10')
  # }
  # p + scale_x_continuous(breaks = breaks_fun, limits = c(0, NA))
  # p <- p + scale_y_break(breaks=c(6, max(DF$HR, na.rm = T)-1), scales = 'free')#breaks_fun
  p <- p + ylab(yLabel)
  p <- p + xlab(xLabel)
  p <- p + ggtitle("")
  p <- p +labs(fill="p.adj")
  p <- p + guides(shape = guide_legend(order = 1,nrow = 1))#,fill = guide_legend(order = 1,nrow = 1)
  

  p <- p + theme_bw(base_family = "Palatino Linotype", base_size = 12)
  p <- p + theme(
    panel.border = element_rect(colour = "black", size=0.2),
    # panel.background = element_rect(fill = "white"),
    # remove the horizontal grid lines
    panel.grid.major.y = element_line(color = "grey",
                                      size = 0.1,
                                      linetype = "dashed") ,
    panel.grid.minor.y = element_line(color = "grey",
                                      size = 0.1,
                                      linetype = "dashed"),
    # panel.grid.minor.y = element_line( size=.1, color="grey", linetype = "dashed" ) ,
    
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.x = element_line( size=.1, color="grey", linetype = "dashed" ), 
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
    # legend.direction='horizontal',
    # legend.box='horizontal',
    legend.justification = 'center', 
    legend.position = 'bottom', 
    # legend.box = 'vertical', 
    legend.box = 'horizontal', 
    legend.box.just = 'center',
    legend.title = element_text(colour="black", size=12, face="plain"))
  p
}


runUniSurvival_Singlevariable <- function(DF,var,varSelection,varSelected,survType,sampleID,varPlotName,printP=FALSE,saveP = FALSE,savePath,onlySign=TRUE,selectFacets=NULL){
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
  
  # DF <- immdata.TcrBcr.full.red
  # var <-"TCR_Richness"
  # varSelection <-mycancertypes.final
  # varSelected <-"tissue"
  # survType <- "OS"
  # sampleID <- "run_accession"
  # varPlotName<-"TCR Richness"
  # onlySign=FALSE
  # selectFacets=mycancertypes.final
  # typeVar <- "cont"
  
  # cat('\n')
  # cat('\n')
  # cat("## ", varPlotName, "{.tabset .tabset-fade .tabset-pills} \n")
  # cat('\n')
  
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
      hrDF.sign <- hrDF.sign[match(selectFacets, hrDF.sign[[varSelected]]),]
    }else{
      hrDF.sign <- res.cat2$hrDF %>% as.data.frame
    }
  }

  # if(isTRUE(saveP)){
  #   res = list(forestCont = res.cont$forestplot,
  #              forestCat = res.cat$forestplot,
  #              km =  res.cat2$kmPlot,
  #              hrDF = res.cat2$hrDF)
  #   # print(paste0(projectRDataPath, paste0("tcga.SURV","_",survType,"_",var),Rdata.suffix))
  #   save(res, file = paste0(projectRDataPath, paste0("tcga.SURV","_",survType,"_",var),Rdata.suffix))
  # }

  if(isTRUE(printP)){
    # cat('\n')
    # cat('\n')
    # cat("### ", "Forest plot - continuous", "\n")
    # cat('\n')

    print(res.cont$forestplot)

    # cat('\n')
    # cat('\n')
    # cat("### ", "Forest plot - categorical", "\n")
    # cat('\n')

    print(res.cat$forestplot)

    # cat('\n')
    # cat('\n')
    # cat("### ", "KM plots", "\n")
    # cat('\n')


    n<-dim(hrDF.sign)[1]/2;
    graphCols <-ceiling(n/floor(sqrt(n)));
    graphRows<-floor(n)

    # p.km<- ggarrange(plotlist = res.cat2$kmPlot[hrDF.sign[[varSelected]]], nrow =  graphRows,ncol=graphCols,common.legend = TRUE, legend="bottom")#

    legend <- get_legend(res.cat2$kmPlot[[1]] + theme(legend.position = "bottom"))
    p.km<- ggarrange(plotlist = res.cat2$kmPlot[hrDF.sign[[varSelected]]], nrow =  graphRows,ncol=graphCols,common.legend = FALSE)#

    print(p.km)
  }else{
    n<-dim(hrDF.sign)[1]/2;
    graphCols <-ceiling(n/floor(sqrt(n)));
    graphRows<-floor(n)
    # print(res.cat2$kmPlot)
    # p.km<- ggarrange(plotlist = res.cat2$kmPlot[hrDF.sign[[varSelected]]], nrow =  graphRows,ncol=graphCols,common.legend = TRUE, legend="bottom")#
    legend <- get_legend(res.cat2$kmPlot[[1]] + theme(legend.position = "bottom"))
    p.km<- ggarrange(plotlist = res.cat2$kmPlot[hrDF.sign[[varSelected]]], nrow =  graphRows,ncol=graphCols,common.legend = FALSE)#

  }
  res = list(forestCont = res.cont$forestplot,
             forestCat = res.cat$forestplot,
             km =  p.km,
             kmleg = legend,
             hrDF = res.cat2$hrDF)
  # res.cat2$hrDF
  res
}

# ## WRAPPER WITH ALL SURV ANALYSIS ##
performSurvivalCalcs <- function(datasetsEvent, dataDF,mainVar, filterVar, survivalVar, sampleVar,covName,covType="cont",
                                 stratifyMethod=NULL,sLevelA = "high",sLevelB = "low",printPlots= FALSE,onlyForest=FALSE, onlyKM = FALSE){
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


getCoxPHres <- function(coxObj,covar,variable,selectedcancer,varType="quant",varLabel=NULL){
  ######################
  ## COX HAZARD RATIO ##
  ######################
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
    distinct(!!SampleCol,.keep_all = TRUE) %>% #
    dplyr::mutate_at(vars(!!endopoint.un,!!event.un),as.numeric) # %>% dplyr::filter(complete.cases(!!event.un,!!endopoint.un))
  #filter(complete.cases(!!endopoint.un,!!event.un,!!predictor.un))
  
  
  if(isTRUE(dichot)){
    # DICHOTOMIZE-- HERE I HAVE TO CHECK IF TISSUE OR DATASET
    if(stratMethod=="median"){
      # Is it a bi-modal distribution?How to select the cutoff?
      cutglobal <- median(DF[[predictor]][!is.na(DF[[predictor]])], na.rm = TRUE)
      
      # Calculate median for cancer
      # cutlocal <- median(DF.sub[[predictor]][!is.na(DF.sub[[predictor]])], na.rm = TRUE)
      cutlocal <- DF.sub %>% dplyr::select(!!predictor.un,!!event.un) %>% dplyr::filter(complete.cases(.)) %>% pull(!!predictor.un) %>% median()
      print(paste0(">>Final total ", nrow(DF.sub)," patients included<<", "\n"))
      print(paste0("LOCAL MEDIAN - ",selectedcancer, " - ",predictor,": ",cutlocal))
    }else if(stratMethod=="mean"){
      cutglobal <- mean(DF[[predictor]][!is.na(DF[[predictor]])], na.rm = TRUE)
      cutlocal <- mean(DF.sub[[predictor]][!is.na(DF.sub[[predictor]])], na.rm = TRUE)
    }else if(stratMethod=="cut"){
      cutglobal <-as.numeric(optCutpointContinuous(DF,endopoint,event,predictor,filterVar=NULL, filterVarSelplotected=NULL,gatherTable=TRUE)[2])
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
      DF.sub2 <- DF.sub %>% dplyr::mutate(!!stratVariable.un := ifelse(!!predictor.un > cutVar,StratLevel1,StratLevel2)) %>% dplyr::mutate(!!stratVariable:=factor(!!stratVariable.un)) %>% dplyr::filter(complete.cases(!!event.un,!!endopoint.un))
      
    }else{
      DF.sub2 <- DF.sub %>% dplyr::mutate(!!stratVariable.un := ifelse(!!predictor.un >= cutVar,StratLevel1,StratLevel2)) %>% dplyr::mutate(!!stratVariable:=factor(!!stratVariable.un)) %>% dplyr::filter(complete.cases(!!event.un,!!endopoint.un))
      
    }
    # DF.sub2[[stratVariable]] <- factor(DF.sub2[[stratVariable]],levels = c())
  }else{
    stratVariable <- predictor
    DF.sub2 <-DF.sub
  }
  
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
  #############################
  ## HazardRatio forest plot ##
  #############################
  
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




fit_all <- function(selectedcancer,tcgaDF,survHRdf,cutPointsTable,covName,
                    stratMethod,signature,stratVariable,SampleCol,endopoint, event,
                    facetVariable,predictor,StratLevel1,StratLevel2,
                    globalLocal,covType="cont",getLegend=FALSE,HRAnnot=FALSE){
  #########################
  ## Kaplan Meier curves ##
  #########################
  
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
  tcgaDF_sign <- dplyr::filter(tcgaDF, !!facetVariable.un %in% selectedcancer) %>% dplyr::filter(complete.cases(!!event.un,!!endopoint.un)) %>% droplevels()
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
  
  # tcgaDF_sign2 <- tcgaDF_sign %>% mutate(!!stratVariable.un := ifelse(!!predictor.un >= cutVar,StratLevel1,StratLevel2)) %>%mutate(!!stratVariable:=factor(!!stratVariable.un))
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
    dplyr::mutate_at(vars(!!endopoint.un,!!event.un),as.numeric)
  
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
      annotateText = c(paste0("HR = ",round(hr,2),"(",round(ci[1],3),"-",round(ci[2],3),")\n","logrankP = ", logrankp)),
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

survComboKM <- function(l1,l2,var){
  # cat('\n')
  # cat('\n')
  # cat("## ", var, "{.tabset .tabset-fade .tabset-pills} \n")
  # cat('\n')
  # 
  
  # l1<-tcga.SURV
  # l2<-aPD1.SURV
  # var<-featuresIn[3:9][1]
  km.tcga <-l1[[var]]
  km.ici <- l2[[var]]
  # names(km.ici) <- mycancertypes.final
  mLeg <- get_legend(km.tcga$SKCM)
  km.tcga<- ggarrange(plotlist = km.tcga,ncol=2, nrow =  2,common.legend = TRUE,legend="bottom")
  km.ici<- ggarrange(plotlist = km.ici,ncol=2, nrow =  2,common.legend = TRUE,legend="bottom")#,legend = "none"
  
  
  
  p<-ggdraw() + 
    draw_plot(km.tcga+ theme(legend.position = "none"),x = 0.0, y = 0, width = .5, height = .95) +
    draw_plot(km.ici+ theme(legend.position = "none"),x = 0.5, y = 0, width = 0.5, height = .95)
  #   # draw_plot(mLeg,x = 0.5, y = 0, width = 0.02, height = .01) +
  #   draw_plot_label(label = c("A: TCGA", "B: ICI datasets"), size = 14,
  #                   x = c(0, .48), y = c(1, 1))
  # print(p)
  p
}



##########################
## PREDICTIONS-MODELS ####
###########################

ROCcalculation_models <- function(modelList, test.data,testX,testY,classes, LevelIndex, responseVar, modelNamesVec, dsids, thres = NULL,rfeE=FALSE,mixedMode = FALSE,trainModelsVec = NULL, rfeModelsVec = NULL){
  
  # modelList = models.final
  # test.data = testDF[[1]]
  # testX =  data.testX
  # testY = testDF[[1]]$response
  # classes = c("CRPR","PD")
  # LevelIndex =1
  # responseVar = "response"
  # modelNamesVec = names(models.final)
  # dsids = 1:length(models.final)
  # rfeE=FALSE
  # mixedMode = TRUE
  # thres = NULL
  # trainModelsVec = c(1:c(length(models.final)-1))
  # rfeModelsVec = c(length(models.final))

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

myBubblePlot <- function(DF, xVar, yVar, sizeVar, colorVar ){
  # DF <- rocs_schemes_merged_DF.sub# %>% dplyr::filter(Model=="TCR Richness") %>% droplevels()
  # xVar <-"Model"
  # yVar <-"scheme"
  # sizeVar <- "Sample_test"
  # colorVar <-"ROC-AUC"
  mycol <- rev(c("#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#ddf1da"))
  xVar.un <- as.name(xVar)
  yVar.un <- as.name(yVar)
  sizeVar.un <- as.name(sizeVar)
  colorVar.un <- as.name(colorVar)
  p <- ggplot(DF, aes(x=!!xVar.un, y=!!yVar.un, size = !!colorVar.un, fill=!!colorVar.un)) +
    viridis::scale_fill_viridis(option = "F",
                                  direction=1,
                                  alpha=0.5,
                                  label = function(x) sprintf("%.2f", x),
                                  limits=c(0.44,.89),
    ) + #breaks=c(0,0.5,0.7,0.8,0.9,1)

    
    scale_fill_gradientn(colours = mycol,breaks=c(0,0.5,0.7,0.8,0.9,1)) +
    
    geom_point(alpha=0.7,shape=21, color="black") +
    geom_text(aes(label = round(!!colorVar.un,3)), color="black", size=4, fontface="bold", family="Palatino Linotype") +
    scale_size(range = c(7,39),limits = c(0.44,.89),breaks=c(0,0.5,0.7,0.8,0.9,1), name="Test sample size (N)", guide = "none") +#c(6,38)
    scale_x_discrete(position = "top") 
  
  # facet_grid(Model ~ ., scales = "free", space = "free") 
  # facet_wrap(Model ~ ., nrow=1) 
  p <- p + theme_bw(base_family = "Palatino Linotype", base_size = 14)
  p <- p + theme(legend.position = "bottom",
                 panel.border = element_blank(),#element_rect(colour = "black", size=0.2),
                 # remove the horizontal grid lines
                 # panel.grid.major.y = element_line(color = "grey",
                 #                                   size = 0.1,
                 #                                   linetype = "dashed") ,
                 # panel.grid.minor.y = element_line(color = "grey",
                 #                                   size = 0.1,
                 #                                   linetype = "dashed"),
                 # 
                 # explicitly set the horizontal lines (or they will disappear too)
                 # panel.grid.major.x = element_line( size=.1, color="grey", linetype = "dashed" ), 
                 # axis.text.x = element_text(angle = 40, hjust = 1,size=10),
                 axis.text.x = element_text(size=13,face="bold"),
                 axis.text.y = element_text(size = 13,face="bold"),#,angle = 90,vjust=0.9
                 axis.ticks.x = element_blank(),
                 axis.ticks.y = element_blank(),
                 # axis.title.x = element_text(size=14, face="bold"),
                 # axis.title.y = element_text(size=14, face="bold"),
                 legend.key = element_rect(size = 2),
                 legend.key.height=unit(0.6,"line"),
                 legend.key.width=unit(2,"line"),
                 legend.text = element_text(size = 12),
                 # legend.background = element_rect(colour = "black"),
                 legend.direction='horizontal',legend.box='horizontal',
                 legend.title = element_text(colour="black", size=14, face="plain"))
  p<-p + labs(x = "", y = "") + guides(fill = guide_colorbar(order = 1,title.vjust = 1,title.hjust = 1))
  p
}


GOsimilarityPlots <- function(goDF, ont="BP", organism="org.Hs.eg.db", simMethod="Rel"){
  simMatrix <- calculateSimMatrix(goDF %>% pull(ID),
                                  orgdb=organism,
                                  ont=ont,
                                  method=simMethod)
  
  scores <- setNames(-log10(goDF %>% pull(qvalue)), goDF %>% pull(ID))
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.Hs.eg.db")
  
  scPlot <-scatterPlot(simMatrix, reducedTerms,size = "size")
  
  scPlot
}


findCommonDiff_GOterms <- function(GO1,GO2, termColumn="Description",GO1name="GO1",GO2name="GO2",
                                   colExclude = c("geneID"),
                                   aggregCols = c("GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "Count"),
                                   keepCols = c(".id", "ID", "Description")){
  # GO1 <- bp_bcr_0.filt.red.DFs$SKCM
  # GO2 <- bp_pdl1.filt.red.DFs$SKCM
  # termColumn="Description"
  # GO1name="BCR"
  # GO2name="PDL1"
  # colExclude = c("")
  # aggregCols = c("GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "Count")
  # keepCols = c(".id", "ID", "Description")
  
  termColumn.un <- as.name(termColumn)
  GO1name.un <- as.name(GO1name)
  GO2name.un <- as.name(GO2name)
  
  go1DF <- as.data.frame(GO1)
  go2DF <- as.data.frame(GO2)
  
  common <-intersect(go1DF[[termColumn]],go2DF[[termColumn]])
  diff1_2 <- setdiff(go1DF[[termColumn]],go2DF[[termColumn]])
  diff2_1 <- setdiff(go2DF[[termColumn]],go1DF[[termColumn]])
  
  DFlist <- list(go1DF %>% dplyr::filter(!!termColumn.un %in% common) %>% droplevels(),
                 go2DF %>% dplyr::filter(!!termColumn.un %in% common) %>% droplevels())
  names(DFlist) <- c(GO1name,GO2name)
  commonDF <-ldply(DFlist) ## NEED to aggregate here, otherwise I have duplicates
  
  if(!nrow(commonDF)==0){
    commonDF2 <-aggregate(as.formula(paste0("cbind(", paste(aggregCols, collapse = ", "),") ~ ",paste(keepCols, collapse = " + "))), 
                          data = commonDF , 
                          FUN = mean, 
                          na.rm = TRUE)
  }else{
    commonDF2 <-commonDF
  }
  
  
  diff1_2DF <-go1DF %>% dplyr::filter(!!termColumn.un %in% diff1_2) %>% droplevels()
  
  diff2_1DF <-go2DF %>% dplyr::filter(!!termColumn.un %in% diff2_1) %>% droplevels()
  
  res = list(common=commonDF,
             diff1_2 = diff1_2DF,
             diff2_1 = diff2_1DF)
  res
}

findCommonDiff_GOSIMterms <- function(GO1,GO2,sim1,sim2, termColumn="parentTerm",GO1name="GO1",GO2name="GO2",
                                      colExclude = c("go", "term",".id"),
                                      aggregCols = c("parentSimScore", "score", "size", "cluster"),
                                      keepCols = c("parentTerm", "parent")){
  # GO1 <- redL_bcr_0.filt$reducedTerms.rel2
  # GO2 <- redL_pdl1.filt$reducedTerms.rel2
  # sim1 <- simL_bcr_0.filt$simMatrix.rel2
  # sim2 <- simL_pdl1.filt$simMatrix.rel2
  # termColumn="parentTerm"
  # GO1name="BCR"
  # GO2name="PDL1"
  # colExclude = c("go", "term",".id")
  # aggregCols = c("parentSimScore", "score", "size", "cluster")
  # keepCols = c("parentTerm", "parent")
  
  
  
  termColumn.un <- as.name(termColumn)
  GO1name.un <- as.name(GO1name)
  GO2name.un <- as.name(GO2name)
  
  go1DF <- as.data.frame(GO1)
  go2DF <- as.data.frame(GO2)
  
  common <-intersect(go1DF[[termColumn]],go2DF[[termColumn]])
  diff1_2 <- setdiff(go1DF[[termColumn]],go2DF[[termColumn]])
  diff2_1 <- setdiff(go2DF[[termColumn]],go1DF[[termColumn]])
  
  DFlist <- list(go1DF %>% dplyr::filter(!!termColumn.un %in% common) %>% droplevels(),
                 go2DF %>% dplyr::filter(!!termColumn.un %in% common) %>% droplevels())
  names(DFlist) <- c(GO1name,GO2name)
  commonDF <-ldply(DFlist) ## NEED to aggregate here, otherwise I have duplicates
  
  if(!nrow(commonDF)==0){
    commonDF2 <-aggregate(as.formula(paste0("cbind(", paste(aggregCols, collapse = ", "),") ~ ",paste(keepCols, collapse = " + "))), 
                          data = commonDF %>% dplyr::select(-all_of(colExclude)), 
                          FUN = mean, 
                          na.rm = TRUE)
  }else{
    commonDF2 <-commonDF
  }
  
  common1 <-go1DF %>% dplyr::filter(!!termColumn.un %in% common) %>% droplevels()
  
  diff1_2DF <-go1DF %>% dplyr::filter(!!termColumn.un %in% diff1_2) %>% droplevels()
  
  diff2_1DF <-go2DF %>% dplyr::filter(!!termColumn.un %in% diff2_1) %>% droplevels()
  
  # Reduce matrices
  sim1.f <- sim1[diff1_2DF$go,diff1_2DF$go]
  sim2.f <- sim2[diff2_1DF$go,diff2_1DF$go]
  simCommon <-sim1[common1$go %>% unique(),common1$go %>% unique()]
  # Gather results
  res = list(common=commonDF2,
             diff1_2 = diff1_2DF,
             diff2_1 = diff2_1DF,
             sim1_2=sim1.f,
             sim2_1=sim2.f,
             sim_common = simCommon)
  res
}


floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)

GO_bubblePlot <- function(DF,sortVar,xVar,yVar,sizeVar,fillVar, aggreg=FALSE,
                          aggregCols=c("parentSimScore", "score", "size", "cluster"),
                          keepCols = c("parentTerm", "parent"),
                          colExclude = c("go", "term"),choosePlot="norm", sizeMin, sizeMax, xMax,seqN, fillMax){
  # DF = bp_pdl1.filt.red.DFs$BLCA
  # sortVar <- "GeneRatio"
  # # factorVar <- "parentTerm"
  # xVar <-"GeneRatio"
  # yVar <-"Description" #  also factor variable
  # sizeVar <- "Count"
  # fillVar <-"p.adjust"
  # aggreg=FALSE
  # aggregCols <- c("parentSimScore", "score", "size", "cluster")
  # keepCols <- c("parentTerm", "parent")
  # colExclude = c("go", "term")
  # choosePlot="norm"
  # sizeMin=10
  # sizeMax=170
  # xMax=ceiling_dec(max(bp_pdl1.filt.red.merged$GeneRatio),2)
  # seqN=6
  # fillMax=max(bp_pdl1.filt.red.merged$p.adjust)
  # maxGR=0.02
  # 
  sortVar.un <- as.name(sortVar)
  xVar.un <- as.name(xVar)
  yVar.un <- as.name(yVar)
  sizeVar.un <- as.name(sizeVar)
  fillVar.un <- as.name(fillVar)
  
  # DF <- bp_bcr_0.filt.red$SKCM %>% as.data.frame() %>% 
  #   dplyr::select(-geneID) %>% 
  #   dplyr::mutate(category=bp_ImmuneRelated.focused.DF$`Parent Description`[match(ID,bp_ImmuneRelated.focused.DF$ID)]) %>%
  #   dplyr::mutate("GeneRatio"=parse_ratio(!!xVar.un))
  
  if(isTRUE(aggreg)){
    DF<-aggregate(as.formula(paste0("cbind(", paste(aggregCols, collapse = ", "),") ~ ",paste(keepCols, collapse = " + "))), 
                  data = DF %>% dplyr::select(-all_of(colExclude)), 
                  FUN = mean, 
                  na.rm = TRUE)
  }
  
  
  if(choosePlot=="sim"){
    p <-DF %>%
      # arrange(desc(!!sortVar.un)) %>%
      arrange(!!sortVar.un) %>%
      mutate("parentTerm" = factor(!!yVar.un, !!yVar.un %>% unique())) %>%
      ggplot(aes(x=!!xVar.un, y=!!yVar.un, size=!!sizeVar.un, fill=!!fillVar.un)) +
      geom_point(alpha=0.7,shape=21, color="black") +
      scale_size(range = c(.1, 20), name="GO similarity score",limits = c(0,24),breaks = c(2,4,6,8,10,12,14,16,20,24)) +
      scale_x_continuous(limits = c(0,24),breaks=seq(0, 24, 2)) +
      # scale_size_continuous(name = "GO similarity score",
      #                       breaks = c(2,4,6,8),
      #                       limits = c(0,9),
      #                       labels = as.character(c(2,4,6,8)),
      #                       range = c(0, 24) ) +
      # scale_fill_viridis(discrete=FALSE, guide=FALSE, option="A") #+
      viridis::scale_fill_viridis(option = "viridis", direction=1, alpha=0.5,breaks = c(2,4,6,8,10,12,14,16),
                                  limits = c(0,17))+
      theme_bw(base_family = "Palatino Linotype", base_size = 15) +
      theme(axis.text.x = element_text(size = 12),
            # axis.text.y = element_text(size = 14),
            axis.text.y = element_text(size = 14, face="bold",
                 colour=),
            axis.title.y=element_text(size = 12),
            axis.title.x=element_text(size = 16),
            text=element_text(family="Palatino Linotype"),
            legend.text=element_text(size = 7),
            legend.title=element_text(size = 10),
            legend.justification = 'center',
            legend.position = 'bottom',
            legend.box = 'vertical',
            legend.box.just = 'center') + guides(size = "none")
  }else if(choosePlot=="norm"){
    # DF.merged <- ldply(DF)
    # min.gr <- min(DF.merged$GeneRatio)
    # max.gr<- max(DF.merged$GeneRatio)
    DF<-DF %>%
      # arrange(desc(!!sortVar.un)) %>%
      arrange(!!sortVar.un) %>%
      mutate("Description" = factor(!!yVar.un, !!yVar.un %>% unique()))
    
    yTextcols <- ifelse(DF %>% pull(category) %in% c("Lymphocyte activation involved in immune response"), '#1E88E5', 
                          ifelse(DF %>% pull(category)  %in% c("Antigen processing and presentation"), '#3C9A83',
                                 ifelse(DF %>% pull(category) %in% c("Immunoglobulin production"), '#C72349',
                                        ifelse(DF %>% pull(category) %in% c("Response to type II interferon"), '#D4841D',
                                               ifelse(DF %>% pull(category)  %in% c("Type II interferon production"), '#9E5F0E',
                                                      ifelse(DF %>% pull(category)  %in% c("Antigen receptor mediated signaling"), '#5EA91F','black'))))))
    # DF[[yVar]] <-gsub("([a-z0-9]* [a-z0-9]*) ", "\\1\n", DF[[yVar]])
    
    DF[[yVar]] <-factor(gsub("([a-z0-9]* [a-z0-9]*) ", "\\1\n", DF[[yVar]]), gsub("([a-z0-9]* [a-z0-9]*) ", "\\1\n", DF[[yVar]]) %>% unique())
    p <-DF %>%
      ggplot(aes(x=!!xVar.un, y=!!yVar.un, size=!!sizeVar.un, fill=!!fillVar.un)) +
      geom_point(alpha=0.7,shape=21, color="black") +
      scale_size(range = c(3,20),limits = c(sizeMin,sizeMax),breaks = seq(sizeMin,sizeMax,length.out=seqN)) + #, name="GO similarity score",limits = c(0,24),breaks = c(2,4,6,8,10,12,14,16,20,24)
      scale_x_continuous(limits = c(0,xMax),breaks=seq.int(from = 0, to=xMax, length.out = 4)%>% round(3)) +
      # scale_size_continuous(name = "GO similarity score",
      #                       breaks = c(2,4,6,8),
      #                       limits = c(0,9),
      #                       labels = as.character(c(2,4,6,8)),
      #                       range = c(0, 24) ) + 
      # scale_fill_viridis(discrete=FALSE, guide=FALSE, option="A") #+
      
      # scale_x_continuous(breaks=seq.int(from = 0, to = maxGR, length.out=4), labels=seq.int(from = 0, to = maxGR, length.out=4))+
      viridis::scale_fill_viridis(option = "viridis", direction=1, alpha=0.5,limits = c(0,fillMax),label = function(x) sprintf("%.1e", x))+ #,breaks = c(2,4,6,8,10,12,14,16),limits = c(0,17)
      theme_bw(base_family = "Palatino Linotype", base_size = 15) +
      theme(axis.text.x = element_text(size = 14),
            # axis.text.y = element_text(size = 14),
            axis.text.y = element_text(size = 15, face="bold",
                                       colour=yTextcols),
            axis.title.y=element_blank(),#element_text(size = 12), 
            axis.title.x=element_text(size = 16),
            text=element_text(family="Palatino Linotype"),
            legend.text=element_text(size = 13),
            legend.title=element_text(size = 15),
            legend.justification = 'center', 
            legend.position = 'bottom', 
            legend.box = 'horizontal', 
            legend.box.just = 'center') + guides(size = guide_legend(override.aes = list(shape=16, color="black"),order = 1,nrow = 2),
                                                 fill = guide_colorbar(order = 0,barwidth = 20,nrow = 1,title.vjust = 1))
    
    
  }
  
  # DF$category <- factor(DF$category, levels=c("Lymphocyte activation involved in immune response", "Antigen processing and presentation", "Immunoglobulin production", "Response to IFNg", "Type II interferon production","Antigen receptor mediated signaling"))
  # myCols <-setNames(c('#1E88E5','#3C9A83','#C72349','#D4841D','#9E5F0E','#5EA91F'),c("Lymphocyte activation involved in immune response", "Antigen processing and presentation", "Immunoglobulin production", "Response to IFNg", "Type II interferon production","Antigen receptor mediated signaling"))
  # 
  # p + ggnewscale::new_scale_color() + geom_line(aes(color=category), data=DF) + scale_color_manual(values=myCols, labels=names(myCols)) + guides(size = guide_legend(override.aes = list(shape=16, color="black")))
  # 
  # p + scale_y_discrete()

  p
}



myTileHeatPlot.2 <- function(DF, xVar, yVar,xLabel, yLabel, labelVar=NULL, colorVar,colorCont=FALSE,facetVar =NULL, removeLegend=TRUE){
  # DF <- df_feats
  # xVar <-"variable"
  # yVar <-"Feature"
  # colorVar <-"value"
  # labelVar <- NULL
  # xLabel = ""
  # yLabel = ""
  # 
  # colorCont=TRUE
  # facetVar =NULL
  # removeLegend=TRUE
  
  
  xVar.un <- as.name(xVar)
  yVar.un <- as.name(yVar)
  colorVar.un <- as.name(colorVar)
  
  if(isTRUE(colorCont)){
    # p <- ggplot(DF, aes(x=!!xVar.un, y=!!yVar.un, fill=!!colorVar.un)) + viridis::scale_fill_viridis(option = "D",
    #                                                                                                  direction=1,
    #                                                                                                  alpha=0.5,
    #                                                                                                  label = function(x) sprintf("%.2f", x)) + geom_tile(alpha=0.7,color="gray")
    # 
    p <- ggplot(DF, aes(x=!!xVar.un, y=!!yVar.un, fill=!!colorVar.un)) + scale_fill_manual(values=c("white","black")) + 
      geom_tile(alpha=0.7,color="gray")
    
  }else{
    p <- ggplot(DF, aes(x=!!xVar.un, y=!!yVar.un, fill=!!colorVar.un)) + scale_fill_manual(values=c("white","black")) + #, + #
      geom_tile(alpha=0.7,color="gray")
    
  }
  if(!is.null(labelVar)){
    labelVar.un <- as.name(labelVar)
    
    p <- p +geom_text(aes(label = round(!!labelVar.un,2)), color="black", size=4, fontface="bold")
    
  }
  
  if(!is.null(facetVar)){
    p <- p + facet_wrap(facetVar,nrow=5,scales='free')
    
  }
  
  p <-p+ scale_x_discrete(position = "top") 
  p <- p + ylab("Boot models")
  p <- p + theme_bw(base_family = "Palatino Linotype", base_size = 12)
  p <- p + theme(legend.position = "bottom",
                 panel.border = element_blank(),#element_rect(colour = "black", size=0.2),
                 # remove the horizontal grid lines
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 axis.text.x = element_text(size=12.5,face="bold", angle=15,hjust=+0.1),
                 axis.text.y = element_text(size=12.5,face="bold"),#element_text(size = 12,face="bold"),#,angle = 90,vjust=0.9
                 axis.ticks.x = element_blank(),
                 axis.ticks.y = element_blank(),
                 # axis.title.x = element_text(size=14, face="bold"),
                 # axis.title.y = element_text(size=14, face="bold"),
                 legend.key = element_rect(size = 2),
                 legend.key.height=unit(0.5,"line"),
                 legend.key.width=unit(1.5,"line"),
                 legend.text = element_text(size = 11),
                 legend.background = element_rect(colour = "black"),
                 legend.direction='horizontal',legend.box='horizontal',
                 legend.title = element_blank())#element_text(colour="black", size=11, face="plain")
  p<-p + labs(x = "", y = "")
  
  if(isTRUE(removeLegend))
    p <- p + theme(legend.position = "none")
  p
}
