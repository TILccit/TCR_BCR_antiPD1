library()

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
  df <- df %>%  dplyr::mutate(prop = round(counts*100/sum(counts), 1),
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
  df <- df %>%  dplyr::mutate(prop = round(counts*100/sum(counts), 1),
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
