
setwd('~/Dropbox/cancer-rna-fish/')

source('plotScripts/RNAFISH/RNAFISH_functions.R')
loadLibraries()

source('plotScripts/RNAFISH/parameters_WM983b_noDrug_20150525.R')
loadLibraries()

dataList <- load_WM983b_noDrug_20150525(percentileForThresholds = 0.02)
plotName <- dataList[[1]]
plotDir <- 'graphs/RNAFISH/'
data <- dataList[[3]]
geneThresholds <- dataList[[4]]

  
  thresholdData <- applyThresholds(data, geneThresholds)
  
  allORs <- data.frame(matrix(ncol=3,nrow=dim(thresholdData[,c(4:22)])[2]^2))
  colnames(allORs)<-c('gene1','gene2','OR')
  k=0
  for (j in 1:(dim(thresholdData)[2]-3)){
    for (i in 1:(dim(thresholdData)[2]-3)){    
      allORs$gene1[i+19*k]<- as.character(colnames(thresholdData)[j+3])
      allORs$gene2[i+19*k] <- as.character(colnames(thresholdData)[i+3])
      crossTable<-table(thresholdData[,c(i+3,j+3)])
      if (dim(crossTable)[1]<2 | dim(crossTable)[2]<2) 
      {allORs$OR[i+19*k] <- NA}
      else {allORs$OR[i+19*k] <- 1/oddsRatio(crossTable,conf.level = 0.95)}
    }
    k=k+1
  }
  temp <- allORs[is.finite(allORs$OR),] 
  print(paste('Maximum odds ratio between resistance markers:',temp[which(temp$OR==max(temp$OR)),]))
  temp2 <- temp[which(log2(temp$OR)>0 & temp$gene1!='VGF' & temp$gene2!='VGF' & temp$gene1!='CCNA2' & temp$gene2!='CCNA2' &
                        temp$gene1!='GAPDH' & temp$gene2!='GAPDH' & temp$gene1!='MITF' & temp$gene2!='MITF' &
                        temp$gene1!='SOX10' & temp$gene2!='SOX10'),]
  print(paste('Minimum odds ratio between resistance markers:',temp2[which(temp2$OR==min(temp2$OR)),]))
  
  
  # which genes do not have any cells above threshold?
  # margin.table(as.matrix(thresholdData),2) # PDGFC, RUNX2, VGF
  
  #change ordering of the variables
  temp <- allORs[is.finite(allORs$OR),]
  allORsOrdered <- temp
  allORsOrdered$gene1 <- factor(allORsOrdered$gene1,c('SOX10','MITF','GAPDH','CCNA2','VGF','FOSL1','RUNX2','CYR61',
                                                      'PDGFC',
                                                      'PDGFRB','NRG1','EGFR','LOXL2','FGFR1','SERPINE1','NGFR','WNT5A',
                                                      'JUN','AXL','VEGFC'))
  allORsOrdered$gene2 <- factor(allORsOrdered$gene2,c('SOX10','MITF','GAPDH','CCNA2','VGF','FOSL1','RUNX2','CYR61',
                                                      'PDGFC',
                                                      'PDGFRB','NRG1','EGFR','LOXL2','FGFR1','SERPINE1','NGFR','WNT5A',
                                                      'JUN','AXL','VEGFC'))
  
  oddsRatioMat <- acast(allORsOrdered, gene1~gene2,value.var='OR')
  lower_tri <- oddsRatioMat
  melted_ORmat <- melt(lower_tri)
  melted_ORmat <- na.omit(melted_ORmat)
  
  myPalette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")), space="Lab")
  
  # we want the yellow color at zero, so that anything blue is odds ratio less than 1
  #maxValue <- log2(max(filter(melted_ORmat,value < Inf)$value)) #max value that is not Inf.
  maxValue <- 6
  values <- seq(-maxValue, maxValue, length=11)
  values <- c(-1.25,-1,-0.75,-0.5,-0.25,values[6:11])
  melted_ORmat$valueCut <- melted_ORmat$value
  melted_ORmat$valueCut[which(melted_ORmat$value>2^6 & melted_ORmat$value<Inf)] <- 2^6
  
  ggplot(melted_ORmat,aes(x=Var1,y=Var2,fill=log2(valueCut)))+
    geom_tile(color='white')+ theme_classic()+
    scale_fill_gradientn(colours = myPalette(10),values=values,rescaler = function(x, ...) x, oob = identity)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    labs(title='log2(odds ratio) heatmap')
  ggsave(paste0(plotDir,plotName,'oddsRatioHeatmap_WithLegend.pdf'),height=8,width=8)
  
  ggplot(melted_ORmat,aes(x=Var1,y=Var2,fill=log2(valueCut)))+
    geom_tile(color='white')+ theme_classic()+
    scale_fill_gradientn(colours = myPalette(10),values=values,rescaler = function(x, ...) x, oob = identity)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    labs(title='log2(odds ratio) heatmap')+
    theme(legend.position='none')
  ggsave(paste0(plotDir,plotName,'oddsRatioHeatmap_NoLegend.pdf'),height=8,width=8)
  
#}
