# load libraries
loadLibraries <- function(){
  require(grid)
  require(gridExtra)
  library(ggplot2)
  library(data.table)
  library(tidyr)
  library(dplyr)
  library(reshape2)
  library(gridExtra)
  library(ggthemes)
  library(MASS)
  library(abd)
  library(RColorBrewer)
  library(mosaic)
}

scattersForJackpotGenes <- function(plotName, plotDir, data){
  noDrugDataMelt <- melt(data, id=c('cellID','Xpos','Ypos'))
  noDrugDataMeltSubset <- filter(noDrugDataMelt, (variable=='EGFR' |
                                                    variable=='WNT5A' | variable=='SERPINE1' |
                                                    variable=='NGFR' | variable=='NRG1' |
                                                    variable=='VEGFC' | variable=='AXL' |
                                                    variable=='LOXL2' | variable=='JUN'))
  noDrugDataMeltSubset$jitter <- jitter(noDrugDataMeltSubset$value,amount=1)
  noDrugDataMeltSubset$jitter[which(noDrugDataMeltSubset$jitter<0)] <- 0
  
  ggplot(noDrugDataMeltSubset, aes(x=jitter,fill=variable))+
    geom_histogram(fill='darkolivegreen1',color='black',size=0.2)+
    geom_rug(aes(color=variable),size=0.3,color='black')+
    theme_classic()+
    facet_wrap(~variable,scales='free')+
    theme(legend.position="none")+
    theme(text=element_text(size=20))+
    theme(axis.line=element_line(size=0.1))
  ggsave(paste0(plotDir, plotName, 'geneHistogramRugs.pdf'),height=12,width=12)
  
  # makes SERPINE1 callout to show zoomed in view.
  temp <- filter(noDrugDataMeltSubset,variable=='SERPINE1')
  ggplot(temp, aes(x=jitter,fill=variable))+
    geom_histogram(fill='darkolivegreen1',color='black',size=0.2)+
    geom_rug(aes(color=variable),size=0.3,color='black')+
    theme_classic()+
    theme(legend.position="none")+
    theme(text=element_text(size=20))+
    ylim(c(0,50))+
    theme(axis.line=element_line(size=0.5))
  ggsave(paste0(plotDir, plotName, 'SERPINE1_jackpotcallout.pdf'))
  
  # Shows two non-jackpot genes
  noDrugDataMeltSubset <- filter(noDrugDataMelt, (variable=='GAPDH' | variable =='MITF' | variable =='SOX10'))
  noDrugDataMeltSubset$jitter <- jitter(noDrugDataMeltSubset$value,amount=1)
  noDrugDataMeltSubset$jitter[which(noDrugDataMeltSubset$jitter<0)] <- 0
  ggplot(noDrugDataMeltSubset, aes(x=jitter,fill=variable))+
    geom_histogram(fill='darkolivegreen1',color='black',size=0.2)+
    geom_rug(aes(color=variable),size=0.3,color='black')+
    theme_classic()+
    facet_wrap(~variable,scales='free')+
    theme(legend.position="none")+
    theme(text=element_text(size=20))+
    theme(axis.line=element_line(size=0.1))
  ggsave(paste0(plotDir, plotName, 'geneHistogramRugs_NonJackpotGenes.pdf'), height=5, width=12)
  
}


twoGeneScatters <- function(plotName, plotDir, data, geneThresholds,gene1, gene2,jitter){
  
  # if we don't have 2 genes defined as inputs to the function, this will just default to EGFR and WNT5A
  if(missing(gene1) | missing(gene2)){
    gene1 <- 'EGFR'
    gene2 <- 'WNT5A'
  }
  dataMelt <- melt(data, id=c('cellID','Xpos','Ypos'))
  
  if(missing(jitter)){
    jitter=1
  }
  
  
  ### Two gene scatter plot for panel C
  meltedPair <- filter(dataMelt, variable==gene1 | variable==gene2)
  meltedPair$jitter <- jitter(meltedPair$value,amount=jitter)
  castPair <- dcast(meltedPair, 'cellID+Xpos+Ypos~variable',value.var='jitter')
  colnames(castPair)[c(4,5)] <- c('gene1','gene2')
  
  myPalette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")), space="Lab")
  ggplot(castPair, aes(x=gene1,y=gene2))+
    stat_bin2d(bins=150)+
    theme_classic()+
    theme(text=element_text(size=20))+
    geom_rug(color='black')+
    xlab(eval(gene2))+ylab(eval(gene1))+
    geom_vline(xintercept=geneThresholds[which(geneThresholds$geneName==gene2),2],linetype='dashed',size=0.5)+
    geom_hline(yintercept=geneThresholds[which(geneThresholds$geneName==gene1),2],linetype='dashed',size=0.5)+
    scale_fill_gradientn(colours = myPalette(10))
  ggsave(paste0(plotDir,plotName, eval(gene1),'vs',eval(gene2),'_withLegend.pdf'), height = 6, width = 6)
  
  ggplot(castPair, aes(x=gene1,y=gene2))+
    stat_bin2d(bins=150)+
    theme_classic()+
    theme(legend.position="none")+
    theme(text=element_text(size=20))+
    xlab(eval(gene2))+ylab(eval(gene1))+
    geom_rug(color='black')+
    geom_vline(xintercept=geneThresholds[which(geneThresholds$geneName==gene2),2],linetype='dashed',size=0.5)+
    geom_hline(yintercept=geneThresholds[which(geneThresholds$geneName==gene1),2],linetype='dashed',size=0.5)+
    scale_fill_gradientn(colours = myPalette(10))
  ggsave(paste0(plotDir,plotName, eval(gene1),'vs',eval(gene2),'_noLegend.pdf'),height=6, width = 6)
  
}

applyThresholds <- function(data, geneThresholds) {
  dataout <- data
  for (i in 1:length(geneThresholds$geneName)) {
    currGene <- as.character(geneThresholds$geneName[i])
    #dataout[[currGene]] <- as.integer(data[[currGene]]>thresholds$threshold[i])
    dataout[[currGene]] <- as.integer(data[[currGene]]>geneThresholds$threshold[i])
  }
  return(dataout)
}

# applyThresholds <- function(data, geneThresholds) {
#   dataout <- data
#   for (i in 1:(dim(data)[2]-3)) {
#     threshold <- filter(geneThresholds, geneName == colnames(data[i+3]))
#     dataout[,i+3] <- as.integer(data[,i+3]>threshold$threshold)
#   }
#   return(dataout)
# }


oddsRatioHeatmap <- function(plotName, plotDir, data, geneThresholds){
  
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
  allORsOrdered <- allORs
  allORsOrdered$gene1 <- factor(allORsOrdered$gene1,c('SOX10','MITF','GAPDH','CCNA2','VGF','FOSL1','RUNX2',
                                                      'PDGFC',
                                                      'PDGFRB','NRG1','EGFR','LOXL2','FGFR1','SERPINE1','NGFR','WNT5A',
                                                      'JUN','AXL','VEGFC'))
  allORsOrdered$gene2 <- factor(allORsOrdered$gene2,c('SOX10','MITF','GAPDH','CCNA2','VGF','FOSL1','RUNX2',
                                                      'PDGFC',
                                                      'PDGFRB','NRG1','EGFR','LOXL2','FGFR1','SERPINE1','NGFR','WNT5A',
                                                      'JUN','AXL','VEGFC'))
  
  oddsRatioMat <- acast(allORsOrdered, gene1~gene2)
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
    scale_fill_gradientn(colours = myPalette(10),values=values)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    labs(title='log2(odds ratio) heatmap')
  ggsave(paste0(plotDir,plotName,'oddsRatioHeatmap_WithLegend.pdf'),height=8,width=8)
  
  ggplot(melted_ORmat,aes(x=Var1,y=Var2,fill=log2(valueCut)))+
    geom_tile(color='white')+ theme_classic()+
    scale_fill_gradientn(colours = myPalette(10),values=values)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    labs(title='log2(odds ratio) heatmap')+
    theme(legend.position='none')
  ggsave(paste0(plotDir,plotName,'oddsRatioHeatmap_NoLegend.pdf'),height=8,width=8)
  
}



numberJackpotsHistogram <- function(plotName, plotDir, data, geneThresholds){
  
  thresholdData <- applyThresholds(data, geneThresholds)
  
  # Remove genes from this analysis that are not resistance markers - GAPDH, CCNA2, SOX10, MITF, VGF
  thresholdData_markersOnly <- thresholdData[,c(1:4,8:16,17:21)]
  numJackpotsPerCell <- data.frame(cellID = thresholdData$cellID,
                                   numberJackpotGenes = margin.table(as.matrix(thresholdData_markersOnly[4:dim(thresholdData_markersOnly)[2]]),1) )
  numCells <- dim(data)[2]
  ggplot(numJackpotsPerCell,aes(x=numberJackpotGenes))+
    geom_bar(aes(y=100*(..count..)/sum(..count..)),binwidth=1,fill='goldenrod2',color='black',size=0.2)+
    theme_classic()+
    stat_bin(aes(y=100*(..count..)/sum(..count..), label=round(100*(..count..)/sum(..count..),2)), 
             geom="text", size=4, binwidth = 1, vjust=-1.5)+
    #stat_bin(aes(y=(100*(..count..)/sum(..count..))+10, label=round(length(cellID)*(100*(..count..)/sum(..count..)),2)), 
    #         geom="text", size=4, binwidth = 1, vjust=-1.5)+
    theme(text=element_text(size=20))+
    scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9),limits=c(0,9))+
    scale_y_continuous(limits=c(0,100))+
    xlab('number of jackpot genes per cell')+
    ylab('percentage of cells')
  ggsave(paste0(plotDir,plotName,'percentJackpotsHistogram.pdf'),height=6,width=7)
  
  ggplot(numJackpotsPerCell,aes(x=numberJackpotGenes))+
    geom_bar(aes(y=(..count..)),binwidth=1,fill='goldenrod2',color='black',size=0.2)+
    theme_classic()+
    stat_bin(aes(y=(..count..), label=round((..count..),4)), 
             geom="text", size=4, binwidth = 1, vjust=-1.5)+
    theme(text=element_text(size=20))+
    scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9),limits=c(0,9))+
    #scale_y_continuous(limits=c(0,100))+
    xlab('number of jackpot genes per cell')+
    ylab('number of cells')
  ggsave(paste0(plotDir,plotName,'numberJackpotsHistogram.pdf'),height=6,width=7)
  
}

gapdhNormalizeData <- function(data,minGapdh){
  dataMelt <- melt(data, id=c('cellID','Xpos','Ypos','GAPDH'))
  dataMelt <- filter(dataMelt, GAPDH>minGapdh) 
  dataMelt$gapdhNorm <- dataMelt$value/dataMelt$GAPDH
  dataCast <- dcast(dataMelt, 'cellID+Xpos+Ypos~variable',value.var='gapdhNorm')
  dataCast$GAPDH <- 1
  dataCast <- dataCast[,c(1:6,22,7:21)]
  return(dataCast)
}


scattersForGAPDHnormData <- function(plotName, plotDir, data){
  noDrugDataMelt <- melt(data, id=c('cellID','Xpos','Ypos'))
  noDrugDataMeltSubset <- filter(noDrugDataMelt, (variable=='EGFR' |
                                                    variable=='WNT5A' | variable=='SERPINE1' |
                                                    variable=='NGFR' | variable=='NRG1' |
                                                    variable=='VEGFC' | variable=='AXL' |
                                                    variable=='LOXL2' | variable=='JUN'))
  
  ggplot(noDrugDataMeltSubset, aes(x=value,fill=variable))+
    geom_histogram(fill='darkolivegreen1',color='black',size=0.2)+
    geom_rug(aes(color=variable),size=0.3,color='black')+
    theme_classic()+
    facet_wrap(~variable,scales='free')+
    theme(legend.position="none")+
    theme(text=element_text(size=20))+
    theme(axis.line=element_line(size=0.1))
  ggsave(paste0(plotDir, plotName, 'geneHistogramRugs.pdf'),height=12,width=12)
  
  # Shows two non-jackpot genes
  noDrugDataMeltSubset <- filter(noDrugDataMelt, (variable=='GAPDH' | variable =='MITF' | variable =='SOX10'))
  ggplot(noDrugDataMeltSubset, aes(x=value,fill=variable))+
    geom_histogram(fill='darkolivegreen1',color='black',size=0.2)+
    geom_rug(aes(color=variable),size=0.3,color='black')+
    theme_classic()+
    facet_wrap(~variable,scales='free')+
    theme(legend.position="none")+
    theme(text=element_text(size=20))+
    theme(axis.line=element_line(size=0.1))
  ggsave(paste0(plotDir, plotName, 'geneHistogramRugs_NonJackpotGenes.pdf'), height=5, width=12)
  
}

histogramsForPrimaryMelanocyteGenes <- function(plotName, plotDir, data){
  noDrugDataMelt <- melt(data, id=c('cellID','Xpos','Ypos'))
  noDrugDataMeltSubset <- filter(noDrugDataMelt, (variable=='FGFR1' |
                                                    variable=='AXL' | variable=='NGFR' |
                                                    variable=='MITF' | variable=='SERPINE1' |
                                                    variable=='GAPDH' | variable=='SOX10' |
                                                    variable=='CCNA2'))
  noDrugDataMeltSubset$jitter <- jitter(noDrugDataMeltSubset$value,amount=1)
  noDrugDataMeltSubset$jitter[which(noDrugDataMeltSubset$jitter<0)] <- 0
  
  ggplot(noDrugDataMeltSubset, aes(x=jitter,fill=variable))+
    geom_histogram(fill='darkolivegreen1',color='black',size=0.2)+
    geom_rug(aes(color=variable),size=0.3,color='black')+
    theme_classic()+
    facet_wrap(~variable,scales='free')+
    theme(legend.position="none")+
    theme(text=element_text(size=20))+
    theme(axis.line=element_line(size=0.1))
  ggsave(paste0(plotDir, plotName, 'geneHistogramRugs.pdf'),height=12,width=12)
}


scattersForJackpotGenesMelanocytes <- function(plotName, plotDir, data){
  dataMelt <- melt(data, id=c('cellID','Xpos','Ypos'))
  dataMeltSubset <- filter(dataMelt, (variable=='SERPINE1' | variable=='NGFR'| variable=='AXL'))
  dataMeltSubset$jitter <- jitter(dataMeltSubset$value,amount=1)
  dataMeltSubset$jitter[which(dataMeltSubset$jitter<0)] <- 0
  
  ggplot(dataMeltSubset, aes(x=jitter,fill=variable))+
    geom_histogram(fill='darkolivegreen1',color='black',size=0.2)+
    geom_rug(aes(color=variable),size=0.3,color='black')+
    theme_classic()+
    facet_wrap(~variable,scales='free')+
    theme(legend.position="none")+
    theme(text=element_text(size=20))+
    theme(axis.line=element_line(size=0.1))
  ggsave(paste0(plotDir, plotName, 'geneHistogramRugs.pdf'),height=3,width=7)
}