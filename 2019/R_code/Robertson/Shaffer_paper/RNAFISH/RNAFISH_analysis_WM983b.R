# Plots jackpot genes in WM983b 
# This is for the supplement

source('plotScripts/RNAFISH/parameters_WM983b_noDrug_20150525.R')

dataList <- load_WM983b_noDrug_20150525(percentileForThresholds = 0.02)
plotName <- dataList[[1]]
plotDir <- dataList[[2]]
noDrugData <- dataList[[3]]
thresholds <- dataList[[4]]

noDrugDataMelt <- melt(noDrugData, id=c('cellID','Xpos','Ypos'))
noDrugDataMeltSubset <- filter(noDrugDataMelt, (variable=='CYR61' |
                                                  variable=='NRG1' | variable=='RUNX2' |
                                                  variable=='NGFR' | variable=='EPHA2' |
                                                  variable=='LOXL2' | variable=='AXL' |
                                                  variable=='FGFR1' | variable=='JUN'))
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