# Make histograms for all genes.

# run all of this from your cancer-rna-fish repo
setwd('~/Dropbox/cancer-rna-fish/')

rm(list = ls())

source('plotScripts/RNAFISH/RNAFISH_functions.R')
loadLibraries()

plotDir <- 'graphs/RNAFISH/histogramsAllGenes/'

source('plotScripts/RNAFISH/parameters_WM9_noDrug_20150618.R')
source('plotScripts/RNAFISH/parameters_WM9_noDrug_20150810.R')

### Data = WM9_noDrug_20150618
dataList <- load_WM9_noDrug_20150618(percentileForThresholds = 0.02)
plotName <- dataList[[1]]
noDrugData <- dataList[[3]]

meltNoDrug <- melt(noDrugData,id=c('cellID','Xpos','Ypos'))
ggplot(meltNoDrug, aes(x=value))+geom_histogram()+
  facet_wrap(~variable, scales='free')+
  theme_classic()+
  geom_rug()
ggsave(paste0(plotDir, plotName, 'histogramAllGenes.pdf'),
       height=10,width=10)


### Data = WM9_noDrug_20150810
dataList <- load_WM9_noDrug_20150810(percentileForThresholds = 0.02)
plotName <- dataList[[1]]
noDrugData <- dataList[[3]]

meltNoDrug <- melt(noDrugData,id=c('cellID','Xpos','Ypos'))
ggplot(meltNoDrug, aes(x=value))+geom_histogram()+
  facet_wrap(~variable, scales='free')+
  theme_classic()+
  geom_rug()
ggsave(paste0(plotDir, plotName, 'histogramAllGenes.pdf'),
       height=10,width=10)

##
data983B <- read.csv('dentistData/WM983b_NoDrug_20150525.txt')
plotName <- 'WM983b_noDrug_20150525'
meltData <- melt(data983B,id=c('cellID','Xpos','Ypos'))
ggplot(meltData, aes(x=value))+geom_histogram()+
  facet_wrap(~variable, scales='free')+
  theme_classic()+
  geom_rug()
ggsave(paste0(plotDir, plotName, 'histogramAllGenes.pdf'),
       height=10,width=10)
