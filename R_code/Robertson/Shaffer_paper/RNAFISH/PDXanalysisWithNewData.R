setwd('~/Dropbox (RajLab)/Papers/cancerpaper_public/paper/')

rm(list = ls())

# load tissue data
data <- read.csv('dentistData/PDX_20160115_WM4335_nonResistant.txt')
source('plotScripts/RNAFISH/RNAFISH_functions.R')
loadLibraries()

# how many cells are high in CYR61 and LOXL2?
filter(data, CYR61>15 & LOXL2>15) # 9 total

# We want to remove cells that don't really have much GAPDH signal 
# (because in these situations there may be a problem with the FISH)
data <- filter(data,GAPDH>=5)
# now we have a total of 6329 cells.

# What is threshold at 0.02?
quantile(data$CYR61,0.98) #9
quantile(data$CYR61,0.99) #12.57
## Arjun suggested ~15 for CYR61
quantile(data$LOXL2,0.98) #19
quantile(data$LOXL2,0.99) #24

# set thresholds here.
CYR61_thresh <- 12
LOXL2_thresh <- 24

### histogram plots
dataMelt <- melt(data, id=c('Xpos','Ypos'))
dataMeltSubset <- filter(dataMelt, (variable=='CYR61' | variable=='LOXL2'))
dataMeltSubset$jitter <- jitter(dataMeltSubset$value,amount=1)
dataMeltSubset$jitter[which(dataMeltSubset$jitter<0)] <- 0

ggplot(dataMeltSubset, aes(x=value,fill=variable))+
  geom_histogram(fill='darkolivegreen1',color='black',size=0.2)+
  geom_rug(aes(x=jitter,color=variable),size=0.3,color='black')+
  theme_classic()+
  facet_wrap(~variable,scales='free')+
  theme(legend.position="none")+
  theme(text=element_text(size=20))+
  theme(axis.line=element_line(size=0.1))+
  
ggsave('graphs/PDXs/PDX_20160115_histograms_CYR61_LOXL2.pdf',height=4,width=6)

gene1 <- 'CYR61'
gene2 <- 'LOXL2'
meltedPair <- filter(dataMelt, variable==gene1 | variable==gene2)
meltedPair$jitter <- jitter(meltedPair$value,amount=1)
castPair <- dcast(meltedPair, 'Xpos+Ypos~variable',value.var='jitter')
colnames(castPair)[c(3,4)] <- c('gene1','gene2')

myPalette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")), space="Lab")
ggplot(castPair, aes(x=gene1,y=gene2))+
  stat_bin2d(bins=150)+
  theme_classic()+
  theme(text=element_text(size=20))+
  geom_rug(color='black')+
  xlab(eval(gene1))+ylab(eval(gene2))+
  geom_vline(xintercept=CYR61_thresh,linetype='dashed',size=0.5)+
  geom_hline(yintercept=LOXL2_thresh,linetype='dashed',size=0.5)+
  scale_fill_gradientn(colours = myPalette(10))
ggsave('graphs/PDXs/PDX_20160115_twoGeneScatter_CYR61_LOXL2_withLegend.pdf', height = 6, width = 7)

ggplot(castPair, aes(x=gene1,y=gene2))+
  stat_bin2d(bins=150)+
  theme_classic()+
  theme(text=element_text(size=20))+
  geom_rug(color='black')+
  xlab(eval(gene1))+ylab(eval(gene2))+
  geom_vline(xintercept=CYR61_thresh,linetype='dashed',size=0.5)+
  geom_hline(yintercept=LOXL2_thresh,linetype='dashed',size=0.5)+
  scale_fill_gradientn(colours = myPalette(10))+
  theme(legend.position="none")
ggsave('graphs/PDXs/PDX_20160115_twoGeneScatter_CYR61_LOXL2.pdf', height = 6, width = 6)


## Contengency table and OR calculations
castPair <- dcast(meltedPair, 'Xpos+Ypos~variable',value.var='value')
thresholdData <- data.frame(CYR61 = as.integer(castPair$CYR61>CYR61_thresh),
                            LOXL2 = as.integer(castPair$LOXL2>LOXL2_thresh))
crossTable<-table(thresholdData)
OR <- 1/mosaic::oddsRatio(crossTable,conf.level = 0.95) 

print(paste('Contengency table: '))
print(eval(crossTable))
print(paste('Odds ratio: ',eval(OR)))

## Gini coeff plots:
library('ineq')
gini <- apply(data[,c(3,5,6)],2,function(x) ineq(x,type='Gini'))
gini <- melt(gini)
gini$gene <- rownames(gini)

ggplot(gini, aes(x=value,label=gene))+
  geom_point(aes(y=1,color=gene))+
  geom_text(y=1,angle=55,hjust=-0.05,vjust=-0.2,aes(color=gene),size=6)+
  theme_classic()+xlim(0,1)+ylim(0.9,1.1)
ggsave('graphs/PDXs/PDX_20160115_gini.pdf',height=5,width=10)


# Total number of cells in this data
print(paste('Total number of cells in the PDX data in Figure 3:', dim(data)[1]))


## Replicate data
data <- read.csv('dentistData/PDX_20150715_WM4299_tissueReplicate.txt')

### histogram plots
dataMelt <- melt(data, id=c('Xpos','Ypos'))
dataMeltSubset <- filter(dataMelt, (variable=='CYR61' | variable=='LOXL2'))
dataMeltSubset$jitter <- jitter(dataMeltSubset$value,amount=1)
dataMeltSubset$jitter[which(dataMeltSubset$jitter<0)] <- 0

ggplot(dataMeltSubset, aes(x=value,fill=variable))+
  geom_histogram(fill='darkolivegreen1',color='black',size=0.2)+
  geom_rug(aes(x=jitter,color=variable),size=0.3,color='black')+
  theme_classic()+
  facet_wrap(~variable,scales='free')+
  theme(legend.position="none")+
  theme(text=element_text(size=20))+
  theme(axis.line=element_line(size=0.1))
  
  ggsave('graphs/PDXs/PDX_20150715_histograms_CYR61_LOXL2.pdf',height=4,width=6)

## Gini coeff plots:
library('ineq')
gini <- apply(data[,c(3,5,6)],2,function(x) ineq(x,type='Gini'))
gini <- melt(gini)
gini$gene <- rownames(gini)

ggplot(gini, aes(x=value,label=gene))+
  geom_point(aes(y=1,color=gene))+
  geom_text(y=1,angle=55,hjust=-0.05,vjust=-0.2,aes(color=gene),size=6)+
  theme_classic()+xlim(0,1)+ylim(0.9,1.1)
ggsave('graphs/PDXs/PDX_20150715_gini.pdf',height=5,width=10)


# Total number of cells in this data
print(paste('Total number of cells in the PDX data in Figure 3:', dim(data)[1]))


# October 2016 data
data1 <- read.csv('dentistData/20161024_WM4325_NGFR_AXL_LOXL2.txt')
data1$tissue <- 'WM4335' # this was mislabeled as WM4325, should be WM4335
data1Melt <- melt(data1, id.vars = c('Xpos','Ypos','tissue'))
data2 <- read.csv('dentistData/20161027_WM3909_NGFR_AXL_LOXL2.txt')
data2$tissue <- 'WM3909'
data2Melt <- melt(data2, id.vars = c('Xpos','Ypos','tissue'))
data3 <- read.csv('dentistData/20161028_WM4325_NRG1_WNT5A_GAPDH.txt')
data3$tissue <- 'WM4335' # this was mislabeled as WM4325, should be WM4335
data3Melt <- melt(data3, id.vars = c('Xpos','Ypos','tissue'))

allData <- rbind(data1Melt,data2Melt,data3Melt)

allData$jitter <- jitter(allData$value,amount=1)
allData$jitter[which(dataMeltSubset$jitter<0)] <- 0

ggplot(allData, aes(x=value,fill=variable))+
  geom_histogram(fill='darkolivegreen1',color='black',size=0.2)+
  geom_rug(aes(x=jitter,color=variable),size=0.3,color='black')+
  theme_classic()+
  facet_wrap(tissue~variable,scales='free')+
  theme(legend.position="none")+
  theme(text=element_text(size=20))+
  theme(axis.line=element_line(size=0.1))
ggsave('graphs/PDXs/histogramsAllGenesNew.pdf',height=15, width=15)
