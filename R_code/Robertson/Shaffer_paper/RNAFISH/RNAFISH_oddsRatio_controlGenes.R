
setwd('~/Dropbox/cancer-rna-fish/')

source('plotScripts/RNAFISH/RNAFISH_functions.R')
source('plotScripts/RNAFISH/parameters_WM9_noDrug_20150810.R')
source('plotScripts/RNAFISH/parameters_WM9_noDrug_20150618.R')

source('plotScripts/RNAFISH/RNAFISH_functions.R')
loadLibraries()

## control data
controlData1 <- read.csv('dentistData/20160502_WM9_noDrug_controlGenes.txt')
controlData1$cellID <- 1:dim(controlData1)[1]
controlData2 <- read.csv('dentistData/20160503_WM9_noDrug_controlGenes.txt')
controlData2$cellID <- 1:dim(controlData2)[1]
#data <- cbind(controlData1,controlData2)

thresholds <- data.frame(geneName = colnames(controlData1)[3:6],threshold = rep(0,4))


## apply thresholding
for (i in 1:4){
  percentThresh <- quantile(controlData1[,c(i+2)],0.98)
  thresholds$threshold[i] <- percentThresh 
}

thresholdData <- applyThresholds(controlData1, thresholds)

# BABAM1-3 vs. KDM5A-4
crossTable<-table(thresholdData[,c(3,4)])
1/mosaic::oddsRatio(crossTable,conf.level = 0.95) #15.97 as GAPDH norm -- 41.74 if not.

# BABAM1-3 vs. LMNA-5
crossTable<-table(thresholdData[,c(3,5)])
1/mosaic::oddsRatio(crossTable,conf.level = 0.95) #22.04 as GAPDH norm -- 43.599 if not.

# KDM5A vs. LMNA
crossTable<-table(thresholdData[,c(4,5)])
1/mosaic::oddsRatio(crossTable,conf.level = 0.95) #10.7 as GAPDH norm -- 33.46 if not.


head(controlData1)
controlData1$KDM5Anorm <- controlData1$KDM5A/controlData1$GAPDH
controlData1$LMNAnorm <- controlData1$LMNA/controlData1$GAPDH
controlData1$BABAM1norm <- controlData1$BABAM1/controlData1$GAPDH

controlData1 <- read.csv('dentistData/20160502_WM9_noDrug_controlGenes.txt')
controlData1$cellID <- 1:dim(controlData1)[1]
data <- controlData1
dataMelt <- melt(data, id=c('Xpos','Ypos','GAPDH','cellID'))
dataMelt <- filter(dataMelt, GAPDH>200) 
dataMelt$gapdhNorm <- dataMelt$value/dataMelt$GAPDH
dataCast <- dcast(dataMelt, 'Xpos+Ypos+cellID~variable',value.var='gapdhNorm')
dataCast$GAPDH <- 1

thresholds <- data.frame(geneName = colnames(dataCast)[4:6],threshold = rep(0,3))

for (i in 1:3){
  percentThresh <- quantile(dataCast[,c(i+3)],0.98)
  thresholds$threshold[i] <- percentThresh 
}

thresholdData <- applyThresholds(dataCast, thresholds)
thresholdData$numHighCells <- rowSums(thresholdData[,c(4:6)])

filter(thresholdData, numHighCells>2)
filter(controlData1, cellID==412)
filter(dataCast, cellID==412)

filter(controlData1,BABAM1==66)

# BABAM1-3 vs. KDM5A-4
crossTable<-table(thresholdData[,c(4,5)])
1/mosaic::oddsRatio(crossTable,conf.level = 0.95) #15.97

# BABAM1-3 vs. LMNA-5
crossTable<-table(thresholdData[,c(4,6)])
1/mosaic::oddsRatio(crossTable,conf.level = 0.95) #22.04

# KDM5A vs. LMNA
crossTable<-table(thresholdData[,c(5,6)])
1/mosaic::oddsRatio(crossTable,conf.level = 0.95) #10.7





## data set 2
controlData2 <- read.csv('dentistData/20160503_WM9_noDrug_controlGenes.txt')
data <- controlData2
dataMelt <- melt(data, id=c('Xpos','Ypos','GAPDH'))
dataMelt <- filter(dataMelt, GAPDH>20) 
dataMelt$gapdhNorm <- dataMelt$value/dataMelt$GAPDH
dataCast <- dcast(dataMelt, 'Xpos+Ypos~variable',value.var='gapdhNorm')
dataCast$GAPDH <- 1

thresholds <- data.frame(geneName = colnames(dataCast)[3:5],threshold = rep(0,3))

for (i in 1:3){
  percentThresh <- quantile(dataCast[,c(i+2)],0.98)
  thresholds$threshold[i] <- percentThresh 
}

thresholdData <- applyThresholds(dataCast, thresholds)

# LMNA vs. KDM5B-4
crossTable<-table(thresholdData[,c(3,4)])
1/mosaic::oddsRatio(crossTable,conf.level = 0.95) #15.97



## Try PCA on the data in figure 3
### Data = WM9_noDrug_20150618
dataList <- load_WM9_noDrug_20150618(percentileForThresholds = 0.02)
plotName <- dataList[[1]]
plotDir <- dataList[[2]]
noDrugData <- dataList[[3]]
thresholds <- dataList[[4]]
gapdhThresholds <- dataList[[5]]
gapdhNormData <- gapdhNormalizeData(noDrugData, 50)

### Data = WM9_noDrug_20150810
dataList <- load_WM9_noDrug_20150810(percentileForThresholds = 0.02)
plotName <- dataList[[1]]
plotDir <- dataList[[2]]
noDrugData <- dataList[[3]]
thresholds <- dataList[[4]]
gapdhThresholds <- dataList[[5]]
gapdhNormData <- gapdhNormalizeData(noDrugData, 50)

noDrugWM9 <- noDrugData[,c(4:22)]
noDrugWM9.log <- log(noDrugWM9)
noDrugWM9.pca <- prcomp(noDrugWM9.log, center=TRUE, scale. = TRUE) # error because of log(0)

noDrugWM9.log <- log(noDrugWM9+0.001)
noDrugWM9.pca <- prcomp(data.frame(noDrugWM9), scale. = TRUE)

plot(noDrugWM9.pca, type = "l")
#library(devtools)
#install_github("ggbiplot", "vqv")

#library(ggbiplot)
g <- ggbiplot(noDrugWM9.pca, obs.scale = 1, var.scale = 1, ellipse = TRUE, 
              circle = TRUE)
# g <- g + scale_color_discrete(name = '')
# g <- g + theme(legend.direction = 'horizontal', 
#                legend.position = 'top')
print(g)

##
noDrugWM9.pca <- princomp(data.frame(noDrugWM9.log), cor=TRUE, score=TRUE)
##
pdf(paste0(plotDir,plotName,'PCA_biplot_log.pdf'))
biplot(noDrugWM9.pca,xlabs=rep(".", dim(noDrugWM9)[1]))
dev.off()

summary(noDrugWM9.pca)

noDrugWM9.pca <- princomp(data.frame(noDrugWM9), cor=TRUE, score=TRUE)
##
pdf(paste0(plotDir,plotName,'PCA_biplot.pdf'))
biplot(noDrugWM9.pca,xlabs=rep(".", dim(noDrugWM9)[1]))
dev.off()

summary(noDrugWM9.pca)

## new data
data3 <- read.csv('dentistData/20160722_WM9_noDrug_JackpotAndControlGenes.txt')
data3 <- data3[,3:6]
noDrugWM9.pca <- princomp(data.frame(data3), cor=TRUE, score=TRUE)
pdf('graphs/RNAFISH/WM9_noDrug_20160722_WM9_noDrug_controlAndJackpot_PCA_biplot.pdf',
    height=8,width=8)
biplot(noDrugWM9.pca, xlabs=rep(".", dim(data3)[1]))
dev.off()
# 56% of variance in component 1
# 28% of variance in component 2
data3[8104,]

#library(vegan)
#mod <- rda(data3, scale = TRUE)
#biplot(mod, scaling = 0.0001)

col4 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","#7FFF7F", 
                           "cyan", "#007FFF", "blue","#00007F"))  
library(corrplot)
M <- cor(data3[,c(2,3,1,4)])
pdf('graphs/RNAFISH/WM9_noDrug_20160722_corrplot.pdf',height=6,width=6)
corrplot(M, method="number",col=col4(20),sig.level=0.01)
dev.off()

## burn in genes
data <-read.csv('dentistData/20160725_WM9noDrugA6G3.txt')

col4 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","#7FFF7F", 
                           "cyan", "#007FFF", "blue","#00007F"))  
library(corrplot)
M <- cor(data[,c(3:6)])
#pdf('graphs/RNAFISH/WM9_noDrug_20160722_corrplot.pdf',height=6,width=6)
corrplot(M, method="number",col=col4(20),sig.level=0.01)

data.pca <- princomp(data.frame(data[,c(3:6)]), cor=TRUE, score=TRUE)
pdf('graphs/RNAFISH/WM9_noDrug_20160725.pdf',height=6,width=6)
biplot(data.pca, xlabs=rep(".", dim(data)[1]))
dev.off()

thresholds <- data.frame(geneName = colnames(data)[3:6],threshold = rep(0,4))

for (i in 1:4){
  percentThresh <- quantile(data[,c(i+2)],0.98)
  thresholds$threshold[i] <- percentThresh 
}

thresholdData <- applyThresholds(data, thresholds)

# EGFR-3 vs. C1S-4
crossTable<-table(thresholdData[,c(3,4)])
1/mosaic::oddsRatio(crossTable,conf.level = 0.95) #15.67

# C1S-4 vs. VCL-5
crossTable<-table(thresholdData[,c(4,5)])
1/mosaic::oddsRatio(crossTable,conf.level = 0.95) #10.4

# EGFR vs. VCL
crossTable<-table(thresholdData[,c(3,5)])
1/mosaic::oddsRatio(crossTable,conf.level = 0.95) #13.68

# EGFR vs GAPDH
crossTable<-table(thresholdData[,c(3,6)])
1/mosaic::oddsRatio(crossTable,conf.level = 0.95) #20

# VCL vs GAPDH
crossTable<-table(thresholdData[,c(5,6)])
1/mosaic::oddsRatio(crossTable,conf.level = 0.95) #50

# C1S vs GAPDH
crossTable<-table(thresholdData[,c(4,6)])
1/mosaic::oddsRatio(crossTable,conf.level = 0.95) #23

##
allORs <- data.frame(matrix(ncol=3,nrow=dim(thresholdData[,c(3:6)])[2]^2))
colnames(allORs)<-c('gene1','gene2','OR')
k=0
for (j in 1:(dim(thresholdData)[2]-2)){
  for (i in 1:(dim(thresholdData)[2]-2)){    
    allORs$gene1[i+4*k]<- as.character(colnames(thresholdData)[j+2])
    allORs$gene2[i+4*k] <- as.character(colnames(thresholdData)[i+2])
    crossTable<-table(thresholdData[,c(i+2,j+2)])
    if (dim(crossTable)[1]<2 | dim(crossTable)[2]<2) 
    {allORs$OR[i+4*k] <- NA}
    else {allORs$OR[i+4*k] <- 1/mosaic::oddsRatio(crossTable,conf.level = 0.95)}
  }
  k=k+1
}
temp <- allORs[is.finite(allORs$OR),] 
oddsRatioMat <- acast(temp, gene1~gene2)
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
  scale_fill_gradientn(colours = myPalette(10),
                       values=values,rescaler = function(x, ...) x, oob = identity)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title='log2(odds ratio) heatmap')
ggsave('graphs/RNAFISH/burnIn_oddsRatioHeatmap_WithLegend.pdf',height=8,width=8)


# Try barplot here
dataEGFR_burnIn <- filter(melted_ORmat,Var1=='EGFR')
ggplot(dataEGFR_burnIn, aes(x=Var2,y=log2(value)))+geom_histogram(stat='identity')+
  theme_classic()+
  ylab('log2(odds ratio) versus EGFR')

# Try dropping GAPDH low cells then ORs normalized by GAPDH?
data <-read.csv('dentistData/20160725_WM9noDrugA6G3.txt')
data.filtered <- filter(data, GAPDH>50)
data.norm <- data.frame(EGFR = data.filtered$EGFR/data.filtered$GAPDH,
                        C1S = data.filtered$C1S/data.filtered$GAPDH,
                        VCL= data.filtered$VCL/data.filtered$GAPDH)

thresholds <- data.frame(geneName = colnames(data.norm),threshold = rep(0,3))
## apply thresholding
for (i in 1:3){
  percentThresh <- quantile(data.norm[,i],0.98)
  thresholds$threshold[i] <- percentThresh 
}
thresholdData <- applyThresholds(data.norm, thresholds)

# EGFR vs. VCL
crossTable<-table(thresholdData[,c(1,3)])
VCL_OR <- 1/mosaic::oddsRatio(crossTable,conf.level = 0.95) #4.49 log2=2.1

# EGFR vs. C1S
crossTable<-table(thresholdData[,c(1,2)])
1/mosaic::oddsRatio(crossTable,conf.level = 0.95) #13 log2 = 3.7
C1S_OR <- 1/mosaic::oddsRatio(crossTable,conf.level = 0.95)

# VCL vs. C1S
crossTable<-table(thresholdData[,c(2,3)])
1/mosaic::oddsRatio(crossTable,conf.level = 0.95) #8


### Try same thing with resistance markees:
dataList <- load_WM9_noDrug_20150810(percentileForThresholds = 0.02)
data <- dataList[[3]]
data.filtered <- filter(data, GAPDH>50)
data.norm <- data.frame(EGFR = data.filtered$EGFR/data.filtered$GAPDH,
                        WNT5A = data.filtered$WNT5A/data.filtered$GAPDH,
                        PDGFRB= data.filtered$PDGFRB/data.filtered$GAPDH,
                        SERPINE1 = data.filtered$SERPINE1/data.filtered$GAPDH,
                        NGFR = data.filtered$NGFR/data.filtered$GAPDH,
                        AXL = data.filtered$AXL/data.filtered$GAPDH,
                        VEGFC = data.filtered$VEGFC/data.filtered$GAPDH)

thresholds <- data.frame(geneName = colnames(data.norm),threshold = rep(0,7))
## apply thresholding
for (i in 1:7){
  percentThresh <- quantile(data.norm[,i],0.98)
  thresholds$threshold[i] <- percentThresh 
}
thresholdData <- applyThresholds(data.norm, thresholds)

# EGFR vs. WNT5A
crossTable<-table(thresholdData[,c(1,2)])
WNT5A_OR <- 1/mosaic::oddsRatio(crossTable,conf.level = 0.95) #4.49 log2=2.1

crossTable<-table(thresholdData[,c(1,3)])
PDGFRB_OR <- 1/mosaic::oddsRatio(crossTable,conf.level = 0.95) 

crossTable<-table(thresholdData[,c(1,4)])
SERPINE1_OR <- 1/mosaic::oddsRatio(crossTable,conf.level = 0.95) 

crossTable<-table(thresholdData[,c(1,5)])
NGFR_OR <- 1/mosaic::oddsRatio(crossTable,conf.level = 0.95) 

crossTable<-table(thresholdData[,c(1,6)])
AXL_OR <- 1/mosaic::oddsRatio(crossTable,conf.level = 0.95) 

crossTable<-table(thresholdData[,c(1,7)])
VEGFC_OR <- 1/mosaic::oddsRatio(crossTable,conf.level = 0.95) 


dataGAPDHnormBurnIn <- data.frame(gene = c('1_VCL','1_C1S','WNT5A','PDGFRB','SERPINE1','NGFR','AXL','VEGFC'),
                                  OR = c(VCL_OR,C1S_OR,WNT5A_OR,PDGFRB_OR,SERPINE1_OR,NGFR_OR,AXL_OR, VEGFC_OR))

ggplot(dataGAPDHnormBurnIn,aes(x=gene,y=log2(OR)))+geom_histogram(stat='identity')


### 

source('plotScripts/RNAFISH/parameters_WM9_noDrug_20150810.R')
dataList <- load_WM9_noDrug_20150810(percentileForThresholds = 0.02)

data <- dataList[[3]]
geneThresholds <- dataList[[4]]

col4 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","#7FFF7F", 
                           "cyan", "#007FFF", "blue","#00007F"))
M <- cor(data[,c(4:22)])
pdf('graphs/RNAFISH/WM9_noDrug_20150810_correlationPlot.pdf',height=10,width=10)
corrplot(M, method="number",col=col4(20),sig.level=0.01)
dev.off()

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
    else {allORs$OR[i+19*k] <- 1/mosaic::oddsRatio(crossTable,conf.level = 0.95)}
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

noDrugEGFRdata <- filter(melted_ORmat,Var1=='EGFR')
noDrugEGFRdata$dataset <- 'old'

dataEGFR_burnIn$dataset <- 'burnIn'
allData <- rbind(noDrugEGFRdata,dataEGFR_burnIn)
allData <- allData[-which(allData$Var2=='EGFR'),]
ggplot(allData, aes(x=Var2,y=log2(value),fill=dataset,color=dataset))+
  geom_histogram(stat='identity',position='jitter')+
  theme_classic()+
  ylab('log2(odds ratio) versus EGFR')
#ggsave('~/Desktop/test.pdf',width=10,height=8)

## pull up sort seq
data1 <- read.csv('RNAseq/countTables/meltedData_20150722_EGFRsortRNArun1.tsv',sep='\t')
data2 <- read.csv('RNAseq/countTables/meltedData_20150831_EGFRsortRNArun2.tsv',sep='\t')
data <- rbind(data1,data2)

# we have to rename the sampleID factors to match our metadata table:
data$sampleID <- plyr::revalue(data$sampleID, c("SortNoDrugMix1"="nodrug_mix_02",
                                                "SortNoDrugMix2"="nodrug_mix_09",
                                                "SortNoDrugMix3"="nodrug_mix_10",
                                                "SortNoDrugEGFR1"="nodrug_EGFR_02",
                                                "SortNoDrugEGFR2"="nodrug_EGFR_09",
                                                "SortNoDrugEGFR3"="nodrug_EGFR_10",
                                                "Sort1WeekDrugMix1"="1week_mix_09",
                                                "week1-mix-2"="1week_mix_10",
                                                "Sort1WeekDrugEGFR1"="1week_EGFR_09",
                                                "week1-EGFRpos-2"="1week_EGFR_10",
                                                "week4-EGFRneg-1"="4weeks_neg_09",
                                                "week4-EGFRneg-2"="4weeks_neg_10",
                                                "week4-EGFRpos-1"="4weeks_EGFR_09",
                                                "week4-EGFRpos-2"="4weeks_EGFR_10",
                                                "week4-mix-1"="4weeks_mix_09",
                                                "week4-mix-2"="4weeks_mix_10",
                                                "week5-EGFRneg"="5weeks_neg",
                                                "week5-mix"="5weeks_mix"))

geneIDtoGeneName <- read.csv('RNAseq/annotations/hg19gene_idToGeneSymbol.tsv',
                             sep = "\t", header = TRUE)

data <- merge(data,geneIDtoGeneName,by='gene_id')

# load up metaData
metaData <- read.csv('RNAseq/metaData/EGFRsortRNAseqMetadata.csv')

dataAndMeta <- left_join(data,metaData,by='sampleID')

# convert counts to rpm
rmCounts <- data.table(dataAndMeta)
setkey(rmCounts, 'sampleID')
totalCountsInSample <- rmCounts[, list(McountsInSample = sum(counts)/10^6), by= 'sampleID']
rmCounts <- rmCounts[totalCountsInSample, ]
rmCounts$rpm <- rmCounts$counts / rmCounts$McountsInSample
dataAllGenes <- rmCounts


# remove genes with max(rpm) < 5
geneExpr <- rmCounts[, list(maxRPM = max(rpm), maxCounts = max(counts)), by = c('GeneSymbol')]
geneExpr <- transform(geneExpr, keepGene = maxRPM > 5)
rmCounts <- merge(rmCounts, subset(geneExpr, keepGene, select=c('GeneSymbol')),by='GeneSymbol')
cdata <- rmCounts
RNAFISHgenes <- filter(cdata, GeneSymbol=='EGFR'|GeneSymbol=='SOX10'|GeneSymbol=='CCNA2'|GeneSymbol=='GAPDH'|
                         GeneSymbol=='WNT5A'|GeneSymbol=='PDGFRB'|GeneSymbol=='SERPINE1'|GeneSymbol=='PDGFC'|
                         GeneSymbol=='VEGFC'|GeneSymbol=='NRG1'|GeneSymbol=='NGFR'|GeneSymbol=='FOSL1'|
                         GeneSymbol=='LOXL2'|GeneSymbol=='AXL'|GeneSymbol=='RUNX2'|GeneSymbol=='MITF'|
                         GeneSymbol=='FGFR1'|GeneSymbol=='JUN'|GeneSymbol=='VGF')

RNAFISHgenes <- filter(RNAFISHgenes, numTime==0)
castRNAFISHgenes <- dcast(RNAFISHgenes,sampleID~GeneSymbol)
castRNAFISHgenes$sampleID <- NULL
castRNAFISHgenes[c(1:3)]/castRNAFISHgenes[c(4:6)]
seq <- colMeans(castRNAFISHgenes[c(1:3)]/castRNAFISHgenes[c(4:6)])
seqMelt <- melt(seq)
seqMelt$GeneSymbol<- rownames(seqMelt)
colnames(seqMelt) <- c('seqValue','GeneSymbol')

# now merge with no drug RNA FISH data
RNAFISH <- noDrugEGFRdata[,c(2,3)]
colnames(RNAFISH) <- c('GeneSymbol','RNAFISHvalue')
plotdata <- left_join(RNAFISH,seqMelt, by='GeneSymbol')

ggplot(plotdata,aes(x=seqValue,y=RNAFISHvalue,label=GeneSymbol))+
  geom_point()+
  geom_text()+
  theme_classic()+
  ylab('OddsRatioWithEGFR by RNA FISH')+
  xlab('EGFR-high/EGFR-mix rpms by RNA-seq')
ggsave('graphs/RNAFISH/plotEGFRenrichmentSeqVsFISH.pdf')



ggplot(plotdata,aes(x=log2(seqValue),y=log2(RNAFISHvalue),label=GeneSymbol))+
  geom_point()+
  geom_text()+
  theme_classic()+
  ylab('log2(OddsRatioWithEGFR) by RNA FISH')+
  xlab('log2(EGFR-high/EGFR-mix rpms) by RNA-seq')
ggsave('graphs/RNAFISH/plotEGFRenrichmentSeqVsFISH_log2.pdf')


## looking at burn in genes
TXNRD1 <- filter(cdata, GeneSymbol=='TXNRD1' & numTime==0)

###
data <- read.csv('dentistData/20160727_WM9noDrug_BurnIn.txt')

# PCA first:
data.pca <- princomp(data.frame(data[,c(3:6)]), cor=TRUE, score=TRUE)
pdf('graphs/RNAFISH/WM9_noDrug_BurnIn_20160727_PCA.pdf',height=6,width=6)
biplot(data.pca, xlabs=rep(".", dim(data)[1]))
dev.off()
# for variance in each component
summary(data.pca)

source('plotScripts/RNAFISH/RNAFISH_functions.R')
thresholds <- data.frame(geneName = colnames(data[,3:6]),threshold = rep(0,4))
## apply thresholding
for (i in 1:4){
  percentThresh <- quantile(data[,i+2],0.98)
  thresholds$threshold[i] <- percentThresh 
}
thresholds$threshold[2] <- 11
thresholdData <- applyThresholds(data[,3:6], thresholds)


allORs <- data.frame(matrix(ncol=3,nrow=dim(thresholdData)[2]^2))
colnames(allORs)<-c('gene1','gene2','OR')
k=0
for (j in 1:(dim(thresholdData)[2])){
  for (i in 1:(dim(thresholdData)[2])){    
    allORs$gene1[i+4*k]<- as.character(colnames(thresholdData)[j])
    allORs$gene2[i+4*k] <- as.character(colnames(thresholdData)[i])
    crossTable<-table(thresholdData[,c(i,j)])
    if (dim(crossTable)[1]<2 | dim(crossTable)[2]<2) 
    {allORs$OR[i+4*k] <- NA}
    else {allORs$OR[i+4*k] <- 1/mosaic::oddsRatio(crossTable,conf.level = 0.95)}
  }
  k=k+1
}
temp <- allORs[is.finite(allORs$OR),] 
oddsRatioMat <- acast(temp, gene1~gene2)
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

#x = factor(x,levels(x)[c(4,5,1:3)])
melted_ORmat$Var1 <- factor(melted_ORmat$Var1,levels(melted_ORmat$Var1)[c(1,3,4,2)])
melted_ORmat$Var2 <- factor(melted_ORmat$Var2,levels(melted_ORmat$Var2)[c(1,3,4,2)])

ggplot(melted_ORmat,aes(x=Var1,y=Var2,fill=log2(valueCut)))+
  geom_tile(color='white')+ theme_classic()+
  scale_fill_gradientn(colours = myPalette(10),
                       values=values,rescaler = function(x, ...) x, oob = identity)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title='log2(odds ratio) heatmap')
ggsave('graphs/RNAFISH/WM9_noDrug_BurnIn_20160727_oddsRatioHeatmap_WithLegend.pdf',height=8,width=8)


col4 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","#7FFF7F", 
                           "cyan", "#007FFF", "blue","#00007F"))
M <- cor(data[,c(3:6)])
pdf('graphs/RNAFISH/WM9_noDrug_BurnIn_20160727_corrPlot.pdf',height=6,width=6)
corrplot(M, method="number",col=col4(20),sig.level=0.01)
dev.off()

# How many cells in this data set?
dim(data)[1]

