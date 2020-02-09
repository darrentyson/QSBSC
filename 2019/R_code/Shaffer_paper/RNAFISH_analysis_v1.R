# RNA FISH figure plots
# this script will generate all the plots needed for the RNA FISH figure

# source the RNAFISH_functions that contains all the requied functions
source('RNAFISH/RNAFISH_functions.R')
loadLibraries()

# source parameter files for the data
source('RNAFISH/parameters_WM9_noDrug_20150618.R')
source('RNAFISH/parameters_WM9_noDrug_20150810.R')
source('RNAFISH/parameters_WM9_48hrs_20150624.R')
source('RNAFISH/parameters_WM9_4weeks_20150701.R')
source('RNAFISH/parameters_WM9_4weeks_20150831.R')
source('RNAFISH/parameters_primaryMelano_20150717.R')

# define your plot directory, this is where all your figures will be saved.
# The path should have a trailing slash to indicate it is a directory. 
# Please do not save these output files within the Git repo
plotDir <- path.expand("~/Desktop/Shaffer_figs/")
if(!file.exists(plotDir)) dir.create(plotDir)

#=============================================================================================================================
### Data = WM9_noDrug_20150618
#=============================================================================================================================
dataList <- load_WM9_noDrug_20150618(percentileForThresholds = 0.02)
plotName <- dataList[[1]]
plotDir <- dataList[[2]]
noDrugData <- dataList[[3]]
thresholds <- dataList[[4]]
gapdhThresholds <- dataList[[5]]
gapdhNormData <- gapdhNormalizeData(noDrugData, 50)
scattersForJackpotGenes(plotName = plotName, plotDir=plotDir,data = noDrugData)
twoGeneScatters(plotName = plotName, plotDir = plotDir, data = noDrugData, gene1='EGFR',gene2='WNT5A',geneThresholds = thresholds)
oddsRatioHeatmap(plotName = plotName, plotDir = plotDir, data = noDrugData, geneThresholds = thresholds)
numberJackpotsHistogram(plotName = plotName, plotDir = plotDir, data = noDrugData, geneThresholds = thresholds)

scattersForGAPDHnormData(plotName=paste0(plotName,'gapdhNorm_'),
                         plotDir, gapdhNormData)
twoGeneScatters(plotName = paste0(plotName,'gapdhNorm_'), plotDir = plotDir,
                data = gapdhNormData, gene1='EGFR',gene2='WNT5A',geneThresholds = gapdhThresholds,
                jitter=0)
oddsRatioHeatmap(plotName = paste0(plotName,'gapdhNorm_'),
                 plotDir=plotDir, data=gapdhNormData, geneThresholds = gapdhThresholds)

#=============================================================================================================================
### Data = WM9_noDrug_20150810
#=============================================================================================================================
dataList <- load_WM9_noDrug_20150810(percentileForThresholds = 0.02)
plotName <- dataList[[1]]
plotDir <- dataList[[2]]
noDrugData <- dataList[[3]]
thresholds <- dataList[[4]]
gapdhThresholds <- dataList[[5]]
gapdhNormData <- gapdhNormalizeData(noDrugData, 50)

#Number of cells in this data set needed for RNA FISH figure legend:
print(paste('Number of cells in WM9_noDrug_20150810:',eval(dim(noDrugData)[1])))

scattersForJackpotGenes(plotName = plotName, plotDir=plotDir,data = noDrugData)
twoGeneScatters(plotName = plotName, plotDir = plotDir, data = noDrugData,geneThresholds = thresholds)
twoGeneScatters(plotName = plotName, plotDir = plotDir, data = noDrugData,geneThresholds = thresholds, gene1='VEGFC',gene2='AXL')
oddsRatioHeatmap(plotName = plotName, plotDir = plotDir, data = noDrugData, geneThresholds = thresholds)
numberJackpotsHistogram(plotName = plotName, plotDir = plotDir, data = noDrugData, geneThresholds = thresholds)

#=============================================================================================================================
#=============================================================================================================================
source('RNAFISH/analysisWM983BCoExpressionAnalysis.R')

# In the text we list the number of cells with greater than 6 jackpot genes and greater than 8.
# Here are those calculations
thresholdData <- applyThresholds(data, geneThresholds)
thresholdData_markersOnly <- thresholdData[,c(1:4,8:16,17:21)]# Remove genes from this analysis that are not resistance markers - GAPDH, CCNA2, SOX10, MITF, VGF
numJackpotsPerCell <- data.frame(cellID = thresholdData$cellID,
                                 numberJackpotGenes = margin.table(as.matrix(thresholdData_markersOnly[4:dim(thresholdData_markersOnly)[2]]),1) )
print(paste('Number of cells with >= 6 jackpot genes:',length(which(numJackpotsPerCell$numberJackpotGenes>=6))))
print(paste('Number of cells with >= 8 jackpot genes:',length(which(numJackpotsPerCell$numberJackpotGenes>=8))))

scattersForGAPDHnormData(plotName = paste0(plotName,'gapdhNorm_'),plotDir, gapdhNormData)
print(paste('Number of cells once we remove cells with GAPDH < 50:', dim(gapdhNormData)[1]))
twoGeneScatters(plotName = paste0(plotName,'gapdhNorm_'), plotDir = plotDir,
                data = gapdhNormData, gene1='EGFR',gene2='WNT5A',geneThresholds = gapdhThresholds,
                jitter=0)
oddsRatioHeatmap(plotName = paste0(plotName,'gapdhNorm_'),
                 plotDir=plotDir, data=gapdhNormData, geneThresholds = gapdhThresholds)

# Contingency table shown in figure is based upon this data.
# Calculate values to show:
thresholdData <- data.frame(AXL = as.integer(noDrugData$AXL>thresholds[13,2]),
                            VEGFC = as.integer(noDrugData$VEGFC>thresholds[12,2]))
crossTable<-table(thresholdData)
OR <- 1/mosaic::oddsRatio(crossTable,conf.level = 0.95)

#=============================================================================================================================
#=============================================================================================================================
print('Contengency table information for callout with heatmap in RNA FISH figure:')
print(paste('Contengency table: '))
print(eval(crossTable))
print(paste('Odds ratio: ',eval(OR)))

# For figure -- need frequency of NGFR, AXL, and EGFR jackpots:
numEGFRjackpots <- sum(noDrugData$EGFR>thresholds[1,2])
numNGFRjackpots <- sum(noDrugData$NGFR>thresholds[9,2])
numAXLjackpots <- sum(noDrugData$AXL>thresholds[13,2])
totalNumCells <- dim(noDrugData)[1]

# print(paste('Total number of cells in the 20160810 no drug data:',totalNumCells))
# print(paste('Total number of EGFR jackpots in the 20160810 no drug data:',numEGFRjackpots))
# print(paste('Total number of NGFR jackpots in the 20160810 no drug data:',numNGFRjackpots))
# print((paste('Total number of AXL jackpots in the 20160810 no drug data:',numAXLjackpots)))
# print(paste('EGFR threshold in the 20160810 no drug data:',thresholds[1,2]))
# print(paste('NGFR threshold in the 20160810 no drug data:',thresholds[9,2]))
# print((paste('AXL threshold in the 20160810 no drug data:',thresholds[13,2])))
#

# ## Also for text
# print('What are the frequencies of jackpots for the genes in panel B of the RNA FISH figure:')
# print(paste('EGFR',sum(as.integer(noDrugData$EGFR>thresholds[1,2]))/dim(noDrugData)[1]))
# print(paste('WNT5A',sum(as.integer(noDrugData$WNT5A>thresholds[5,2]))/dim(noDrugData)[1])) #2
# print(paste('SEPRINE1',sum(as.integer(noDrugData$SERPINE1>thresholds[8,2]))/dim(noDrugData)[1])) #3
# print(paste('NGFR',sum(as.integer(noDrugData$NGFR>thresholds[9,2]))/dim(noDrugData)[1])) #4
# print(paste('NRG1',sum(as.integer(noDrugData$NRG1>thresholds[10,2]))/dim(noDrugData)[1])) #5
# print(paste('VEGFC',sum(as.integer(noDrugData$VEGFC>thresholds[12,2]))/dim(noDrugData)[1])) #6
# print(paste('JUN',sum(as.integer(noDrugData$JUN>thresholds[18,2]))/dim(noDrugData)[1])) #7
# print(paste('LOXL2',sum(as.integer(noDrugData$LOXL2>thresholds[15,2]))/dim(noDrugData)[1])) #8
# print(paste('AXL',sum(as.integer(noDrugData$AXL>thresholds[13,2]))/dim(noDrugData)[1])) #9


# ## does ORs plot change much with changes in threshold??
# dataList <- load_WM9_noDrug_20150810(percentileForThresholds = 0.02)
# plotName <- dataList[[1]]
# plotDir <- dataList[[2]]
# noDrugData <- dataList[[3]]
# thresholds <- dataList[[4]]
# gapdhThresholds <- dataList[[5]]
# gapdhNormData <- gapdhNormalizeData(noDrugData, 50)
# plotName <- paste0(plotName, '02thresh_')
# oddsRatioHeatmap(plotName = plotName, plotDir = plotDir, data = noDrugData, geneThresholds = thresholds)
#
# plotName_half <- paste0(plotName, 'halfThresh_')
# thresholdsHalf <- thresholds
# thresholdsHalf[,2] <- thresholds$threshold/2
# oddsRatioHeatmap(plotName = plotName_half, plotDir = plotDir, data = noDrugData, geneThresholds = thresholdsHalf)
#
# plotName_2X <- paste0(plotName, '2XThresh_')
# thresholds2X <- thresholds
# thresholds2X[,2] <- thresholds$threshold*2
# oddsRatioHeatmap(plotName = plotName_2X, plotDir = plotDir, data = noDrugData, geneThresholds = thresholds2X)
#
# ### Data = WM9_48hrs_20150624
# dataList <- load_WM9_48hrs_20150624(percentileForThresholds = 0.02)
# plotName <- dataList[[1]]
# plotDir <- dataList[[2]]
# noDrugData <- dataList[[3]]
# thresholds <- dataList[[4]]
#

scattersForJackpotGenes(plotName = plotName, plotDir=plotDir,data = noDrugData)
twoGeneScatters(plotName = plotName, plotDir = plotDir, data = noDrugData,geneThresholds = thresholds)
# this makes the figure 3a
twoGeneScatters(plotName = plotName, plotDir = plotDir, data = noDrugData,geneThresholds = thresholds, gene1='AXL',gene2='VEGFC')
oddsRatioHeatmap(plotName = plotName, plotDir = plotDir, data = noDrugData, geneThresholds = thresholds)
numberJackpotsHistogram(plotName = plotName, plotDir = plotDir, data = noDrugData, geneThresholds = thresholds)

#=============================================================================================================================
#=============================================================================================================================

