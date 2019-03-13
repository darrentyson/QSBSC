# New burn-in data

rm(list = ls())
library(ggplot2)

setwd('~/Dropbox/cancer-rna-fish/')

dataNoDrug <- read.csv(file = 'extractedData/dentistData/burnInGenes/20161023_WM9_noDrug_burnIns.txt')
dataNoDrug$condition <- 'noDrug'
# How many cells?
dim(dataNoDrug)

dataResistant <- read.csv(file = 'extractedData/dentistData/burnInGenes/20161026_WM9_ResistantB2_2_BurnIn_Pool1.txt')
dataResistant$condition <- 'resistant line'
#How many cells?
dim(dataResistant)

allData <- rbind(dataNoDrug, dataResistant)
allDataMelt <- melt(allData,id.vars = c('Xpos','Ypos','condition'))

ggplot(allDataMelt, aes(x=condition, y=value,color=condition,fill=condition))+
  geom_boxplot(color='black')+facet_wrap(~variable,scales='free')+theme_classic()
ggsave('graphs/RNAFISH/boxplotsNewBurnInData_APCDD1_MAP1B.pdf')


## replicates for new burn-in genes
dataNoDrug <- read.csv(file = 'extractedData/dentistData/burnInGenes/20161103_noDrug_pool1burnIn_replicates.txt')
dataNoDrug$condition <- 'noDrug'
dim(dataNoDrug)

dataResistant <- read.csv(file = 'extractedData/dentistData/burnInGenes/20161108_WM9_resistant_burnIn_replicates.txt')
dataResistant$condition <- 'resistant line'
dim(dataResistant)

allData <- rbind(dataNoDrug, dataResistant)
allDataMelt <- melt(allData,id.vars = c('Xpos','Ypos','condition'))

ggplot(allDataMelt, aes(x=condition, y=value,color=condition,fill=condition))+
  geom_boxplot(color='black')+facet_wrap(~variable,scales='free')+theme_classic()
ggsave('graphs/RNAFISH/boxplotsNewBurnInData_APCDD1_MAP1B_replicates.pdf')

# how many cells in this data set?
dim(allDataMelt)


# TXNRD1 in an untreated sample
dataNoDrug <- read.csv('extractedData/dentistData/20160727_WM9noDrug_BurnIn.txt')
dataNoDrug$condition <- 'noDrug'

# TXNRD1 in a sample with drug for 4 weeks:
dataResistant <- read.csv(file = 'extractedData/dentistData/burnInGenes/20161004_WM9_week4_TXNRD1.txt')
dataResistant$condition <- 'week_4_data'

allData <- rbind(dataNoDrug, dataResistant)
allDataMelt <- melt(allData,id.vars = c('Xpos','Ypos','condition'))

ggplot(allDataMelt, aes(x=condition, y=value,color=condition,fill=condition))+
  geom_boxplot(color='black')+facet_wrap(~variable,scales='free')+theme_classic()
ggsave('graphs/RNAFISH/boxplotsNewBurnInData_TXNRD1_EGFR_LOXL2.pdf')

# TXNRD1 in a purified resistant sample:
dataResistant <- read.csv(file = 'extractedData/dentistData/burnInGenes/20161109_WM9_resistant_TXNRD1.txt')
dataResistant$condition <- 'resistant_data'

allData <- rbind(dataNoDrug, dataResistant)
allDataMelt <- melt(allData,id.vars = c('Xpos','Ypos','condition'))

ggplot(allDataMelt, aes(x=condition, y=value,color=condition,fill=condition))+
  geom_boxplot(color='black')+facet_wrap(~variable,scales='free')+theme_classic()
ggsave('graphs/RNAFISH/boxplotsNewBurnInData_resistant_TXNRD1_EGFR_LOXL2.pdf')

#ggplot(filter(allDataMelt,variable=='TXNRD1'), aes(x=condition, y=value,color=condition,fill=condition))+
#  geom_boxplot(color='black')+theme_classic()
#ggsave('boxplotsNewBurnInData_resistant_TXNRD1only.pdf')


