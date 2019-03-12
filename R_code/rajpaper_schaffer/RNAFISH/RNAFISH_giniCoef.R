# calculate Gini coef

# run all of this from your cancer-rna-fish repo
setwd('~/Dropbox/cancer-rna-fish/')

rm(list = ls())

source('plotScripts/RNAFISH/RNAFISH_functions.R')
loadLibraries()
library('ineq')

# source parameter files for the data
source('plotScripts/RNAFISH/parameters_WM9_noDrug_20150618.R')
source('plotScripts/RNAFISH/parameters_WM9_noDrug_20150810.R')
source('plotScripts/RNAFISH/parameters_primaryMelano_20150717.R')


### Data = WM9_noDrug_20150618
dataList <- load_WM9_noDrug_20150618(percentileForThresholds = 0.02)
plotName1 <- dataList[[1]]
plotDir1 <- dataList[[2]]
noDrugData1 <- dataList[[3]]

Gini_20150618 <- apply(noDrugData1[,4:22],2,function(x) ineq(x,type='Gini'))

### Data = WM9_noDrug_20150810
dataList <- load_WM9_noDrug_20150810(percentileForThresholds = 0.02)
plotName2 <- dataList[[1]]
plotDir2 <- dataList[[2]]
noDrugData2 <- dataList[[3]]

Gini_20150810 <- apply(noDrugData2[,4:22],2,function(x) ineq(x,type='Gini'))

### We need a good way to display this information.
Gini <- data.frame(genes = names(Gini_20150618), Gini_set1 = Gini_20150618, Gini_set2 = Gini_20150810)
Gini2 <- Gini[c(1:6,8:10,12:15,18),] # Here we remove PDGFC, RUNX2, and VGF because they don't have any expression.

ggplot(Gini2, aes(x=Gini_set1,label=genes))+
  geom_point(aes(y=1,color=genes))+
  geom_text(y=1,angle=55,hjust=-0.05,vjust=-0.2,aes(color=genes),size=6)+
  theme_classic()+xlim(0,1)+ylim(0.9,1.1)
#ggsave(paste0(plotDir1,plotName1,'Gini.pdf'),height=5,width=10)

ggplot(Gini2, aes(x=Gini_set2,label=genes))+
  geom_point(aes(y=1,color=genes))+
  geom_text(y=1,angle=55,hjust=-0.05,vjust=-0.2,aes(color=genes),size=6)+
  theme_classic()+xlim(0,1)+ylim(0.9,1.1)
#ggsave(paste0(plotDir2,plotName2,'Gini.pdf'),height=5,width=10)

ggplot(Gini2, aes(x=Gini_set1,y=Gini_set2,label=genes))+
  geom_point(aes(color=genes))+geom_text(aes(color=genes),hjust=-0.05,vjust=-0.2)+
  theme_classic()+ylab('Gini_WM9_noDrug_20150810')+xlab('Gini_WM9_noDrug_20150618')+xlim(0,1.1)+ylim(0,1.1)+
  geom_abline(intercept=0,slope=1,color='gray')
ggsave(paste0('~/Dropbox/cancer-rna-fish/graphs/RNAFISH/giniReplicates/WM9_noDrug.pdf'),height=6,width=7)

## control data
controlData1 <- read.csv('extractedData/dentistData/20160502_WM9_noDrug_controlGenes.txt')
controlData2 <- read.csv('extractedData/dentistData/20160503_WM9_noDrug_controlGenes.txt')

control1_gini <- apply(controlData1[,3:6],2,function(x) ineq(x,type='Gini'))
control2_gini <- apply(controlData2[,3:5],2,function(x) ineq(x,type='Gini'))
controls <- rbind(cbind(melt(control1_gini),melt(control1_gini)),cbind(melt(control2_gini),melt(control2_gini)))
controls$genes <- row.names(controls)
controls <- controls[,c(3,1,2)]
colnames(controls) <- c('genes','Gini_set1','Gini_set2')

Gini3 <- rbind(Gini2,controls)
ggplot(Gini3, aes(x=Gini_set2,label=genes))+
  geom_point(aes(y=1,color=genes))+
  geom_text(y=1,angle=55,hjust=-0.05,vjust=-0.2,aes(color=genes),size=4)+
  theme_classic()+ylim(0.9,1.1)+xlim(0,1)
ggsave(paste0('~/Dropbox/cancer-rna-fish/graphs/RNAFISH/Gini_ExtraControlGenes.pdf'),height=3,width=10)

## Add in the countries
countries <- read.csv('extractedData/giniData/rawdata_2172.csv', header=FALSE)
countriesFinal <- countries[c(2,8,10,15,20,26,43,66,98,119,126,140),]

ggplot(countriesFinal, aes(x=V3*0.01,label=V2))+
  geom_point(aes(y=1),color='blue')+
  geom_text(y=1,angle=55,hjust=1.2,vjust=0,size=3)+
  theme_classic()+ylim(0.9,1.1)+xlim(0.15,0.9)
ggsave(paste0('~/Dropbox/cancer-rna-fish/graphs/RNAFISH/Gini_countries.pdf'),height=3,width=10)

ggplot(countriesFinal, aes(x=V3*0.01))+
  geom_point(data=Gini3,aes(x=Gini_set1,y=1),color='red')+
  geom_text(data=Gini3,aes(x=Gini_set1,y=1,label=genes),angle=55,hjust=-0.05,vjust=-0.2,size=3)+
  theme_classic()+ylim(0.9,1.1)+xlim(0.15,0.9)
ggsave(paste0('~/Dropbox/cancer-rna-fish/graphs/RNAFISH/Gini_genesFormatedLikeCountries.pdf'),height=3,width=10)

ggplot(countriesFinal, aes(x=V3*0.01))+
  geom_point(aes(y=1),color='blue')+
  geom_text(aes(label=V2,y=1),angle=55,hjust=1.2,vjust=0,size=3)+
  geom_point(data=Gini3,aes(x=Gini_set1,y=1),color='red')+
  geom_text(data=Gini3,aes(x=Gini_set1,y=1,label=genes),angle=55,hjust=0.05,vjust=-0.2,size=3)+
  theme_classic()+ylim(0.9,1.1)+xlim(0.15,0.9)
ggsave(paste0('~/Dropbox/cancer-rna-fish/graphs/RNAFISH/Gini_countriesAndGenes.pdf'),height=3,width=10)


## let's put primary melanocytes on the same plot as WM9 data.
dataList <- load_primaryMelano_20150717()
primaryMelanoData <- dataList[[3]]
# This removes - EGFR, WNT5A, PDGFRB, PDGFC, NRG1, FOSL1,LOXL2, RUNX2,JUN,VGF,VEGFC
primaryMelanoGini <- apply(primaryMelanoData[,c(5:7,11,12,16,17,20)],2,function(x) ineq(x,type='Gini'))
primaryMelanoGini <- melt(primaryMelanoGini)
primaryMelanoGini$genes <- rownames(primaryMelanoGini)
primaryMelanoGini$dataset <- 'primary melanocytes'

WM9_Gini <- apply(noDrugData2[,c(4:9,11:13,15:18,20:21)],2,function(x) ineq(x,type='Gini')) 
# We remove the genes that are not really expressed in this data set -- PDGFC, FOSL1, RUNX2, VGF
WM9_Gini  <- melt(WM9_Gini)
WM9_Gini$genes <- rownames(WM9_Gini)
WM9_Gini  <- WM9_Gini[,c(2,1)]

colnames(controls) <- c('genes','value','value2')
WM9_Gini  <- rbind(WM9_Gini ,controls[,c(1:2)])
WM9_Gini$dataset <- 'WM9_20150810' 

## 983b data
source('plotScripts/RNAFISH/parameters_WM983b_noDrug_20150525.R')
dataList <- load_WM983b_noDrug_20150525(percentileForThresholds = 0.02)
WM983bData <- dataList[[3]]
# Only want the genes that are expressed. -- removing EGFR,PDGFRB,PDGFC,WNT5A,JUN,VEGFC  
WM983b_Gini <- apply(WM983bData[,c(5:8,11:13,16:22)],2,function(x) ineq(x,type='Gini')) 
WM983b_Gini <- melt(WM983b_Gini)
WM983b_Gini$genes <- rownames(WM983b_Gini)
WM983b_Gini$dataset <- 'WM983b_20150525'

## 1205Lu data
Lu1205_hyb1 <- read.csv('extractedData/dentistData/20160507_1205Lu_hyb1.txt')
Lu1205_hyb3 <- read.csv('extractedData/dentistData/20160507_1205Lu_hyb3.txt')
Lu1205_hyb3 <- Lu1205_hyb3[,c(1:5)]
Lu_1 <- apply(Lu1205_hyb1[,3:6],2,function(x) ineq(x,type='Gini'))
Lu_3 <- apply(Lu1205_hyb3[,3:5],2,function(x) ineq(x,type='Gini'))
Lu1205_Gini <- rbind(melt(Lu_1),melt(Lu_3))
Lu1205_Gini$genes <- rownames(Lu1205_Gini)
Lu1205_Gini$dataset <- '1205Lu_20160507'

# replicates 1205Lu
Lu1205_hyb1 <- read.csv('extractedData/dentistData/20160611_1205Lu_hyb1.txt')
Lu1205_hyb3 <- read.csv('extractedData/dentistData/20160611_1205Lu_hyb3.txt')
Lu1205_hyb3 <- Lu1205_hyb3[,c(1:5)]
Lu_1 <- apply(Lu1205_hyb1[,3:6],2,function(x) ineq(x,type='Gini'))
Lu_3 <- apply(Lu1205_hyb3[,3:5],2,function(x) ineq(x,type='Gini'))
Lu1205_Gini2 <- rbind(melt(Lu_1),melt(Lu_3))
Lu1205_Gini2$genes2 <- rownames(Lu1205_Gini2)
Lu1205_Gini2$dataset <- '1205Lu_20160611'
colnames(Lu1205_Gini2) <- c('value2','genes2','dataset2')

LuReps <- cbind(Lu1205_Gini,Lu1205_Gini2)
ggplot(LuReps, aes(x=value,y=value2,label=genes))+
  geom_point(aes(color=genes))+geom_text(aes(color=genes),hjust=-0.05,vjust=-0.2)+
  theme_classic()+ylab(LuReps$dataset2[1])+xlab(LuReps$dataset[1])+xlim(0,1.1)+ylim(0,1.1)+
  geom_abline(intercept=0,slope=1,color='gray')
ggsave('graphs/RNAFISH/giniReplicates/Lu1205.pdf',height=6,width=7)

#SK-MEL-28 data
SKMEL_hyb1 <- read.csv('extractedData/dentistData/20160513_SKMEL28_hyb1.txt')
SKMEL_hyb3 <- read.csv('extractedData/dentistData/20160513_SKMEL28_hyb3.txt')
SKMEL_1 <- apply(SKMEL_hyb1[,4:6],2,function(x) ineq(x,type='Gini'))
SKMEL_3 <- apply(SKMEL_hyb3[,3:5],2,function(x) ineq(x,type='Gini'))
SKMEL_Gini <- rbind(melt(SKMEL_1),melt(SKMEL_3))
SKMEL_Gini$genes <- rownames(SKMEL_Gini)
SKMEL_Gini$dataset <- 'SKMEL28_20160513'

# replicates SKMEL28
SKMEL_hyb1 <- read.csv('extractedData/dentistData/20160614_SKMEL28_hyb1.txt')
SKMEL_hyb3 <- read.csv('extractedData/dentistData/20160614_SKMEL28_hyb3.txt')
SKMEL_hyb3 <- SKMEL_hyb3[,c(1:5)]
SKMEL_1 <- apply(SKMEL_hyb1[,4:6],2,function(x) ineq(x,type='Gini'))
SKMEL_3 <- apply(SKMEL_hyb3[,3:5],2,function(x) ineq(x,type='Gini'))
SKMEL_Gini2 <- rbind(melt(SKMEL_1),melt(SKMEL_3))
SKMEL_Gini2$genes <- rownames(SKMEL_Gini2)
SKMEL_Gini2$dataset <- 'SKMEL28_20160614'
colnames(SKMEL_Gini2) <- c('value2','genes2','dataset2')

SKMELReps <- cbind(SKMEL_Gini,SKMEL_Gini2)
ggplot(SKMELReps, aes(x=value,y=value2,label=genes))+
  geom_point(aes(color=genes))+geom_text(aes(color=genes),hjust=-0.05,vjust=-0.2)+
  theme_classic()+ylab(SKMELReps$dataset2[1])+xlab(SKMELReps$dataset[1])+xlim(0,1.1)+ylim(0,1.1)+
  geom_abline(intercept=0,slope=1,color='gray')
ggsave('graphs/RNAFISH/giniReplicates/SKMEL28.pdf',height=6,width=7)


## PC9 data
PC9_pool1_1 <- read.csv('extractedData/dentistData/20160522_PC9_pool1.txt') 
PC9data <- PC9_pool1_1
PC9_pool2_1 <- read.csv('extractedData/dentistData/20160522_PC9_pool2.txt') 
PC9_pool1_2 <- read.csv('extractedData/dentistData/20160602_PC9_pool1.txt') 
PC9_pool2_2 <- read.csv('extractedData/dentistData/20160602_PC9_pool2.txt') 
PC9_pool1_1 <- apply(PC9_pool1_1[,3:6],2,function(x) ineq(x,type='Gini'))
PC9_pool2_1 <- apply(PC9_pool2_1[,4:5],2,function(x) ineq(x,type='Gini'))
PC9_pool1_2 <- apply(PC9_pool1_2[,3:6],2,function(x) ineq(x,type='Gini'))
PC9_pool2_2 <- apply(PC9_pool2_2[,4:5],2,function(x) ineq(x,type='Gini'))
PC9_1 <- rbind(melt(PC9_pool1_1),melt(PC9_pool2_1))
PC9_2 <- rbind(melt(PC9_pool1_2),melt(PC9_pool2_2))
PC9_1$genes <- rownames(PC9_1)
PC9_1$dataset <- 'PC9_20160522'
PC9_2$genes <- rownames(PC9_2)
PC9_2$dataset <- 'PC9_20160602'
colnames(PC9_2) <- c('value2','genes2','dataset2')
PC9Reps <- cbind(PC9_1,PC9_2)
ggplot(PC9Reps, aes(x=value,y=value2,label=genes))+
  geom_point(aes(color=genes))+geom_text(aes(color=genes),hjust=-0.05,vjust=-0.2)+
  theme_classic()+ylab(PC9Reps$dataset2[1])+xlab(PC9Reps$dataset[1])+xlim(0,1.1)+ylim(0,1.1)+
  geom_abline(intercept=0,slope=1,color='gray')
ggsave('graphs/RNAFISH/giniReplicates/PC9.pdf',height=6,width=7)


## MDA-MB-231
MDAMB231_pool1 <- read.csv('extractedData/dentistData/20160517_MDAMB231_pool1.txt')
MDAMB231_pool2 <- read.csv('extractedData/dentistData/20160517_MDAMB231_pool2.txt')
MDAgini_1 <- apply(MDAMB231_pool1[,3:6],2,function(x) ineq(x,type='Gini'))
MDAgini_2 <- apply(MDAMB231_pool2[,3:6],2,function(x) ineq(x,type='Gini'))
MDAgini <- rbind(melt(MDAgini_1),melt(MDAgini_2))
MDAgini$genes <- rownames(MDAgini)
MDAgini$dataset <- 'MDAMB231_20160517'

# replicates MDAMB231
MDAMB231_pool1 <- read.csv('extractedData/dentistData/20160606_MDAMB231_pool1.txt')
MDAMB231_pool2 <- read.csv('extractedData/dentistData/20160606_MDAMB231_pool2.txt')
MDAgini_1 <- apply(MDAMB231_pool1[,3:6],2,function(x) ineq(x,type='Gini'))
MDAgini_2 <- apply(MDAMB231_pool2[,3:6],2,function(x) ineq(x,type='Gini'))
MDAgini2 <- rbind(melt(MDAgini_1),melt(MDAgini_2))
MDAgini2$genes <- rownames(MDAgini2)
MDAgini2$dataset <- 'SKMEL28_20160606'
colnames(MDAgini2) <- c('value2','genes2','dataset2')

MDAReps <- cbind(MDAgini,MDAgini2)
ggplot(MDAReps, aes(x=value,y=value2,label=genes))+
  geom_point(aes(color=genes))+geom_text(aes(color=genes),hjust=-0.05,vjust=-0.2)+
  theme_classic()+ylab(MDAReps$dataset2[1])+xlab(MDAReps$dataset[1])+xlim(0,1.1)+ylim(0,1.1)+
  geom_abline(intercept=0,slope=1,color='gray')
ggsave('graphs/RNAFISH/giniReplicates/MDAMB231.pdf',height=6,width=7)

#Hela cells
Hela_pool1 <- read.csv('extractedData/dentistData/20160527_Hela_pool1.txt')
Hela_pool2 <- read.csv('extractedData/dentistData/20160527_Hela_pool2.txt')
Hela_1 <- apply(Hela_pool1[,3:6],2,function(x) ineq(x,type='Gini'))
Hela_2 <- apply(Hela_pool2[,3:6],2,function(x) ineq(x,type='Gini'))
Helagini <- rbind(melt(Hela_1),melt(Hela_2))
Helagini$genes <- rownames(Helagini)
Helagini$dataset <- 'Hela_20160527'

# replicates HeLa
Hela_pool1 <- read.csv('extractedData/dentistData/20160712_Hela_pool1.txt')
Hela_pool2 <- read.csv('extractedData/dentistData/20160712_Hela_pool2.txt')
Hela_1 <- apply(Hela_pool1[,3:6],2,function(x) ineq(x,type='Gini'))
Hela_2 <- apply(Hela_pool2[,3:6],2,function(x) ineq(x,type='Gini'))
Helagini2 <- rbind(melt(Hela_1),melt(Hela_2))
Helagini2$genes <- rownames(Helagini2)
Helagini2$dataset <- 'Hela_20160712'
colnames(Helagini2) <- c('value2','genes2','dataset2')

HelaReps <- cbind(Helagini,Helagini2)
ggplot(HelaReps, aes(x=value,y=value2,label=genes))+
  geom_point(aes(color=genes))+geom_text(aes(color=genes),hjust=-0.05,vjust=-0.2)+
  theme_classic()+ylab(HelaReps$dataset2[1])+xlab(HelaReps$dataset[1])+xlim(0,1.1)+ylim(0,1.1)+
  geom_abline(intercept=0,slope=1,color='gray')
ggsave('graphs/RNAFISH/giniReplicates/HeLa.pdf',height=6,width=7)

#SH-SY-5Y cells
SHSY5Y_pool1 <- read.csv('extractedData/dentistData/20160524_SH_SY_5Y_pool1.txt')
SHSY5Y_pool2 <- read.csv('extractedData/dentistData/20160524_SH_SY_5Y_pool2.txt')
SHSY5Ygini_1 <- apply(SHSY5Y_pool1[,c(3:4,6)],2,function(x) ineq(x,type='Gini'))
SHSY5Ygini_2 <- apply(SHSY5Y_pool2[,3:6],2,function(x) ineq(x,type='Gini'))
SHSY5Ygini <- rbind(melt(SHSY5Ygini_1),melt(SHSY5Ygini_2))
SHSY5Ygini$genes <- rownames(SHSY5Ygini)
SHSY5Ygini$dataset <- 'SHSY5Y_20160524'

# replicates SHSY5Y
SHSY5Y_pool1 <- read.csv('extractedData/dentistData/20160616_SHSY5Y_pool1.txt')
SHSY5Y_pool2 <- read.csv('extractedData/dentistData/20160616_SHSY5Y_pool2.txt')
SHSY5Ygini_1 <- apply(SHSY5Y_pool1[,c(3:4,6)],2,function(x) ineq(x,type='Gini'))
SHSY5Ygini_2 <- apply(SHSY5Y_pool2[,3:6],2,function(x) ineq(x,type='Gini'))
SHSY5Ygini2 <- rbind(melt(SHSY5Ygini_1),melt(SHSY5Ygini_2))
SHSY5Ygini2$genes <- rownames(SHSY5Ygini2)
SHSY5Ygini2$dataset <- 'SKMEL28_20160606'
colnames(SHSY5Ygini2) <- c('value2','genes2','dataset2')

SHSY5YReps <- cbind(SHSY5Ygini,SHSY5Ygini2)
ggplot(SHSY5YReps, aes(x=value,y=value2,label=genes))+
  geom_point(aes(color=genes))+geom_text(aes(color=genes),hjust=-0.05,vjust=-0.2)+
  theme_classic()+ylab(SHSY5YReps$dataset2[1])+xlab(SHSY5YReps$dataset[1])+xlim(0,1.1)+ylim(0,1.1)+
  geom_abline(intercept=0,slope=1,color='gray')
ggsave('graphs/RNAFISH/giniReplicates/SHSY5Y.pdf',height=6,width=7)

#apply(noDrugData1,2,function(x) max(x))
data <- rbind(WM9_Gini,primaryMelanoGini,WM983b_Gini,Lu1205_Gini,SKMEL_Gini, PC9_1, MDAgini, SHSY5Ygini,Helagini)
data$dataset <- as.factor(data$dataset)
data$dataset = factor(data$dataset,levels(data$dataset)[c(2,3,4,6,5,1,7,9,8)])

ggplot(data, aes(x=value,y=as.factor(dataset), label=genes))+
  geom_point(aes(color=genes))+
  geom_text(angle=55,hjust=-0.05,vjust=-0.2,aes(color=genes),size=4)+
  theme_classic()+xlim(0,1)
ggsave('graphs/RNAFISH/Gini_panel_forFigure.pdf',width=11,height=9)



# Number of cells in each data set

#wm9 no drug
dim(noDrugData2)

#wm983b data
dim(WM983bData)

#SKMEL28
dim(SKMEL_hyb1)

#1205Lu
dim(Lu1205_hyb1)

# primary melanocytes
dim(primaryMelanoData)

# SH-SY5Y
dim(SHSY5Y_pool1)

# PC9
dim(PC9data)

#MD-MB-231
dim(MDAMB231_pool1)

# HeLa
dim(Hela_pool1)
