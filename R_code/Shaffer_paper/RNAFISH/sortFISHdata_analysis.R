
# run all of this from your cancer-rna-fish repo
setwd('~/Dropbox/cancer-rna-fish/')

rm(list = ls())

source('plotScripts/RNAFISH/RNAFISH_functions.R')
source('plotScripts/RNAFISH/parameters_WM9_noDrug_20150810.R')
loadLibraries()
WM9_noDrug_list <- load_WM9_noDrug_20150810(0.02)
thresholds_WM9noDrug <- data.frame(WM9_noDrug_list[4])
thresholds_WM9noDrug[17,2] <- 80 # Change FGFR1 because 24 is really low considering how much these cells have
thresholds_WM9noDrug[13,2] <- 30

# source parameter files for the data
source('plotScripts/RNAFISH/parameters_WM9_sortEGFRwell_20150725.R')
source('plotScripts/RNAFISH/parameters_WM9_sortMixwell_20150725.R')

## Data = WM9_sortEGFRwell_20150725
dataList <- load_WM9_sortEGFRwell_20150725()
plotName <- dataList[[1]]
plotDir <- 'graphs/suppplementGraphs/'
data_EGFR <- dataList[[3]]
thresholds <- dataList[[4]]

numberJackpotsHistogram(plotName = plotName,plotDir = plotDir,data = data_EGFR, geneThresholds = thresholds)

## Data = WM9_sortEGFRwell_20150725
dataList <- load_WM9_sortMixwell_20150725()
plotName <- dataList[[1]]
plotDir <- 'graphs/suppplementGraphs/'
data_mix <- dataList[[3]]
thresholds <- dataList[[4]]

numberJackpotsHistogram(plotName = plotName,plotDir = plotDir,data = data_mix, geneThresholds = thresholds)

data_EGFR_cast <- melt(data_EGFR, id = c('cellID','Xpos','Ypos'))
data_EGFR_cast$dataset <- 'EGFRhigh'
data_mix_cast <- melt(data_mix, id = c('cellID','Xpos','Ypos'))
data_mix_cast$dataset <- 'mix'
data <- rbind(data_EGFR_cast, data_mix_cast)
data$jitter <- jitter(data$value,amount=1)
data$jitter[which(data$jitter<0)] <- 0

#####
gene='EGFR'
percentEGFR <- 100*dim(filter(data,variable==eval(gene) & dataset =='EGFRhigh' & value>thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName==eval(gene))]))[1]/dim(filter(data,variable==eval(gene) & dataset =='EGFRhigh'))[1]
percentMix <- 100*dim(filter(data,variable==eval(gene) & dataset =='mix' & value>thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName==eval(gene))]))[1]/dim(filter(data,variable==eval(gene) & dataset =='mix'))[1]
print(paste('Percentage', eval(gene), 'jackpots in EGFR high ', percentEGFR))
print(paste('Percentage', eval(gene), 'jackpots in mix ',percentMix))

ggplot(filter(data,variable=='EGFR'),aes(x=jitter,fill=dataset))+
  geom_histogram(color='black')+facet_wrap(~dataset,ncol=1)+
  geom_rug()+theme_classic()+
  annotate("text",x=c(30,75),y=500,label=c(percentEGFR,percentMix))+
  geom_vline(xintercept=thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName=='EGFR')])
ggsave(filename = 'graphs/suppplementGraphs/RNAFISH_sortEGFR.pdf', height = 9, width =7)

# fold change here:
temp <- filter(data,variable=='EGFR')
high <- filter(temp, dataset=='EGFRhigh')
mix <- filter(temp, dataset=='mix')
mean(high$value)/mean(mix$value)

#####
gene='WNT5A'
percentEGFR <- 100*dim(filter(data,variable==eval(gene) & dataset =='EGFRhigh' & value>thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName==eval(gene))]))[1]/dim(filter(data,variable==eval(gene) & dataset =='EGFRhigh'))[1]
percentMix <- 100*dim(filter(data,variable==eval(gene) & dataset =='mix' & value>thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName==eval(gene))]))[1]/dim(filter(data,variable==eval(gene) & dataset =='mix'))[1]
print(paste('Percentage', eval(gene), 'jackpots in EGFR high ', percentEGFR))
print(paste('Percentage', eval(gene), 'jackpots in mix ',percentMix))

ggplot(filter(data,variable=='WNT5A'),aes(x=jitter,fill=dataset))+
  geom_histogram(color='black')+facet_wrap(~dataset,ncol=1)+
  geom_rug()+theme_classic()+
  annotate("text",x=c(50,100),y=500,label=c(percentEGFR,percentMix))+
  geom_vline(xintercept=thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName=='WNT5A')])
ggsave(filename = 'graphs/suppplementGraphs/RNAFISH_sortWNT5A.pdf', height = 9, width =7)

#####
gene='SERPINE1'
percentEGFR <- 100*dim(filter(data,variable==eval(gene) & dataset =='EGFRhigh' & value>thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName==eval(gene))]))[1]/dim(filter(data,variable==eval(gene) & dataset =='EGFRhigh'))[1]
percentMix <- 100*dim(filter(data,variable==eval(gene) & dataset =='mix' & value>thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName==eval(gene))]))[1]/dim(filter(data,variable==eval(gene) & dataset =='mix'))[1]
print(paste('Percentage', eval(gene), 'jackpots in EGFR high ', percentEGFR))
print(paste('Percentage', eval(gene), 'jackpots in mix ',percentMix))

ggplot(filter(data,variable=='SERPINE1'),aes(x=jitter,fill=dataset))+
  geom_histogram(color='black')+facet_wrap(~dataset,ncol=1)+
  geom_rug()+theme_classic()+
  annotate("text",x=c(50,110),y=500,label=c(percentEGFR,percentMix))+
  geom_vline(xintercept=thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName=='SERPINE1')])
ggsave(filename = 'graphs/suppplementGraphs/RNAFISH_sortSERPINE1.pdf', height = 9, width =7)

######
gene='PDGFRB'
percentEGFR <- 100*dim(filter(data,variable==eval(gene) & dataset =='EGFRhigh' & value>thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName==eval(gene))]))[1]/dim(filter(data,variable==eval(gene) & dataset =='EGFRhigh'))[1]
percentMix <- 100*dim(filter(data,variable==eval(gene) & dataset =='mix' & value>thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName==eval(gene))]))[1]/dim(filter(data,variable==eval(gene) & dataset =='mix'))[1]
print(paste('Percentage', eval(gene), 'jackpots in EGFR high ', percentEGFR))
print(paste('Percentage', eval(gene), 'jackpots in mix ',percentMix))

ggplot(filter(data,variable=='PDGFRB'),aes(x=jitter,fill=dataset))+
  geom_histogram(color='black')+facet_wrap(~dataset,ncol=1)+
  geom_rug()+theme_classic()+
  annotate("text",x=c(15,25),y=500,label=c(percentEGFR,percentMix))+
  geom_vline(xintercept=thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName=='PDGFRB')])
ggsave(filename = 'graphs/suppplementGraphs/RNAFISH_sortPDGFRB.pdf', height = 9, width =7)

####
gene='NRG1'
percentEGFR <- 100*dim(filter(data,variable==eval(gene) & dataset =='EGFRhigh' & value>thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName==eval(gene))]))[1]/dim(filter(data,variable==eval(gene) & dataset =='EGFRhigh'))[1]
percentMix <- 100*dim(filter(data,variable==eval(gene) & dataset =='mix' & value>thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName==eval(gene))]))[1]/dim(filter(data,variable==eval(gene) & dataset =='mix'))[1]
print(paste('Percentage', eval(gene), 'jackpots in EGFR high ', percentEGFR))
print(paste('Percentage', eval(gene), 'jackpots in mix ',percentMix))

ggplot(filter(data,variable=='NRG1'),aes(x=jitter,fill=dataset))+
  geom_histogram(color='black')+facet_wrap(~dataset,ncol=1)+
  geom_rug()+theme_classic()+
  annotate("text",x=c(20,60),y=500,label=c(percentEGFR,percentMix))+
  geom_vline(xintercept=thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName=='NRG1')])
ggsave(filename = 'graphs/suppplementGraphs/RNAFISH_sortNRG1.pdf', height = 9, width =7)

####
gene='FGFR1'
percentEGFR <- 100*dim(filter(data,variable==eval(gene) & dataset =='EGFRhigh' & value>thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName==eval(gene))]))[1]/dim(filter(data,variable==eval(gene) & dataset =='EGFRhigh'))[1]
percentMix <- 100*dim(filter(data,variable==eval(gene) & dataset =='mix' & value>thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName==eval(gene))]))[1]/dim(filter(data,variable==eval(gene) & dataset =='mix'))[1]
print(paste('Percentage', eval(gene), 'jackpots in EGFR high ', percentEGFR))
print(paste('Percentage', eval(gene), 'jackpots in mix ',percentMix))

ggplot(filter(data,variable=='FGFR1'),aes(x=jitter,fill=dataset))+
  geom_histogram(color='black')+facet_wrap(~dataset,ncol=1)+
  geom_rug()+theme_classic()+
  annotate("text",x=c(100,200),y=500,label=c(percentEGFR,percentMix))+
  geom_vline(xintercept=thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName=='FGFR1')])
ggsave(filename = 'graphs/suppplementGraphs/RNAFISH_sortFGFR1.pdf', height = 9, width =7)

###
gene='AXL'
percentEGFR <- 100*dim(filter(data,variable==eval(gene) & dataset =='EGFRhigh' & value>thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName==eval(gene))]))[1]/dim(filter(data,variable==eval(gene) & dataset =='EGFRhigh'))[1]
percentMix <- 100*dim(filter(data,variable==eval(gene) & dataset =='mix' & value>thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName==eval(gene))]))[1]/dim(filter(data,variable==eval(gene) & dataset =='mix'))[1]
print(paste('Percentage', eval(gene), 'jackpots in EGFR high ', percentEGFR))
print(paste('Percentage', eval(gene), 'jackpots in mix ',percentMix))

ggplot(filter(data,variable==gene),aes(x=jitter,fill=dataset))+
  geom_histogram(color='black')+facet_wrap(~dataset,ncol=1)+
  geom_rug()+theme_classic()+
  annotate("text",x=c(50,200),y=200,label=c(percentEGFR,percentMix))+
  geom_vline(xintercept=thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName==gene)])
ggsave(filename = 'graphs/suppplementGraphs/RNAFISH_sortAXL.pdf', height = 9, width =7)


#### Try MITF....
gene='MITF'
percentEGFR <- 100*dim(filter(data,variable==eval(gene) & dataset =='EGFRhigh' & value>thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName==eval(gene))]))[1]/dim(filter(data,variable==eval(gene) & dataset =='EGFRhigh'))[1]
percentMix <- 100*dim(filter(data,variable==eval(gene) & dataset =='mix' & value>thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName==eval(gene))]))[1]/dim(filter(data,variable==eval(gene) & dataset =='mix'))[1]
print(paste('Percentage', eval(gene), 'jackpots in EGFR high ', percentEGFR))
print(paste('Percentage', eval(gene), 'jackpots in mix ',percentMix))
medianEGFR <- median(filter(data,variable==eval(gene) & dataset =='EGFRhigh')$value)
medianMix <- median(filter(data,variable==eval(gene) & dataset =='mix')$value)
print(paste('Median', eval(gene), 'jackpots in EGFR high ', medianEGFR))
print(paste('Median', eval(gene), 'jackpots in mix ', medianMix))
meanEGFR <- mean(filter(data,variable==eval(gene) & dataset =='EGFRhigh')$value)
meanMix <- mean(filter(data,variable==eval(gene) & dataset =='mix')$value)
print(paste('Mean', eval(gene), 'jackpots in EGFR high ', meanEGFR))
print(paste('Mean', eval(gene), 'jackpots in mix ', meanMix))

ggplot(filter(data,variable==gene),aes(x=jitter,fill=dataset))+
  geom_histogram(color='black')+facet_wrap(~dataset,ncol=1)+
  geom_rug()+theme_classic()+
  annotate("text",x=c(100,200),y=500,label=c(medianEGFR,medianMix))
  #geom_vline(xintercept=thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName==gene)])
ggsave(filename = 'graphs/suppplementGraphs/RNAFISH_sortMITF.pdf', height = 9, width =7)


###
gene='SOX10'
percentEGFR <- 100*dim(filter(data,variable==eval(gene) & dataset =='EGFRhigh' & value>thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName==eval(gene))]))[1]/dim(filter(data,variable==eval(gene) & dataset =='EGFRhigh'))[1]
percentMix <- 100*dim(filter(data,variable==eval(gene) & dataset =='mix' & value>thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName==eval(gene))]))[1]/dim(filter(data,variable==eval(gene) & dataset =='mix'))[1]
print(paste('Percentage', eval(gene), 'jackpots in EGFR high ', percentEGFR))
print(paste('Percentage', eval(gene), 'jackpots in mix ',percentMix))
medianEGFR <- median(filter(data,variable==eval(gene) & dataset =='EGFRhigh')$value)
medianMix <- median(filter(data,variable==eval(gene) & dataset =='mix')$value)
print(paste('Median', eval(gene), 'jackpots in EGFR high ', medianEGFR))
print(paste('Median', eval(gene), 'jackpots in mix ', medianMix))

ggplot(filter(data,variable==gene),aes(x=jitter,fill=dataset))+
  geom_histogram(color='black')+facet_wrap(~dataset,ncol=1)+
  geom_rug()+theme_classic()+
  annotate("text",x=c(300,500),y=200,label=c(medianEGFR,medianMix))
  #geom_vline(xintercept=thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName==gene)])
ggsave(filename = 'graphs/suppplementGraphs/RNAFISH_sortSOX10.pdf', height = 9, width =7)

###
gene='GAPDH'
percentEGFR <- 100*dim(filter(data,variable==eval(gene) & dataset =='EGFRhigh' & value>thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName==eval(gene))]))[1]/dim(filter(data,variable==eval(gene) & dataset =='EGFRhigh'))[1]
percentMix <- 100*dim(filter(data,variable==eval(gene) & dataset =='mix' & value>thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName==eval(gene))]))[1]/dim(filter(data,variable==eval(gene) & dataset =='mix'))[1]
print(paste('Percentage', eval(gene), 'jackpots in EGFR high ', percentEGFR))
print(paste('Percentage', eval(gene), 'jackpots in mix ',percentMix))
medianEGFR <- median(filter(data,variable==eval(gene) & dataset =='EGFRhigh')$value)
medianMix <- median(filter(data,variable==eval(gene) & dataset =='mix')$value)
print(paste('Median', eval(gene), 'jackpots in EGFR high ', medianEGFR))
print(paste('Median', eval(gene), 'jackpots in mix ', medianMix))

ggplot(filter(data,variable==gene),aes(x=jitter,fill=dataset))+
  geom_histogram(color='black')+facet_wrap(~dataset,ncol=1)+
  geom_rug()+theme_classic()+
  annotate("text",x=c(1000,1500),y=200,label=c(medianEGFR,medianMix))
  #geom_vline(xintercept=thresholds_WM9noDrug$threshold[which(thresholds_WM9noDrug$geneName==gene)])
ggsave(filename = 'graphs/suppplementGraphs/RNAFISH_sortGAPDH.pdf', height = 9, width =7)

EGFRhigh <- filter(data, dataset=='EGFRhigh' & variable=='EGFR')
EGFRmix <- filter(data, dataset=='mix' & variable=='EGFR')

t.test(EGFRhigh$value,EGFRmix$value)

