# 

source('plotScripts/RNAFISH/RNAFISH_functions.R')
loadLibraries()

data <- read.csv('dentistData/WM9_resistant_20150802.txt')
print(paste('Number of cells in resistance data:',dim(data)[1]))
data2 <- data[,c(1:4,8,9,11,12,13,15,16,21)]
data2$dataset <- 'resistant'
dataMelt <- melt(data2, id.vars = c('cellID','Xpos','Ypos','dataset'))

data_noDrug <- read.csv('dentistData/WM9_noDrug_20150810.txt')
data_noDrug2 <- data_noDrug[,c(1:4,8,9,11,12,13,15,16,21)]
data_noDrug2$dataset <- 'no drug'
data_noDrugMelt <- melt(data_noDrug2, id.vars = c('cellID','Xpos','Ypos','dataset'))
print(paste('Number of cells in no drug data:',dim(data_noDrug)[1]))

allData <- rbind(dataMelt, data_noDrugMelt)

ggplot(allData,aes(x=dataset,y=value,fill=dataset))+
  geom_boxplot()+
  facet_wrap(~variable, scales='free')+
  theme_classic()
ggsave(filename = 'graphs/suppplementGraphs/noDrugVsResistant_boxplot.pdf',height=8, width=8)

# for asterisks on these plots
group_by(allData, dataset, variable) %>% dplyr::summarise(median(value), quantile(value,0.98))


