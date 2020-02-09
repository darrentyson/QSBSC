library(ggplot2) #Load in Seurat package
library(data.table)
library(gridExtra)
# setwd('C:/Users/Dora Obodo/VANDERBILT/1st Year/SPRING 2019/MM2 Quant Sys Bio SC/QSBSC/R_code/Obodo')
dataCellTag = fread('morrislab.celltag.dataset.csv', header= TRUE)
# tag_metadata.HF1

#Create stacked bar plots from metadata

celltypeplot = ggplot(data = dataCellTag, aes(x = Timepoint, fill = factor(CellType.Monocle))) + 
        geom_bar(position = 'fill') + labs(title = 'Stacked Bar plot of cell types', 
                                           caption = 'This is a plot of the labelled cell types for each cell in the data. ')
        scale_y_continuous(labels = scales::percent)



stateplot = ggplot(data = dataCellTag, aes(x = Timepoint, fill = factor(State.Monocle))) + 
  geom_bar(position = 'fill') + labs(title = 'Stacked Bar plot of cell states', 
                                     caption = 'This is a plot of the labelled cell states for each cell in the data. ')

tsneclusters = ggplot(dataCellTag, aes(x = tSNE_1, y = tSNE_2, color = factor(Cluster.Seurat))) + geom_point(size = 0.75) +
  labs(title = 't-SNE plot of genetically distinct cell clusters from t-sne', 
       caption = 'This is a plot of the distribution of genetically distinct clusters in the data ')


clusterplot = ggplot(data = dataCellTag, aes(x = Timepoint, fill = factor(Cluster.Seurat))) + 
  geom_bar(position = 'fill') + labs(title = 'Stacked Bar plot of cell clusters from t-sne', 
                                     caption = 'This is a plot of the distribution of states for each cell in the data. \nUnfortunately, they did not include metadata for which cells \n belong to which of the 4 transition states ')



plots.list = list(celltypeplot,stateplot, tsneclusters, clusterplot)
plots = marrangeGrob(plots.list, nrow = 2, ncol = 1)


ggsave("created_plots_from_Biddy.pdf", plots,  width = 21, height = 29.7, units = "cm")





