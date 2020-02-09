# Biddy data analysi/
setwd('C:/Users/Dora Obodo/VANDERBILT/1st Year/SPRING 2019/MM2 Quant Sys Bio SC/QSBSC/R_code/Biddy_paper')
library(readxl)
filepaths <- list.files('data', full.names=TRUE)

d <- lapply(filepaths, function(x) read_xlsx(x, skip=2))
names(d) <- paste('SuppData',2:8,sep='')
