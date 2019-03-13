# Biddy data analysis
library(readxl)
filepaths <- list.files('data', full.names=TRUE)

d <- lapply(filepaths, function(x) read_xlsx(x, skip=2))
names(d) <- paste('SuppData',2:8,sep='')
