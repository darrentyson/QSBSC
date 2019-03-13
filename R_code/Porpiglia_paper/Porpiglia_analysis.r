# attempting to load FlowRepository data for Porpiglia paper

library(FlowSOM)
library(flowCore)
library(openCyto)
library(MEM)
library(tidyverse)


# a <- lapply(pf.FCS[[1]], function(x)  as.data.frame(x@exprs))
# n <- basename(names(a))
# 
# for(i in n) a[[grep(i,names(a))]] <-  
#     cbind( a[[grep(i,names(a))]],
#            data.frame(sample=rep(i,nrow(a[[grep(i,names(a))]]))))
# 
# a <- do.call(rbind,a)
# rownames(a) <- NULL
# 
# s <- a[grep('weeks_1',a$sample),][1:100000,]
# 
# #xyplot(`DNA1(Ir191)Dd` ~ `Cell_length` | `sample`, data=a, pch='.')
# 
# markers <- lapply(pf.FCS, function(x)  unique(unlist(sapply(x, function(z) z@parameters@data['desc']))))

f <- list.files('ZY3K_filtered_files',full.names=TRUE)
d <-  lapply(f, read.FCS)

n <- unique(unlist(lapply(d, colnames)))
