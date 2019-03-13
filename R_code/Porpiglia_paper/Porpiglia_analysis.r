# attempting to load FlowRepository data for Porpiglia paper

library(FlowSOM)
library(flowCore)
library(openCyto)
library(MEM)
library(tidyverse)

# get file names of FCS files
f <- list.files('ZY3K_filtered_files',full.names=TRUE)

# read files into list
d <-  lapply(f, read.FCS)

# colnames are unique metals detected in different Cytof channels
channels <- unique(unlist(lapply(d, colnames)))
# extract marker names from descriptions
markers <- unique(unlist(sapply(d, function(z) z@parameters@data['desc'])))

# keep only expression level data
d <- lapply(d, function(x) as.data.frame(x@exprs))

# extract sample name and replicate number from file names
# first remove extraneous
samp_names <- gsub('EP20160824-InjuryTimeCourse-','',basename(f))
samp_names <- gsub('_alpha7__CD9_.fcs','',samp_names)
# extract replicate number from samp_names
rep_num <- sapply(samp_names, function(x) strsplit(x,'_')[[1]][2])
# remove replicate number from samp_names
samp_names <- sapply(samp_names, function(x) strsplit(x,'_')[[1]][1])


# add columns for sample name and replicate
d <- lapply(seq_along(d), 
    function(x) cbind(
        d[[x]],
        data.frame( sample=rep(samp_names[x],nrow(d[[x]])),
                    rep=rep(rep_num[x],nrow(d[[x]])))))

# use same number of cells from each sample and replicate
# lowest number of cells in any sample is 1854
# set a seed if you want this to always produce the same result
set.seed(1854)
d <- lapply(d, function(x) sample_n(x,1854))

# combine all data into single data.frame (tibble)
d <- do.call(rbind, d)
colnames(d) <- c(markers,'sample','rep')

# remove event length and beadDist columns
d <- d[,!colnames(d) %in% c('beadDist','Event_length')]
# put sample and rep columns first 
d <- d[,c('sample','rep',head(colnames(d),ncol(d)-2))]

### Data should now be ready for applying dimensionality reduction, clustering, etc
