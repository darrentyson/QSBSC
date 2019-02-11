# execute this code to install all necessary packages for Dr. Irish's lectures
source("https://bioconductor.org/biocLite.R")
biocLite("flowCore", ask=FALSE)
biocLite("FlowSOM", ask=FALSE)

getPackageIfNeeded <- function(pkg) {
  if (!require(pkg, character.only=TRUE))
    install.packages(pkgs=pkg, dependencies=TRUE)
}

pkgs <- c('ggplots','ggplots2','hexbin','viridis','ggExtra','tidyverse','devtools')

sapply(pkgs,getPackageIfNeeded)
