# thresholding script...

thresholdPercentile <- function(data, percentile) {
  dataout <- data
  thresholds <- data.frame(geneName=colnames(data)[4:22],threshold=rep(0,19))
  for (i in 1:(dim(data)[2]-3)) {
    thresholds$threshold[i] <- quantile(data[,i+3],1-percentile)
    
    #dataout[,i+3] <- as.integer(data[,i+3]>threshold)
  }
  return(dataout)
}

getWM9thresholds <- function(){
  geneThresholds <- data.frame(geneName=c('EGFR','SOX10','CCNA2','GAPDH','WNT5A','PDGFRB',
                                          'PDGFC','SERPINE1','NGFR','NRG1','FOSL1',
                                          'VEGFC','AXL','MITF','LOXL2','RUNX2','FGFR1',
                                          'JUN','VGF'),threshold=c(10, #EGFR
                                                                   150, #Sox10
                                                                   75, #CCNA2
                                                                   750, #gapdh
                                                                   15, #WNT5A
                                                                   5, #PDGFRB
                                                                   10, #PDGFC
                                                                   30, #serpine1
                                                                   100, #ngfr
                                                                   15, #nrg1
                                                                   5, #fosl1
                                                                   10, #vegfc
                                                                   50, #axl
                                                                   120, #mitf
                                                                   10, #loxl2
                                                                   10, #runx2
                                                                   15, #fgfr1
                                                                   20, #jun
                                                                   200)) #vgf
  return(geneThresholds)
}

thresholdManualWM9NoDrug <- function(data) {
  geneThresholds <- getWM9thresholds()
  dataout <- data
  for (i in 1:(dim(data)[2]-3)) {
    threshold <- filter(geneThresholds, geneName == colnames(data[i+3]))
    dataout[,i+3] <- as.integer(data[,i+3]>threshold$threshold)
  }
  return(dataout)
  
}
  

thresholdManualPrimaryMelanocytes <- function(data) {
  geneThresholds <- data.frame(geneName=c('EGFR','SOX10','CCNA2','GAPDH','WNT5A','PDGFRB',
                                          'PDGFC','SERPINE1','NGFR','NRG1','FOSL1',
                                          'VEGFC','AXL','MITF','LOXL2','RUNX2','FGFR1',
                                          'JUN','VGF'),threshold=c(3, #EGFR
                                                                   200, #Sox10
                                                                   15, #CCNA2
                                                                   500, #gapdh
                                                                   15, #WNT5A
                                                                   10, #PDGFRB
                                                                   10, #PDGFC
                                                                   7, #serpine1
                                                                   25, #ngfr
                                                                   6, #nrg1
                                                                   12, #fosl1
                                                                   8, #vegfc
                                                                   20, #axl
                                                                   75, #mitf
                                                                   13, #loxl2
                                                                   12, #runx2
                                                                   30, #fgfr1
                                                                   4, #jun
                                                                   5)) #vgf
  dataout <- data
  for (i in 1:(dim(data)[2]-3)) {
    threshold <- filter(geneThresholds, geneName == colnames(data[i+3]))
    dataout[,i+3] <- as.integer(data[,i+3]>threshold$threshold)
  }
  return(dataout)
  
}
thresholdManualColonyData <- function(data) {
  geneThresholds <- data.frame(geneName=c('EGFR','SOX10','CCNA2','GAPDH','WNT5A','PDGFRB',
                                          'PDGFC','SERPINE1','NGFR','NRG1','FOSL1',
                                          'VEGFC','AXL','MITF','LOXL2','RUNX2','FGFR1',
                                          'JUN','VGF'),threshold=c(10, #EGFR
                                                                   150, #Sox10
                                                                   6, #CCNA2
                                                                   450, #gapdh
                                                                   10, #WNT5A
                                                                   10, #PDGFRB
                                                                   11, #PDGFC
                                                                   25, #serpine1
                                                                   7, #ngfr - look at data
                                                                   10, #nrg1
                                                                   11, #fosl1
                                                                   11, #vegfc
                                                                   11, #axl
                                                                   100, #mitf - look at data
                                                                   10, #loxl2
                                                                   10, #runx2
                                                                   15, #fgfr1
                                                                   15, #jun
                                                                   25)) #vgf - look at data
  dataout <- data
  for (i in 1:(dim(data)[2]-3)) {
    threshold <- filter(geneThresholds, geneName == colnames(data[i+3]))
    dataout[,i+3] <- as.integer(data[,i+3]>threshold$threshold)
  }
  return(dataout)
  
}


thresholdManualWM94weeks <- function(data) {
  geneThresholds <- data.frame(geneName=c('EGFR','SOX10','CCNA2','GAPDH','WNT5A','PDGFRB',
                                          'PDGFC','SERPINE1','NGFR','NRG1','FOSL1',
                                          'VEGFC','AXL','MITF','LOXL2','RUNX2','FGFR1',
                                          'JUN','VGF'),threshold=c(10, #EGFR
                                                                   500, #Sox10
                                                                   75, #CCNA2
                                                                   750, #gapdh
                                                                   10, #WNT5A
                                                                   10, #PDGFRB
                                                                   10, #PDGFC
                                                                   50, #serpine1
                                                                   70, #ngfr
                                                                   15, #nrg1
                                                                   5, #fosl1
                                                                   20, #vegfc
                                                                   50, #axl
                                                                   120, #mitf
                                                                   30, #loxl2
                                                                   10, #runx2
                                                                   10, #fgfr1
                                                                   20, #jun
                                                                   20)) #vgf
  dataout <- data
  for (i in 1:(dim(data)[2]-3)) {
    threshold <- filter(geneThresholds, geneName == colnames(data[i+3]))
    dataout[,i+3] <- as.integer(data[,i+3]>threshold$threshold)
  }
  return(dataout)
  
}
