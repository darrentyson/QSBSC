# parameters file for WM9 sort mix well

load_WM9_sortMixwell_20150725 <- function(){
  # define directories
  plotDir <- plotDir
  
  # this will be the begining of the name for every plot generated by this script
  plotName <- 'WM9_sortMixwell_20150725'
  
  # Load data 
  data <- read.csv('dentistData/WM9_sortMixWell_20150725.txt')
  
  # Load thresholds 
  thresholds <- data.frame(geneName=c('EGFR','SOX10','CCNA2','GAPDH','WNT5A','PDGFRB',
                                      'PDGFC','SERPINE1','NGFR','NRG1','FOSL1',
                                      'VEGFC','AXL','MITF','LOXL2','RUNX2','FGFR1',
                                      'JUN','VGF'),threshold=c(10, #EGFR
                                                               150, #Sox10
                                                               25, #CCNA2 - no expression
                                                               750, #gapdh
                                                               15, #WNT5A
                                                               10, #PDGFRB
                                                               10, #PDGFC - nothing.
                                                               21, #serpine1
                                                               106, #ngfr
                                                               10, #nrg1
                                                               10, #fosl1 - no signal
                                                               10, #vegfc
                                                               30, #axl
                                                               120, #mitf
                                                               15, #loxl2 - no signal
                                                               10, #runx2 - no signal
                                                               80, #fgfr1
                                                               33, #jun
                                                               90)) #vgf - no signal
  
  WM9_sortMixwell_20150725 <- list(plotName,plotDir,data,thresholds)
  return(WM9_sortMixwell_20150725)
}


