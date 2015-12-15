### Preliminaries 
library(dplyr)

#----------------------------------------------------
### Binning Functions
#----------------------------------------------------

## Standard 1d Rectangular Binning
StandRectBin1d <- function(xs, origin, width){
  binxs <- origin + width*( ceiling((xs-origin)/width) - 0.5 )
  binxs[xs == origin] <- origin + width/2
  return(binxs)
}


## Random 1d Rectangular Binning
# Establish bounding bin centers for each observation
# Then use Unif(0,1) draws compared to assignment probs to allocate
# Reassignment for values below first center (bx1) or
#     above highest center (bxJ)
RandRectBin1d <- function(xs, origin, width){
  bx1 <- origin + width/2
  bxJ <- origin + width*(floor((max(xs)-bx1)/width) + .5)
  lbs <- bx1 + width*floor((xs-bx1)/width)
  ubs <- bx1 + width*ceiling((xs-bx1)/width)
  # initially assign all values to upper bound
  binxs <- ubs
  # then use runif to reassign based on distance from 
  plower <- (ubs - xs)/width
  lowerindex <- (plower > runif(length(xs), 0, 1))
  binxs[lowerindex] <- lbs[lowerindex]
  binxs[xs < bx1] <- bx1
  binxs[xs > bxJ] <- bxJ
  return(binxs)
}


### Create binning function that does not calculate loss
# to be used to track computation time for binning

## 2d Rectangular Binning (for either Standard or Random)
# standard is straight forward extension of 1d binning
# random binning needs post processing to calculate minimum spatial loss
RectBin2d <- function(xs,ys, originx, originy, widthx, widthy, type="standard"){
  if(type=="standard"){
    tempdat <- data.frame(xs = xs, ys=ys,
                          binxs = StandRectBin1d(xs,originx,widthx),
                          binys = StandRectBin1d(ys,originy,widthy))
  }
  if(type=="random"){
    tempdat<- data.frame(xs = xs, ys=ys,
                         binxs =  RandRectBin1d(xs,originx,widthx),
                         binys =  RandRectBin1d(ys,originy,widthy))
  }
  outdat <- tempdat %>% 
    group_by(binxs,binys) %>% 
    summarise(binfreq = length(xs))
  return(outdat[,1:3])
}


## Frequency Binning  
# allows for standard or quantile binning of either  bin counts or log bin counts (4 combinations)
# input requires binned data output, number of freq breaks and type of freq binning
# output of frequency bin values, labels and loss are attached the original then returned
freqBin <- function(binout, binType="standard", ncolor, logCount=FALSE){ 
  cs <- binout$binfreq
  # add frequency bin centers and bin labels of the form "(a,b]" to binout
  # binning depends on type / log counts
  if (logCount)  cs <- log(binout$binfreq)
  if(binType=="standard"){
    #width <- ceiling((max(cs)-min(cs))/ncolor)
    width <- (max(cs)-min(cs))/ncolor 
    binout$freqgroup <- round(StandRectBin1d(cs, min(cs) , width),5)
    binout$freqlabel <- paste("(",round(binout$freqgroup - width/2,1),
                                   ",",round(binout$freqgroup + width/2,1),"]",sep="")
    #close interval for smallest counts
    closeidx <- binout$freqlabel == min(binout$freqlabel)
    binout$freqlabel[closeidx] <- paste("[",round(min(binout$freqgroup)- width/2,1),
                                             ",",round(min(binout$freqgroup) + width/2,1),"]",sep="")
  } 
  if(binType=="quantile"){
    binout$freqgroup <- as.numeric(round(QuantBin1d(cs, ncolor),5))
    quantbounds <- unique(quantile(cs, (0:ncolor)/ncolor))
    if (length(quantbounds)-1 < ncolor) warning("two or more quantiles of the data have same value due to many bins with equivalent frequencies, color will be rendered equivalently for bins with matching quantiles")
    binout$freqlabel <- cut(cs, breaks=quantbounds, include.lowest=TRUE)
  }
 return(binout)
}

