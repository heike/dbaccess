### Add 2010 to 2013 to baseball data
setwd("C:\\Users\\maurerkt\\Documents\\GitHub\\dbaccess\\versionForJCGS\\data")

source("RectangularBinningFunctions.R")

load(file="pitching.Rdata")

pitch2013 <- read.csv("Pitching2013.csv")[,c("G","SO")]
head(pitch2013)

binout1 <- RectBin2d(xs=pitch2013$G,ys=pitch2013$SO, 
                     originx = .5, originy = -.5, 
                     widthx = 1, widthy = 1, type="standard")
names(binout1) <- names(d1)
head(binout1)

# It turns out we already have the 2013 updated data version, we never updated the year in the text of the paper