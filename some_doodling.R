library(geoR); library(MASS); library(akima); library(fields); library(RSAGA)
#Adapting data from matrix form, to original + square 50x50
original_mapping = grid.to.xyz(t(volcano))

#Constructing square mapping with adapted axes.
square_mapping = interp(x = original_mapping$x, y = original_mapping$y, z = original_mapping$z,
                        xo = seq(1,max(original_mapping$x),length.out = 50),
                        yo = seq(1,max(original_mapping$y),length.out = 50))
square_mapping$x = seq(50); square_mapping$y = seq(50)

prior ranges; prior variance; prior field; numeric approximation integrat;
exponential covariance func; matern?; squared exp?; different polynomial trends?;

krigerestimate of posterior

different designs
multivariate normal for everything!