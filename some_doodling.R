library(geoR); library(MASS); library(akima); library(fields); library(RSAGA)
#Adapting data from matrix form, to original + square 50x50
original_mapping = grid.to.xyz(t(volcano))

#Constructing square mapping with adapted axes.
square_mapping = interp(x = original_mapping$x, y = original_mapping$y, z = original_mapping$z,
                        xo = seq(1,max(original_mapping$x),length.out = 50),
                        yo = seq(1,max(original_mapping$y),length.out = 50))
square_mapping$x = seq(50); square_mapping$y = seq(50)

#Constructing isotropic correlation matrix. Now using Exponential
corr_matrix <- function(pred_points, sample_points, range){
  ncolumns = dim(pred_points)[1]; nrows = dim(sample_points)[1]
  
  corr_matrix = matrix(numeric(nrows*ncolumns),nrow=nrows)
  
  for (i in 1:ncolumns){
    corr_matrix[,i] = sqrt( rowSums((sample_points - matrix(numeric(2*nrows) + 1,ncol=2)%*%diag(pred_points[i,]))^2) )
  }
  return( Exponential(corr_matrix, range = range) )
}

prior_field = expand.grid(seq(10), seq(10))

#prior ranges; prior variance; prior field; numeric approximation integrat;
#exponential covariance func; matern?; squared exp?; different polynomial trends?;

#krigerestimate of posterior

#different designs
#multivariate normal for everything!