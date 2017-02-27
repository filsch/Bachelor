library(geoR); library(MASS); library(akima); library(fields); library(RSAGA)
set.seed(4545)
#Adapting data from matrix form
original_mapping = grid.to.xyz(t(volcano))

#Constructing square mapping with adapted axes to a 50x50 
#square_mapping = expand.grid(x = seq(1, 50), y= seq(1,50))
square_mapping = expand.grid(y = seq(1:50), x = seq(1:50))
square_mapping$z = as.vector(interp(x = original_mapping$x, y = original_mapping$y, z = original_mapping$z,
                                    xo = seq(1,max(original_mapping$x),length.out = 50),
                                    yo = seq(1,max(original_mapping$y),length.out = 50))$z)


#Constructing isotropic correlation matrix. Now using Exponential
corr_matrix <- function(pred_points, sample_points, range){
  ncolumns = dim(pred_points)[1]; nrows = dim(sample_points)[1]
  
  corr_matrix = matrix(numeric(nrows*ncolumns),nrow=nrows)
  zero_matrix = matrix(numeric(2*ncolumns),nrow=2)
  for (i in 1:ncolumns){
    corr_matrix[,i] = sqrt( colSums((t(sample_points) - (zero_matrix + pred_points[i,]))^2))
  }
  return( Exponential(corr_matrix, range = range) )
}

#--------------------------------------------------------------------------------------------------------

#Assigning prior values and estimating beta
alpha = 5; beta = 0.5; measurement_sigma = 0.1; 
#"Uniform"
prior_field = expand.grid(sigma2=qgamma(0.09*seq(1:10),alpha,beta), range = seq(1,10))
prior_field$prob = numeric(100) + 0.01
#Centre mean, 0.4*std indrement out from centre (but why? because)
values = round(seq(-5,4)*0.4*sqrt(alpha)/beta + alpha/beta)
lower = c(0, values[1:length(values) - 1] + diff(values)/2); upper = c(values[2:length(values)] - diff(values)/2, Inf)
prior_field = expand.grid(sigma2=values, range = seq(1,10))
prior_field$prob = rep(0.1*(pgamma(upper, alpha,beta) - pgamma(lower,alpha,beta)), times = 10)


#--------------------------------------------------------------------------------------------------------

#Generating sample from grid by design
#Regular, 64 samples
grid_samples = expand.grid(x = seq(1, 50, 8), y= seq(1,50,8))
#grid_samples$z = square_mapping$z[matrix(cbind(grid_samples$x, grid_samples$y), nrow=49,ncol=2)]
#(To check the values intersect(which(square_mapping$y==1),which(square_mapping$x==34)))
grid_samples$z = square_mapping$z[(grid_samples$x-1)*50 + grid_samples$y] + rnorm(0, measurement_sigma)


#--------------------------------------------------------------------------------------------------------
# Fitting trend w.r.t. covariance matrix and sample positions
input_correlation = cbind(grid_samples$x,grid_samples$y)

glm_object = myglm_lite('quadratic', grid_samples, corr_matrix(input_correlation, input_correlation, 5))

#--------------------------------------------------------------------------------------------------------

#Predicting the data
prediction_grid = expand.grid(y = seq(1,50), x = seq(1,50))




#--------------------------------------------------------------------------------------------------------
#KRIGING:
#Computing the predictor:
posterior_kriger = krige.conv(coords=cbind(grid_samples$x,grid_samples$y), data=grid_samples$z, locations = prediction_grid,
                              krige=krige.control(type.krige= "OK", trend.d = "2nd", trend.l = "2nd", cov.model="exponential", cov.pars=c(sigma2, range)));
#Plotting some results
par(mfrow=c(2,2))
plot_values = list(z=matrix(posterior_kriger$predict,nrow=50,ncol=50,byrow=FALSE), x = seq(0,50), y = seq(0,50))
image.plot(plot_values, main="Kriging predictions, no noise in datamodel", zlim=c(80,200))

plot_values$z = matrix(sqrt(posterior_kriger$krige.var),nrow=50,ncol=50,byrow=FALSE) 
image.plot(plot_values, main="Kriging std.dev., no noise on data Y", zlim=c(0,5))

plot_values$z = square_mapping$z - matrix(posterior_kriger$predict,nrow=50,ncol=50,byrow=FALSE) 

image.plot(plot_values, main="Kriging vs actual", zlim=c(-30, 30))
plot(grid_samples$x,grid_samples$y, main="Observed points")
