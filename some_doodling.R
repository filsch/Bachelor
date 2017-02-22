library(geoR); library(MASS); library(akima); library(fields); library(RSAGA)
set.seed(4545)
#Adapting data from matrix form
original_mapping = grid.to.xyz(t(volcano))

#Constructing square mapping with adapted axes to a 50x50 
#square_mapping = expand.grid(x = seq(1, 50), y= seq(1,50))
square_mapping = interp(x = original_mapping$x, y = original_mapping$y, z = original_mapping$z,
                        xo = seq(1,max(original_mapping$x),length.out = 50),
                        yo = seq(1,max(original_mapping$y),length.out = 50))
square_mapping$x = seq(0,50); square_mapping$y = seq(0,50)

#Constructing isotropic correlation matrix. Now using Exponential
corr_matrix <- function(pred_points, sample_points, range){
  ncolumns = dim(pred_points)[1]; nrows = dim(sample_points)[1]
  
  corr_matrix = matrix(numeric(nrows*ncolumns),nrow=nrows)
  
  for (i in 1:ncolumns){
    corr_matrix[,i] = sqrt( rowSums((sample_points - matrix(numeric(2*nrows) + 1,ncol=2)%*%diag(pred_points[i,]))^2) )
  }
  return( Exponential(corr_matrix, range = range) )
}

#Generating sample from grid by design
#Regular, 64 samples
grid_samples = expand.grid(x = seq(1, 50, 8), y= seq(1,50,8))
grid_samples$z = square_mapping$z[matrix(cbind(grid_samples$x, grid_samples$y), nrow=49,ncol=2)]

#Assigning prior values and estimating beta
alpha = 5; beta = 0.5;
#"Uniform"
prior_field = expand.grid(sigma2=qgamma(0.09*seq(1:10),alpha,beta), range = seq(1,10))
prior_field$prob = numeric(100) + 0.01
#Centre mean, 0.4*std indrement out from centre
values = seq(-5,4)*0.4*sqrt(alpha)/beta + alpha/beta
lower = c(0, values[2:length(values)] - 0.2*sqrt(alpha)/beta); upper = c(values[1:length(values)-1]  + 0.2*sqrt(alpha)/beta,Inf)
prior_field = expand.grid(sigma2=values, range = seq(1,10))
prior_field$prob = rep(0.1*(pgamma(upper, alpha,beta) - pgamma(lower,alpha,beta)), times = 10)

# GLM must be done manually, to account for covariance (How to in GLM?) 
# beta = glm(grid_samples$z ~ grid_samples$x*grid_samples$y + I(grid_samples$x*grid_samples$x) + I(grid_samples$y*grid_samples$y))

prediction_grid = expand.grid(seq(50),seq(50))

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


#prior ranges; prior variance; prior field; numeric approximation integrat;
#exponential covariance func; matern?; squared exp?; different polynomial trends?;

#krigerestimate of posterior

#different designs
#multivariate normal for everything!