#Auxilliary functions for bachelor-thesis
library(fields); library(akima); library(geoR); library(MASS); library(spatial); library(RSAGA)

#Constructing square map with adapted axes 
squareMap <- function(map, grid_size){
  square_map = expand.grid(y = seq(1:grid_size), x = seq(1:grid_size))
  square_map$z = as.vector(interp(x = map$x, y = map$y, z = map$z,
                                      xo = seq(1,max(map$x),length.out = grid_size),
                                      yo = seq(1,max(map$y),length.out = grid_size))$z)
  return(square_map)
}

#Constructs an isotropic correlation matrix based on the Exponential and Matern correlation function 
correlationMatrix <- function(sample_points, prediction_points, range, correlation_function = 0, v=1.5){
  ncolumns = dim(prediction_points)[1]; nrows = dim(sample_points)[1]
  
  corr_matrix = matrix(numeric(nrows*ncolumns),nrow=nrows, ncol=ncolumns)
  zero_matrix = matrix(numeric(2*nrows),nrow=2)
  
  #Calculating distances, one column of the matrix at the time
  for (i in 1:ncolumns){
    corr_matrix[,i] = sqrt( colSums((t(sample_points) - (zero_matrix + as.numeric(prediction_points[i,])) )^2))
  }
  if (correlation_function == "exponential"){
    return( Exponential(corr_matrix, range = range) )
  }
  else if(correlation_function == "matern"){
    return(Matern(corr_matrix, range = range, smoothness=v))
  } else {cat('Choose either exponential or matern','\n'); return()}
}

#Samples n samples with position from the grid
#(To check the values: z[ intersect(which(square_map$y==a),which(square_map$x==b))) ]
gridSampler <- function(n, map, design, noise = 0){
  if (design == 'regular'){
    xmax = max(map$x);   ymax = max(map$y);   xmin = min(map$x);   ymin = min(map$y)
    
    #Designed to have a perfect regular grid, disregarding the actual n
    xincrement = round( (xmax - xmin)/(round(sqrt(n)) - 1) )
    yincrement = round( (ymax - ymin)/(round(sqrt(n)) - 1) )
    
    samples = expand.grid(x = seq(xmin, xmax, xincrement), y = seq(ymin, ymax, yincrement))
    samples$z = map$z[(samples$x-1)*xmax + samples$y]
  } 
  else if (design == 'random'){
    xmax = max(map$x);   ymax = max(map$y);   xmin = min(map$x);   ymin = min(map$y)
    nr_samples = s
  }
  
  samples$z = samples$z + rnorm(length(samples$z),0, noise)    
  return(samples)
}

#Constructs a "priorfield", in this instance using a uniform prior for range and gamma for GRF variance
priorField <- function(size_sigma, size_range, alpha=5, beta=0.5){
  #Centre mean, with some (theoretical) std deviation increments from the theoretical mean
  increment = 2/floor(size_sigma/2)
  values = round(seq(-floor(size_sigma/2),-floor(size_sigma/2) + size_sigma - 1)*increment*sqrt(alpha)/beta + alpha/beta)
  lower = c(0, values[1:length(values) - 1] + diff(values)/2); upper = c(values[2:length(values)] - diff(values)/2, Inf)
  
  prior_field = expand.grid(sigma2=values, range = seq(1, size_range))
  prior_field$prob = rep(1/size_range * (pgamma(upper, alpha,beta) - pgamma(lower,alpha,beta)), times = size_range)
  return(prior_field)
}

#Performs GLS estimation of the samples
glmLite <- function(trend, data, correlation_matrix){
  if (trend == 'cubic'){
    formula = z ~ x*y + I(x^2) + I(y^2) + I(y^3) + I(y^3)
  }
  else if (trend == 'quadratic'){
    formula = z ~ x*y + I(x^2) + I(y^2)
  }
  else if (trend == 'linear'){
    formula = z ~ x + y 
  }
  else if (trend == 'simple'){
    formula = z ~ 1
  }
  else {cat('Choose either simple, linear or quadratic trend','\n'); return()}
  mf = model.frame(formula = formula, data = data);
  X  = model.matrix(attr(mf, "terms"), data = mf);
  Y  = model.response(mf);
  
  #Performing GLS 
  W = solve(correlation_matrix);
  covariance = solve(t(X)%*%W%*%X);
  beta = covariance%*%t(X)%*%W%*%Y;
  return(list(coefficients = beta, covar = covariance, formula = formula, X = X, terms = attr(mf, "terms")))
}

#Constructs posterior distribution
posteriorDistribution <- function(size_predictionx, size_predictiony, correlation_function = 'exponential', glm_est, GRFsigma2 = 1){
  
  prediction_grid = expand.grid(y = seq(1,size_predictiony), x = seq(1,size_predictionx));
  prediction_grid$z = numeric(size_predictionx*size_predictiony)
  
  mf = model.frame(formula = glm_est$formula, data = prediction_grid)
  Xprediction  = model.matrix(attr(mf, "terms"), data = mf)
  Xobserved = glm_est$X
  
  prediction_grid$z <- NULL
  prediction_pos = cbind(prediction_grid$x, prediction_grid$y)
  
  correlation_prediction_observed = correlationMatrix(cbind(samples$x,samples$y), prediction_pos, range = 5, correlation_function);
  correlation_prediction = correlationMatrix(prediction_pos, prediction_pos, range = 5, correlation_function);
  correlation_observed = correlationMatrix(cbind(samples$x,samples$y),cbind(samples$x,samples$y), range = 5, correlation_function);
  
  fitted_predictions = Xprediction%*%glm_est$coefficients;
  covariance_predictions_observed = Xprediction%*%glm_est$covar%*%t(Xobserved);
  
  #Constructing the different parts of the variance matrix of the multivariate normal
  sigma21 = t(correlation_prediction_observed * GRFsigma2) + covariance_predictions_observed
  sigma22 = correlation_observed * GRFsigma2;
  sigma11 = correlation_prediction * GRFsigma2;
  
  #Fitting predictions and variance 
  variance = prediction_grid
  prediction_grid$z = fitted_predictions + (sigma21)%*%solve(sigma22)%*%(samples$z - Xobserved%*%glm_object$coefficients)
  variance$z = diag( sigma11 + (sigma21)%*%solve(sigma22)%*%(t(sigma21)) )
  return(list(prediction = prediction_grid, variance = variance))
}