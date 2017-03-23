#Auxilliary functions for bachelor-thesis, relevant for stats
library(fields); library(akima); library(geoR); library(MASS); library(spatial); library(RSAGA); library(lattice)

#----------------------------------------------------------------------------------------

#Constructs an isotropic correlation matrix based on the Exponential and Matern correlation function 
correlationMatrix <- function(sample_points, prediction_points, range, correlation_function = 0, v=1.5){
  
  corr_matrix = constructDistanceMatrix(sample_points, prediction_points)
  
  if (correlation_function == "exponential"){
    return( Exponential(corr_matrix, range = range) )
  }
  else if(correlation_function == "matern"){
    return(Matern(corr_matrix, range = range, smoothness=v))
  } else {cat('Choose either exponential or matern','\n'); return()}
}

#----------------------------------------------------------------------------------------

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

#----------------------------------------------------------------------------------------

#Performs GLS estimation of the samples
glmLite <- function(trend, data, correlation_matrix, sigma2 = 1){
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
  W = solve(correlation_matrix*sigma2);
  covariance = solve(t(X)%*%W%*%X);
  beta = covariance%*%t(X)%*%W%*%Y;
  return(list(coefficients = beta, covar = covariance, formula = formula, X = X, terms = attr(mf, "terms")))
}

#----------------------------------------------------------------------------------------

#Constructs posterior distribution
posteriorDistribution <- function(size_predictionx, size_predictiony, samples, correlation_function = 'exponential', glm_est, GRFsigma2 = 1){
  
  prediction_grid = expand.grid(x = seq(1,size_predictiony), y = seq(1,size_predictionx));
  prediction_grid$z = numeric(size_predictionx*size_predictiony)
  
  mf = model.frame(formula = glm_est$formula, data = prediction_grid)
  Xprediction  = model.matrix(attr(mf, "terms"), data = mf)
  Xobserved = glm_est$X
  
  prediction_grid$z <- NULL
  prediction_pos = cbind(prediction_grid$x, prediction_grid$y)
  sampling_pos = cbind(samples$x,samples$y)
  
  correlation_prediction_observed = correlationMatrix(sampling_pos, prediction_pos, range = 5, correlation_function);
  
  correlation_prediction = correlationMatrix(prediction_pos, prediction_pos, range = 5, correlation_function);
  
  correlation_observed = correlationMatrix(sampling_pos, sampling_pos, range = 5, correlation_function);
  
  fitted_predictions = Xprediction%*%glm_est$coefficients;
  covariance_predictions_observed = Xprediction%*%glm_est$covar%*%t(Xobserved);
  
  #Constructing the different parts of the variance matrix of the multivariate normal
  sigma12 = correlation_prediction_observed * GRFsigma2 + covariance_predictions_observed
  sigma22 = correlation_observed * GRFsigma2;
  sigma11 = correlation_prediction * GRFsigma2;
  sigma12sigma22 = sigma12%*%solve(sigma22)
  
  #Fitting predictions and variance 
  variance = prediction_grid
  prediction_grid$z = fitted_predictions + sigma12sigma22%*%(samples$z - Xobserved%*%glm_object$coefficients)
  variance$z = diag( sigma11 + sigma12sigma22%*%(t(sigma12)) )
  
  fit = prediction_grid; fit$z = fitted_predictions
  beforeFit = prediction_grid; beforeFit$z = sigma12sigma22%*%(samples$z - Xobserved%*%glm_object$coefficients)
  
  othervariance = variance
  variance = xyz.to.grid(variance)
  prediction_grid = xyz.to.grid(prediction_grid)
  fit = xyz.to.grid(fit)
  beforeFit = xyz.to.grid(beforeFit)
  
  return(list(prediction = prediction_grid, variance = variance, fit = fit, beforeFit = beforeFit,  othervariance = othervariance))
}
