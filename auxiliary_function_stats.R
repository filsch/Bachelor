#Auxilliary functions for bachelor-thesis, relevant for stats
library(fields); library(akima); library(geoR); library(MASS); library(spatial); library(RSAGA); library(lattice)


#Constructs an isotropic correlation matrix based on the Exponential and Matern correlation function 
correlationMatrix <- function(corr_matrix, range, correlation_function = 0, v=1.5){
  if (correlation_function == "exponential"){
    return( Exponential(corr_matrix, range = range) )
  }
  else if(correlation_function == "matern"){
    return(Matern(corr_matrix, range = range, smoothness=v))
  } else {cat('Choose either exponential or matern','\n'); return()}
}


#Constructs a "priorfield", in this instance using a uniform prior for range and gamma for GRF variance
priorField <- function(size_sigma, size_tau, sigma2_prior = "gamma", rau_prior = "exponential", sigma2_param = ..., tau_param = ...){
  #Centre mean, with some (theoretical) std deviation increments around the mean set by parameters
  if (sigma2_prior == "gamma"){
    shape = sigma2_param$alpha 
    rate = sigma2_param$beta
    if( shape <= 0 || rate <= 0){
      cat('You need positive parameters for Gamma distribution \n')
      return(NULL)
    }
    half_adjustment = 1/(2*size_sigma)
    discretization = seq(0,1,length.out=size_sigma + 1)
    field_values_sigma2 = discretization[1:size_sigma] + half_adjustment
    field_values_sigma2 = qgamma(field_values_sigma2, shape=shape, rate=rate)
    probability_sigma = diff(discretization)
    
    #-----
    #curve(dgamma(x,shape=shape,rate=rate),from=0,to=7)
    #points(field_values_sigma2,dgamma(field_values_sigma2,shape=shape,rate=rate),col="red",pch=4)
  }
  if (tau_prior == "exponential"){
    lambda = tau_param$lambda
    if(lambda <= 0){
      cat('You need positive parameters for Exponential distribution \n')
      return(NULL)
    }
    half_adjustment = 1/(2*size_tau)
    discretization = seq(0,1,length.out=size_tau + 1)
    field_values_tau = discretization[1:size_tau] + half_adjustment
    field_values_tau = qexp(field_values_tau, rate=lambda)
    probability_tau = diff(discretization)
  }
  
  prior_field = expand.grid(sigma2=sigma2, range = seq(1, size_range))
  prior_field$prob = rep(1/size_range * (pgamma(upper, alpha,beta) - pgamma(lower,alpha,beta)), times = size_range)
  return(prior_field)
}


#Performs GLS estimation of the samples
glmLite <- function(trend, onlyFormula = FALSE, data, covariance_matrix){
  if (trend == 'cubic'){
    formula = z ~ x*y + I(x^2) + I(y^2) + I(x^3) + I(y^3)
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
  if (onlyFormula == TRUE){
    return(formula)
  }
  mf <- model.frame(formula = formula, data = data);
  X <- model.matrix(attr(mf, "terms"), data = mf);
  Y <- model.response(mf);
  
  #Performing GLS 
  W <- chol2inv(chol(covariance_matrix));
  covariance <- chol2inv(chol(crossprod(X, W %*% X)));
  beta <- tcrossprod(covariance, X) %*% W %*% Y; 
  return(list(coefficients = beta, covar = covariance, formula = formula, X = X, terms = attr(mf, "terms")))
}


#Constructs posterior distribution
posteriorDistribution <- function(correlation_function = 'exponential', glm_object, GRFsigma2 = 1,
                                  sampling_noise, dist_smpl_pred, dist_pred_pred, dist_smpl_smpl,
                                  Xobs, Xpred, samples){
  
  c_p_s <- correlationMatrix(dist_smpl_pred, range = range, correlation_function);
  
  c_p <- correlationMatrix(dist_pred_pred, range = range, correlation_function);

  c_s <- correlationMatrix(dist_smpl_smpl, range = range, correlation_function);
  
  fitted_predictions <- Xpred%*%glm_object$coefficients;
  
  XpCOVtXo <- tcrossprod(Xpred %*% glm_object$covar, Xobs); #tcrossprod same as X %*% W %*% t(X), but slightly faster
  XpCOVtXp <- tcrossprod(Xpred %*% glm_object$covar, Xpred);
  XoCOVtXo <- tcrossprod(Xobs  %*% glm_object$covar, Xobs);
  
  #Constructing the different parts of the variance matrix of the conditional multivariate normal
  sigma12 <- XpCOVtXo + c_p_s * GRFsigma2
  sigma22 <- XoCOVtXo + c_s * GRFsigma2 + sampling_noise;
  sigma11 <- XpCOVtXp + c_p * GRFsigma2;
  sigma22inv <- chol2inv(chol(sigma22)) #Inverse by Cholesky, a constant factor faster
  sigma12sigma22 <- sigma12 %*% sigma22inv   
  
  #Fitting predictions and variance 
  #variance = prediction_grid
  prediction <- fitted_predictions + sigma12sigma22 %*% (samples - Xobs%*%glm_object$coefficients)
  variance <- diag( sigma11 - tcrossprod(sigma12sigma22, sigma12) )
  
  #fit = prediction_grid; fit$z = fitted_predictions
  #beforeFit = prediction_grid; beforeFit$z = sigma12sigma22%*%(samples$z - Xobs%*%glm_object$coefficients)
  
  #Saving as matrix for "tighter" plot
  #variance = xyz.to.grid(variance)
  #prediction_grid = xyz.to.grid(prediction_grid)
  #fit = xyz.to.grid(fit)
  #beforeFit = xyz.to.grid(beforeFit)
  
  return(list(prediction = prediction, variance = variance))#, fit = fit, beforeFit = beforeFit))
}
