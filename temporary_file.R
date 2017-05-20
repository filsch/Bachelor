for (range in prior_range){
  #New correlation_matrix for each range
  correlation_matrix <- correlationMatrix(corr_matrix=s_dist, range = range, correlation_function = "exponential")
  
  #Obtaining GLS fit with updated correlation_matrix
  glm_object <- glmLite(trend = trend, data = samples, correlation_matrix = correlation_matrix)
  covar <- glm_object$covar
  
  #Predicting the data with glm estimates of trend
  c_p_s <- correlationMatrix(s_p_dist, range = range, correlation_function);
  c_p <- correlationMatrix(p_dist, range = range, correlation_function);
  c_s <- correlationMatrix(s_dist, range = range, correlation_function);
  
  fitted_predictions <- Xpred%*%glm_object$coefficients;
  
  XpCOVtXo <- tcrossprod(Xpred %*% glm_object$covar, Xobs); #tcrossprod same as X %*% W %*% t(X), but slightly faster
  XpCOVtXp <- tcrossprod(Xpred %*% glm_object$covar, Xpred);
  XoCOVtXo <- tcrossprod(Xobs  %*% glm_object$covar, Xobs);
  corr12 <- XpCOVtXo + c_p_s
  corr11 <- XpCOVtXp + c_p
  fitted_observed <- Xobs%*%glm_object$coefficients
  
  for (sigma2 in prior_sigma2){
    cat('Iteration for prior pair: (',sigma2,',',range,') \n')
    
    #Obtaining joint probability for this particular (sigma2, range)
    probability <- 0.01#prior$prob[intersect(which(prior$sigma2 == sigma2), which(prior$range == range))]
    
    #Constructing the different parts of the variance matrix of the conditional multivariate normal
    corr22inv <- chol2inv(chol( (XoCOVtXo + c_s)*sigma2 + sampling_noise)) #Inverse by Cholesky, a constant factor faster
    corr12corr22inv <- corr12 %*% corr22inv   
    
    #Fitting predictions and variance 
    prediction <- fitted_predictions + sigma2*corr12corr22inv %*% (samples - fitted_observed)
    variance <- sigma2 * diag( corr11 - sigma2*tcrossprod(corr12corr22inv, corr12) )
    
    #Summing the evaluated values
    posterior_distribution$prediction$z <- posterior_distribution$prediction$z + temp_posterior$prediction*probability
    posterior_distribution$variance$z <- posterior_distribution$variance$z + temp_posterior$variance*probability
  }
}

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
