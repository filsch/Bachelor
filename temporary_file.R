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