#Auxilliary functions for bachelor-thesis, relevant for stats
library(fields); library(akima); library(geoR); library(MASS); library(spatial); library(RSAGA); library(lattice)


#Constructs an isotropic correlation matrix based on the Exponential and Matern correlation function 
correlationMatrix <- function(corr_matrix, range, correlation_function = 0, v=1.5){
  if (correlation_function == "exponential"){
    return(matern(corr_matrix, phi = range, kappa=0.5) )
  }
  else if(correlation_function == "matern"){
    return(matern(corr_matrix, phi = range, kappa=v))
  } else {cat('Choose either exponential or matern','\n'); return()}
}


#Constructs a "priorfield", in this instance using a uniform prior for range and gamma for GRF variance
prior <- function(size, prior = "gamma", param = ...){
  #Centre mean, with some (theoretical) std deviation increments around the mean set by parameters
  if (prior == "gamma"){
    shape = param$alpha 
    rate = param$beta
    if( shape <= 0 || rate <= 0){
      cat('You need positive parameters for Gamma distribution \n')
      return(NULL)
    }
    half_adjustment = 1/(2*size)
    discretization = seq(0,1,length.out=size + 1)
    field_values = discretization[1:size] + half_adjustment
    field_values = qgamma(field_values, shape=shape, rate=rate)
    probability = 1/size
    
    #-----
    #curve(dgamma(x,shape=shape,rate=rate),from=0,to=7)
    #points(field_values_sigma2,dgamma(field_values_sigma2,shape=shape,rate=rate),col="red",pch=4)
  }
  if (prior == "exponential"){
    lambda = param$lambda
    if(lambda <= 0){
      cat('You need positive parameters for Exponential distribution \n')
      return(NULL)
    }
    half_adjustment = 1/(2*size)
    discretization = seq(0,1,length.out=size + 1)
    field_values = discretization[1:size] + half_adjustment
    field_values = qexp(field_values, rate=lambda)
    probability = diff(discretization)
  }
  
  return( list(values = field_values, probability = probability) )
}


#Performs GLS estimation of the samples
glmLite <- function(trend, onlyFormula = FALSE, data, covariance_matrix, intercept = TRUE){
  if (trend == 'cubic'){
    formula = "z ~ x*y + I(x^2) + I(y^2) + I(x^3) + I(y^3)"
  }
  else if (trend == 'quadratic'){
    formula = "z ~ x*y + I(x^2) + I(y^2)"
  }
  else if (trend == 'linear'){
    formula = "z ~ x + y"
  }
  else if (trend == 'linear interaction'){
    formula = "z ~ x*y"
  }
  else if (trend == 'simple'){
    formula = "z ~ 1"
  }
  else {cat('Choose either simple, linear or quadratic trend','\n'); return()}
  if (!intercept){
    formula = formula(paste(formula,"-1"))
  }
  if (onlyFormula){
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
posteriorDistributionIntegration <- function(samples, correlation_function = 'exponential', rho_alpha, rho_beta, tau_lambda,
                                  sampling_noise, prediction_xyz, prior_tau, prior_rho,trend,intercept, variance_only=FALSE,dof=1){
  formula = glmLite(trend = trend, onlyFormula = TRUE, intercept = intercept)
  mf <- model.frame(formula = formula, data = prediction_xyz)
  Xpred  = model.matrix(attr(mf, "terms"), data = mf)
  p_pos = cbind(prediction_xyz$x, prediction_xyz$y)

  dimension_samples = dim(samples)[1];
  
  mf <- model.frame(formula = formula, data = samples)
  Xobs  = model.matrix(attr(mf, "terms"), data = mf)
  
  s_pos = cbind(samples$x, samples$y)
  
  p_dist = rdist( p_pos, p_pos )
  s_dist = rdist( s_pos, s_pos )
  s_p_dist = rdist(p_pos, s_pos)
  
  #Setting data that are independent of range and sigma outside the integration
  data_noise = sampling_noise * diag(dimension_samples)
  #Keeping GLS estimates fixed for easing numerical integration
  rho_hat = rho_alpha / rho_beta;
  tau_hat = 1/tau_lambda;
  correlation_matrix <- correlationMatrix(corr_matrix=s_dist, range = tau_hat, correlation_function = covariance_function,v=dof)
  covariance_samples = rho_hat*correlation_matrix + data_noise
  
  glm_object <- glmLite(trend = trend, data = samples, covariance_matrix = covariance_samples, intercept = intercept)
  covar <- glm_object$covar
  
  #Fitting the quadratic and linear forms with GLS estimates
  fitted_predictions <- Xpred%*%glm_object$coefficients;
  fitted_observed <- Xobs%*%glm_object$coefficients;
  
  #Fitted the expectation of [z | tau, rho]
  diff_zbeta = samples$z - fitted_observed
  
  XpCOVtXo <- tcrossprod(Xpred %*% covar, Xobs); #tcrossprod same as X %*% W %*% t(X), but slightly faster
  XpCOVtXp <- tcrossprod(Xpred %*% covar, Xpred);
  XoCOVtXo <- tcrossprod(Xobs  %*% covar, Xobs);
  if(!variance_only){
  j = 0  
  for (tau in prior_tau$values){
    #Predicting the data with glm estimates of trend
    c_p_s <- correlationMatrix(s_p_dist, range = tau, correlation_function = correlation_function, v=dof);
    c_p <- correlationMatrix(p_dist, range = tau, correlation_function = correlation_function,v=dof);
    c_s <- correlationMatrix(s_dist, range = tau, correlation_function = correlation_function,v=dof);
    for (rho in prior_rho$values){
      j = j + 1
      #cat('Iteration number ', j, ' for prior pair: ( Rho: ',rho,', tau: ',tau,') \n')
      
      #Constructing the different parts of the variance matrix of the posterior conditional multivariate normal
      cov12 <- XpCOVtXo + c_p_s*rho
      cov22inv <- chol2inv(chol(XoCOVtXo + c_s*rho + data_noise)) #Inverse by Cholesky, a constant factor faster
      cov12cov22inv <- cov12 %*% cov22inv   
      
      #Fitting posterior predictions and variance 
      prediction <- fitted_predictions + cov12cov22inv %*% (diff_zbeta)
      variance <- XpCOVtXp + c_p*rho - tcrossprod(cov12cov22inv, cov12)
      
      #Constructing values for distribution of Z(s) | tau, rho, expected already fitted
      log_z_evaluated[j] = -1/2 *(log(det(cov22inv)) + crossprod( (diff_zbeta), cov22inv%*%(diff_zbeta) ) )
      
      
      #Storing the evaluated values
      posterior_expected[[j]] <- prediction
      posterior_variance[[j]] <- variance
      }
  }
  
  expected = prediction_xyz
  expected$z = expected$z*0
  variance = prediction_xyz
  variance$z = variance$z*0
  z_evaluated = exp(log_z_evaluated)
  norm_z = sum(z_evaluated)
  covmatrix = posterior_variance[[1]]*0
  
  for (i in 1:j){
    expected$z = expected$z + posterior_expected[[i]]*z_evaluated[i]
    variance$z = variance$z + diag(posterior_variance[[i]])*z_evaluated[i]
    covmatrix = covmatrix + (posterior_variance[[i]])*z_evaluated[i]
  }
  expected$z = expected$z / norm_z
  variance$z = variance$z / norm_z #Making it standard deviation
  covmatrix = covmatrix / norm_z
  
  fit = prediction_xyz
  fit$z=fitted_predictions
  return(list(expected = expected, sd = variance, fit = fit, covmatrix = covmatrix))
  } else {
    j = 0  
    for (tau in prior_tau$values){
      #Predicting the data with glm estimates of trend
      c_p_s <- correlationMatrix(s_p_dist, range = tau, correlation_function = correlation_function,v=dof);
      c_p <- correlationMatrix(p_dist, range = tau, correlation_function = correlation_function,v=dof);
      c_s <- correlationMatrix(s_dist, range = tau, correlation_function = correlation_function,v=dof);
      for (rho in prior_rho$values){
        j = j + 1
        #cat('Iteration number ', j, ' for prior pair: ( Rho: ',rho,', tau: ',tau,') \n')
        
        #Constructing the different parts of the variance matrix of the posterior conditional multivariate normal
        cov12 <- XpCOVtXo + c_p_s*rho
        cov22inv <- chol2inv(chol(XoCOVtXo + c_s*rho + data_noise)) #Inverse by Cholesky, a constant factor faster
        cov12cov22inv <- cov12 %*% cov22inv   
        
        #Fitting posterior predictions and variance 
        variance <- diag( XpCOVtXp + c_p*rho - tcrossprod(cov12cov22inv, cov12) )
        
        #Constructing values for distribution of Z(s) | tau, rho, expected already fitted
        log_z_evaluated[j] = -1/2 *(log(det(cov22inv)) + crossprod( (diff_zbeta), cov22inv%*%(diff_zbeta) ) )
        
        
        #Storing the evaluated values
        posterior_variance[[j]] <- variance
      }
    }
    
    variance = prediction_xyz
    variance$z = variance$z*0
    z_evaluated = exp(log_z_evaluated)
    norm_z = sum(z_evaluated)
    
    for (i in 1:j){
      variance$z = variance$z + posterior_variance[[i]]*z_evaluated[i]
    }
    variance$z = sqrt((variance$z/norm_z + sampling_noise)) #Making it standard deviation
    
    return(list(sd = variance))
    }
}

posteriorGRF <- function(added_design,samples, correlation_function = 'exponential', rho_alpha, rho_beta, tau_lambda,
                         sampling_noise,trend,intercept,prediction_grid,dof=1){
  
  formula = glmLite(trend = trend, onlyFormula = TRUE, intercept = intercept)
  mf <- model.frame(formula = formula, data = added_design)
  Xpred  = model.matrix(attr(mf, "terms"), data = mf)
  p_pos = cbind(added_design$x, added_design$y)
  
  dimension_samples = dim(samples)[1];
  
  mf <- model.frame(formula = formula, data = samples)
  Xobs  = model.matrix(attr(mf, "terms"), data = mf)
  
  s_pos = cbind(samples$x, samples$y)
  
  p_dist = rdist( p_pos, p_pos )
  s_dist = rdist( s_pos, s_pos )
  s_p_dist = rdist(p_pos, s_pos )
  
  #Setting data that are independent of range and sigma outside the integration
  data_noise = sampling_noise * diag(dimension_samples)
  #Keeping GLS estimates fixed for easing numerical integration
  rho_hat = rho_alpha / rho_beta;
  tau_hat = 1/tau_lambda;
  correlation_matrix <- correlationMatrix(corr_matrix=s_dist, range = tau_hat, correlation_function = correlation_function, v=dof)
  covariance_samples = rho_hat*correlation_matrix + data_noise
  
  glm_object <- glmLite(trend = trend, data = samples, covariance_matrix = covariance_samples, intercept = intercept)
  covar <- glm_object$covar
  
  #Fitting the quadratic and linear forms with GLS estimates
  fitted_predictions <- Xpred%*%glm_object$coefficients;
  fitted_observed <- Xobs%*%glm_object$coefficients;
  
  diff_zbeta = samples$z - fitted_observed
  
  XpCOVtXo <- tcrossprod(Xpred %*% covar, Xobs);
  XpCOVtXp <- tcrossprod(Xpred %*% covar, Xpred);
  XoCOVtXo <- tcrossprod(Xobs  %*% covar, Xobs);
  
  c_p_s <- correlationMatrix(s_p_dist, range = tau_hat, correlation_function = correlation_function, v=dof);
  c_p <- correlationMatrix(p_dist, range = tau_hat, correlation_function = correlation_function,v=dof);
  c_s <- correlationMatrix(s_dist, range = tau_hat, correlation_function = correlation_function,v=dof);

    #Constructing the different parts of the variance matrix of the posterior conditional multivariate normal
    cov12 <- XpCOVtXo + c_p_s*rho_hat
    cov22inv <- chol2inv(chol(XoCOVtXo + c_s*rho_hat + data_noise)) #Inverse by Cholesky, a constant factor faster
    cov12cov22inv <- cov12 %*% cov22inv   
    
    #Fitting posterior predictions and variance 
    expected <- fitted_predictions + cov12cov22inv %*% (diff_zbeta)
    variance <- XpCOVtXp + c_p*rho_hat - tcrossprod(cov12cov22inv, cov12)
    variance_inv = chol2inv(chol(variance))
  
    return(exp(-1/2 * (log(det(variance)) + crossprod( (added_design$z - expected), variance_inv%*%(added_design$z - expected) ) ) ) )
}
