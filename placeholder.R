system.time(
  for (range in prior_range){
    #New correlation_matrix for each range
    correlation_matrix <- correlationMatrix(corr_matrix=s_dist, range = range, correlation_function = "exponential")
    
    #Obtaining GLS fit with updated correlation_matrix
    glm_object <- glmLite(trend = trend, data = samples, correlation_matrix = correlation_matrix)
    covar <- glm_object$covar
    
    for (sigma2 in prior_sigma2){
      cat('Iteration for prior pair: (',sigma2,',',range,') \n')
      
      #Constructing the only dependency on sigma2 in glm_object
      glm_object$covar <- covar*sigma2
      
      #Obtaining joint probability for this particular (sigma2, range)
      probability <- 0.01#prior$prob[intersect(which(prior$sigma2 == sigma2), which(prior$range == range))]
      
      #Predicting the data with glm estimates of trend
      temp_posterior <- posteriorDistribution(glm_object=glm_object, GRFsigma2 = sigma2, sampling_noise = sampling_noise,
                                              dist_smpl_pred = s_p_dist, dist_pred_pred = p_dist,
                                              dist_smpl_smpl = s_dist, Xobs = Xobs, Xpred = Xpred, samples = samples$z)
      
      #Summing the evaluated values
      posterior_distribution$prediction$z <- posterior_distribution$prediction$z + temp_posterior$prediction*probability
      posterior_distribution$variance$z <- posterior_distribution$variance$z + temp_posterior$variance*probability
    }
  }
)