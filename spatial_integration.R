#Source auxilliary functions
source('~/Documents/Bachelor/auxiliary_function_procedure.R')
source('~/Documents/Bachelor/auxiliary_function_stats.R')
set.seed(4545)

#TODO:
# - Varians blir stor ved numerisk integrasjon
# - Begynne ?? notere i latex

#Adapting data from matrix form
original_mapping = grid.to.xyz(t(volcano))

n = 100; #smoothness of mapping and prediction.
n_samples = ; #how many samples to take
size_range = 10 #size of range discretization
size_sigma = 10 #size of sigma2 discretization
sampling_noise = 1; #sd of data collecting (on Z(s))
trend = 'quadratic' #What trend to use for the regression of GRF


#Constructs map as a matrix
map = constructSquareMap(original_mapping, n)

#--------------------------------------------------------------------------------------------------------

#Constructing priorfield
prior = priorField(size_sigma = size_sigma, size_range = size_range) 

#--------------------------------------------------------------------------------------------------------

#Generating sample from grid by design
samples = gridSampler(n = n_samples, map=map, design="regular", noise = sampling_noise)
#--------------------------------------------------------------------------------------------------------

# Fitting trend w.r.t. covariance matrix and sample position
correlation_sample = cbind(samples$x,samples$y)

#Constructing distance matrix that is only dependent on samples
distance_matrix = constructDistanceMatrix(correlation_sample,correlation_sample)

#Extracting unique range and sigma2 for easy iterations
prior_range = unique(prior$range)
prior_sigma2 = unique(prior$sigma2)

#Initializing posterior_distribution with empty data
prediction_grid = expand.grid(x = seq(1,n), y = seq(1,n));
prediction_grid$z = numeric(n*n)
posterior_distribution = list(prediction = prediction_grid,
                              variance = prediction_grid)


formula = glmLite(trend = trend, onlyFormula = TRUE)

mf = model.frame(formula = formula, data = prediction_grid)
Xpred  = model.matrix(attr(mf, "terms"), data = mf)

mf = model.frame(formula = formula, data = samples)
Xobs  = model.matrix(attr(mf, "terms"), data = mf)

prediction_pos = cbind(prediction_grid$x, prediction_grid$y)

sampling_pos = cbind(samples$x,samples$y)
dist_smpl_pred = constructDistanceMatrix(sampling_pos,prediction_pos)
dist_pred_pred = constructDistanceMatrix(prediction_pos,prediction_pos)
dist_smpl_smpl = constructDistanceMatrix(sampling_pos,sampling_pos)



#Numerically integration over the prior domain
for (range in prior_range){
  #New correlation_matrix for each range
  correlation_matrix = correlationMatrix(corr_matrix=distance_matrix, range = range, correlation_function = "exponential")
  #Obtaining GLS fit with updated correlation_matrix
  glm_object = glmLite(trend=trend, data = samples, correlation_matrix = correlation_matrix)
  covar = glm_object$covar
  
  for (sigma2 in prior_sigma2){
    cat('Iteration for prior pair: (',sigma2,',',range,') \n')
    
    #Constructing the only dependency on sigma2 in glm_object
    glm_object$covar = covar/sigma2
    
    #Obtaining joint probability for this particular (sigma2, range)
    probability = prior$prob[intersect(which(prior$sigma2 == sigma2), which(prior$range == range))]
    
    #Predicting the data with glm estimates of trend
    temp_posterior = posteriorDistribution(glm_object=glm_object, GRFsigma2 = sigma2, sampling_noise = samples$noise,
                                           dist_smpl_pred = dist_smpl_pred, dist_pred_pred = dist_pred_pred,
                                           dist_smpl_smpl = dist_smpl_smpl, Xobs = Xobs, Xpred = Xpred)
    
    posterior_distribution$prediction$z = posterior_distribution$prediction$z + temp_posterior$prediction*probability
    posterior_distribution$variance$z = posterior_distribution$variance$z + temp_posterior$variance*probability
  }
}
points = matrix(numeric(n*n), nrow=n,ncol=n)
points[cbind(samples$x,samples$y)] = 1
zlim = c(min(c(posterior_distribution$prediction$z,original_mapping$z)), max(max(posterior_distribution$prediction$z, original_mapping$z)))
par(mfrow=c(2,2))
image.plot(xyz.to.grid(original_mapping), zlim=zlim, "Original data")
image.plot(xyz.to.grid(posterior_distribution$prediction),zlim=zlim, main="Prediction")
image.plot(xyz.to.grid(posterior_distribution$variance), main="Estimated variance")
image(points, main="Observed points", col=c("white", "red"))
