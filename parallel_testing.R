#parallel testing
#Source auxilliary functions
source('~/Documents/Bachelor/auxiliary_function_procedure.R')
source('~/Documents/Bachelor/auxiliary_function_stats.R')
library(foreach)
set.seed(4545)

#TODO:
# - Varians blir stor ved numerisk integrasjon
# - Begynne ?? notere i latex

#Adapting data from matrix form
original_mapping = grid.to.xyz(t(volcano))

n = 30; #smoothness of mapping and prediction.
n_samples = 8^2; #how many samples to take
size_range = 10 #size of range discretization
size_sigma = 10 #size of sigma2 discretization
sampling_noise = 5; #sd of data collecting (on Z(s))
trend = 'cubic' #What trend to use for the regression of GRF


#Constructs map as a matrix
map = reshapeMap(map=original_mapping, grid_size = n, type="square")

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

sampling_noise = samples$noise*diag(length(samples$x))

#Numerically integration over the prior domain
correlation = foreach(i=1:size_range) %dopar% {correlationMatrix(corr_matrix=distance_matrix, range = prior_range[i], correlation_function = "exponential") 
correlationMatrix(corr_matrix=distance_matrix, range = prior_range[i], correlation_function = "exponential")
correlationMatrix(corr_matrix=distance_matrix, range = prior_range[i], correlation_function = "exponential")
correlationMatrix(corr_matrix=dist_smpl_smpl, range = prior_range[i], correlation_function = "exponential")}


glm_objects = list()
covar = list()
system.time(
foreach(i=1:size_range) %dopar% {glmLite(trend=trend, data = samples, correlation_matrix = correlationMatrix(corr_matrix=distance_matrix, range = prior_range[i], correlation_function = "exponential"))}
)
system.time(
for (i in 1:size_range){
  glm_objects[[i]] = glmLite(trend=trend, data = samples, correlation_matrix = correlationMatrix(corr_matrix=distance_matrix, range = prior_range[i], correlation_function = "exponential"))
}
)

system.time(
  for (range in prior_range){

    for (sigma2 in prior_sigma2){
      cat('Iteration for prior pair: (',sigma2,',',range,') \n')
      
      fitted_predictions <- Xpred%*%glm_object$coefficients;
      
      XpBETAtXo <- tcrossprod(Xpred %*% glm_object$covar, Xobs);
      XpBETAtXp <- tcrossprod(Xpred %*% glm_object$covar, Xpred);
      XoBETAtXo <- tcrossprod(Xobs  %*% glm_object$covar, Xobs);
      
      #Constructing the different parts of the variance matrix of the conditional multivariate normal
      sigma12 <- XpBETAtXo + correlation_prediction_observed * GRFsigma2
      sigma22 <- XoBETAtXo + correlation_observed * GRFsigma2 + sampling_noise;
      sigma11 <- XpBETAtXp + correlation_prediction * GRFsigma2;
      sigma22inv <- chol2inv(chol(sigma22))
      sigma12sigma22 <- sigma12 %*% sigma22inv   
      
      #Fitting predictions and variance 
      #variance = prediction_grid
      prediction <- fitted_predictions + sigma12sigma22 %*% (samples$z - Xobs%*%glm_object$coefficients)
      variance <- diag( sigma11 - tcrossprod(sigma12sigma22, sigma12) )
    }
  }
)
points = matrix(numeric(n*n), nrow=n,ncol=n)
points[cbind(samples$x,samples$y)] = 1
zlim = c(min(c(posterior_distribution$prediction$z,original_mapping$z)), max(max(posterior_distribution$prediction$z, original_mapping$z)))
par(mfrow=c(2,2))
image.plot(xyz.to.grid(original_mapping), zlim=zlim, main="Original data")
image.plot(xyz.to.grid(posterior_distribution$prediction),zlim=zlim, main="Prediction")
image.plot(xyz.to.grid(posterior_distribution$variance), main="Estimated variance")
image(points, main="Observed points", col=c("white", "red"))
