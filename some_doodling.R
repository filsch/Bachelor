#Source auxilliary functions
source('~/Documents/Bachelor/auxiliary_function_procedure.R')
source('~/Documents/Bachelor/auxiliary_function_stats.R')
set.seed(4545)

#TODO:
# - Lage flere funksjoner for oversiktlighet
# - Varians prediksjon er feil -> GLM + Ikke brukt sigma i GLM
# - Begynne ?? notere i latex
# - Finne ut hvordan behandle matrisedimensjoner. Kan alle arve av alle?
# - Legge inn default verdier for stoerrelse. Burde ikke reshape flere ganger. Stick to the size

#Adapting data from matrix form
original_mapping = grid.to.xyz(t(volcano))

n = 50
#Constructs map as a matrix
map = constructSquareMap(original_mapping, n)

#--------------------------------------------------------------------------------------------------------

#Constructing priorfield
size_range = 10
size_sigma = 10
prior = priorField(size_sigma = size_sigma, size_range = size_range) 

  #--------------------------------------------------------------------------------------------------------

#Generating sample from grid by design
samples = gridSampler(n = 50, map=map, design="regular", noise = 1)
#--------------------------------------------------------------------------------------------------------

# Fitting trend w.r.t. covariance matrix and sample position
correlation_sample = cbind(samples$x,samples$y)

#Constructing distance matrix that is only dependent on samples
distance_matrix = constructDistanceMatrix(correlation_sample,correlation_sample)

#Extracting unique range and sigma2 for easy iterations
prior_range = unique(prior$range)
prior_sigma2 = unique(prior$sigma2)

#Initializing posterior_distribution with empty data
posterior_distribution = list(prediction = matrix(numeric(n*n),nrow=n,ncol=n),
                              variance = matrix(numeric(n*n),nrow=n,ncol=n))

#Numerically integration over the prior domain
for (range in prior_range){
  #New correlation_matrix for each range
  correlation_matrix = correlationMatrix(corr_matrix=distance_matrix, range = range, correlation_function = "exponential")

  for (sigma2 in prior_sigma2){
    #Obtaining joint probability for this particular (sigma2, range)
    cat('Iteration for prior pair: (',sigma2,',',range,') \n')
    probability = prior$prob[intersect(which(prior$sigma2 == sigma2), which(prior$range == range))]
    
    #Obtaining GLS fit with updated correlation_matrix
    glm_object = glmLite('quadratic', samples, correlation_matrix, sigma2 = sigma2)
    
    #Predicting the data with glm estimates of trend
    temp_posterior = posteriorDistribution(size_predictionx = n, size_predictiony = n, samples = samples, glm_object = glm_object, GRFsigma2 = sigma2)
    
    posterior_distribution$prediction = posterior_distribution$prediction + temp_posterior$prediction*probability
    posterior_distribution$variance = posterior_distribution$variance + temp_posterior$variance*probability
  }
}
image.plot(xyz.to.grid(original_mapping), zlim=c(50,200))
#image.plot(posterior_distribution$fit)
#image.plot(posterior_distribution$beforeFit)
par(mfrow=c(1,2))
image.plot(posterior_distribution$prediction,zlim=c(50,200))
image.plot(posterior_distribution$variance)
par(mfrow=c(1,1))
plot(samples$x,samples$y, main="Observed points")

#--------------------------------------------------------------------------------------------------------

sigma2=1
prediction_grid$z <- NULL
prediction_grid = expand.grid(x = seq(1,50), y = seq(1,50))

#KRIGING:
#Computing the predictor:
posterior_kriger = krige.conv(coords=cbind(samples$x,samples$y), data=samples$z, locations = prediction_pos,
                              krige=krige.control(type.krige= "OK", trend.d = "2nd", trend.l = "2nd", cov.model="exponential", cov.pars=c(1, 5)));
#Plotting some results
par(mfrow=c(2,2))
plot_values = list(z=matrix(posterior_kriger$predict,nrow=50,ncol=50,byrow=FALSE), x = seq(0,50), y = seq(0,50))
image.plot(plot_values, main="Kriging predictions, no noise in datamodel", zlim=c(80,200))

plot_values$z = matrix(sqrt(posterior_kriger$krige.var),nrow=50,ncol=50,byrow=FALSE) 
image.plot(plot_values, main="Kriging std.dev., no noise on data Y", zlim=c(0,5))

plot_values$z = map$z - matrix(posterior_kriger$predict,nrow=50,ncol=50,byrow=FALSE) 

image.plot(plot_values, main="Kriging vs actual", zlim=c(-30, 30))
plot(samples$x,samples$y, main="Observed points")
