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

par(mfrow=c(3,1))
#Adapting data from matrix form
original_mapping = grid.to.xyz(t(volcano))
image.plot(original_mapping)

#Constructs map as a matrix
map = constructSquareMap(original_mapping, 50)

#--------------------------------------------------------------------------------------------------------

#Constructing priorfield
prior = priorField(10,10) 

#--------------------------------------------------------------------------------------------------------

#Generating sample from grid by design
samples = gridSampler(50, map,"regular", noise = 1)
#--------------------------------------------------------------------------------------------------------

# Fitting trend w.r.t. covariance matrix and sample position
correlation_sample = cbind(samples$x,samples$y)

glm_object = glmLite('quadratic', samples,
                     correlationMatrix(correlation_sample, correlation_sample,
                                       range = 5, correlation_function = 'exponential'))

#--------------------------------------------------------------------------------------------------------

#Predicting the data with glm estimates of trend
posterior_distribution = posteriorDistribution(50, 50, samples = samples, glm_est=glm_object)
image.plot(posterior_distribution$fit)
image.plot(posterior_distribution$beforeFit)
#image.plot(posterior_distribution$prediction)
image.plot(posterior_distribution$variance)
image.plot(posterior_distribution$othervariance)
plot(samples$x,samples$y, main="Observed points")
par(mfrow=c(1,1))

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
