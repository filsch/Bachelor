source('~/Documents/Bachelor/auxiliary_function_procedure.R')
source('~/Documents/Bachelor/auxiliary_function_stats.R')


tempdata = read.table('~/Documents/Bachelor/tempdata.txt', header = TRUE)
sail_lines = read.table('~/Documents/Bachelor/sail_lines.txt')

names(tempdata)[1] <- 'x'; names(tempdata)[2] <- 'y'; names(tempdata)[3] <- 'z'
sail_lines = list(run1 = data.frame(x = tempdata$x[sail_lines$V1], y = tempdata$y[sail_lines$V1]),
                  run2 = data.frame(x = tempdata$x[sail_lines$V2], y = tempdata$y[sail_lines$V2]),
                  run3 = data.frame(x = tempdata$x[sail_lines$V3], y = tempdata$y[sail_lines$V3]),
                  run4 = data.frame(x = tempdata$x[sail_lines$V4], y = tempdata$y[sail_lines$V4]),
                  run5 = data.frame(x = tempdata$x[sail_lines$V5], y = tempdata$y[sail_lines$V5]),
                  run6 = data.frame(x = tempdata$x[sail_lines$V6], y = tempdata$y[sail_lines$V6]),
                  run7 = data.frame(x = tempdata$x[sail_lines$V7], y = tempdata$y[sail_lines$V7]),
                  run8 = data.frame(x = tempdata$x[sail_lines$V8], y = tempdata$y[sail_lines$V8]))
tempdata = data.frame(x = c(1, tempdata$x), y = c(1, tempdata$y), z = c(5.513, tempdata$z))

#--------------
#par(mfrow=c(4,2))
new_tempdata = xyz.to.grid(tempdata)
for (i in 1:8){
  #i = attr(sail_lines,"names")[i]
  runs = new_tempdata*0
  runs[cbind(sail_lines[[i]]$x,sail_lines[[i]]$y)] = runs[cbind(sail_lines[[i]]$x,sail_lines[[i]]$y)] + 1
  image.plot(runs, col=c('black','lightgreen'))
}
#---------------

#Setting initial parameters
nx = 10; ny = 10; trend = 'quadratic'; total_grid_size = 100

#Constructing grid for prediction
prediction_grid = reshapeMap(map=tempdata, grid_size = total_grid_size, type="reduced")  
prediction_xyz = grid.to.xyz(prediction_grid)
original_grid = xyz.to.grid(map=tempdata)
original_xyz = tempdata
posterior_distribution = list(prediction = prediction_xyz, variance = prediction_xyz)
posterior_distribution$prediction$z = posterior_distribution$prediction$z*0
posterior_distribution$variance$z = posterior_distribution$variance$z*0


#Adapting axes on sail_lines. Sampling as well
sail_lines_xyz = list()
sail_lines_grid = list()
for (i in 1:8){
  sail_lines_xyz[[i]] = adaptLines(sail_lines[[i]], prediction_grid, original_grid)
  sail_lines_grid[[i]] = prediction_grid*0
  sail_lines_grid[[i]][cbind(sail_lines_xyz[[i]]$x, sail_lines_xyz[[i]]$y)] = sail_lines_grid[[i]][cbind(sail_lines_xyz[[i]]$x, sail_lines_xyz[[i]]$y)] + 1
  #image.plot(sail_lines_grid[[i]], nlevel=2, col=c('black','blue'))
  
  sail_lines_xyz[[i]]$z = prediction_grid[cbind(sail_lines_xyz[[i]]$x, sail_lines_xyz[[i]]$y)]
}

#Sampling some initial data. Gjøre satelitt dataene over hele eller bare på punkter? Her hele.
#Computing variogram, assuming isotropic relation
#Could add nuggeteffect for variations at miniscale
satelite_samples = gridSampler(nx = nx, ny = ny, map=prediction_grid, design="regular", noise=0)
satelite_samples$z = prediction_grid[cbind(x=satelite_samples$x, y=satelite_samples$y)]

v = variog(coords = cbind(x = satelite_samples$x, y = satelite_samples$y), data = satelite_samples$z, trend = "1st")
vfit = variofit(v, ini.cov.pars = ini, cov.model = "exponential", weights = "cressie")
#Denne trenger en matrise på ini, så den får valgt den som er best. Hvorfor har initialbetingelsene så mye å si?

plot(v, type="b", main = "Variogram for satelite data")
image.plot(prediction_grid, main = "Interpolated data for reduced grid")

#image.plot(satelite_samples, main = "Sampling positions", nlevel=2, col=c("black","green"))
#----
#Her legger vi til funksjoner for å tilpasse rangeverdier
prior_sigma2 = seq(10); prior_range = seq(10)
prior = list(sigma2 = prior_sigma2, range = prior_range, prob = rep(1/100,100))
#----

#Obtaining design matrices for all locations and sampling locations


i=1
#for (i in 1:8){
samples = sail_lines_xyz[[i]]
formula = glmLite(trend = trend, onlyFormula = TRUE)

mf <- model.frame(formula = formula, data = prediction_xyz)
Xpred  = model.matrix(attr(mf, "terms"), data = mf)

#Square sampling, no lines
#samples = gridSampler(nx = nx, ny = ny, map=prediction_grid, design="regular", noise=0)
#samples$z = prediction_grid[cbind(x=satelite_samples$x, y=satelite_samples$y)]
  
mf <- model.frame(formula = formula, data = samples)
Xobs  = model.matrix(attr(mf, "terms"), data = mf)


p_pos = cbind(prediction_xyz$x, prediction_xyz$y)
s_pos = cbind(samples$x, samples$y)

p_dist = constructDistanceMatrix( p_pos, p_pos )
s_dist = constructDistanceMatrix( s_pos, s_pos )
s_p_dist = constructDistanceMatrix( s_pos, p_pos )


#test denne i morgen
noise = 0
len = length(prior_sigma2)
sampling_noise = list()
for (i in 1:len){
  sampling_noise$prior_sigma2[[i]] = noise / prior_sigma2[[i]]
}

#Numerically integration over the prior domain
system.time(
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
    corr22 <- XoCOVtXo + c_s #+ sampling_noise;
    corr11 <- XpCOVtXp + c_p
    fitted_observed <- Xobs%*%glm_object$coefficients
    
    for (sigma2 in prior_sigma2){
      cat('Iteration for prior pair: (',sigma2,',',range,') \n')
      
      #Obtaining joint probability for this particular (sigma2, range)
      probability <- 0.01#prior$prob[intersect(which(prior$sigma2 == sigma2), which(prior$range == range))]
      
      #Constructing the different parts of the variance matrix of the conditional multivariate normal
      corr22inv <- chol2inv(chol(corr22)) #Inverse by Cholesky, a constant factor faster
      corr12corr22inv <- corr12 %*% corr22inv   
      
      #Fitting predictions and variance 
      prediction <- fitted_predictions + corr12corr22inv %*% (samples - fitted_observed)
      variance <- sigma2 * diag( corr11 - tcrossprod(corr12corr22inv, corr12) )
      
      #Summing the evaluated values
      posterior_distribution$prediction$z <- posterior_distribution$prediction$z + temp_posterior$prediction*probability
      posterior_distribution$variance$z <- posterior_distribution$variance$z + temp_posterior$variance*probability
    }
  }
)

