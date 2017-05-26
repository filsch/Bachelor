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

#new_tempdata = xyz.to.grid(tempdata)
#for (i in 1:8){
#  #i = attr(sail_lines,"names")[i]
#  runs = new_tempdata*0
#  runs[cbind(sail_lines[[i]]$x,sail_lines[[i]]$y)] = runs[cbind(sail_lines[[i]]$x,sail_lines[[i]]$y)] + 1
#}

#Setting initial parameters
nx = 10; ny = 10; trend = 'linear interaction'; intercept = TRUE; total_grid_size = 50

#Constructing grid for prediction
prediction_grid = reshapeMap(map=tempdata, grid_size = total_grid_size, type="reduced")  
prediction_xyz = grid.to.xyz(t(prediction_grid)); 
prediction_xyz$y = rev(prediction_xyz$y)

satelite_xyz = grid.to.xyz(t(reshapeMap(map=tempdata, grid_size = total_grid_size, type="reduced")))
satelite_xyz$y = rev(satelite_xyz$y)#*100/16
satelite_xyz$x = satelite_xyz$x#*100/16

original_grid = xyz.to.grid(map=tempdata)
original_xyz = tempdata

posterior_expected = list()
posterior_variance = list()
log_z_evaluated = numeric()


#Adapting axes on sail_lines. Sampling as well
sail_lines_xyz = list()
sail_lines_grid = list()
for (i in 1:8){
  sail_lines_xyz[[i]] = adaptLines(sail_lines[[i]], prediction_grid, original_grid)
  sail_lines_xyz[[i]]$z = prediction_grid[cbind(sail_lines_xyz[[i]]$x, sail_lines_xyz[[i]]$y)]
  sail_lines_grid[[i]] = prediction_grid*0
  sail_lines_grid[[i]][cbind(sail_lines_xyz[[i]]$x, sail_lines_xyz[[i]]$y)] = sail_lines_grid[[i]][cbind(sail_lines_xyz[[i]]$x, sail_lines_xyz[[i]]$y)] + 1
}

#Sampling some initial data. Gj??re satelitt dataene over hele eller bare p?? punkter? Her hele.
#Computing variogram, assuming isotropic relation
#Could add nuggeteffect for variations at miniscale

v = variog(coords = cbind(x = satelite_xyz$x, y = satelite_xyz$y), data = satelite_xyz$z, trend = "1st") #+ rnorm(n=length(prediction_xyz$z),mean=0,sd=0.5)
plot(v, type="b", main = "Variogram for satelite data")

vfit1 = variofit(v, ini.cov.pars = c(mean(v$v),5), cov.model = "exponential", weights = "equal", fix.nugget = TRUE)
#vfit2 = variofit(v, ini.cov.pars = c(mean(v$v),5), cov.model = "matern", kappa=3/2, weights = "cressie", fix.nugget = TRUE)
#vfit3 = variofit(v, ini.cov.pars = c(mean(v$v),10), cov.model = "matern", kappa=5/2, weights = "cressie", fix.nugget = TRUE)
#lines(1:70,vfit$cov.pars[1]-vfit$cov.pars[1]*matern(u=1:70,phi=vfit$cov.pars[2],kappa=1/2), col="red",lwd=2)
#lines(1:70,vfit2$cov.pars[1]-vfit2$cov.pars[1]*matern(u=1:70,phi=vfit2$cov.pars[2],kappa=3/2), col="darkgreen",lwd=2)
#lines(1:70,vfit3$cov.pars[1]-vfit3$cov.pars[1]*matern(u=1:70,phi=vfit3$cov.pars[2],kappa=5/2), col="skyblue",lwd=2)
#legend(x="bottomright",c("Variogram","Expoential","Matern 3/2","Matern 5/2"), lty=c(1,1,1,1), lwd=c(1,2,2,2),col=c("black","red","darkgreen","skyblue"))
rho_beta = vfit1$cov.pars[1]/var(v$v); rho_alpha = rho_beta*vfit1$cov.pars[1]; tau_lambda = 1/vfit1$cov.pars[2]

#Constructing prior
prior_rho = prior(10,"gamma",list(alpha=rho_alpha,beta=rho_beta))
prior_tau = prior(10,"exponential",list(lambda=tau_lambda))

prior_field = expand.grid(rho=prior_rho$values, tau=prior_tau$values)
prior_field$prob = prior_rho$probability * prior_tau$probability

#Sampling noise
sigma2 = 1.5

formula = glmLite(trend = trend, onlyFormula = TRUE, intercept = intercept)
mf <- model.frame(formula = formula, data = prediction_xyz)
Xpred  = model.matrix(attr(mf, "terms"), data = mf)
p_pos = cbind(prediction_xyz$x, prediction_xyz$y)

system.time(
for (k in 1:8){
#samples = gridSampler(nx=10,ny=10,map = prediction_grid, design="regular", noise=sigma2)
samples =  sail_lines_xyz[[k]]
samples$z = samples$z + rnorm(dim(samples)[1], mean = 0, sd=sqrt(sigma2))
dimension = dim(samples)[1];

mf <- model.frame(formula = formula, data = samples)
Xobs  = model.matrix(attr(mf, "terms"), data = mf)

s_pos = cbind(samples$x, samples$y)

p_dist = constructDistanceMatrix( p_pos, p_pos )
s_dist = constructDistanceMatrix( s_pos, s_pos )
s_p_dist = constructDistanceMatrix( s_pos, p_pos )

#Setting data that are independent of range and sigma outside the integration
  data_noise = sigma2 * diag(dimension)
  #Keeping GLS estimates fixed for easing numerical integration
  rho_hat = rho_alpha / rho_beta;
  tau_hat = 1/tau_lambda;
  correlation_matrix <- correlationMatrix(corr_matrix=s_dist, range = tau_hat, correlation_function = "exponential")
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
  
  #Obtaining joint probability, is constant
  probability <- prior_field$prob[1]

j = 0  
  for (tau in prior_tau$values){
    #Predicting the data with glm estimates of trend
    c_p_s <- correlationMatrix(s_p_dist, range = tau, correlation_function = "exponential");
    c_p <- correlationMatrix(p_dist, range = tau, correlation_function = "exponential");
    c_s <- correlationMatrix(s_dist, range = tau, correlation_function = "exponential");
    for (rho in prior_rho$values){
      j = j + 1
      #cat('Iteration number ', j, ' for prior pair: ( Rho: ',rho,', tau: ',tau,') \n')
      
      #Constructing the different parts of the variance matrix of the posterior conditional multivariate normal
      cov12 <- XpCOVtXo + c_p_s*rho
      cov22inv <- chol2inv(chol(XoCOVtXo + c_s*rho + data_noise)) #Inverse by Cholesky, a constant factor faster
      cov12cov22inv <- cov12 %*% cov22inv   
      
      #Fitting posterior predictions and variance 
      prediction <- fitted_predictions + cov12cov22inv %*% (diff_zbeta)
      variance <- diag( XpCOVtXp + c_p*rho - tcrossprod(cov12cov22inv, cov12) )
      
      #Constructing values for distribution of Z(s) | tau, rho, expected already fitted
      log_z_evaluated[j] = -1/2 *(log(det(cov22inv)) + crossprod( (diff_zbeta), cov22inv%*%(diff_zbeta) ) )
      
      
      #Storing the evaluated values
      posterior_expected[[j]] <- prediction
      posterior_variance[[j]] <- variance
      
      #posterior_distribution$prediction$z <- posterior_distribution$prediction$z + prediction*probability
      #posterior_distribution$variance$z <- posterior_distribution$variance$z + variance*probability
    }
  }
 
expected = prediction_xyz
expected$z = expected$z*0
variance = prediction_xyz
variance$z = variance$z*0
z_evaluated = exp(log_z_evaluated)
norm_z = sum(z_evaluated)

for (i in 1:j){
  expected$z = expected$z + posterior_expected[[i]]*z_evaluated[i]
  variance$z = variance$z + posterior_variance[[i]]*z_evaluated[i]
}
expected$z = expected$z / norm_z
variance$z = sqrt(variance$z / norm_z) #Making it standard deviation

#test = prediction_grid*0
#test[cbind(samples$x,samples$y)] = 1
#image.plot(test)
fitted_regression = prediction_xyz
fitted_regression$z = fitted_regression$z*0 + fitted_predictions

par(mfrow=c(1,3))
#image.plot(samples, main="Samples", zlim=c(0,15))

image.plot(expected, main="Expected", zlim=c(3,9))

lines(sail_lines_xyz[[k]]$x, sail_lines_xyz[[k]]$y,col="black",lwd=2)

#lines(sail_lines_xyz[[k+4]]$x, sail_lines_xyz[[k+4]]$y,col="black",lwd=2)
#image.plot(fitted_regression, main="Fitted regression", zlim=c(0,15))
#lines(samples$x,samples$y,col="black", lwd=2)
fitted_regression$z = expected$z - fitted_regression$z
#image.plot(fitted_regression, main="expected - fitted")
#lines(samples$x,samples$y,col="black", lwd=2)

#image.plot(prediction_xyz, main=paste("Original k =",k,sep=" "), zlim=c(3,9));

fitted_regression$z = prediction_xyz$z - expected$z

image.plot(fitted_regression, main="Error")

image.plot(variance, main=paste("Std. dev. k =",k,sep=" "));
#lines(samples$x,samples$y,col="black", lwd=2)
}
)

#Sammenligne med fikserte parametere -> optimerte 
