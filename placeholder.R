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
nx = 10; ny = 10; trend = 'linear interaction'; intercept = TRUE; total_grid_size = 50;
covariance_function = "matern"; n_tau=10; n_rho=10; sigma2 = 0.64; dof=1.5;

#Constructing grid for prediction
prediction_grid = reshapeMap(map=tempdata, grid_size = total_grid_size, type="reduced")  
prediction_xyz = grid.to.xyz(t(prediction_grid)); 
prediction_xyz$y = rev(prediction_xyz$y)
prediction_xyz$x = prediction_xyz$x + 1
prediction_xyz$y = prediction_xyz$y + 1

satelite_xyz = grid.to.xyz(t(reshapeMap(map=tempdata, grid_size = total_grid_size, type="reduced")))
satelite_xyz$y = rev(satelite_xyz$y)
satelite_xyz$x = satelite_xyz$x

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

#Computing variogram, assuming isotropic relation
#Could add nuggeteffect for variations at miniscale

v = variog(coords = cbind(x = satelite_xyz$x, y = satelite_xyz$y), data = satelite_xyz$z, trend = "1st") #+ rnorm(n=length(prediction_xyz$z),mean=0,sd=0.5)
plot(v, type="b", main = "Variogram for satelite data")

vfit1 = variofit(v, ini.cov.pars = c(mean(v$v),5), cov.model = covariance_function, weights = "equal", fix.nugget = TRUE, kappa = dof)
rho_beta = vfit1$cov.pars[1]/var(v$v); rho_alpha = rho_beta*vfit1$cov.pars[1]; tau_lambda = 1/vfit1$cov.pars[2]

#Constructing prior
prior_rho = prior(10,"gamma",list(alpha=rho_alpha,beta=rho_beta))
prior_tau = prior(10,"exponential",list(lambda=tau_lambda))

prior_field = expand.grid(rho=prior_rho$values, tau=prior_tau$values)
prior_field$prob = prior_rho$probability * prior_tau$probability

variance=data.frame("Design 1"=numeric(8),"Design 2"=numeric(8),
                    "Design 3"=numeric(8),"Design 4"=numeric(8),
                    "Design 5"=numeric(8),"Design 6"=numeric(8),
                    "Design 7"=numeric(8),"Design 8"=numeric(8))
system.time(
  for(k in 1:8){
    samples =  sail_lines_xyz[[k]]
    samples$z = samples$z + rnorm(dim(samples)[1], mean = 0, sd=sqrt(sigma2))
    
    results = posteriorDistributionIntegration(samples=samples, correlation_function = covariance_function, rho_alpha=rho_alpha,
                                               rho_beta=rho_beta, tau_lambda=tau_lambda,
                                               sampling_noise=sigma2, prediction_xyz=prediction_xyz, prior_tau=prior_tau, 
                                               prior_rho=prior_rho,trend=trend,intercept=intercept, dof=dof)
    
    
    par(mfrow=c(1,3))
    image.plot(results$expected, main=paste("Expected, design ",k))
    points(samples$x, samples$y,col="black",lwd=2,pch="x")
    test = prediction_xyz
    test$z = results$expected$z - prediction_xyz$z
    image.plot(results$sd, main=paste("Std. dev. design ",k))
    image.plot(test, main=paste("Residuals, design ",k))
    
    sum_posterior = 0
    posterior = 0
    for (i in ((1:8)[1:8 != k])){
      added_design = sail_lines_xyz[[i]];
      sum_posterior = 0
      
      total = rbind(samples,added_design)
      total$z <- NULL
      total = total[!duplicated(total),]
      added_design = total[(dim(samples)[1] + 1):dim(total)[1],]
      std=0
      indexes = numeric(0)
      for (j in 1:dim(added_design)[1]){
        ind = which(results$expected$x==added_design$x[j])
        new = ind[which(results$expected$y[results$expected$x==added_design$x[j]] == added_design$y[j])]
        indexes=c(indexes,new)
      }
      
      mu = results$expected$z[indexes]
      C = chol(results$covmatrix[indexes,indexes] + diag(length(indexes))*sigma2)
      Sigma_inv = chol2inv(C)
      det_sigma = det(Sigma_inv)
      for (j in 1:1){
        added_design$z = C%*%rnorm(length(indexes), mean=0, sd=1) + mu
        proposed = rbind(samples, added_design)
        cat(i,",",j,"\n")
        temp = posteriorDistributionIntegration(samples=proposed, correlation_function = covariance_function, rho_alpha=rho_alpha,
                                                rho_beta=rho_beta, tau_lambda=tau_lambda,
                                                sampling_noise=sigma2, prediction_xyz=prediction_xyz, prior_tau=prior_tau, 
                                                prior_rho=prior_rho,trend=trend,intercept=intercept, variance_only=TRUE, dof=dof)
        posterior = exp(-1/2 * (log(det_sigma) + crossprod( (added_design$z - mu), Sigma_inv%*%(added_design$z - mu) ) ) ) 
        cat("Posterior: ", posterior, " k = ",k, " i = ", i, "\n")
        std = std + sum(temp$sd$z)*posterior
        sum_posterior = sum_posterior + posterior
      }
      #image.plot(temp$sd, main=paste("Std. deviation",i))
      #points(samples$x,samples$y,pch="x", cex=1.2, col="black")
      #points(added_design$x,added_design$y, pch="x", cex=1, col="red")
      variance[[k]][i] = std/(sum_posterior*1)
    }
  })

m = 8
new_samples = sail_lines_xyz[[m]] + rnorm(dim(sail_lines_xyz[[m]])[1], mean=0,sd=sqrt(sigma2))
new_samples = rbind(samples,new_samples)

new_results = posteriorDistributionIntegration(samples=new_samples, correlation_function = covariance_function, rho_alpha=rho_alpha,
                                               rho_beta=rho_beta, tau_lambda=tau_lambda,
                                               sampling_noise=sigma2, prediction_xyz=prediction_xyz, prior_tau=prior_tau, 
                                               prior_rho=prior_rho,trend=trend,intercept=intercept)

par(mfrow=c(1,3))
image.plot(new_results$expected, main="Posterior expected")
points(new_samples$x, new_samples$y,col="black",lwd=2,pch="x")
residuals = new_results$expected
residuals$z = residuals$z - prediction_xyz$z
image.plot(residuals, main="Posterior residuals")
image.plot(new_results$sd, main="Posterior std. dev.");

