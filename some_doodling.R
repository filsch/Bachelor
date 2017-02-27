library(geoR); library(MASS); library(akima); library(fields); library(RSAGA)
set.seed(4545)

#TODO:
# - Lage flere funksjoner for oversiktlighet
# - Gjennomsøke myglm_lite
# - Begynne å notere i latex
# - Koble prior med prediksjonen
# - Undersøke hvorfor 'linear' og 'quadratic' gir tilsynelatende motsatt resultater





#Adapting data from matrix form
original_mapping = grid.to.xyz(t(volcano))
image.plot(original_mapping)

map = squareMap(original_mapping, 50)

#--------------------------------------------------------------------------------------------------------

#Assigning prior values and estimating beta
alpha = 5; beta = 0.5; measurement_sigma = 0.1;
#Centre mean, 0.4*std indrement out from centre (but why? because)
values = round(seq(-5,4)*0.4*sqrt(alpha)/beta + alpha/beta)
lower = c(0, values[1:length(values) - 1] + diff(values)/2); upper = c(values[2:length(values)] - diff(values)/2, Inf)
prior_field = expand.grid(sigma2=values, range = seq(1,10))
prior_field$prob = rep(0.1*(pgamma(upper, alpha,beta) - pgamma(lower,alpha,beta)), times = 10)


#--------------------------------------------------------------------------------------------------------

#Generating sample from grid by design
samples = gridSampler(50, map,"regular")

#--------------------------------------------------------------------------------------------------------
# Fitting trend w.r.t. covariance matrix and sample position
correlation_sample = cbind(grid_samples$x,grid_samples$y)

glm_object = myglm_lite('linear', grid_samples, corr_matrix(correlation_sample, correlation_sample, 5))

#--------------------------------------------------------------------------------------------------------

#Predicting the data with glm estimates of trend
prediction_grid = expand.grid(x = seq(1,50), y = seq(1,50)); prediction_grid$z = numeric(2500)

mf = model.frame(formula = glm_object$formula, data = prediction_grid)
Xs0  = model.matrix(attr(mf, "terms"), data = mf)
Xs = glm_object$X

fitted_glm = Xs0%*%glm_object$coefficients
prediction_pos= cbind(prediction_grid$x, prediction_grid$y)

correlation_sample_and_pred = corr_matrix(sample_correlation_input, prediction_pos, range = 5);
correlation_predictions = corr_matrix(prediction_pos, prediction_pos, range = 5);
correlation_samples = corr_matrix(sample_correlation_input,sample_correlation_input, range = 5);

covar_beta_Xs0_Xs = Xs0%*%glm_object$covar%*%t(Xs);

#Using 1 as placeholder for variance sigmacorrelation_samplescorrelation_samples
sigma21 = correlation_sample_and_pred*1 + covar_beta_Xs0_Xs;
sigma22 = correlation_samples*1;
sigma11 = correlation_predictions*1;

#En eller annen mindre feil her - noe unøyaktig
prediction_grid$z = fitted_glm - (sigma21)%*%solve(sigma22)%*%(grid_samples$z - Xs%*%glm_object$coefficients)
image.plot(prediction_grid)


#--------------------------------------------------------------------------------------------------------
#KRIGING:
#Computing the predictor:
posterior_kriger = krige.conv(coords=cbind(grid_samples$x,grid_samples$y), data=grid_samples$z, locations = prediction_grid,
                              krige=krige.control(type.krige= "OK", trend.d = "2nd", trend.l = "2nd", cov.model="exponential", cov.pars=c(sigma2, range)));
#Plotting some results
par(mfrow=c(2,2))
plot_values = list(z=matrix(posterior_kriger$predict,nrow=50,ncol=50,byrow=FALSE), x = seq(0,50), y = seq(0,50))
image.plot(plot_values, main="Kriging predictions, no noise in datamodel", zlim=c(80,200))

plot_values$z = matrix(sqrt(posterior_kriger$krige.var),nrow=50,ncol=50,byrow=FALSE) 
image.plot(plot_values, main="Kriging std.dev., no noise on data Y", zlim=c(0,5))

plot_values$z = square_mapping$z - matrix(posterior_kriger$predict,nrow=50,ncol=50,byrow=FALSE) 

image.plot(plot_values, main="Kriging vs actual", zlim=c(-30, 30))
plot(grid_samples$x,grid_samples$y, main="Observed points")
