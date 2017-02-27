#Auxilliary functions for bachelor-thesis

#Constructing square map with adapted axes 
squareMap <- function(map, grid_size){
  square_map = expand.grid(y = seq(1:grid_size), x = seq(1:grid_size))
  square_map$z = as.vector(interp(x = map$x, y = map$y, z = map$z,
                                      xo = seq(1,max(map$x),length.out = grid_size),
                                      yo = seq(1,max(map$y),length.out = grid_size))$z)
  return(square_map)
}

#Constructs an isotropic correlation matrix based on the Exponential correlation function 
correlationMatrix <- function(sample_points, prediction_points, range){
  ncolumns = dim(prediction_points)[1]; nrows = dim(sample_points)[1]
  
  corr_matrix = matrix(numeric(nrows*ncolumns),nrow=nrows, ncol=ncolumns)
  zero_matrix = matrix(numeric(2*nrows),nrow=2)
  
  #Calculating distances, one column of the matrix at the time
  for (i in 1:ncolumns){
    corr_matrix[,i] = sqrt( colSums((t(sample_points) - (zero_matrix + prediction_points[i,]))^2))
  }
  return( Exponential(corr_matrix, range = range) )
}

#Samples n samples with position from the grid
#(To check the values: z[ intersect(which(square_map$y==a),which(square_map$x==b))) ]
gridSampler <- function(n, map, design, measurement_sigma = 0){
  if (design == 'regular'){
    #Designed to have a perfect regular grid, disregarding the actual n
    xincrement = round( (max(map$x) - min(map$x))/round(sqrt(n)))
    yincrement = round( (max(map$y) - min(map$y))/round(sqrt(n)))
    
    samples = expand.grid(x = seq(min(map$x), max(map$x), xincrement), y= seq(min(map$y),max(map$y), yincrement))
    samples$z = map$z[(samples$x-1)*50 + samples$y]
    
  }
  samples$z = samples$z + rnorm(length(samples$z),0, measurement_sigma)    
  return(samples)
}
#(To check the values intersect(which(square_map$y==1),which(square_map$x==34)))



#Performs GLS estimation of the samples
glmLite <- function(trend, data){
  if (trend == 'quadratic'){
    formula = z ~ x*y + I(x^2) + I(y^2)
  }
  else if (trend == 'linear'){
    formula = z ~ x + y 
  }
  else if (trend == 'simple'){
    formula = z ~ 1
  }
  else {cat('Choose either simple, linear or quadratic trend','\n'); return()}
  mf = model.frame(formula = formula, data = data);
  X  = model.matrix(attr(mf, "terms"), data = mf);
  Y  = model.response(mf);
  
  #Performing GLS 
  W = solve(corr_matrix);
  covariance = solve(t(X)%*%W%*%X);
  beta = covariance%*%t(X)%*%W%*%Y;
  return(list(coefficients = beta, covar = covariance, formula = formula, X = X, terms = attr(mf, "terms")))
}