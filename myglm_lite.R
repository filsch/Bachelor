myglm_lite <- function(trend, data, W){
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
  mf = model.frame(formula = formula, data = data)
  X  = model.matrix(attr(mf, "terms"), data = mf)
  Y  = model.response(mf)
  temp = t(X)%*%W
  covariance = solve( temp%*%X)
  beta = covariance%*%temp%*%Y
  return(list(coefficients = beta, covar = covariance))
}