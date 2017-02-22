data = read.table("https://www.math.ntnu.no/emner/TMA4315/2016h/Assignment3/wikimountains.txt",header=TRUE)
#library(plot3D)

myglm = function(formula, data = list(), family, ...){
  # Extract model matrix & responses
  mf = model.frame(formula = formula, data = data)
  X  = model.matrix(attr(mf, "terms"), data = mf)
  y  = model.response(mf)
  terms = attr(mf, "terms")
  if(family == "gaussian"){
    if (attr(terms,"intercept") == 1){
      rss_null <- t(y-mean(y))%*%(y-mean(y))
      p = ncol(X) - 1
    } else{
      rss_null <- t(y)%*%y
      p = ncol(X)
    }

  beta_hat = solve(t(X)%*%X)%*%t(X)%*%y;      residuals = y - X%*%beta_hat;       rss_beta = t(y - X%*%beta_hat)%*%(y - X%*%beta_hat)
  n = nrow(X);                                sigma2 = as.numeric(rss_beta/n);
  beta_hat_cov = solve((t(X)%*%X)) * sigma2;  fitted = as.vector(X %*% beta_hat)
  R_squared = 1 - rss_beta/rss_null;          adj_R_squared = R_squared - (1 - R_squared) * (p)/(n - p)

  est = list(terms = terms, model = mf, coefficients = beta_hat, beta_cov = beta_hat_cov,
             residuals = residuals, n=n, p = p, sigma2=sigma2,
             rss_beta = rss_beta, rss_null = rss_null, fitted = fitted, y = y,
             R_squared = R_squared, adj_R_squared = adj_R_squared, x = X,
             family = "gaussian")
  }
  else if(family == "poisson"){
    offset = model.offset(mf)
    if (is.null(offset) == 1){
      offset = 0
    }
    loglik = function(offset, beta, X, y) {
        return (sum(dpois(y, lambda = exp(offset + X%*%beta),  log=TRUE)))
    }
    beta = optim(par = numeric(dim(X)[2]), fn = loglik, method="BFGS", hessian = TRUE, control = list(fnscale=-1), X = X, offset = offset, y = y)
    if (attr(terms,"intercept") == 1){
    betanull = optim(par = numeric(1), fn = loglik, method="BFGS", hessian = TRUE, control = list(fnscale=-1), X = 1, offset = offset, y = y)
    } else {
      betanull = list(par=0)
    }
    ls = (y * log(y) - y)
    ln = (y * (offset + as.numeric(betanull$par)) - exp(offset + as.numeric(betanull$par)))
    lp = (y * (offset + X %*% beta$par) - exp(offset + X %*% beta$par))

    null_deviances = 2*sum(ls - ln);   residuals = 2*sum(ls - lp)
    res_df = dim(X)[1] - length(beta$par)
    AIC = 2*length(beta$par) - 2*sum(lp - lfactorial(y))
    deviances = sign(y - exp(offset + X %*% beta$par)) * sqrt(2*(ls - lp))

    if (attr(terms,"intercept") == 1){
      null_df = dim(X)[1] - 1
    } else {
      null_df = dim(X)[1]
    }
    yhat = offset + X%*%beta$par
    yres = y - exp(yhat)
      est = list(terms = terms, mf=mf, y = y, x = X, model = mf, offset = offset,
               coefficients = beta$par, beta_cov = solve(-beta$hessian),
               res_df = res_df, null_deviances = null_deviances,
               residuals = residuals, AIC = AIC, family=family,
               deviances = deviances, null_df = null_df, data = data,
               yhat = yhat, yres = yres, iterations = "Not used")
  }
  else if(family == "binomial" || family == "geometric"){
    if (family == "geometric"){
      n = numeric(length(y)) + 1
    } else {
    n = matrix(rowSums(y))
    ratios = matrix(y[,1]/rowSums(y))
    y = matrix(y[,1])
    }
  
    beta = IRLS(n, y, X)
    
    if (attr(terms,"intercept") == 1){
      betanull = IRLS(n, y, numeric(length(y)) + 1)
      null_df = dim(X)[1] - 1
      
    } else {
      betanull = list(par=0)
      null_df = dim(X)[1]
    }
    
    ls = y*log(y) + (n - y)*log(n - y) - n*log(n);        ls[is.na(ls)] = 0
    ln = y*as.numeric(betanull$par) - n*log(1 + exp(as.numeric(betanull$par)))
    lp = y*X %*% beta$par - n * log(1 + exp(X %*% beta$par))

    deviances = sign(ratios - exp(X %*% beta$par) / (1 + exp(X %*% beta$par)) ) * sqrt(2*(ls-lp))
    residuals = 2*sum(ls - lp)
    res_df = dim(X)[1] - length(beta$par)
    null_deviances = 2*sum(ls - ln)
    AIC = 2*length(beta$par) - 2*sum(lp + lfactorial(n) - lfactorial(y) - lfactorial(n-y))

    est = list(terms = terms, mf = mf, y = y, ratios = ratios,
               x = X, model = mf, coefficients = beta$par,
               beta_cov = beta$covmatrix, res_df = res_df,
               null_deviances = null_deviances,
               residuals = residuals, AIC = AIC, family = family,
               deviances = deviances, null_df = null_df, data = data,
               iterations = beta$iterations)

  }

  est$call = match.call()
  est$formula = formula
  class(est) = 'myglm'
  return(est)
}

print.myglm = function(object, ...){
  cat('Call: \n')
  print(object$call)
  cat('\nCoefficients: \n')
  if (object$family == "poisson"){
    frame = data.frame(t(object$coefficients), row.names=''); colnames(frame) <- attr(object$x,"dimnames")[[2]]
    print(frame, right = FALSE, digits = 4)
    cat('\n')
    cat('Degrees of Freedom: ', object$null_df, ' Total (i.e. Null);  ', object$res_df, 'Residual \n')
    cat('Null deviance: ', object$null_deviances, '\n')
    cat('Residual deviance: ', object$residuals, '\t AIC: ', object$AIC)
    }
  else if (object$family == "gaussian") {
    cat(x$coefficients[,])
  }
  else if (object$family == "binomial"){
    frame = data.frame(t(object$coefficients), row.names=''); colnames(frame) <- attr(object$x,"dimnames")[[2]]
    print(frame, right = FALSE, digits = 4)
    cat('\n')
    cat('Degrees of Freedom: ', object$null_df, ' Total (i.e. Null);  ', object$res_df, 'Residual \n')
    cat('Null deviance: ', object$null_deviances, '\n')
    cat('Residual deviance: ', object$residuals, '\t AIC: ', object$AIC)
  }
}

summary.myglm = function(object, ...){
  if (object$family == "gaussian"){
  coeff = object$coefficients; resid = object$residuals; p = object$p; n = object$n
  std_error = z_value = z_score = numeric(length(coeff))
  for (i in 1:length(coeff)){
    z_value[i]   = coeff[i] / sqrt(object$beta_cov[i,i])
    std_error[i] = sqrt(object$beta_cov[i,i])
  }
  z_score   = 1 - pnorm(abs(z_value), mean = 0, sd = 1)
  chi_stat = ((object$rss_null - object$rss_beta)/p) / (object$rss_beta/(n - p))
  chi_p_value = 1 - pchisq(chi_stat, p)
  chi_p_value[chi_p_value < 2e-16] <- '< 2e-16';   z_score[z_score < 2e-16] <- '< 2e-16'

  options(scipen=0)
  cat("Call: \n")
  print(object$call)
  cat('\n Residuals: \n')
  frame = data.frame(min(resid), quantile(resid, 0.25), median(resid),
                     quantile(resid, 0.75), max(resid), row.names = '')
  colnames(frame) <- c('Min','1stQ','Median','3rdQ','Max')
  print(frame,digits = 3)
  cat("\n Coefficients: \n")
  frame = data.frame(coeff, std_error, z_value, z_score)
  colnames(frame) <- c('Estimate','Std. Error','Z-value','Z-score')
  print(frame)
  cat('--- \n')
  cat("Residual standard error: ", sqrt(object$sigma2), " on ", n-p, " degrees of freedom \n")
  cat("Multiple R-squared: ", round(object$R_squared,5), "\t Adjusted R-squared:  ", round(object$adj_R_squared,5), "\n")
  cat("'Chi2'-statistic: ", round(chi_stat,3), " on ", p, "DF,  p-value: ", chi_p_value)

  }
  else if (object$family == "poisson" || object$family == "binomial") {

    coefficients = as.array(object$coefficients); deviances = object$deviances
    std_error = z_value = z_score = numeric(length(coefficients))

    std_error = sqrt(diag(object$beta_cov))
    z_value   = coefficients / std_error

    z_score   = 1 - pnorm(abs(z_value), mean = 0, sd = 1); z_score[z_score == 0] = '< 2e-16'
    row.names(coefficients) = attr(object$x,"dimnames")[[2]]

    cat("Call: \n")
    print(object$call)
    cat('\n Deviance Residuals: \n')
    frame = data.frame(min(deviances), quantile(deviances, 0.25), median(deviances),
                       quantile(deviances, 0.75), max(deviances), row.names = '')
    colnames(frame) <- c('Min','1stQ','Median','3rdQ','Max')
    print(frame, digits = 5)
    cat("\n Coefficients: \n")
    frame = data.frame(coefficients, std_error, z_value, z_score)
    colnames(frame) <- c('Estimate','Std. Error','Z-value','Z-score')
    print(frame, digits=5)
    cat('--- \n')
    cat("Null deviance: ", object$null_deviances, " on ", object$null_df, " degrees of freedom \n")
    cat("Residual deviance: ", object$residuals, " on ", object$res_df, " degrees of freedom \n")
    cat("AIC: ", object$AIC, '\n \n')
    cat('Number of IRLS iterations:', object$iterations)
  }
}

plot.myglm = function(object, ...){
  if (object$family == "gaussian"){
    plot(object$fitted, object$y, xlab='Fitted', ylab='Observed', main='Observed vs fitted values')
  }
  else if (object$family == "poisson"){
    colors <- 1:(max(object$y)+1)
plot(x=object$yhat, y=object$yres, col=colors[object$y+1], main="Residuals v. Fitted")
  }
  else if (object$family == "binomial"){
    par = object$coefficients
    mesh = mesh(seq(min(object$mf$height),max(object$mf$height),500),
                seq(min(object$mf$prominence),max(object$mf$prominence),500))
    height = mesh$x
    prominence = mesh$y
    intercept_bool <- attr(object$terms, "intercept")
    if(intercept_bool)
      intercept = par[1] 
      
    fitted_values = exp(intercept + height*par[intercept_bool + 1] + prominence*par[intercept_bool + 2]) /
      (1 + exp(intercept + height*par[intercept_bool + 1] + prominence*par[intercept_bool + 2]))
    surf3D(x = height, y = prominence, z = fitted_values, colkey=FALSE, bty="b2", main="Regressed prob. of successful attempt",
           phi = 20, theta = 130, zlim=c(0,1), xlab='Height',ylab='Prominence',zlab='Probability of success')
  }
}

#HER ER DET EN FEIL. PROEV ANOVA UTEN INTERCEPT
anova.myglm = function(object, ...){
  if (object$family == "gaussian"){
  comp = attr(object$terms, "term.labels")
  dof = Sum_Sq = numeric(length(comp)+1)
  # Name of response
  response = deparse(object$terms[[2]])
  # Fit the sequence of models
  txtFormula = paste(response, "~", sep = "")
  model = list()
  extra <- 0
  for(numComp in 1:length(comp)){
    if(numComp == 1){
      txtFormula = paste(txtFormula, comp[numComp])
    }
    else{
      txtFormula = paste(txtFormula, comp[numComp], sep = "+")
    }
    formula = formula(txtFormula)
    model[[numComp]] = myglm(formula=formula, data = object$model,family="gaussian")

    intercept_bool <- attr(object$terms, "intercept")
    num_levels <- length(levels(model[[numComp]]$model[[comp[numComp]]]))

    if (num_levels == 0){
      dof[numComp] = 1
    } else {
      dof[numComp] = num_levels - 1

      if (extra == 0 && intercept_bool == 0){
        extra = numComp
      }

    }

    if (numComp == 1){
      Sum_Sq[numComp] = object$rss_null - model[[numComp]]$rss_beta
    } else {
      Sum_Sq[numComp] = model[[numComp - 1]]$rss_beta - model[[numComp]]$rss_beta
    }

    factors_combi <- attr(model[[numComp]]$terms,"factors")[,numComp]
    if(sum(factors_combi) > 1){
      factors_combi_bool = as.logical(attr(model[[numComp]]$terms,"factors")[-1,numComp])
      multi_dof <- 1
      for (i in 1:length(factors_combi_bool)){
        if (factors_combi_bool[i] != 0){
          multi_dof = multi_dof * dof[i]
        }
      }
      dof[numComp] = multi_dof
    }
  }

  if (attr(object$terms, "intercept") == 0){
    dof[extra] = dof[extra] + 1
  }
  Sum_Sq[length(Sum_Sq)] = object$rss_null - sum(Sum_Sq);           dof[length(dof)] = nrow(object$x) - ncol(object$x);
  Mean_Sq = Sum_Sq/dof;                                             Chi2_value = Mean_Sq[1:length(Mean_Sq)-1] / Mean_Sq[length(Mean_Sq)];
  Chi2_score = 1 - pchisq(Chi2_value,dof[length(dof) - 1])

  Chi2_score[Chi2_score < 2e-16] <- '< 2e-16'

  frame = data.frame(dof, Sum_Sq, Mean_Sq, c(round(Chi2_value,5), ''), c(Chi2_score,''))
  colnames(frame) <- c('Df','Sum Sq','Mean Sq','Chi2-value','Chi2-score')
  rownames(frame) <- c(comp, "Residuals")
  cat('Analysis of Variance Table\n')
  cat(c('Response: ', response, '\n'), sep='')
  print(frame)
  cat('--- \n')}
  else if (object$family == "poisson" || object$family == "binomial"){
    comp = attr(object$terms, "term.labels")
    dof_resid = deviances = numeric(length(comp)+1)
    dev = dof = numeric(length(comp))

    # Name of response
    response = deparse(object$terms[[2]])
    # Fit the sequence of models
    txtFormula = paste(response, "~", sep = "")
    model = list()
    extra <- 0
    offset <- attr(object$terms,"offset")
    if (!is.null(offset)){
      txtFormula = paste(txtFormula, rownames(attr(object$terms,"factors"))[offset])
    }
    for(numComp in 1:length(comp)){
      if(numComp == 1){
        if(!is.null(offset)){
          txtFormula = paste(txtFormula, comp[numComp], sep = "+")
        } else {
          txtFormula = paste(txtFormula, comp[numComp])
        }
      }
      else{
        txtFormula = paste(txtFormula, comp[numComp], sep = "+")
      }
      formula = formula(txtFormula)

      if (object$family == "poisson" ){
        model[[numComp]] = myglm(formula=formula, data = object$data, family="poisson")
        link = "log"
      } else if (object$family == "binomial"){
        model[[numComp]] = myglm(formula=formula, data = object$data, family="binomial")
        link = "logit"
      }
      intercept_bool <- attr(object$terms, "intercept")
      num_levels <- length(levels(model[[numComp]]$model[[comp[numComp]]]))

      if (num_levels == 0){
        dof[numComp] = 1
      } else {
        dof[numComp] = num_levels - 1
        if (extra == 0 && intercept_bool == 0){
          extra = numComp
        }

      }
      deviances[1] = object$null_deviances
      deviances[numComp + 1] = model[[numComp]]$residuals

      factors_combi <- attr(model[[numComp]]$terms,"factors")[,numComp]
      if(sum(factors_combi) > 1){
        factors_combi_bool = as.logical(attr(model[[numComp]]$terms,"factors")[-1,numComp])
        if (!is.null(offset)){
          factors_combi_bool = factors_combi_bool[2:length(factors_combi_bool)]
        }
        multi_dof <- 1
        print(factors_combi_bool)
        for (i in 1:length(factors_combi_bool)){
          if (factors_combi_bool[i] != 0){
            multi_dof = multi_dof * dof[i]
          }
        }
        dof[numComp] = multi_dof
      }
    }

    if (attr(object$terms, "intercept") == 0){
      dof[extra] = dof[extra] + 1
    }
    dof_resid[1] = object$null_df
    for (i in 1:length(dof)){
      dev[i] = deviances[i] - deviances[i+1]
      dof_resid[1 + i] = dof_resid[i] - dof[i]
    }

    frame = data.frame(c('', dof), c('', round(dev)), dof_resid, deviances)
    colnames(frame) <- c('Df', 'Deviance', 'Resid. Df', 'Resid. Dev')
    rownames(frame) <- c('NULL', comp)

    cat('Analysis of Deviance Table \n \n')
    cat('Model: ', object$family, ' link: ',link, '\n \n')
    cat('Response: ', rownames(attr(object$terms,"factors"))[1], '\n \n')
    cat('Terms added sequentially (first to last) \n \n \n')
    print(frame, digits = 1)
  }
}

IRLS = function(n, y, X){
  iter = 0
  z <- log(y + 0.5) - log(n - y + 0.5)
  beta <- solve(t(X) %*% X) %*% t(X) %*% z
  beta_prev = 10000000
  epsilon = 0.001
  while(sum(abs(beta - beta_prev)) > epsilon) {
    iter = iter + 1
    eta <- X %*% beta;              mu <- (exp(eta) / (1 + exp(eta)))*n
    W <- diag( as.numeric(mu*(n - mu)/n) );     z <- eta + n*(y - mu)/(mu*(n - mu))
    beta_prev <- beta;              beta <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% z
  }
  beta_covmatrix = solve(t(X) %*% W %*% X)
  return (list(par = as.numeric(beta), covmatrix = beta_covmatrix, iterations = iter))
}
