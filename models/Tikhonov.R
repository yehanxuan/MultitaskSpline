# Quadratic regularisation 
# A FPCA-based algorithm. 

Tikhonov_Estimate = function(X, Y, rho, beta_true = NULL){
  # Y can be multiple columns. 
  X_c = t(X) - colMeans(X)
  X_c = t(X_c)
  C = t(X_c)%*%X_c/nrow(X)
  Y_c = t(Y) - colMeans(Y)
  Y_c = t(Y_c)
  # Here we need to divide a Mont here due to matrix computation
  # Note the inverse operator, not the inverse matrix
  g = t(X_c)%*%Y_c/nrow(X)
  btilde = solve((C+rho*diag(ncol(X)))/ncol(X), g)
  # The paper does not give, We can also have the estimation of mean 
  alpha = t(colMeans(Y)- colMeans(X)%*%btilde/ncol(X))
  return(list( "alpha" = alpha ,"beta" = btilde))
}

Tikhonov_Selection = function(X, Y, RhoSeq, nFold = 10, Select_Method = "CV",cvMembership = NULL , beta_true = NULL){
  # Y can have multiple columns
  X_c = t(X) - colMeans(X)
  X_c = t(X_c)
  AveMat = matrix(1/nrow(X_c), nrow(X_c), nrow(X_c))
  C = t(X_c)%*%X_c/nrow(X_c)
  Y_c = t(Y) - colMeans(Y)
  Y_c = t(Y_c)
  g = t(X_c)%*%Y_c/nrow(X_c)
  
  nSample = nrow(X_c)
  if (is.null(cvMembership)){
    cvMembership =  getCVPartition(nSample, nFold)
  }
  
  if (Select_Method == "CV") {
    ErrorSeq = rep(1e+10, length(RhoSeq))
    for (i in 1:length(RhoSeq)) {
      testerror = rep(1e+7, nFold)
      rho = RhoSeq[i]
      for (cf in 1:nFold) {
        testIndex = ( cvMembership == cf)
        Xtest = X[testIndex, ]
        Ytest = Y[testIndex, , drop = F]
        Xtrain = X[!testIndex, ]
        Ytrain = Y[!testIndex, , drop = F]
        fit_train = Tikhonov_Estimate(Xtrain, Ytrain, rho)
        betaTrain = fit_train$betahat 
        alphaTrain = fit_train$alphahat
        Ytesthat = Xtest%*%betaTrain/ncol(Xtest) + rep(1, nrow(Xtest)) %*% t(alphaTrain)
        testerror[cf] = sum((Ytesthat - Ytest)^2)/nrow(Ytest)
      }
      ErrorSeq[i] = mean(testerror)
    }
    index = which.min(ErrorSeq) 
    opt_lambda = RhoSeq[index]
    Estimate = Tikhonov_Estimate(X, Y, opt_lambda, beta_true)
    beta = Estimate$beta
    alpha = Estimate$alpha
  }
  return(list("alpha" = alpha, "beta" = beta, "ErrorSeq" = ErrorSeq, "opt_lambda" = opt_lambda))
}
  

Tikhonov_Selection_Comp = function(X, Y, RhoSeq, nFold = 10, Select_Method = "CV", cvMembership = NULL, beta_true = NULL) {
  X_c = t(X) - colMeans(X)
  X_c = t(X_c)
  C = t(X_c)%*%X_c/nrow(X)
  Y_c = t(Y) - colMeans(Y)
  Y_c = t(Y_c)
  g = t(X_c)%*%Y_c/nrow(X)
  nSample = nrow(X)
  
  ErrorMat = matrix(1e10, ncol(Y), length(RhoSeq))
  opt_lambda = rep(0, ncol(Y))
  btildeMat = matrix(0, ncol(X), ncol(Y))
  alphaMat = matrix(0, ncol(Y), 1)
  if (Select_Method == "CV"){
    if (is.null(cvMembership)){
      cvMembership =  getCVPartition(nSample, nFold)
    }
    for (i in 1:ncol(Y)){
      YVec = Y[ , i, drop = F]
      for (j in 1:length(RhoSeq)){
        testerror = rep(1e7, nFold)
        rho = RhoSeq[j]
        for (cf in 1:nFold){
          testIndex = (cvMembership == cf) 
          Xtest = X[testIndex, ]
          Ytest = YVec[testIndex, , drop = F]
          Xtrain = X[!testIndex, ]
          Ytrain = YVec[!testIndex, , drop = F]
          fit_train = Tikhonov_Estimate(Xtrain, Ytrain, rho)
          betaTrain = fit_train$beta
          alphaTrain = fit_train$alpha
          Ytesthat = as.numeric(alphaTrain) + Xtest%*%betaTrain/ncol(Xtest)
          testerror[cf] = mean((Ytest -Ytesthat)^2)
        }
        ErrorMat[i, j] = mean(testerror)
      }
      index = which.min(ErrorMat[i, ])
      opt_lambda[i] = RhoSeq[index]
      Estimate = Tikhonov_Estimate(X, YVec, opt_lambda[i])
      btildeMat[ ,i] = Estimate$beta
      alphaMat[i, ] = Estimate$alpha
    }
  } else if (Select_Method == "GCV") {
    
  }
  return(list("alpha" = alphaMat, "beta" = btildeMat, "ErrorMat" = ErrorMat, "opt_lambda" = opt_lambda))
}

Tikhonov_Estimate_MultiCov = function(X_approx, Y, opt_lambda){
  Mont = ncol(X_approx[[1]])
  betaMat = matrix(0, Mont, ncol(Y))
  alphaMat = matrix(0, ncol(Y), 1)
  for (i in 1:ncol(Y)) {
    YVec = Y[, i, drop = F]
    Xmat = X_approx[[i]]
    fit = Tikhonov_Estimate(Xmat, YVec, opt_lambda[i])
    betaMat[, i] = fit$beta
    alphaMat[i, ] = fit$alpha
  }
  return(list("alpha" = alphaMat, "beta" = betaMat))
}

Tikhonov_MultiCov = function(X_approx, Y, RhoSeq, Select_Method = "CV", nFold = 10, cvMembership = NULL, beta_true = NULL) {
  Mont = ncol(X_approx[[1]])
  betaMat = matrix(0, Mont, ncol(Y))
  alphaMat = matrix(0, ncol(Y), 1)
  opt_lambda = rep(0, ncol(Y))
  for (i in 1:ncol(Y)) {
    YVec = Y[, i, drop = F]
    Xmat = X_approx[[i]]
    fit = Tikhonov_Selection_Comp(Xmat, YVec, RhoSeq, nFold, Select_Method, cvMembership)
    betaMat[, i] = fit$beta
    alphaMat[i, ] = fit$alpha
    opt_lambda[i] = fit$opt_lambda
  }
  nSample = unlist(lapply(X_approx, nrow))
  Yhat = matrix(0, nSample[1], ncol(Y))
  for (i in 1:length(X_approx)) {
    Yhat[, i] = X_approx[[i]] %*% betaMat[, i]/ncol(X_approx[[i]]) + as.numeric(alphaMat[i, ])
  }
  return(list("beta" = betaMat, "alpha" = alphaMat, "opt_lambda" = opt_lambda, "Yhat" = Yhat))
}
