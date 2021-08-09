# RKHS functional linear regression yuan10
library(pracma)
library(rdist)
library(doSNOW)
library(microbenchmark)
library(ggplot2)
library(geigen)

# Compute reproducing kernel 
ComputeKernel = function(Mont){
  grid = seq(1/Mont, 1, 1/Mont)
  DistMat = rdist::pdist(grid)
  K = pracma::bernoulli(2,grid)%*%t(pracma::bernoulli(2,grid))/4 - pracma::bernoulli(4, DistMat)/24
  return(K)
}

# RKHS estimate for single task 
RKHS_Estimate = function(X, Y, lambda, beta_true = NULL, CovMat = NULL) {
  # Y can only be a vector here (single task).
  Mont = ncol(X)
  Y = as.matrix(Y)
  nSample = nrow(Y)
  
  X_mean = matrix(colMeans(X), ncol(X), 1)  
  X_c = t(X) - colMeans(X)
  X_c = t(X_c)
  Y_mean = mean(Y)
  grid = seq(1/Mont, 1, 1/Mont)
  # basis for H0
  Xi = cbind(grid^(0), grid^(1))
  # Compute matrix T
  MatT = X_c%*%Xi/Mont
  K = ComputeKernel(Mont)
  Sigma = X_c%*%K%*%t(X_c)/(Mont^2)
  W = Sigma + lambda*diag(nSample)    #Multiple a term nsample as in paper! No multiple, consider the parameter is too small
  
  invW = solve(W)
  invMid = solve(t(MatT)%*%invW%*%MatT)
  d1 = invMid%*%t(MatT)%*%invW
  tmp = invW%*%MatT
  c1 = invW - tmp%*%invMid%*%t(tmp) 
  
  Mat_tmp = (Xi%*%d1 + (K%*%t(X_c)/Mont)%*%c1) # Modify X_c
  betahat = Mat_tmp%*%Y
  hatMat = X_c%*%Mat_tmp/Mont        ## Modify the X to X_c
  hatMat = hatMat + matrix(1/nSample, nSample, nSample)
  alphahat = Y_mean - t(X_mean)%*%betahat/Mont
  
  if (!is.null(beta_true)) {
    Est_error = sum( (betahat - beta_true)^2 )/Mont 
    Pred_error = PredErr_true(CovMat, betahat, beta_true)
  } else {
    Est_error = NULL
    Pred_error = NULL
    #cat("True beta required!")
  }
  # Compute the GCV of lambda
  gcv = GCV(hatMat, Y)
  
  return(list("alpha" = alphahat, "beta" = betahat, "Est_error" = Est_error, "Pred_error" = Pred_error, "GCV" = gcv))
}

RKHS_Selection = function(X, Y, lambdaSeq, Select_Method = "CV", nFold = 10, cvMembership = NULL, beta_true = NULL, CovMat = NULL){
  if (Select_Method == "CV") {
    nSample = nrow(X)
    if (is.null(cvMembership)){
      cvMembership =  getCVPartition(nSample, nFold)
    }
    ErrorVec = rep(1e+10, length(lambdaSeq))
    for (i in 1:length(lambdaSeq)) {
      testerror = rep(1e7, nFold) 
      lambda = lambdaSeq[i]
      for (cf in 1:nFold) {
        testIndex = ( cvMembership == cf)
        Xtest = X[testIndex, ]
        Ytest = Y[testIndex, , drop = F]
        Xtrain = X[!testIndex, ]
        Ytrain = Y[!testIndex, , drop = F]
        fit_train = RKHS_Estimate(Xtrain, Ytrain, lambda)
        betaTrain = fit_train$beta
        alphaTrain = fit_train$alpha
        YtestHat = as.numeric(alphaTrain) + Xtest%*%betaTrain/ncol(Xtest)
        testerror[cf] =  sum( (YtestHat - Ytest)^2)/nrow(Ytest)
      }
      ErrorVec[i] = mean(testerror)
    }
    index = which.min(ErrorVec)
    opt_lambda = lambdaSeq[index]
    if (!is.null(beta_true)) {
      fit = RKHS_Estimate(X, Y, opt_lambda, beta_true, CovMat)
      opt_Est = fit$Est_error
      opt_Pred = fit$Pred_error
    } else {
      fit = RKHS_Estimate(X, Y, opt_lambda)
      opt_Est = fit$Est_error
      opt_Pred = fit$Pred_error
    }
    opt_beta = fit$beta
    opt_alpha = fit$alpha
  }
  
  return(list("lambda" = opt_lambda, "alpha" = opt_alpha, "beta" = opt_beta, "Est_error" = opt_Est, "Pred_error" = opt_Pred))
}


RKHS_Wrap = function(Xmat, Ymat, lambdaSeq, Select_Method = "CV", nFold = 10, cvMembership = NULL, beta_true = NULL, CovMat = NULL) {
  betaHatMat = matrix(0, ncol(Xmat), ncol(Ymat))
  alphaMat = matrix(0, ncol(Ymat), 1)
  EstError = rep(0, ncol(Ymat))
  PredError = rep(0, ncol(Ymat))
  opt_lambda = rep(0, ncol(Ymat))
  nSample = nrow(Xmat)
  Mont = ncol(Xmat)
  if (is.null(cvMembership)){
    cvMembership = getCVPartition(nSample, nFold)
  }
  for (i in 1:ncol(Ymat)) {
    if (!is.null(beta_true)){
      fit = RKHS_Selection(Xmat, Ymat[, i, drop = F], lambdaSeq, Select_Method, nFold, 
                           cvMembership, beta_true[, i, drop = F], CovMat)
      EstError[i] = fit$Est_error
      PredError[i] = fit$Pred_error
      #PredError[i] = PredErr_true(CovMat, fit$beta,  beta_true[ ,i, drop = F])
      opt_lambda[i] = fit$lambda
      betaHatMat[, i] = fit$beta
      alphaMat[i, ] = fit$alpha
    } else {
      fit = RKHS_Selection(Xmat, Ymat[, i, drop = F], lambdaSeq, Select_Method, nFold,
                           cvMembership)
      opt_lambda[i] = fit$lambda
      betaHatMat[, i] = fit$beta 
      alphaMat[i, ] = fit$alpha
    }
  }
  if (!is.null(beta_true)){
    MSE = EstError
    PredE = PredError
  } else {
    MSE = NULL
    PredE = NULL
  }
  lambda = opt_lambda
  beta = betaHatMat
  alpha = alphaMat
  return(list("MSE" = MSE, "PredE" = PredE, "lambda" = lambda, "beta" = beta, "alpha" = alpha))
}


RKHS_MultiCov = function(X_approx, Y, lambda){
  Mont = ncol(X_approx[[1]])
  betaHatMat = matrix(0, Mont, ncol(Y))
  alphaMat = matrix(0, ncol(Y), 1)
  gcvVec = rep(0, ncol(Y))
  for (i in 1:ncol(Y)){
    YVec = Y[, i, drop = F]
    Xmat = X_approx[[i]]
    fit  = RKHS_Estimate(Xmat, YVec, lambda)
    betaHatMat[, i] = fit$beta
    alphaMat[i, ] = fit$alpha
    gcv[i] = fit$GCV
  }
  return(list("beta" = betaHatMat, "alpha" = alphaMat, "GCV" = gcv))
}


RKHS_Selection_MultiCov = function(X_approx, Y, lambdaSeq, Select_Method = "CV",
                                   nFold = 10, cvMembership = NULL) {
  Mont = ncol(X_approx[[1]])
  betaHatMat = matrix(0, Mont, ncol(Y))
  alphaMat = matrix(0, ncol(Y), 1)
  opt_lambda = rep(0, ncol(Y))
  for (i in 1:ncol(Y)){
    Xmat = X_approx[[i]]
    YVec = Y[, i, drop = F]
    select = RKHS_Selection(Xmat, YVec, Select_Method, nFold, cvMembership)
    opt_lambda[i] = select$lambda
    betaHatMat[, i] = select$beta
    alphaMat[i, ] = select$alpha
  }
  lambda = opt_lambda
  beta = betaHatMat
  alpha = alphaMat
  
  nSample = unlist(lapply(X_approx, nrow))
  
  Yhat = matrix(0, nSample[1], ncol(Y))
  for (i in 1:length(X_approx)) {
    Yhat[, i] = X_approx[[i]]%*%beta[, i]/ncol(X_approx[[i]]) + as.numeric(alpha[i,])
  }
  return(list("lambda" = lambda, "beta" = beta, "alpha" = alpha, "Yhat" = Yhat))
}


