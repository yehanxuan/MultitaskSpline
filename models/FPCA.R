library(mvtnorm)

# Estimate by FPCA method, use threshold like in Hall06
FPCA_Estimate = function(X, Y, t, beta_true = NULL) {
  # t: threshold
  # Y can be only 1 column or multiple columns
  X_c = t(X) - colMeans(X)
  X_c = t(X_c)
  C = t(X_c)%*%X_c/nrow(X)
  # Singular value decomposition of covariance matrix
  decomp = eigen(C)
  Phi = decomp$vectors
  # normalize 
  #tmp = Phi%*%t(Phi)/ncol(X_c)  #shao hou xiugai 
  Phi = sqrt(nrow(Phi))*Phi
  K = (decomp$values)/nrow(Phi)
  Y_c = t(Y) - colMeans(Y)
  Y_c = t(Y_c)
  g = t(X_c)%*%Y_c/nrow(X)     # T \times M
  g_coeff = (t(Phi)%*%g)/nrow(g)
  
  Ind = K > t
  Ktrun = K[Ind == TRUE]
  m = length(Ktrun)
  if ((m == 0)||(m == 1) ){
    m = 1
    Ktrun = K[1]
    KtrunInv = 1/Ktrun
    bhat = Phi[ ,1:m, drop = F]%*%KtrunInv%*%(g_coeff[1:m, , drop = F])
  } else {
    KtrunInv = 1/Ktrun
    bhat = Phi[ ,1:m, drop = F]%*%diag(KtrunInv)%*%(g_coeff[1:m, , drop = F])
  }  
  
  alphahat = t( colMeans(Y) - colMeans(X)%*%bhat/ncol(X))
  return(list("alpha" = alphahat, "beta" = bhat, "p" = m))
}

FPCA_Selection = function(X, Y, tSeq, Kfold, Select_Method = "CV", cvMembership = NULL, beta_true = NULL){
  nSample = nrow(X)
  if (is.null(cvMembership)){
    cvMembership = getCVPartition(nSample, nFold)
  }
  
  ErrorSeq = rep(1e10, length(tSeq))
  if (method == "CV"){
    for (i in 1:length(tSeq)){
      t = tSeq[i]
      testerror = rep(1e7, nFold)
      for (cf in 1:nFold){
        testIndex = (cvMembership == cf)
        Xtest = X[testIndex, ]
        Ytest = Y[testIndex, , drop = F]
        Xtrain = X[!testIndex, ]
        Ytrain = Y[!testIndex, , drop = F]
        fit_train = FPCA_Estimate(Xtrain, Ytrain, t)
        betaTrain = fit_train$beta
        alphaTrain = fit_train$alpha
        Ytesthat = Xtest%*%betaTrain/ncol(Xtest) + rep(1, nrow(Xtest)) %*% t(alphaTrain)
        testerror[cf] = sum((Ytesthat - Ytest)^2)/nrow(Ytest)
      }
      ErrorSeq[i] = mean(testerror)
    }
    index = which.min(ErrorSeq)
  } else if (method == "AIC") {
    AICVec = c()
    for (i in 1:length(tSeq)){
      fit = FPCA_Estimate(X, Y, tSeq[i])
      eta = rep(1, nrow(X)) %*% t(fit$alpha) + X %*% fit$beta/ncol(X)
      AIC = sum(-Y*eta + eta^2/2) + 2*select$p
      AICVec = c(AICVec, AIC) 
    }
    index = which.min(AICVec)
  }
  topt = tSeq[index]
  Estimate = FPCA_Estimate(X, Y, topt)
  alpha_opt = Estimate$alpha
  beta_opt = Estimate$beta
  opt_error = ErrorSeq[index]
  return(list("alpha" = alpha_opt, "beta" = beta_opt, "ErrorSeq" = ErrorSeq, "topt" = topt))
}

# FPCA method. Picking the number of components. It is the same with using threshold essentially
# Single task FPCA
FPCA_Estimate_Comp = function(X, YVec, p){
  X_c = t(X) - colMeans(X)
  X_c = t(X_c)
  C = t(X_c)%*%X_c/nrow(X)
  decomp = eigen(C)
  Phi = decomp$vectors
  Phi = sqrt(nrow(Phi)) * Phi
  eigenvalues = (decomp$values)/nrow(Phi)
  Z_c = (X_c %*% Phi)/nrow(Phi)
  YVec_c = t(YVec) - colMeans(YVec)
  YVec_c = t(YVec_c)
  
  Z_trun = Z_c[, 1:p, drop = F]
  b_est =  solve(t(Z_trun)%*%Z_trun, t(Z_trun)%*%YVec_c)
  alpha = mean(YVec) - colMeans(Z_trun)%*%b_est
  eta = as.numeric(alpha) + Z_trun%*%b_est 
  beta = Phi[, 1:p]%*%b_est
  return(list("beta" = beta, "alpha" = alpha))
}

FPCA_Selection_Comp = function(X, Y, pSeq, Select_Method = "CV", nFold = 10, cvMembership = NULL){
  Y = as.matrix(Y)
  X_c = t(X) - colMeans(X)
  X_c = t(X_c)
  C = t(X_c)%*%X_c/nrow(X)
  decomp = eigen(C)
  Phi = decomp$vectors
  Phi = sqrt(nrow(Phi)) * Phi
  eigenvalues = (decomp$values)/nrow(Phi)
  Z_c = (X_c %*% Phi)/nrow(Phi)
  Y_c = t(Y) - colMeans(Y)
  Y_c = t(Y_c)
  g = t(X_c)%*%Y_c/nrow(X)
  
  nSample = nrow(X)
  Mont = ncol(X)
  betaMat = matrix(0, Mont, ncol(Y))
  alphaMat = matrix(0, ncol(Y), 1)
  pSelect = rep(0, ncol(Y))
  
  if (is.null(cvMembership)){
    cvMembership =  getCVPartition(nSample, nFold)
  }
  
  if (Select_Method == "CV"){
    ErrorMat = matrix(0, ncol(Y), length(pSeq))
    for (i in 1:ncol(Y)) {
      YVec = Y[ , i, drop = F]
      YVec_c = YVec - colMeans(YVec)
      for (j in 1:length(pSeq)) {
        p = pSeq[j]
        testerror = rep(1e7, nFold)
        for (cf in 1:nFold) {
          testIndex = ( cvMembership == cf) 
          Xtest = X[testIndex, ]
          Ytest = YVec[testIndex, , drop = F]
          Xtrain = X[!testIndex, ]
          Ytrain = YVec[!testIndex, , drop = F]
          fit_train = FPCA_Estimate_Comp(Xtrain, Ytrain, p)
          betaTrain = fit_train$beta
          alphaTrain = fit_train$alpha
          Ytesthat = as.numeric(alphaTrain) + Xtest%*%betaTrain/ncol(Xtest)
          testerror[cf] = sum((Ytest - Ytesthat)^2)/nrow(Ytest)
        }
        ErrorMat[i, j] = mean(testerror)
      }
      index = which.min(ErrorMat[i, ])
      pSelect[i] = pSeq[index]
      Estimate = FPCA_Estimate_Comp(X, YVec, pSeq[index])
      betaMat[, i] = Estimate$beta
      alphaMat[i, ] = Estimate$alpha
    }
  } else if (Select_Method == "AIC") {
    bSelect = list()
    for (i in 1:ncol(Y)){
      YVec = Y[ , i, drop = F]
      YVec_c = YVec - mean(YVec)
      AICVec  = c()
      bmatAIC = list()
      alphaAIC = list()
      for (j in 1:length(pSeq)){
        p = pSeq[j]
        Z_trun = Z_c[ , 1:p, drop = F]
        b_est = solve(t(Z_trun)%*%Z_trun, t(Z_trun)%*%YVec_c)
        alpha = mean(YVec) - colMeans(Z_trun)%*%b_est
        eta = as.numeric(alpha) + Z_trun%*%b_est
        AIC = sum(-YVec*eta + eta^2/2) + 2*p
        AICVec = c(AICVec, AIC)
        bmatAIC[[j]] = b_est
        alphaAIC[[j]] = as.numeric(alpha)
      }
      pSelect[i] = pSeq[which.min(AICVec)]
      bSelect[[i]] = bmatAIC[[which.min(AICVec)]]
      betaMat[,i] = Phi[ , 1:pSelect[i] ] %*% bSelect[[i]]
      alphaMat[i, ] = alphaAIC[[which.min(AICVec)]]
    }
  }
  return(list("beta" = betaMat, "alpha" = alphaMat, "pSelect" = pSelect))
}


FPCA_MultiCov = function(X_approx, Y, pSeq, Select_Method = "CV", nFold = 10, cvMembership = NULL) {
  # X_approx should be a list. Y can be matrix
  Mont = ncol(X_approx[[1]])
  betaMat = matrix(0, Mont, ncol(Y))
  alphaMat = matrix(0, ncol(Y), 1)
  pSelect = rep(0, ncol(Y))
  for (i in 1:ncol(Y)){
    YVec = Y[, i, drop = F]
    Xmat = X_approx[[i]] 
    fit = FPCA_Selection_Comp(Xmat, YVec, pSeq, Select_Method, nFold, cvMembership)
    betaMat[ , i] = fit$beta
    alphaMat[i, ] = fit$alpha
    pSelect[i] = fit$pSelect
  }
  nSample = unlist(lapply(X_approx, nrow))
  Yhat = matrix(0, nSample[1], ncol(Y))
  for (i in 1:length(X_approx)){
    Yhat[, i] = X_approx[[i]] %*% betaMat[, i]/ncol(X_approx[[i]]) + as.numeric(alphaMat[i,]) 
  }
  return(list("beta" = betaMat, "alpha" = alphaMat, "pSelect" = pSelect, "Yhat" = Yhat))
}
