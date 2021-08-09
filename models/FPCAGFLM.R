FPCA_GFLM = function(X, Y, pSeq, TypeVec, beta_true = NULL, CovMat = NULL) {
  X_c = t(X) - colMeans(X)
  X_c = t(X_c)
  C = t(X_c)%*%X_c/nrow(X)
  decomp = eigen(C)
  Phi = decomp$vectors
  Phi = sqrt(nrow(Phi))*Phi  
  eigenvalues = (decomp$values)/nrow(Phi)
  Z_c = (X_c %*% Phi)/nrow(Phi)
  
  Mont = ncol(X)
  betaMat = matrix(0, Mont, ncol(Y))
  bSelect = list()
  pSelect = rep(0, ncol(Y))
  for (i in 1:ncol(Y)){
    YVec = Y[ , i, drop = F]
    AICVec  = c()
    bmatAIC = list()
    for(j in 1:length(pSeq)){
      p = pSeq[j]
      Z_trun = Z_c[ , 1:p]
      Z_trun = cbind(1, Z_trun)
      if (TypeVec[i] == 0){
        b_est = solve(t(Z_trun)%*%Z_trun, t(Z_trun)%*%YVec)
        eta = Z_trun%*%b_est
        AIC = sum(-YVec*eta + eta^2/2) + 2*p
      } else if (TypeVec[i] == 1){
        tol = 1
        b_init = rep(1, ncol(Z_trun))
        mu = as.vector(exp(Z_trun%*%b_init))
        P = mu/(1+mu)
        while (tol > 1e-3){
          W = P*(1-P)
          fprime = t(Z_trun)%*%(-YVec + P)
          Hessian = t(Z_trun)%*%(W*Z_trun)
          b_est = b_init - 0.1*solve(Hessian, fprime)
          tol = sqrt( sum( (b_est - b_init)^2 ) )
          mu = as.vector(exp(Z_trun%*%b_est))
          P = mu/(1+mu)
          b_init = b_est
        }
        eta = Z_trun%*%b_est
        AIC = sum(-YVec*eta + log(1+exp(eta))) + 2*p
      } else if (TypeVec[i] == 2){
        tol = 1
        b_init = rep(1, ncol(Z_trun))
        mu = as.vector(exp(Z_trun%*%b_init))
        while (tol > 1e-3){
          fprime = t(Z_trun) %*% (-YVec + mu)
          Hessian = t(Z_trun) %*% (mu * Z_trun)
          b_est = b_init - 0.1*solve(Hessian, fprime)
          tol = sqrt( sum( (b_est - b_init)^2 ) )
          mu = as.vector(exp(Z_trun%*%b_est))
          b_init = b_est
        }
        eta = Z_trun%*%b_est
        AIC = sum(-YVec*eta + exp(eta)) + 2*p
      }
      AICVec = c(AICVec, AIC)
      bmatAIC[[j]] = b_est
    }
    pSelect[i] = pSeq[which.min(AICVec)]
    bSelect[[i]] = bmatAIC[[which.min(AICVec)]]
    betaMat[, i] = cbind(1, Phi[, 1:pSelect[i]])%*%bSelect[[i]]
  }
  
  if (!is.null(beta_true)){
    MSE = colSums((betaMat - beta_true)^2)/Mont
    if (!is.null(CovMat)){ 
      PredE = diag(PredErr_true(CovMat, betaMat, beta_true))
    }
  }
  return(list("beta" = betaMat, "MSE" = MSE, "PredE" = PredE, "pSelect" = pSelect))
}

oneReplicate_GFLMFPCA = function(seedJ){
  set.seed(seedJ + repID * 300)
  source("OneRep/oneRep-GLMFPCA.R")
  fit = FPCA_GFLM(Xmat, Ymat, pSeq, TypeVec, beta_true, CovMat)
  #fit = GLMFPCA(Xmat, Ymat, pSeq, Mont, TypeVec, beta_true, CovMat)
  MSE = fit$MSE
  PredE = fit$PredE
  return(list("MSE" = MSE, "PredE" = PredE, "lambda" = NULL, "rank" = NULL))
}


oneReplicateWrap_GFLMFPCA = function(seedJ){
  try({
    eval = oneReplicate_GLMFPCA(seedJ)
  })
  return(eval)
}



