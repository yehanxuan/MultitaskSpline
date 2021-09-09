# Compute Graph Laplacian 
library(matrixcalc)
Graph_Laplacian = function(S, h){
  size = nrow(S)
  K = matrix(0, size, size)
  for (i in 1:size){
    for (j in 1:i){
      # take u = 2
      ## Modify the kernel here to see any difference 09/08/2021
      # K[i, j] =  1/(2*pi*(h^2))*exp(-sum((S[i, ] - S[j, ])^2)/h^2)
      K[i, j] = 1/(2*pi*(h^3)) * exp(-sum((S[i, ] - S[j, ])^2)/(h^2)  )
      K[j, i] = K[i, j] 
    }
  }
  d = rowMeans(K)
  W = d^(-1/2) * K %*% diag(d^(-1/2))
  D = rowSums(W)
  Gamma_un = diag(D) - W
  Gamma_rw = diag(size) - (1/D) * W
  #return(Gamma_rw)
  return(Gamma_un)
}

# Single predictor 
Graph_Spline = function(X, Y, S, mOrder, nKnots, lambda, eta, h, beta_true = NULL, CovMat = NULL){
  Mont = ncol(X)
  grid = seq(1/Mont, 1, 1/Mont)
  tmin = 0
  tmax = 1
  splineBasis = new(orthoSpline, tmin, tmax, mOrder, nKnots)
  basisMat = splineBasis$evalSpline(grid)
  
  Zmat = X %*% t(basisMat)/Mont 
  Zmean = colMeans(Zmat)
  Zmat_c = t(t(Zmat) - Zmean)
  Ymean = colMeans(Y)
  Y_c = t(t(Y) - Ymean) 
  Omega = splineBasis$get_Omega()   # This Omega corresponds to Gamma in the article 
  Omega = Omega/max(abs(Omega))
  K = splineBasis$getDoF()  # spline degree of freedom
  M = ncol(as.matrix(Y))
  nSample = nrow(Zmat)
  
  # Laplacian 
  Gamma = Graph_Laplacian(S, h)
  Cov = t(Zmat_c)%*%Y_c/nSample
  Cor = t(Zmat_c)%*%Zmat_c/nSample
  
  vecCov = vec(Cov)
  p1 = kronecker(diag(1, M), Cor)
  p2 = lambda*kronecker(diag(1,M), Omega)
  #p3 = eta*kronecker(Gamma, diag(1, K))
  p3 = eta*kronecker(Gamma, Cor) # In single task 09/08/2021
  p4 = lambda*eta*kronecker(Gamma, Omega)
  Core = p1 + p2 + p3 + p4
  vecB = solve(Core, vecCov)
  BHat = matrix(vecB, nrow = K)
  alphaHat = Ymean - t(BHat)%*%Zmean
  
  vecZ_c = kronecker(diag(1, M), Zmat_c)
  HatS = vecZ_c %*% solve(Core, t(vecZ_c))
  GCV = (1/(nSample*M))*( norm(Y_c - Zmat_c%*%BHat, "F")^2 )/ ( (1 - sum(diag(HatS))/(nSample*M))^2 )
  
  betaHat = t(basisMat) %*% BHat
  if (!is.null(beta_true)){
    EstError = colSums((betaHat - beta_true)^2)/Mont
    PredError = diag(PredErr_true(CovMat, betaHat, beta_true))
  } else {
    EstError = NULL
    PredError = NULL
  }
  return(list("alpha" = alphaHat, "beta" = betaHat, "EstError" = EstError, "PredError" = PredError, "GCV" = GCV))
}

CV_Graph_Spline = function(X, Y, S, mOrder, nKnots, lambdaSeq, etaSeq, hSeq, nFold = 10, cvMembership = NULL, beta_true = NULL) {
  nSample = nrow(X)
  if (is.null(cvMembership)){
    cvMembership =  getCVPartition(nSample, nFold)
  }
  L1 = length(lambdaSeq)
  L2 = length(etaSeq)
  L3 = length(hSeq)
  ErrorArray = array(1, dim = c(L1, L2, L3))
  for (i in 1:L1){
    for (j in 1:L2){
      for (k in 1:L3){
        testerror = rep(1e+7, nFold)
        for (cf in 1:nFold){
          testIndex = ( cvMembership == cf)
          Xtest = X[testIndex, ]
          Ytest = Y[testIndex, , drop = F]
          Xtrain = X[!testIndex, ]
          Ytrain = Y[!testIndex, , drop = F]
          fit_train = Graph_Spline(Xtrain, Ytrain, S, mOrder, nKnots, lambdaSeq[i],
                                   etaSeq[j], hSeq[k])
          YtestHat = Xtest%*%fit_train$beta/ncol(Xtest) + rep(1, nrow(Xtest)) %*% t(fit_train$alpha)
          testerror[cf] = sum( (Ytest - YtestHat)^2 )/nrow(Ytest)
        }
        ErrorArray[i, j, k] = mean(testerror)
      }
    }
  }
  index = which(ErrorArray == min(ErrorArray), arr.ind = TRUE)[1, ]
  index1 = index[1]
  index2 = index[2]
  index3 = index[3]
  opt_lambda = lambdaSeq[index1]
  opt_eta = etaSeq[index2]
  opt_h = hSeq[index3]
  opt_error = ErrorArray[index1, index2, index3]
  return(list("ErrorArray" = ErrorArray, "opt_error" = opt_error, "opt_lambda" = opt_lambda, "opt_eta" = opt_eta, "opt_h" = opt_h))
}



# Multi-predictor 
Graph_Spline_MultiCov = function(X_approx, Y, S, mOrder, nKnots, lambda, eta, h){
  # h:Bandwidth
  Mont = ncol(X_approx[[1]])
  tSeq = seq(0, 1, length.out = Mont)
  tmin = 0
  tmax = 1
  splineBasis = new(orthoSpline, tmin, tmax, mOrder, nKnots)
  basisMat = splineBasis$evalSpline(tSeq)
  ZList = list()
  for (i in 1:length(X_approx)){
    ZList[[i]] = X_approx[[i]] %*% t(basisMat)/ncol(X_approx[[i]])
  }
  Zmean_List = lapply(ZList, colMeans)
  ZList_c = list()
  for (i in 1:length(ZList)){
    ZList_c[[i]] = t(t(ZList[[i]]) - Zmean_List[[i]])
  }
  
  Ymean = colMeans(Y)
  Y_c = t(t(Y) - Ymean)
  
  Omega = splineBasis$get_Omega()   # This Omega corresponds to Gamma in the article 
  Omega = Omega/max(abs(Omega))
  K = splineBasis$getDoF()  # spline degree of freedom
  M = ncol(as.matrix(Y_c))
  nSample = unlist(lapply(ZList, nrow))
  
  Gamma = Graph_Laplacian(S,h)
  CovList = list()
  for (i in 1:length(ZList)){
    CovList[[i]] = (1/M)*t(ZList_c[[i]])%*%Y_c[ ,i]/nSample[i]
    # Divided by M here
  }
  
  CorList = list()
  for (i in 1:length(ZList)){
    CorList[[i]] = (1/M)*t(ZList_c[[i]])%*%ZList_c[[i]]/nSample[i]
  }
  
  CorSum = Reduce("+", CorList) # 1/(MN) \sum t(Z)Z
  
  p1 = bdiag(CorList)
  p2 = lambda * kronecker(diag(1,M), Omega)
  #p3 = eta*kronecker(Gamma, diag(1, K))
  p3 = eta*kronecker(Gamma, CorSum)
  p4 = lambda*eta*kronecker(Gamma, Omega)
  Core = p1 + p2 + p3 + p4
  
  # Vectorization 
  vecCov = unlist(CovList)
  vecB = solve(Core, vecCov)
  BHat = matrix(vecB, nrow = K)
  alphaHat = c()
  for (i in 1:M){
    alphaHat[i] =  Ymean[i] - Zmean_List[[i]]%*%BHat[,i]
  }
  betaHat = t(basisMat) %*% BHat
  Yhat = matrix(0, nSample[1], M)
  for (i in 1:length(ZList)) {
    Yhat[ ,i]  = ZList[[i]]%*%BHat[ ,i] + as.numeric( alphaHat[i])
  }
  return(list("alpha" = alphaHat, "beta" = betaHat, "B" = BHat, "Yhat" = Yhat ))
}

CV_Graph_Spline_MultiCov = function(X_approx, Y, S, mOrder, nKnots, lambdaSeq, etaSeq, hSeq, nFold = 10, cvMembership = NULL){
  nSample = nrow(Y)
  M = ncol(Y)
  if (is.null(cvMembership)){
    cvMembership =  getCVPartition(nSample, nFold)
  }
  L1 = length(lambdaSeq)
  L2 = length(etaSeq)
  L3 = length(hSeq)
  ErrorArray = array(1e10, dim = c(L1, L2, L3))
  for (i in 1:L1){
    for (j in 1:L2){
      for (k in 1:L3){
        
        testerror = rep(0, nFold)
        for (cf in 1:nFold){
          testIndex = ( cvMembership == cf)
          Xtest_List = list()
          for (t in 1:length(X_approx)){
            Xtest_List[[t]] = X_approx[[t]][testIndex, ]
          }
          Ytest = Y[testIndex, , drop = F]
          for (t in 1:length(X_approx)){
            Xtrain_List[[t]] = X_approx[[t]][!testIndex, ]
          }
          Ytrain = Y[!testIndex, , drop = F]
          
          fit_train = Graph_Spline_MultiCov(Xtrain_List, Ytrain, S, mOrder, nKnots, lambdaSeq[i],
                                            etaSeq[j], hSeq[k])
          for (v in 1:length(Xtest_List)) {
            
            testerror[cf] = testerror[cf] + 
              mean( (Ytest[ , v, drop = F] - Xtest_List[[v]]%*%fit_train$beta[, v]/ncol(Xtest_List[[v]]) - fit_train$alpha[v])^2)
          }
        }
        ErrorArray[i, j, k] = mean(testerror)/M
      }
    }
  }
  index = which(ErrorArray == min(ErrorArray), arr.ind = TRUE)[1, ]
  index1 = index[1]
  index2 = index[2]
  index3 = index[3]
  opt_lambda = lambdaSeq[index1]
  opt_eta = etaSeq[index2]
  opt_h = hSeq[index3]
  opt_error = ErrorArray[index1, index2, index3]
  return(list("ErrorArray" = ErrorArray, "opt_error" = opt_error, "opt_lambda" = opt_lambda, "opt_eta" = opt_eta, "opt_h" = opt_h))
}


Graph_Spline_MultiCov_SelectKnots = function(X_approx, Y, S, mOrder, K.set, etaSeq, hSeq, nFold = 10, cvMembership = NULL) {
  nSample = nrow(Y)
  M = ncol(Y)
  if (is.null(cvMembership)){
    cvMembership =  getCVPartition(nSample, nFold)
  }
  L1 = length(K.set)
  L2 = length(etaSeq)
  L3 = length(hSeq)
  ErrorArray = array(1e10, dim = c(L1, L2, L3))
  for (i in 1:L1) {
    K = K.set[i]
    for (j in 1:L2){
      eta = etaSeq[j]
      for (k in 1:L3) {
        h = hSeq[k]
        testerror = rep(0, nFold)
        for (cf in 1:nFold){
          testIndex = (cvMembership == cf)
          for (t in 1:length(X_approx)){
            Xtest_List[[t]] = X_approx[[t]][testIndex, ]
          }
          Ytest = Y[testIndex, , drop = F]
          for (t in 1:length(X_approx)){
            Xtrain_List[[t]] = X_approx[[t]][!testIndex, ]
          }
          Ytrain = Y[!testIndex, , drop = F]
          fit_train = Graph_Spline_MultiCov(Xtrain_List, Ytrain, S, mOrder, K, lambda = 0, eta, h)
          
          for (v in 1:length(Xtest_List)) {
            
            testerror[cf] = testerror[cf] + 
              mean( (Ytest[ , v, drop = F] - Xtest_List[[v]]%*%fit_train$beta[, v]/ncol(Xtest_List[[v]]) - fit_train$alpha[v])^2)
          }
        }
        ErrorArray[i, j, k] = mean(testerror)/M
      }
    }
  }
  index = which(ErrorArray == min(ErrorArray), arr.ind = TRUE)[1, ]
  index1 = index[1]
  index2 = index[2]
  index3 = index[3]
  opt_K = K.set[index1]
  opt_eta = etaSeq[index2]
  opt_h = hSeq[index3]
  opt_error = ErrorArray[index1, index2, index3]
  return(list("ErrorArray" = ErrorArray, "opt_error" = opt_error, "opt_K" = opt_K, "opt_eta" = opt_eta, "opt_h" = opt_h))
}




Graph_Spline_Wrap = function(X, Y, S, mOrder, nKnots, lambdaSeq, etaSeq, hSeq, nFold = 10,
                             beta_true = NULL, CovMat = NULL) {
  select = CV_Graph_Spline(X, Y, S, mOrder, nKnots, lambdaSeq, etaSeq, hSeq, nFold, beta_true = beta_true)
  lambda = select$opt_lambda
  eta = select$opt_eta
  h = select$opt_h
  fit = Graph_Spline(X, Y, S, mOrder, nKnots, lambda, eta, h, beta_true = beta_true, CovMat = CovMat)
  MSE = fit$EstError
  PredE = diag( PredErr_true(CovMat, fit$beta, beta_true))
  return(list("beta"= fit$beta, "MSE" = MSE, "PredE" = PredE, "optlambda" = lambda, "opteta" = eta, "opth" = h))
}


oneReplicateWrap_Graph = function(seedJ){
  #try({
  eval = oneReplicate_Graph(seedJ)
  #})
  return(eval)
}

oneReplicate_Graph = function(seedJ){
  set.seed(seedJ + repID * 300)
  source("./oneRep/oneRep-Graph.R")
  if (method == "FPCA") {
    fit = FPCA_Selection_Comp(Xmat, Ymat, lambdaSeq, Select_Method, nFold)
    betahat = fit$beta
    MSE = colSums((betahat - beta_true)^2)/Mont 
    PredE = diag(PredErr_true(CovMat, betahat, beta_true))
    lambda = fit$pSelect
    eta = NULL
    h = NULL
  } else if (method == "Tikhonov") {
    fit = Tikhonov_Selection_Comp(Xmat, Ymat, lambdaSeq, nFold, Select_Method)
    betahat = fit$beta
    MSE = colSums((betahat - beta_true)^2)/Mont
    PredE = diag(PredErr_true(CovMat, betahat, beta_true))
    lambda = fit$opt_lambda
    eta = NULL
    h = NULL
  } else if (method == "RKHS") {
    fit = RKHS_Wrap(Xmat, Ymat, lambdaSeq, Select_Method, nFold, beta_true = beta_true, CovMat = CovMat)
    MSE = fit$MSE
    PredE = fit$PredE
    lambda = fit$lambda
    betahat = fit$beta
    eta = NULL
    h = NULL
  } else if (method == "PSpline") {
    fit = Graph_Spline_Wrap(Xmat, Ymat, S, simu_order, simu_knots, lambdaSeq, etaSeq, hSeq,
                      nFold, beta_true = beta_true, CovMat = CovMat)
    MSE = fit$MSE
    PredE = fit$PredE
    lambda = fit$optlambda
    eta = fit$opteta
    h = fit$opth
  } 
  return(list("MSE" = MSE, "PredE" = PredE, "lambda" = lambda, "eta"= eta, "h" =h))
}


oneReplicate_MultiCov_Graph = function(seedJ){
  set.seed(seedJ + repID * 300)
  source("./oneRep/oneRep-Wind.R")
  if (method == "PSpline") {
    select = CV_Graph_Spline_MultiCov(Xtrain_List, Ytrain, S, order, nknots, lambdaSeq, etaSeq, hSeq, nFold)
    fit = Graph_Spline_MultiCov(Xtrain_List, Ytrain, S, order, nknots, lambda = select$opt_lambda,
                                eta = select$opt_eta, h =select$opt_h)
    YtestHat = matrix(0, nrow(Xtest_List[[1]]), ncol(Ymat))
    for (i in 1:length(Xtest_List)) {
      YtestHat[, i] =  Xtest_List[[i]]%*% fit$beta[, i]/ncol(Xtest_List[[i]]) + as.numeric( fit$alpha[i])
    }
    
    lambda = select$opt_lambda
    eta = select$opt_eta
    h = select$opt_h
    MSP = mean( colMeans( (YtestHat - Ytest)^2 ) ) 
  } else if (method == "RKHS") {
    fit = RKHS_Selection_MultiCov(Xtrain_List, Ytrain, lambdaSeq)
    YtestHat = matrix(0, nrow(Xtest_List[[1]]), ncol(Ymat))
    for (i in 1:length(Xtest_List)) {
      YtestHat[, i] = Xtest_List[[i]]%*%fit$beta[, i]/ncol(Xtest_List[[i]]) + as.numeric(fit$alpha[i, ])
    }
    lambda = fit$lambda
    eta = NULL
    h = NULL
    MSP = mean(colMeans( (YtestHat - Ytest)^2 )  )
  } else if (method == "FPCA") {
    fit = FPCA_MultiCov(Xtrain_List, Ytrain, pSeq, Select_Method, nFold)
    YtestHat = matrix(0, nrow(Xtest_List[[1]]), ncol(Ymat) )
    for (i in 1:length(Xtest_List)){
      YtestHat[, i] = Xtest_List[[i]] %*% fit$beta[, i]/ncol(Xtest_List[[i]]) + as.numeric(fit$alpha[i,]) 
    }
    lambda = fit$pSelect
    eta = NULL
    h = NULL
    MSP = mean(colMeans( (YtestHat - Ytest)^2 )  )
  } else if (method == "Tikhonov") {
    
    fit = Tikhonov_MultiCov(Xtrain_List, Ytrain, RhoSeq, Select_Method = "CV", nFold)
    YtestHat = matrix(0, nrow(Xtest_List[[1]]), ncol(Ymat) )
    for (i in 1:length(Xtest_List)) {
      YtestHat[, i] = Xtest_List[[i]]%*%fit$beta[, i]/ncol(Xtest_List[[i]]) + as.numeric(fit$alpha[i,])
    }
    lambda = fit$opt_lambda
    eta = NULL
    h = NULL
    MSP = mean( colMeans( (YtestHat - Ytest)^2) )
  } else if (method == "SelectKnots") {
    Select = Graph_Spline_MultiCov_SelectKnots(Xtrain_List, Ytrain, S, order, K.set, etaSeq, hSeq, 
                                               nFold)
    lambda = Select$opt_K
    # Here the output lambda is optimal K 
    eta = Select$opt_eta
    h = Select$opt_h
    fit = Graph_Spline_MultiCov(Xtrain_List, Ytrain, S, order, Select$opt_K, lambda = 0,
                                Select$opt_eta, Select$opt_h)
    YtestHat = matrix(0, nrow(Xtest_List[[1]]), ncol(Ymat) )
    for (i in 1:length(Xtest_List)) {
      YtestHat[, i] =  Xtest_List[[i]]%*% fit$beta[, i]/ncol(Xtest_List[[i]]) + as.numeric( fit$alpha[i])
    }
    MSP = mean( colMeans( (YtestHat - Ytest)^2) )
  }
  return(list("MSP" = MSP, "lambda" = lambda, "eta"= eta, "h" =h))
}

oneReplicateWrap_MultiCov_Graph = function(seedJ){
  eval = oneReplicate_MultiCov_Graph(seedJ)
  return(eval)
}





