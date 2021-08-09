ReduceRankSpline = function(X, Y, mOrder, nKnots, rank, lambda, beta_true = NULL, CovMat = NULL, InitList = NULL){
  Mont = ncol(X)
  grid = seq(1/Mont, 1, 1/Mont)
  tmin = 0
  tmax = 1
  splineBasis = new(orthoSpline, tmin, tmax, mOrder, nKnots)
  basisMat = splineBasis$evalSpline(grid)
  
  Zmean = t( colMeans(X) %*% t(basisMat)/Mont)
  Zmat_c  = t(t(X) - colMeans(X)) %*% t(basisMat)/Mont
  Zmat = Xmat %*% t(basisMat)/Mont
  Ymean = colMeans(Y)
  Y_c = t(t(Y) - Ymean)
  ## Use optimization algorithm
  Omega = splineBasis$get_Omega()
  Omega = Omega/max(abs(Omega))
  
  K = splineBasis$getDoF()
  M = ncol(as.matrix(Y))
  nSample = nrow(Zmat)
  
  Cor = t(Zmat_c)%*%Y_c
  Cov = t(Zmat_c)%*%Zmat_c/nrow(Zmat_c)
  gcv = NULL
  
  
  if (is.null(InitList)){
    UInit = matrix(rnorm(rank*M), M, rank)
    UInit = qr.Q(qr(UInit))
    
    VInit = matrix(rnorm(rank*rank), rank, rank)
    VInit = qr.Q(qr(VInit))
    
    
    AInit = UInit%*%t(VInit)
    DInit = solve(Cov + lambda*Omega, Cor%*%AInit/nrow(Zmat))
    InitList = list(AInit, DInit)
  }else {
    AInit = InitList[[1]]
    DInit = InitList[[2]]
  }
  
  Ydiff = Y_c - Zmat_c %*% DInit %*% t(AInit)
  obj = norm(Ydiff, "F")^2/nrow(Y) + lambda*sum( diag(t(DInit)%*%Omega%*%DInit) ) # Add square here 
  tol = 1
  iter = 0
  while (tol > 1e-5) {
    SVD = svd(t(Cor)%*%DInit)
    U = SVD$u
    V = SVD$v
    A = U%*%t(V)
    D = solve(Cov + lambda*Omega, Cor%*%A/nrow(Zmat_c))
    Ydiff = Y_c - Zmat_c %*% D %*% t(A)
    tmp = norm(Ydiff, "F")^2/nrow(Y) + lambda*sum( diag(t(D)%*%Omega%*%D) )
    obj = cbind(obj, tmp)
    
    tol = norm(D - DInit, "F")
    DInit = D
    iter = iter + 1
  }
  SFinal = list(A, D)
  CHat = D%*%t(A)
  alphaHat = Ymean - t(CHat)%*%Zmean
  betaHat =  t(basisMat) %*% CHat
  
  if (!is.null(beta_true)){
    EstError = colSums((betaHat - beta_true)^2)/Mont
    X_c = t(X) - colMeans(X)
    X_c = t(X_c)
    PredError = diag(PredErr_true(CovMat, betaHat, beta_true))
  } else {
    EstError = NULL
    PredError = NULL
  }
  return(list("alpha"= alphaHat, "beta" = betaHat, "SFinal" = SFinal, "EstError" = EstError, "PredError" = PredError, "obj" = obj, "GCV" = gcv))
}


Spline_Selection = function(X, Y, mOrder, nKnots, lambdaSeq, rankSeq, nFold = 10, cvMembership = NULL, beta_true = NULL, InitList = NULL){
  nSample = nrow(X)
  if (is.null(cvMembership)){
    cvMembership =  getCVPartition(nSample, nFold)
  }
  
  L1 = length(rankSeq)
  L2 = length(lambdaSeq)
  ErrorMat = matrix(1e10, L1, L2)
  for (i in 1:L1) {
    for (j in 1:L2){
      testerror = rep(1e7, nFold)
      for (cf in 1:nFold){
        testIndex = ( cvMembership == cf)
        Xtest = X[testIndex, ]
        Ytest = Y[testIndex, , drop = F]
        Xtrain = X[!testIndex, ]
        Ytrain = Y[!testIndex, , drop = F]
        fit_train = ReduceRankSpline(Xtrain, Ytrain, mOrder, nKnots, rankSeq[i], lambdaSeq[j])
        betaTrain = fit_train$beta
        alphaTrain = fit_train$alpha
        YtestHat = Xtest %*% betaTrain/ncol(Xtest) + rep(1, nrow(Xtest))%*%t(alphaTrain)
        testerror[cf] = sum((YtestHat - Ytest)^2)/nrow(Ytest)
      }
      ErrorMat[i, j] = mean(testerror)
    }
  }
  index = which(ErrorMat == min(ErrorMat), arr.ind = TRUE)
  index1 = index[1]
  index2 = index[2]
  opt_rank = rankSeq[index1]
  opt_lambda = lambdaSeq[index2]
  opt_error = ErrorMat[index1, index2]
  return(list( "ErrorMat" = ErrorMat, "opt_lambda" = opt_lambda, "opt_rank" = opt_rank))
}

Spline_Select_Knots = function(X, Y, mOrder, r.set, K.set, nFold = 10, 
                               cvMembership = NULL, beta_true = NULL, InitList = NULL) {
  # r.set: set of rank
  # K.set: set of knots to be selected
  tmin = 0
  tmax = 1
  nSample = nrow(Xmat)
  Mont = ncol(Xmat)
  
  if (is.null(cvMembership)){
    cvMembership = getCVPartition(nSample, nFold)
  }
  L1 = length(r.set)
  L2 = length(K.set)
  ErrorMat = matrix(1e+10, L1, L2)
  for (i in 1:L1){
    r = r.set[i]
    for (j in 1:L2){
      K = K.set[j]
      testerror = rep(1e7, nFold)
      for (cf in 1:nFold){
        testIndex = (cvMembership == cf)
        Xtest = Xmat[testIndex, ]
        Ytest = Ymat[testIndex, , drop = F]
        Xtrain = Xmat[!testIndex, ]
        Ytrain = Ymat[!testIndex, , drop = F]
        fit_train = ReduceRankSpline(Xtrain, Ytrain, mOrder, K, r, lambda = 0,
                                     beta_true, CovMat, InitList)
        betaTrain = fit_train$beta
        alphaTrain = fit_train$alpha
        YtestHat = Xtest%*%betaTrain/ncol(Xtest) + rep(1, nrow(Xtest))%*%t(alphaTrain)
        testerror[cf] = sum( (YtestHat - Ytest)^2 )/nrow(Ytest)
      }
      ErrorMat[i, j] = mean(testerror)
    }
  }
  index = which(ErrorMat == min(ErrorMat), arr.ind = TRUE)
  index1 = index[1]
  index2 = index[2]
  opt_rank = r.set[index1]
  opt_K = K.set[index2]
  opt_error = ErrorMat[index1, index2]
  return(list("ErrorMat" = ErrorMat, "opt_lambda" = opt_K, "opt_rank" = opt_rank, "opt_error" = opt_error ))
}



# Wrap function for regression setting 
Multi_Regression = function(Xmat, Ymat, rankSeq, lambdaSeq, method, nFold = 10, Select_Method = "CV",
                 beta_true, CovMat, simu_knots = NULL,simu_order = NULL, K.set = NULL ) {
  if (method == "RKHS") {
    EstError = rep(0,ncol(Ymat))
    PredError = rep(0, ncol(Ymat))
    lambda = rep(0, ncol(Ymat))
    for (i in 1:ncol(Ymat)) {
      fit = RKHS_Selection(Xmat, Ymat[, i, drop = F], lambdaSeq, Select_Method,
                           nFold, beta_true = beta_true[ ,i, drop = F], CovMat = CovMat)
      EstError[i] = fit$Est_error
      PredError[i] = fit$Pred_error
      lambda[i] = fit$lambda
    }
    MSE = EstError
    PredE = PredError
    lambda = lambda 
    rank = NULL
  } else if (method == "Tikhonov") {
    fit = Tikhonov_Selection_Comp(Xmat, Ymat, lambdaSeq, nFold, Select_Method, beta_true = beta_true)
    betaHat = fit$beta
    MSE = colSums((betaHat - beta_true)^2)/ncol(Xmat) 
    PredE = diag(PredErr_true(CovMat, betaHat, beta_true))
    lambda = fit$opt_lambda
    rank = NULL
  } else if (method == "FPCA") {
    # task by task 
    # lambdaSeq should be sequence of integer 
    fit = FPCA_Selection_Comp(Xmat, Ymat, pSeq = lambdaSeq, Select_Method, 
                              nFold)
    betaHat = fit$beta
    MSE = colSums((betaHat - beta_true)^2)/ncol(Xmat)
    PredE = diag(PredErr_true(CovMat, betaHat, beta_true))
    lambda = fit$pSelect
    rank = NULL
  } else if (method == "PSpline"){
    Select = Spline_Selection(Xmat, Ymat, simu_order, simu_knots, 
                              lambdaSeq, rankSeq, nFold, beta_true = beta_true)
    fit = ReduceRankSpline(Xmat, Ymat, simu_order, simu_knots, Select$opt_rank, Select$opt_lambda,
                           beta_true = beta_true, CovMat = CovMat)
    MSE = fit$EstError
    PredE = diag( PredErr_true(CovMat, fit$beta, beta_true))
    lambda = Select$opt_lambda 
    rank = Select$opt_rank
  } else if (method == "SelectKnots") {
    Select = Spline_Select_Knots(Xmat, Ymat, simu_order, rankSeq, K.set, 
                                 nFold)
    opt_K = Select$opt_lambda # optimal K 
    lambda = opt_K
    rank = Select$opt_rank
    fit = ReduceRankSpline(Xmat, Ymat, simu_order, opt_K, rank, lambda = 0, 
                           beta_true = beta_true, CovMat = CovMat)
    betaHat = fit$beta
    MSE = colSums((betaHat - beta_true)^2)/ncol(Xmat)
    PredE = diag( PredErr_true(CovMat, fit$beta, beta_true))
  }
  return(list("MSE" = MSE, "PredE" = PredE, "optlambda" = lambda, "optrank" = rank))
}



oneReplicateWrap_Reg = function(seedJ){
  try({ 
    eval = oneReplicate_Reg(seedJ)
  })
  return(eval)
}


oneReplicate_Reg = function(seedJ){
  set.seed(seedJ + repID * 300)
  source("./OneRep/oneRep-Reg.R")
  if (!is.null(K.set)){
    fit = Multi_Regression(Xmat, Ymat, rankSeq, lambdaSeq, method, nFold, Select_Method,
                           beta_true = beta_true, CovMat, simu_knots, simu_order, K.set)
  } else {
    fit = Multi_Regression(Xmat, Ymat, rankSeq, lambdaSeq, method, nFold, Select_Method,
                           beta_true = beta_true, CovMat, simu_knots, simu_order)
  }
  
  MSE = fit$MSE
  PredE = fit$PredE
  return(list("MSE" = MSE, "PredE" = PredE, "lambda" = fit$optlambda, "rank" = fit$optrank))
}





