# GFLM Spline 
GLMReduceRankSolve = function(Ymat, Zmat, Intercept, InitA, InitD, lambda, optParam, cvParam, Omega, TypeVec){
  GLMReduceSolveObj = new(algorithmGLMReduce, Ymat, Zmat)
  GLMReduceSolveObj$initAD(Intercept, InitA, InitD)
  if(!is.null(optParam))
    GLMReduceSolveObj$setParameters(optParam)
  if(!is.null(lambda))
    GLMReduceSolveObj$setPenalty(lambda)
  if(!is.null(cvParam)){
    GLMReduceSolveObj$setCVFold(cvParam$cvMembership)
    GLMReduceSolveObj$activateCV(cvParam$cf)
  }
  
  GLMReduceSolveObj$setDataType(TypeVec) ####
  GLMReduceSolveObj$set_penaltyMatrix(Omega)
  GLMReduceSolveObj$run()
  alpha_opt = GLMReduceSolveObj$get_alpha()
  Aopt = GLMReduceSolveObj$get_A()
  Dopt = GLMReduceSolveObj$get_D()
  
  if (!is.null(cvParam)){
    oobError = GLMReduceSolveObj$outOfBagError()
  }
  rm(GLMReduceSolveObj)
  return(list(alpha = alpha_opt, Aopt = Aopt, Dopt = Dopt, oobError = oobError))
}


GLMReduceRank_Selection = function(Ymat, Zmat, M, K, lambdaSeq, rankSeq, optParam, nFold = 10, Omega, TypeVec){
  nSample = nrow(Ymat)
  cvMembership = getCVPartition(nSample, nFold)
  
  L1 = length(rankSeq)
  L2 = length(lambdaSeq)
  
  ErrorMat = matrix(0, L1, L2)
  alphaList = vector("list", length = L1)
  AInitList = vector("list", length = L1)
  DInitList = vector("list", length = L1)
  
  
  for (i in 1:L1){
    Init_alpha = matrix(rnorm(M), M, 1)
    InitA = matrix(rnorm(rankSeq[i]*M), M, rankSeq[i])
    InitD = matrix(rnorm(rankSeq[i]*K), K, rankSeq[i])
    alphaList[[i]] = Init_alpha
    AInitList[[i]] = InitA
    DInitList[[i]] = InitD
    for (j in 1:L2){
      testerror = rep(0, nFold)
      for (cf in 1:nFold){
        cvParam = list(cvMembership = cvMembership, cf = cf)
        fit = GLMReduceRankSolve(Ymat, Zmat, Init_alpha, InitA, InitD, lambdaSeq[j], optParam, cvParam, Omega, TypeVec)
        testerror[cf] = mean(fit$oobError)
      }
      ErrorMat[i, j] = mean(testerror)
    }
  }
  
  index = which(ErrorMat== min(ErrorMat), arr.ind = TRUE)
  index1 = index[1]
  index2 = index[2]
  opt_rank = rankSeq[index1]
  opt_lambda = lambdaSeq[index2]
  opt_error = ErrorMat[index1, index2]
  
  alpha_opt = alphaList[[index1]]
  Aopt = AInitList[[index1]]
  Dopt = DInitList[[index1]]
  fit_final = GLMReduceRankSolve(Ymat, Zmat, alpha_opt, Aopt, Dopt, opt_lambda, optParam, cvParam, Omega, TypeVec)
  alpha_final = fit_final$alpha
  Afinal = fit_final$Aopt
  Dfinal = fit_final$Dopt
  
  return(list( "opt_lambda" = opt_lambda, "opt_rank" = opt_rank,"alpha_opt" = alpha_final, "Aopt" = Afinal, "Dopt" = Dfinal,"ErrorMat" = ErrorMat))
}

Multi_Mixed_ReduceRank = function(Zmat, Ymat, beta_true, CovMat, Mont, order, nknots, nFold, rankSeq, lambdaSeq, TypeVec){
  # GLMReduceOptParams = list(iterMax = 1e4, epsilon = 1e-7, verbose = 1)
  GLMReduceOptParams = list(iterMax = 5e3, epsilon = 1e-4, alphaD = 0.618, betaD = 0.8, gammaD=0.27,verbose = 1 ,lr = 0.5, momentum = 0.618)
  # gamma \in (0, 0.5)
  M = ncol(Ymat)
  K = ncol(Zmat)
  select = GLMReduceRank_Selection(Ymat, Zmat, M, K, lambdaSeq, rankSeq, GLMReduceOptParams, nFold, Omega, TypeVec)
  optrank = select$opt_rank
  alphaInit = select$alpha_opt
  AInit = select$Aopt
  DInit = select$Dopt
  nSample = nrow(Zmat)
  cf = -1
  cvMembership = getCVPartition(nSample, nFold)
  cvParam = list(cvMembership = cvMembership, cf = cf)
  fit = GLMReduceRankSolve(Ymat, Zmat, alphaInit,  AInit, DInit, select$opt_lambda, GLMReduceOptParams, cvParam, Omega, TypeVec)
  alpha_opt = fit$alpha
  Aopt = fit$Aopt
  Dopt = fit$Dopt
  CHat = Dopt%*%t(Aopt)
  beta = t(basisMat) %*% CHat
  MSE = colSums((beta - beta_true)^2)/Mont
  PredE = diag( PredErr_true(CovMat, beta, beta_true))
  
  return(list("alpha" = alpha_opt, "beta" = beta, "MSE" = MSE, "PredE" = PredE, "optlambda" = select$opt_lambda, "optrank" = optrank ))
}

oneReplicateWrap_GFLM = function(seedJ){
  try({
    eval = oneReplicate_GFLM(seedJ)
  })
  return(eval)
}


oneReplicate_GFLM = function(seedJ){
  set.seed(seedJ + repID * 300)
  source("OneRep/oneRep-GLMLowRank.R")
  fit = Multi_Mixed_ReduceRank(Zmat, Ymat, beta_true, CovMat, Mont, simu_order, simu_knots, nFold,
                               rankSeq, lambdaSeq, TypeVec)
  MSE = fit$MSE
  PredE = fit$PredE
  return(list("MSE" = MSE, "PredE" = PredE, "lambda" = fit$optlambda, "rank" = fit$optrank))
}

