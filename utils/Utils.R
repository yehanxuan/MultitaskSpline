# GCV function 
GCV = function(hatMat, Y){
  Yhat = hatMat%*%Y
  gcv = 1/length(Y)*sum((Yhat - Y)^2)/(1- sum(diag(hatMat))/length(Y))^2
  return(gcv)
}


#  GSVD Simple function 
gsvd_simple = function(A, B){
  gsvdO = gsvd(A, B)
  Z = gsvd.R(gsvdO) %*% t(gsvdO$Q)
  result = list(alpha = gsvdO$alpha,
                beta = gsvdO$beta,
                U = gsvdO$U,
                Z = Z)
  return(result)
}


# GCV fast 
gcv_fast = function(Y_mat , alpha, beta, U, lambdaVec){
  p = length(alpha)
  gcvValue = rep(0, length(lambdaVec))
  
  for (i in 1:length(lambdaVec)){
    wt = (alpha^2)/(alpha^2 + lambdaVec[i]*(beta^2))
    SY = t(Y_mat)%*%U[ , 1:p]
    SY = wt * t(SY)
    SY= U[, 1:p]%*%SY 
    resid = Y_mat - SY
    df = sum(wt)/length(Y_mat)
    gcvValue[i] = 1/length(Y_mat) * sum(resid^2)/((1 - df)^2)
  }
  return(gcvValue)
}


# Prediction error using true covariance 
PredErr_true = function(C, betahat,beta0){
  error =  (t(betahat - beta0)%*%C%*%(betahat - beta0)/(nrow(C)^2))
  return(error)
}