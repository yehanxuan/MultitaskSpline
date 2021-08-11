################
### Real data###
################
### Wind Energy#
################

rm(list = ls())
library(readr)
library(dplyr)
library(data.table)
library(lubridate)
library(FDABasics)
library(Matrix)
library(ggplot2)
library(grid)

args = commandArgs(trailingOnly = TRUE)
method = args[1]

source("./models/MultiRegression.R")
source("./models/FPCA.R")
source("./models/Tikhonov.R")
source("./models/RKHS.R")
source("./models/GraphicModel.R")
source("./utils/Utils.R")
source("./utils/crossvalidation.R")

Wind_location = read_csv("./Dataset/WindEnergy/Wind_farm_location.csv")
ID = Wind_location['ID']

XData = read_csv("./Dataset/WindEnergy/X_List.csv")
XData
YData = read_csv("./Dataset/WindEnergy/YData.csv" )
YData
colnames(YData) <- colnames(XData)

Inter = intersect(unlist(ID), colnames(XData) )
S = Wind_location[Wind_location$ID %in% Inter, ]['coordinates']
S = apply(S,1, as.character)
S =  do.call( rbind, lapply(strsplit(S, split = ','), as.numeric) )

XData = XData[, c('timestamp', Inter)]
YData = YData[, c('timestamp', Inter)]
Ymat = as.matrix(YData[ , -1])
colnames(Ymat) <- NULL

Xframes = XData[ ,-1] %>% 
  group_by(hour = format(XData$timestamp,format='%Y-%m-%d %H')) 
Xframes
X = group_split(Xframes)

#X_List = list()
#frames = data.frame()
#for (i in 1:length(X)){
#  if ( i%%6 == 0){
#    frames = rbind(frames, X[[i]][, -21])
#    X_List[[i/6]] = frames
#    frames = NULL
#  } else if ( i%%6 != 0 ){
#    frames =rbind(frames, X[[i]][, -21])
#  }
#}

##### Normalize on [0, 1] grid, suppose [0, 1] represents 6 hours.
nGrid = 72
tmin = 0
tmax = 1
tVec = seq(tmin, tmax, length.out = 72)
Mont = 1000
tSeq = seq(tmin, tmax, length.out =  Mont)
order = 4
nknots = 20 ## larger knots here 

#splineObj = new(orthoSpline, tmin, tmax, order, nknots)
#X_approx = list()
#for (j in 1:20){
#  mat = c()
#  for (i in 1:length(X_List)){
#    fitSeq = fitMeanCurveSingle(tVec, X_List[[i]] %>% pull(j), splineObj, lambda = 0)
#    ySeq = fitSeq$meanFunction(tSeq)
#    mat = rbind(mat, ySeq)
#  }
#  rownames(mat) = NULL
#  X_approx[[j]] = mat
#}
#save(X_approx, file = "./Dataset/X_approx.RData")

load("./Dataset/X_approx.RData")

nSample = nrow(X_approx[[1]])
Kfold = 10
cvMembership =  getCVPartition(nSample, 5) 
cf = 1 
testIndex = (cvMembership == cf)
Xtest_List = list()
Xtrain_List = list()
for (i in 1:length(X_approx)){
  Xtest_List[[i]] = X_approx[[i]][testIndex, ]
  Xtrain_List[[i]] = X_approx[[i]][!testIndex, ]
}
Ytest = Ymat[testIndex, , drop = F]
Ytrain = Ymat[!testIndex, , drop = F]



result = list()
if (method == "PSpline"){
  lambdaSeq = exp(seq(-15, -2, length.out = 5))
  etaSeq = exp(seq(-10, 1, length.out = 5))
  hSeq = exp(seq(-7, -1, length.out = 5))
  select = CV_Graph_Spline_MultiCov(Xtrain_List, Ytrain, S, order, nknots, lambdaSeq, etaSeq, hSeq, Kfold)
  
  Wind_Reduce = Graph_Spline_MultiCov(Xtrain_List, Ytrain, S, order, nknots, lambda = select$opt_lambda,
                                     eta = select$opt_eta, h =select$opt_h)
  #Wind_Reduce = Graph_Spline_MultiCov(Xtest_List, Ytest, S, order, nknots, lambda = select$opt_lambda,
  #                                    eta = select$opt_eta, h =select$opt_h)
  Yhat_Spline = matrix(0, nrow(Xtest_List[[1]]), ncol(Ymat) )
  for (i in 1:length(Xtest_List)) {
    Yhat_Spline[, i] = Xtest_List[[i]]%*%Wind_Reduce$beta[, i]/ncol(Xtest_List[[i]]) + as.numeric(Wind_Reducde$alpha[i])
  }
  MSP = mean( colMeans( (Yhat_Spline - Ytest)^2)) 
  
  result = c(Wind_Reduce, list(MSP = MSP)) 
} else if (method == "RKHS"){
  lambdaSeq = exp(seq(-12, 0, length.out = 7))
  Wind_RKHS = RKHS_Selection_MultiCov(Xtrain_List, Ytrain, lambdaSeq)
  Yhat_RKHS = matrix(0, nrow(Xtest_List[[1]]), ncol(Y))
  for (i in 1:length(Xtest_List)) {
    Yhat_RKHS[, i] = Xtest_List[[i]]%*%Wind_RKHS$beta[, i]/ncol(Xtest_List[[i]]) + as.numeric(Wind_RKHS$alpha[i, ])
  }
  
  #Fit = RKHS_Estimate(Xtest_List, Ytest, Wind_RKHS$lambda)
  result = c(Wind_RKHS, list( MSP = mean( colMeans( (Yhat_RKHS - Ytest)^2) ) ))
} else if (method == "FPCA"){
  pSeq = seq(1, 10, length.out = 5)
  Wind_FPCA = FPCA_MultiCov(Xtrain_List, Ytrain, pSeq, Kfold, Select_Method = "CV")
  
  Yhat_FPCA = matrix(0, nrow(Xtest_List[[1]]), ncol(Y) )
  for (i in 1:length(Xtest_List)){
    Yhat_FPCA[, i] = Xtest_List[[i]] %*% Wind_FPCA$beta[, i]/ncol(Xtest_List[[i]]) + as.numeric(Wind_FPCA$alpha[i,]) 
  }
  MSP = mean(colMeans( (Yhat_FPCA - Ytest)^2))
  
  result = c(Wind_FPCA, list( MSP =  MSP))
} else if (method == "Tikhonov"){
  RhoSeq = exp(seq(-15, 0, length.out = 1))
  Wind_Tikhonov = Tikhonov_MultiCov(Xtrain_List, Ytrain, RhoSeq, Kfold, Select_Method = "CV")
  
  Yhat_Tikhonov = matrix(0, nrow(Xtest_List[[1]]), ncol(Y) )
  for (i in 1:length(Xtest_List)) {
    Yhat_Tikhonov[, i] = Xtest_List[[i]]%*%Wind_Tikhonov$beta[, i]/ncol(Xtest_List[[i]]) + as.numeric(Wind_Tikhonov$alpha[i,])
  }
  MSP = mean(colMeans( (Yhat_Tikhonov - Ytest)^2))
  result = c(Wind_Tikhonov, list(MSP = MSP) ) 
}


Final = result 

savePath = paste0("./WindGraph/method-", method, ".RData")


save(Final, file = savePath)



