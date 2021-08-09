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

X_List = list()
frames = data.frame()
for (i in 1:length(X)){
  if ( i%%6 == 0){
    frames = rbind(frames, X[[i]][, -21])
    X_List[[i/6]] = frames
    frames = NULL
  } else if ( i%%6 != 0 ){
    frames =rbind(frames, X[[i]][, -21])
  }
}

##### Normalize on [0, 1] grid, suppose [0, 1] represents 6 hours.
nGrid = 72
tmin = 0
tmax = 1
tVec = seq(tmin, tmax, length.out = 72)
Mont = 1000
tSeq = seq(tmin, tmax, length.out =  Mont)
order = 4
nknots = 20 ## larger knots here 

splineObj = new(orthoSpline, tmin, tmax, order, nknots)

X_approx = list()
for (j in 1:20){
  mat = c()
  for (i in 1:length(X_List)){
    fitSeq = fitMeanCurveSingle(tVec, X_List[[i]] %>% pull(j), splineObj, lambda = 0)
    ySeq = fitSeq$meanFunction(tSeq)
    mat = rbind(mat, ySeq)
  }
  rownames(mat) = NULL
  X_approx[[j]] = mat
}

Kfold = 10
result = list()
if (method == "GraphMulti"){
  lambdaSeq = exp(seq(-15, -2, length.out = 5))
  etaSeq = exp(seq(-10, 1, length.out = 5))
  hSeq = exp(seq(-7, -1, length.out = 5))
  select = CVGraph_Multi(X_approx, Ymat, S, order, nknots, lambdaSeq, etaSeq, hSeq, Kfold)
  Wind_Reduce = ReduceRankGraph_MultiCov(X_approx, Ymat, S, order, nknots, lambda = select$opt_lambda,
                                         eta = select$opt_eta, h =select$opt_h)
  mean( colMeans( (Wind_Reduce$Yhat - Ymat)^2)) 
  
  result = c(Wind_Reduce, list(MSE = mean( colMeans( (Wind_Reduce$Yhat - Ymat)^2)))) 
} else if (method == "RKHS"){
  lambdaSeq = exp(seq(-12, 0, length.out = 7))
  Wind_RKHS = RKHS_Selection_MultiCov(X_approx, Ymat, lambdaSeq)
  
  result = c(Wind_RKHS, list( MSE = mean( colMeans( (Wind_RKHS$Yhat - Ymat)^2) ) ))
} else if (method == "FPCA"){
  pSeq = seq(1, 20, by = 1)
  Wind_FPCA = FPCA_MultiCov(X_approx, Ymat, pSeq, Kfold, Select_Method = "AIC")
  mean(colMeans( (Wind_FPCA$Yhat - Ymat)^2))
  
  result = c(Wind_FPCA, list( MSE = mean(colMeans( (Wind_FPCA$Yhat - Ymat)^2)) ))
} else if (method == "Tikhonov"){
  RhoSeq = exp(seq(-15, 0, length.out = 1))
  Wind_Tikhonov = Tikhonov_MultiCov(X_approx, Ymat, RhoSeq, Kfold, Select_Method = "CV")
  mean(colMeans( (Wind_Tikhonov$Yhat - Ymat)^2))
  result = c(Wind_Tikhonov, list( MSE = mean(colMeans( (Wind_Tikhonov$Yhat - Ymat)^2)) ) )
}


Final = result 

savePath = paste0("./WindGraph/method-", method, ".RData")


save(Final, file = savePath)



