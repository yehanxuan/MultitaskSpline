### Main Spectral data ####
rm(list = ls())
library(FDABasics)
source("./models/MultiRegression.R")
source("./models/FPCA.R")
source("./models/Tikhonov.R")
source("./models/RKHS.R")

args = commandArgs(trailingOnly = TRUE)
method = args[1]

### only sample the first 1000 data
Y = spcTable[1:1000, ][, c("logg", "teff", "feh")]
X_raw = spcData[1:1000, ]

## The wavelength point of each light curve  
wavelength = 10^seq(3.57, 3.9, length.out=2000)
logwavelength = log(wavelength)
grids = ( logwavelength - min(logwavelength))/(max(logwavelength) - min(logwavelength))

X = X_raw 
nSample = dim(X)[1]
cvMembership = getCVPartition(nSample, 5)
cf = 1
test_index =  which(cvMembership == cf)
# Standarize 
Y_new = apply(Y, 2, function(x) (x - mean(x))/sd(x) )
Y = Y_new

X_test = X[test_index,]
X_train = X[-test_index, ]
Y_test = Y[test_index, ]
Y_train = Y[-test_index, ]


Kfold = 10
nKnots = 30
result = list()
if (method == "Reduce"){
  mOrder = 4
  Mont = dim(X)[2]
  lambdaSeq = exp(seq(-45, -17, length.out = 5)) 
  RankSeq = c(1:3)
  Select = CVSelection_Multi(X_train, Y_train, Mont, Kfold, mOrder, nKnots, lambdaSeq,
                             RankSeq, method = "Reduce")
  fit = ReduceRankSPline(X_train, Y_train, Mont, mOrder, nKnots, Select$opt_rank,
                         Select$opt_lambda)
  Ytest_reduce = (X_test %*% fit$beta)/Mont + matrix(1, nrow(X_test), 1)%*%t(fit$alpha)
  colnames(Ytest_reduce) <- colnames(Y_test)
  MSP = colSums( (Y_test - Ytest_reduce)^2 )/nrow(Y_test) 
  result = c(fit, list(MSP = MSP))
} else if (method == "FPCA"){
  pSeq = seq(1, 15, by = 1)
  Mont = ncol(X)
  fit = FPCA_Graph(X_train, Y_train, pSeq, Kfold, Select_Method = "CV")
  Ytest_FPCA = (X_test %*% fit$beta)/ncol(X_test) + matrix(1, nrow(X_test), 1)%*%t(fit$alpha)
  colnames(Ytest_FPCA) <- colnames(Y_test)
  MSP = colSums((Y_test - Ytest_FPCA)^2)/nrow(Y_test)
  result = c(fit, list(MSP = MSP))
} else if (method == "Tikhonov"){
  RhoSeq =  exp(seq(-5, 2 , length.out = 7)) 
  fit_Tikhonov = Tikhonov_Graph(X_train, Y_train, RhoSeq, Kfold, Select_Method = "CV")
  Ytest_Tikhonov = (X_test%*% fit_Tikhonov$beta)/ncol(X_test) + 
    matrix(1, nrow(X_test) ,1)%*%t(fit_Tikhonov$alpha)
  colnames(Ytest_Tikhonov) <- colnames(Y_test)
  MSP = colSums((Y_test - Ytest_Tikhonov)^2)/nrow(Y_test)
  result = c(fit_Tikhonov, list(MSP = MSP))
} else if (method == "RKHS"){
  lambdaSeq = exp(seq(-45, -17, length.out = 5))
  fit_RKHS = RKHS_Wrap(X_train, Y_train, lambdaSeq)
  Ytest_RKHS = (X_test %*% fit_RKHS$beta)/ncol(X_test) + 
    matrix(1, nrow(X_test) ,1)%*%t(fit_RKHS$alpha)
  colnames(Ytest_RKHS) <- colnames(Y_test)
  MSP =  colSums((Y_test - Ytest_RKHS)^2)/nrow(Y_test)
  result = c(fit_RKHS, list(MSP = MSP))
}

Final = result 
savePath = paste0("./SpectraData/method-", method, "-", nKnots, ".RData")

save(Final, file = savePath)


