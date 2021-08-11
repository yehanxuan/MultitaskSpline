rm(list = ls())
library(snowfall)

args <- commandArgs(trailingOnly = TRUE)
settingsID = as.numeric(args[1]) ## simulation settings 1/2
simu_n = as.numeric(args[2])
method = args[3]
repID = as.numeric(args[4])
simu_knots = as.numeric(args[5])
sigma_e = as.numeric(args[6])

method = "PSpline"
settingsID = 3
simu_n = 100
repID = 2
simu_knots = 15
sigma_e = 0.5

source("./Dataset/Generate_data.R")
source("./loadAll.R")
Rcpp::sourceCpp("models/gaussianOj.cpp")
library(FDABasics)
library(rOptManifold)
source("models/RKHS.R")
source("models/FPCA.R")
source("models/Tikhonov.R")
source("models/MultiRegression.R")
source("./models/GraphicModel.R")

source(paste0("./simuSetting/simuSetting-", settingsID, ".R"))
source(paste0("./OneRep/oneRep-Graph.R"))
source("utils/resultGLM.R")

savepath = paste0("./data/setting-", settingsID,"-", simu_n, "-", method, "-",simu_knots, "-", repID, "-", sigma_e, ".RData")

#if (file.exists(savepath) ) stop("file exists")

nCPUS = 2
maxIter = 2

result1 = list()
result2 = list()
result3 = list()
result4 = list()
result5 = list()

for (i in 1:ceiling(maxIter/nCPUS)){
  print(i)
  sfInit(parallel = TRUE, cpus = nCPUS)
  sfExportAll()
  sfLibrary(FDABasics)
  sfLibrary(fAdditiveModel)
  sfLibrary(rOptManifold)
  sfLibrary(matrixcalc)
  sBegin = (i-1)*nCPUS +1 
  sEnd = min(i*nCPUS, maxIter)
  seedSeq = seq(sBegin, sEnd, by = 1)
  tmp = sfClusterApplyLB(seedSeq, oneReplicateWrap_Graph)
  sfStop()
  tmp1 = lapply(tmp, function(x) x[[1]] )
  tmp2 = lapply(tmp, function(x) x[[2]] )
  tmp3 = lapply(tmp, function(x) x[[3]] )
  tmp4 = lapply(tmp, function(x) x[[4]] ) 
  tmp5 = lapply(tmp, function(x) x[[5]] )
  result1 = c(result1, tmp1)
  result2 = c(result2, tmp2)
  result3 = c(result3, tmp3)
  result4 = c(result4, tmp4)
  result5 = c(result5, tmp5)
}


Final = list("Est" = result1, "Pred" = result2, "lambda" = result3, "eta" = result4, "h" = result5)
save(Final, file = savepath)

Est = compute(result1)
Est

Pred = compute(result2)
Pred

lambda = result3
lambda

eta = result4
eta

h = result5
h
