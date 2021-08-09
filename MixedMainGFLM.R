rm(list = ls())
library(snowfall)

args <- commandArgs(trailingOnly = TRUE)
settingsID = as.numeric(args[1]) ## simulation settings 1/2
simu_n = as.numeric(args[2])
method = args[3] # model fitting method
repID = as.numeric(args[4])
simu_knots = as.numeric(args[5])
DataType = args[6]
sigma_e = as.numeric(args[7])

settingsID = 2
simu_n = 200
method ="GLMFPCA"
repID = 2
simu_knots = 20
DataType = "GLMexact"
sigma_e = 0.5 


source("./loadAll.R")
source("./Dataset/Generate_data.R")
source("./models/SplineGFLM.R")
source("./models/FPCAGFLM.R")

source(paste0("./simuSetting/simuSetting-", settingsID, ".R"))
source(paste0("./oneRep/oneRep-", method, ".R"))

savepath = paste0("./data/setting-", settingsID,"-",simu_n, "-", method, "-", simu_knots, "-",repID, "-",
                  DataType, "-", sigma_e, ".RData")

nCPUS = 4
maxIter = 4


result1 = list()
result2 = list()
result3 = list()
result4 = list()

if (method == "GLMLowRank"){
  for (i in 1:ceiling(maxIter/nCPUS)){
    print(i)
    sfInit(parallel = TRUE, cpus = nCPUS)
    sfExportAll()
    sfLibrary(FDABasics)
    sfLibrary(fAdditiveModel)
    sfLibrary(rOptManifold)
    sBegin = (i-1)*nCPUS +1 
    sEnd = min(i*nCPUS, maxIter)
    seedSeq = seq(sBegin, sEnd, by = 1)
    tmp = sfClusterApplyLB(seedSeq, oneReplicateWrap_GFLM)
    sfStop()
    tmp1 = lapply(tmp, function(x) x[[1]] )
    tmp2 = lapply(tmp, function(x) x[[2]] )
    tmp3 = lapply(tmp, function(x) x[[3]] )
    tmp4 = lapply(tmp, function(x) x[[4]] )
    result1 = c(result1, tmp1)
    result2 = c(result2, tmp2)
    result3 = c(result3, tmp3)
    result4 = c(result4, tmp4)
  }
} else if (method == "GLMFPCA"){
  for (i in 1:ceiling(maxIter/nCPUS)){
    print(i)
    sfInit(parallel = TRUE, cpus = nCPUS)
    sfExportAll()
    sfLibrary(FDABasics)
    sfLibrary(fAdditiveModel)
    sfLibrary(rOptManifold)
    sBegin = (i-1)*nCPUS +1 
    sEnd = min(i*nCPUS, maxIter)
    seedSeq = seq(sBegin, sEnd, by = 1)
    tmp = sfClusterApplyLB(seedSeq, oneReplicateWrap_GFLMFPCA)
    sfStop()
    tmp1 = lapply(tmp, function(x) x[[1]] )
    tmp2 = lapply(tmp, function(x) x[[2]] )
    tmp3 = lapply(tmp, function(x) x[[3]] )
    tmp4 = lapply(tmp, function(x) x[[4]] )
    result1 = c(result1, tmp1)
    result2 = c(result2, tmp2)
    result3 = c(result3, tmp3)
    result4 = c(result4, tmp4)
  }
}


Final = list("Est" = result1, "Pred" = result2, "lambda" = result3, "rank" = result4)
save(Final, file = savepath)


Est = compute(result1)
Est

Pred = compute(result2)
Pred

lambda = result3
lambda

rank = result4






