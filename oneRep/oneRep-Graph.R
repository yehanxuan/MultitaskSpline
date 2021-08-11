simuData = Generate_GraphData(Mont, simu_n, M, sigma_e)
Ymat = simuData$Y
Xmat = simuData$X
CovMat = simuData$CovMat
beta_true = simuData$beta_true
S = simuData$location
