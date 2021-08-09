simuData = Generate_GraphData(Mont, simu_n, M, sigma_e)
Y_mat = simuData$Y
X_mat = simuData$X
CovMat = simuData$CovMat
beta_true = simuData$beta_true
S = simuData$location
