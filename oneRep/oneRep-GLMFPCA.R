## Generate Mixed data
nu = 1.5
if (DataType == "GLMexact"){
    simuData = GenMulti_Mixed(Mont, simu_n, M, sigma_e, nu, TypeVec)
} else if (DataType == "GLMHybrid"){
    simuData = GenMulti_Mixed_Hybrid(Mont, simu_n, M, sigma_e, nu, TypeVec)
}

grid = seq(1/Mont, 1, 1/Mont)
Ymat = simuData$Y
Xmat = simuData$X
CovMat = simuData$C

beta_true = simuData$beta

