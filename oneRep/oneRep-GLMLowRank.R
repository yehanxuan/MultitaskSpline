splineObj = new(orthoSpline, simu_tmin, simu_tmax, simu_order, simu_knots)


Omega = splineObj$get_Omega()
Omega = Omega/ mean(abs(Omega))
simu_K = splineObj$getDoF()


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
basisMat = splineObj$evalSpline(grid)
Zmat = Xmat%*%t(basisMat)/Mont

#AInit = matrix(rnorm(Rank*M), M, Rank)
#DInit = matrix(rnorm(Rank*simu_K), simu_K, Rank)

