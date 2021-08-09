
if (DataType == "Reg"){
    K.set = seq(5, 20, by = 2)
    nu = 2
    simuData = Generate_MultiData(Mont, simu_n, M ,sigma_e, nu)
    Ymat = simuData$Y
    Xmat = simuData$X
    CovMat = simuData$CovMat
    beta_true = simuData$beta_true    
} else if (DataType == "misalign"){
    k0 = 10
    simuData = Generate_cai2011(Mont, simu_n, M, sigma_e, CovType = "misalign", k0)
    Ymat = simuData$Y
    Xmat = simuData$X
    CovMat = simuData$CovMat
    beta_true = simuData$beta_true
} else if (DataType == "Haar"){
    simuData = Generate_cai2011(Mont, simu_n, M, sigma_e, CovType = "Haar")
    Ymat = simuData$Y
    Xmat = simuData$X
    CovMat = simuData$CovMat
    beta_true = simuData$beta_true
}
