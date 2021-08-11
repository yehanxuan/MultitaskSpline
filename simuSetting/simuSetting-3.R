Mont = 1000
M = 20
simu_tmin = 0
simu_tmax = 1
simu_order = 4
nFold = 10
Select_Method = "CV"
if (method == "PSpline"){
    lambdaSeq = exp(seq(-5, 2, length.out = 3))
} else if (method == "Tikhonov"){
    lambdaSeq = exp(seq(-5,4, length.out = 10)) 
} else if (method == "FPCA"){
    lambdaSeq = seq(1, 8, by = 1)
} else if (method == "RKHS"){
    lambdaSeq = exp(seq(-7, 1, length.out = 7))
}

etaSeq = exp(seq(-5, 2, length.out = 3))
hSeq = exp(seq(-5, 2, length.out = 4))
