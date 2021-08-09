Mont = 1000
M = 6 
nFold = 10
TypeVec = c(0, 1, 2, 0, 1, 2)


simu_tmin = 0
simu_tmax = 1
simu_order = 4

if ( simu_knots <= 10){
    lambdaSeq = exp(seq(-8, 1.5, length.out = 10))  # K =5, 10 
} else {
    lambdaSeq = exp(seq(-1, 4, length.out = 7))
}

rankSeq = rev(2:5)

if (method == "GLMFPCA"){
    pSeq = seq(3, 20, by = 1)
}
