
Mont = 1000
M = 5

nFold = 10

## Spline 
simu_tmin = 0
simu_tmax = 1
simu_order = 4

## Tunning parameter
rankSeq = rev(2:5)

if (method == "RKHS"){
    lambdaSeq  = exp(seq(-13,0, length.out = 10)) 
} else if (method == "Tikhonov"){
    lambdaSeq = exp(seq(-2,3.5, length.out = 10)) 
} else if (method == "FPCA"){
    lambdaSeq = exp(seq(log(0.001), log(0.5), length.out = 10))
} else if (method == "Reduce"){
    if (simu_knots <= 7){
        lambdaSeq = exp(seq(-10,0, length.out = 7))  # K <= 7
    } else if ( (simu_knots > 7) && (simu_knots <= 10) ){
        lambdaSeq = exp(seq(-4,0, length.out = 7))  # K = 8, 9
    } else if ((simu_knots > 10) && (simu_knots <= 20)){
        lambdaSeq = exp(seq(0,3.4, length.out = 15))   # K = 20
    } else {
        lambdaSeq = exp(seq(2,5, length.out = 7)) 
    }
}


