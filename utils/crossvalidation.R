getCVPartition = function(nSample, nFold){
    nPerFold = ceiling(nSample / nFold)
    cvMembership = rep(1:nFold, nPerFold)
    shuffleV = sample(1:nSample, size = nSample, replace = FALSE)
    cvMembership = cvMembership[shuffleV]
    return(cvMembership)
}






