
#### Generate Multitask regression data ####
Generate_MultiData = function(Mont, nsample, M, sigma_e, nu){
    tSeq = seq(0, 1, length.out = Mont) # dense sequence on [0,1]
    deltaT = 1/Mont
    #psi_1 = sqrt(2)*sin(pi*tSeq)
    psi_1 = sqrt(2)*cos(2*pi*tSeq) #rep(1, nSeq) #\psi_1(t) = 1 is a constant function
    psi_2 = sqrt(2)*cos(3*pi*tSeq)
   # psi_2 = sqrt(2)*sin(3*pi*tSeq) #\psi_2(t) = t^2     # psi becomes the orthonormal basis
    psi_3 = sqrt(2)*cos(4*pi*tSeq)
    #psi_3 = sqrt(2)*sin(5*pi*tSeq)

    PSI = cbind(psi_1, psi_2, psi_3)
    PSI_SVD = svd(PSI)
    PSI = PSI_SVD$u %*% diag(PSI_SVD$d) # orthonormalize
    eigen_beta = 5*diag(seq(1, 3)^-1.5)
    A = matrix(rnorm(3*M), M ,3)
    A = qr.Q(qr(A))
    beta_true = PSI %*% eigen_beta %*% t(A) # true regression function
    
    ## Generate Data
    # phi_0 = rep(1,Mont)
    # phi_1 = sqrt(2)*sin(2*pi*tSeq)
    # phi_2 = sqrt(2)*cos(2*pi*tSeq)
    # phi_3 = sqrt(2)*sin(4*pi*tSeq)
    # phi_4 = sqrt(2)*cos(4*pi*tSeq)
    # phi_5 = sqrt(2)*sin(6*pi*tSeq)
    # phi_mat = cbind(phi_0,phi_1,phi_2,phi_3,phi_4,phi_5)
    # nu = 2
    # eigen_x = diag(seq(1,6)^-(nu/2))
    # score_mat = matrix(rnorm(6*nsample), 6, nsample)
    #X_mat = t(phi_mat %*% eigen_x%*%score_mat)
    #C = phi_mat%*%eigen_x%*%eigen_x%*%t(phi_mat)
    k = 50
    Phi = matrix(1, Mont, k)
    for ( i in 2:k ){
      Phi[ ,i] = sqrt(2)*cos((i-1)*pi*tSeq)
    }
    weight = c(1:k)
    # nu = 1.5
    eigen_x =diag( (-1)^(weight+1)*(weight^(-nu/2)) )
    # like Z in Cai
    Score_mat = matrix(rnorm(k*nsample), k, nsample)
    X_mat = t(Phi%*%eigen_x%*%Score_mat)
    Y_mat = X_mat %*% beta_true  * deltaT + matrix(rnorm(nsample*M, sd = sigma_e) , nsample,M)
    C = Phi%*%eigen_x%*%eigen_x%*%t(Phi)
    
    
    return(list("beta_true" = beta_true, "X" = X_mat, "Y" = Y_mat, "CovMat" = C))
}

# Generate data represented by spline basis 
Generate_MultiData_Spline = function(Mont, nSample, M, r = 3, sigma_e, nu) {
  tSeq = seq(0, 1, length.out = Mont) 
  deltaT = 1/Mont
  
  A = matrix(rnorm(r*M), M, r)
  A = qr.Q(qr(A))
  
  k = 50
  Phi = matrix(1, Mont, k)
  for ( i in 2:k ){
    Phi[ ,i] = sqrt(2)*cos((i-1)*pi*tSeq)
  }
  weight = c(1:k)
  eigen_x =diag( (-1)^(weight+1)*(weight^(-nu/2)) )
  Score_mat = matrix(rnorm(k*nsample), k, nsample)
  X_mat = t(Phi%*%eigen_x%*%Score_mat)
}


#### Generate data for Mixed generalized linear model #### Exact case
GenMulti_Mixed = function(Mont, nsample, M, sigma_e, nu,  TypeVec){
    tSeq = seq(0, 1, length.out = Mont) # dense sequence on [0,1]
    deltaT = 1/Mont
    psi_1 = sqrt(2)*cos(2*pi*tSeq) #rep(1, nSeq) #\psi_1(t) = 1 is a constant function
    psi_2 = sqrt(2)*sin(2*pi*tSeq) #\psi_2(t) = t^2
    psi_3 = sqrt(2)*cos(4*pi*tSeq)
    PSI = cbind(psi_1, psi_2, psi_3)
    PSI_SVD = svd(PSI)
    PSI = PSI_SVD$u %*% diag(PSI_SVD$d) # orthonormalize
    eigen_beta = 5*diag(seq(1, 3)^-2)
    A = matrix(rnorm(3*M), M ,3)
    A = qr.Q(qr(A))
    beta_true = PSI %*% eigen_beta %*% t(A) # true regression function
     ## Generate Data
    # phi_0 = rep(1,Mont)
    # phi_1 = sqrt(2)*sin(2*pi*tSeq) 
    # phi_2 = sqrt(2)*cos(2*pi*tSeq) 
    # phi_3 = sqrt(2)*sin(4*pi*tSeq) 
    # phi_4 = sqrt(2)*cos(4*pi*tSeq) 
    # phi_5 = sqrt(2)*sin(6*pi*tSeq)
    # phi_mat = cbind(phi_0,phi_1,phi_2,phi_3,phi_4,phi_5)
    # nu = 2
    # eigen_x = diag(seq(1,6)^-(nu/2))
    # score_mat = matrix(rnorm(6*nsample), 6, nsample)
    #    X_mat = t(phi_mat %*% eigen_x%*%score_mat)
    
    #    lambda = exp(X_mat%*%beta_true/Mont)
    #    P = lambda/(1+lambda)
    k = 50 
    Phi = matrix(1, Mont, k)
    for ( i in 2:k ){
        Phi[ ,i] = sqrt(2)*cos((i-1)*pi*tSeq)
    }
    weight = c(1:k)
    # nu = 1.5
    eigen_x =diag( (-1)^(weight+1)*(weight^(-nu/2)) )
    # like Z in Cai
    Score_mat = matrix(rnorm(k*nsample), k, nsample)
    
    
    Y_mat = matrix(0, nsample, M)
    for (i in 1:length(TypeVec)){
        X_mat = t(Phi %*% eigen_x%*%Score_mat)
        lambda = exp(X_mat%*%beta_true/Mont)
        P = lambda/(1+lambda)
        if (TypeVec[i] == 0 ){
            Y_mat[ ,i] = X_mat %*% beta_true[ ,i]*deltaT + rnorm(nsample, sd = sigma_e)
        } else if (TypeVec[i] == 1){
            Y_mat[ ,i] = 1*(runif(nsample) < P[ ,i])
        } else if (TypeVec[i] == 2){
            Y_mat[ ,i] = rpois(nsample, lambda[, i])
        }
    }
    C = Phi%*%eigen_x%*%eigen_x%*%t(Phi)
    return(list("X" = X_mat, "Y" = Y_mat, "beta" = beta_true, "CovMat" = C, "TypeVec" = TypeVec) )
}

#### Generate data for Mixed generalized linear model #### Hybrid case
GenMulti_Mixed_Hybrid = function(Mont, nsample, M, sigma_e, nu, TypeVec){
    tSeq = seq(0, 1, length.out = Mont) 
    deltaT = 1/Mont
    psi_1 = sqrt(2)*cos(2*pi*tSeq) #rep(1, nSeq) #\psi_1(t) = 1 is a constant function
    psi_2 = sqrt(2)*sin(2*pi*tSeq) #\psi_2(t) = t^2
    psi_3 = sqrt(2)*cos(4*pi*tSeq)
    psi_4 = sqrt(2)*sin(4*pi*tSeq)
    psi_5 = sqrt(2)*cos(6*pi*tSeq)
    psi_6 = sqrt(2)*sin(6*pi*tSeq)
    PSI = cbind(psi_1, psi_2, psi_3, psi_4, psi_5, psi_6)
    PSI_SVD = svd(PSI)
    PSI = PSI_SVD$u %*% diag(PSI_SVD$d)
    eigen_beta = diag(c(5.0, 2.0, 1.0, 0.25, 0.015, 0.001))
    A = matrix(rnorm(6*M), M ,6)
    A = qr.Q(qr(A))
    beta_true = PSI %*% eigen_beta %*% t(A) # true regression function
    
    k = 50 
    Phi = matrix(1, Mont, k)
    for ( i in 2:k ){
        Phi[ ,i] = sqrt(2)*cos((i-1)*pi*tSeq)
    }
    weight = c(1:k)
    # nu = 1.5
    eigen_x =diag( (-1)^(weight+1)*(weight^(-nu/2)) )
    # like Z in Cai
    Score_mat = matrix(rnorm(k*nsample), k, nsample)
    Y_mat = matrix(0, nsample, M)
    
    for (i in 1:length(TypeVec)){
        X_mat = t(Phi %*% eigen_x%*%Score_mat)
        lambda = exp(X_mat%*%beta_true/Mont)
        P = lambda/(1+lambda)
        if (TypeVec[i] == 0 ){
            Y_mat[ ,i] = X_mat %*% beta_true[ ,i]*deltaT + rnorm(nsample, sd = sigma_e)
        } else if (TypeVec[i] == 1){
            Y_mat[ ,i] = 1*(runif(nsample) < P[ ,i])
        } else if (TypeVec[i] == 2){
            Y_mat[ ,i] = rpois(nsample, lambda[, i])
        }
    }
    C = Phi%*%eigen_x%*%eigen_x%*%t(Phi)
    return(list("X" = X_mat, "Y" = Y_mat, "beta" = beta_true, "CovMat" = C, "TypeVec" = TypeVec) )
}


#### Generate graph regularized data ####
Generate_GraphData = function(Mont, nsample, M, sigma_e){
    tSeq = seq(0, 1, length.out = Mont)
    deltaT = 1/Mont
    s1 = runif(M)
    s2 = runif(M)
    k = 10
    Psi = matrix(0, Mont, k)
    for (i in 1:k){
        Psi[ , i] = sqrt(2)*cos(i*pi*tSeq)
    }
    eigen_beta = seq(1,k)^(-2)
    #Psi%*%eigen_beta
    # Location vector
    #L = (s1 + 1) * (s2^2 -s2 + 3)
    L = (s1 + 1) * (s2^2 -s2 + 1)
    beta_true = Psi%*%eigen_beta%*%t(L)
    
    
    kx = 50
    Phi = matrix(1, Mont, kx)
    for (i in 2:kx){
        Phi[ ,i] = sqrt(2)*cos((i-1)*pi*tSeq)
    }
    weight = c(1:kx)
    nu = 1.5
    eigen_x =diag( (-1)^(weight+1)*(weight^(-nu/2)) )
    Score_mat = matrix(rnorm(kx*nsample), kx, nsample)
    X_mat = t(Phi%*%eigen_x%*%Score_mat)
    Y_mat = X_mat%*%beta_true * deltaT + matrix(rnorm(nsample*M, sd = sigma_e) , nsample,M)
    C = Phi%*%eigen_x%*%eigen_x%*%t(Phi)
    S = cbind(s1, s2)
    return(list("beta_true" = beta_true, "X" = X_mat, "Y" = Y_mat, "CovMat" = C, "location" = S))
}


GenerateData = function(Mont, k, nu, nsample, sigma, method){   ### Tony Cai
    # Assume the interval is [0, 1]
    grid = seq(1/Mont, 1, 1/Mont)
    Phi = matrix(1, Mont, k)
    for ( i in 2:k ){
        Phi[ ,i] = sqrt(2)*cos((i-1)*pi*grid)
    }
    # generate zeta
    if (method == "close"){
        zeta = rep(1,k)
        for (i in 2:k){
            if ((i >= 2)&(i <= 4)){
                zeta[i] = 0.2*((-1)^(i+1))*(1 - 0.0001*i)
            } else {
                zeta[i] = 0.2*((-1)^(i+1))*( (5*floor(i/5))^(-nu/2) - 0.0001*(i %% 5))
            }
        } } else if (method == "well"){
            weight = c(1:k)
            zeta = (-1)^(weight+1)*(weight^(-nu/2))
        }
    # generate beta0
    weight = c(1:k)
    weight1 = 4*((-1)^(weight+1)) * (weight^(-2))
    beta0 = Phi%*%weight1
    # generate Z, X  
    X = matrix(0, Mont, nsample)
    for (i in 1:nsample){
        Z = runif(k, min = -sqrt(3), max = sqrt(3))
        X[, i] = Phi%*%(zeta*Z)
    }
    # generate Y, we use Monte Carlo computation here
    Y = t(X)%*%beta0/Mont + rnorm(nsample, 0 ,sd = sigma)
    
    C= Phi%*%((zeta^2)*t(Phi))
    return(list("X" = t(X), "Y" = Y, "beta0" = beta0, "CovMat" = C))
}


Generate_cai2011 = function(Mont, nsample, M, sigma_e, CovType, k0 = 5){
  tSeq = seq(0, 1, length.out = Mont) 
  deltaT = 1/Mont
  psi_1 = sqrt(2)*cos(2*pi*tSeq) 
  # psi_2 = sqrt(2)*cos(3*pi*tSeq) 
  psi_2 = sqrt(2)*cos(3*pi*tSeq) 
  psi_3 = sqrt(2)*cos(4*pi*tSeq)
  PSI = cbind(psi_1, psi_2, psi_3)
  PSI_SVD = svd(PSI)
  PSI = PSI_SVD$u %*% diag(PSI_SVD$d) # orthonormalize
  eigen_beta = 5*diag(seq(1, 3)^-2)
  A = matrix(rnorm(3*M), M ,3)
  A = qr.Q(qr(A))
  beta_true = PSI %*% eigen_beta %*% t(A) # true regression function
  
  k = 50
  Phi = matrix(1, Mont, k)
  if (CovType == "misalign"){
    for ( i in 2:k ){
      Phi[ ,i] = sqrt(2)*cos((i-1)*pi*tSeq)
    }
    weight = c(1:k)
    eigen_x = diag( (-1)^(weight+1) * ( abs(weight - k0) + 1)^{-1} )
  } else if (CovType == "Haar"){
    weight = c(1:k)
    eigen_x =diag( (-1)^(weight+1)*(weight^(-1)) )
    #sSeq = floor( log(weight, 2))
    #lSeq = (weight + 1) - 2^(sSeq)
    for (i in 2:k){
      s = floor( log(i-1, 2))
      l = i -1 +1 - 2^(s)
      Phi[, i] = mapply(Haar, t = tSeq, s, l)
    }
  }
  
  Score_mat = matrix(rnorm(k*nsample), k, nsample)
  X_mat = t(Phi%*%eigen_x%*%Score_mat)
  Y_mat = X_mat %*% beta_true *deltaT + matrix(rnorm(nsample*M, sd = sigma_e) , nsample,M)
  C = Phi%*%eigen_x%*%eigen_x%*%t(Phi)
  return(list("beta_true" = beta_true, "X" = X_mat, "Y" = Y_mat, "CovMat" = C))
}


## Haar function ##
h_1_0 = function(t){
  if ((t >= 0) & (t <= 1/2)){return(1)}
  else if ((t >= 1/2) & (t <= 1)){return(-1)}
  else{return(0)}
}


Haar = function(t,s,l){
  return( 2^(s/2) * h_1_0( 2^(s) * t-l+1) )
}

# }
