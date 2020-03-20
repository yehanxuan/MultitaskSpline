library(FDABasics)
library(rOptManifold)
Rcpp::sourceCpp("./code/gaussianOj.cpp")

set.seed(1000)
nSeq = 1000 # grid sequence length
tSeq = seq(0, 1, length.out = nSeq) # dense sequence on [0,1]
deltaT = 1/nSeq

N = 1000 # number of samples
M = 4 # number of tasks
sigma_e = 0.001 # sigma for noise


psi_1 = sin(2*pi*tSeq) #rep(1, nSeq) #\psi_1(t) = 1 is a constant function
psi_2 = cos(2*pi*tSeq) #\psi_2(t) = t^2
PSI = cbind(psi_1, psi_2)
PSI_SVD = svd(PSI)
PSI = PSI_SVD$u %*% diag(PSI_SVD$d) # orthonormalize

A = matrix(rnorm(2*M), M ,2)
A = qr.Q(qr(A))

beta_true = PSI %*% t(A) # true regression function

## Generate Data
phi_0 = rep(1,nSeq) * 5
phi_1 = sin(2*pi*tSeq) * 5
phi_2 = cos(2*pi*tSeq) * 4
phi_3 = sin(4*pi*tSeq) * 3
phi_4 = cos(4*pi*tSeq) * 2
phi_5 = sin(6*pi*tSeq) * 1
phi_mat = cbind(phi_0,phi_1,phi_2,phi_3,phi_4,phi_5)


set.seed(102)
score_mat = matrix(rnorm(6*N), 6, N)
X_mat = t(phi_mat %*% score_mat)
Y_mat = X_mat %*% beta_true  * deltaT + matrix(rnorm(N*M, sd = sigma_e) , N,M)

# parameters
tmin = 0
tmax = 1
order = 4
nknots = 4
splineBasis = new(orthoSpline, tmin, tmax, order, nknots)
basisMat = splineBasis$evalSpline(tSeq)

# Least Square Estimation
Z_mat = X_mat %*% t(basisMat)*deltaT
tmp = t(Z_mat) %*% Z_mat 
betaHat = solve(tmp, t(Z_mat) %*% Y_mat)

UInit = svd(betaHat)$u[,1:2]
VInit = svd(betaHat)$v[,1:2]
WInit = diag(svd(betaHat)$d[1:2])
InitList = list(UInit, WInit, VInit)

k = 2
#par(mfrow=c(2,1), mar=c(2,4,0,0))
plot(tSeq, beta_true[,k], type = "l")
lines(tSeq, t(basisMat) %*% betaHat[,k], col = "red", type = "l")

Omega = splineBasis$get_Omega()
Omega = Omega/max(abs(Omega))

# Manifold Estimation
optObj = new(gaussianObj, Y_mat, Z_mat)
optObj$set_penaltyMatrix(Omega)
optObj$set_tuningParameter(0)
optObj$objF(InitList)
# Random Initial Point
Rank = 3 #rank
K = splineBasis$getDoF() # spline degree of freedom
UInit = matrix(rnorm(Rank*K), K, Rank)
UInit = qr.Q(qr(UInit))
VInit = matrix(rnorm(Rank*M), M,Rank)
VInit = qr.Q(qr(VInit))
WInit = diag(rep(1,Rank))
InitList = list(UInit, WInit, VInit)
optObj$objF(InitList)
#optObj$gradF(InitList)

controlList = list(alpha = 0.1, tol = 1e-8, 
                   sigma = 1e-3, beta = 0.618,
                   iterMax = 1e4, verbose = 1)


problem = new(manifoldOpt) # manifold optimization method
problem$select_algorithm("sd") # steepest descent method
problem$set_eigenRegUnsym(K, M, Rank) # U*W*V^T \in R^{K\times M}, rank Rank
# The manifold loss version
problem$setObjective(optObj$objF)
problem$setGradient(optObj$gradF)
problem$initial_point(InitList) # Initial Point
problem$update_control(controlList) # Control 
problem$solve() # Execute 
SFinal = problem$get_optimizer() # Get optimal Solution

#SFinal = InitList
optObj$objF(SFinal)

CHat = SFinal[[1]] %*% SFinal[[2]] %*% t(SFinal[[3]])
k = 2
plot(tSeq, beta_true[,k], type = "l")
lines(tSeq, t(basisMat) %*% CHat[,k], col = "red", type = "l")

plot(X_mat %*% t(basisMat) %*% CHat[,k] * deltaT, Y_mat[,k], pch = 20)
