###--------------------------------------------------------------------------###
###-- Sparse TVP VECMs with an application to modeling electricity prices ---###
###------------------ Hauzenberger, Pfarrhofer, & Rossini -------------------###
###----------------- International Journal of Forecasting -------------------###
###--------------------------------------------------------------------------###
###-------------------------- TVP-VECM estimation ---------------------------###
###--------------------------------------------------------------------------###
rm(list=ls())

###--------------------------------------------------------------------------###
###-------------- Load packages and load additional functions ---------------###
###--------------------------------------------------------------------------###
# Packages
require(MASS)
require(Matrix) 
require(mvtnorm)
require(zoo)
require(Rcpp)
require(stochvol)
require(glasso)

# Functions
source("aux_funcs/aux-file.R")
sourceCpp("aux_funcs/kf.cpp")

# Directories
wdir <- ""        # Working directory  
resdir <- "res/"  # Results directory
dir.create(resdir)

###--------------------------------------------------------------------------###
###------------------ Set up to produce a single forecast -------------------###
###--------------------------------------------------------------------------###
# Estimation setup
roll.window <- 365 # Rolling window of 365 observations
end.in <- 2019 + 364/365 # Last sample period
p <- 2         # No. of lags
p.tilda <- p+1 # No. of lags in levels (p+1)
model <- "VECM-TVP-iSV-t"
type <- substr(model,1,4)
cons.model <- (substr(model, 6,8) == "TIV")
sv <- substr(model, 10,14) # triangular (cholesky) SV ("iSV")

# MCMC setup
nburn <- 100
nsave <- 100
thin  <- 3
ntot <- nburn + thin*nsave
save.set <- seq(nburn+thin, ntot, thin)
jrep <- 0 

# Data setup
load("ELP-data/DE_elp_hourly.rda")
Yraw <- price.data # Yraw is days x 24 (hours)
# Average over night hours
Yraw <- cbind(Yraw[,c(8:18)], rowMeans(Yraw[,c(1:7, 19:24)]))
colnames(Yraw) <- c(paste0(rep("h",11), 8:18), "night")

Yraw.out <- window(Yraw, start = (end.in + 1/365), end = (end.in + 1/365), extend = TRUE)
Yraw     <- window(Yraw, (end.in - (roll.window-1)/365),  end = end.in)

# For estimation, we standardized the data 
Yraw.mu <- apply(Yraw, 2, mean); Yraw.sd <- apply(Yraw, 2, sd) # Can be used to re-scale forecasts
Yraw    <- apply(Yraw, 2, function(x){(x-mean(x))/sd(x)})
    
Z <- mlag(Yraw,1) #y_t-1
Y.lag <- mlag(Yraw,p+1) # lags in levels: y_t-1, ..., y_(t-p-1)
Y.lag <- Y.lag[(p+2):nrow(Y.lag),]
Z <- Z[(p+2):nrow(Z),]

# Bring it in VECM-form: Delta y_t = c + a_t b' y_t-1 + A_1t Delta y_t-1 + ... A_pt Delta y_t-p + e_t
Y <- diff(Yraw) # Delta y_t
Xlag <- mlag(Y,p) #Delta y_(t-1):(t-p)
X <- Xlag[(p+1):nrow(Xlag),]
Y <- Y[(p+1):nrow(Y),]
Yraw <- Yraw[(p+2):nrow(Yraw),]
X.levels <- Y.lag

# Get key dimensions of the data
M <- ncol(Y)
r <- M # Maximum number of cointegration relationships
K <- ncol(X) 
M.beta <- M 
KK <- K + r   
K.var <- K + M.beta 
T <- nrow(Y)
k <- K*M
v <- (M*(M-1))/2
id_free_cov <- lower.tri(diag(M))

###--------------------------------------------------------------------------###
###----------- Prior setup, starting values and storage matrices ------------###
###--------------------------------------------------------------------------###
# Priors for SV-t
sv_priors <- specify_priors(
  mu = sv_normal(mean = 0, sd = 1), # prior on unconditional mean in the state equation
  phi = sv_beta(shape1 = 25, shape2 = 1.5), #informative prior to push the model towards a random walk in the state equation (for comparability)
  sigma2 = sv_gamma(shape = 0.5, rate = 0.5), # Gamma prior on the state innovation variance
  nu = sv_exponential(0.1),
  rho = sv_constant(0))
  
# Initialization of SV processes
svdraw <- list(mu = 0, phi = 0.99, sigma = 0.01, nu = 10000, rho = 0, beta = NA, latent0 = 0)

# Variances 
sv.idio <- list()
for (jj in 1:M) sv.idio[[jj]] <- svdraw
ht <- Hvar <- matrix(0.1,T,M) #idiosyncratic variances
pars_var <- matrix(0,4,M)
tau.fat <- matrix(1,T,M)

# HS prior on VAR coefficients
tau.A <- tau.o <- 1
zeta.A <- zeta.o <- 1
nu.A <- nu.o <- rep(1,k+r*M)
lambda.A <- lambda.o <- rep(1,k+r*M)
  
# HS prior on covariances
tau.cov <- tau.zeta <- 1
zeta.cov <- zeta.zeta <- 1
nu.cov <- nu.zeta <- rep(1,v)
lambda.cov <- lambda.zeta <- rep(1,v)

# Prior mean on VAR coefficients
pr.mean <- 0
A.prior <- matrix(0,2*KK,M)
A.prior[(r+1):(r+M),] <- diag(M)*pr.mean

# Prior mean on cointegration relationships
H.mean <- matrix(0,M.beta,r)
diag(H.mean) <- 1
Pmat <- diag(1, M.beta)

# Starting values and OLS quantities (if possible)
beta_draw  <- H.mean
bb.sqrt <- diag(1, r)

D <- D_draw<- Z%*%beta_draw
X.full <- cbind(D,X)

# Stack coefficients and regressors
A.full <- solve(crossprod(X.full))%*%crossprod(X.full,Y)
alpha <- t(A.full[1:r,,drop=FALSE]) 
A_draw <- A.full[(r+1):(K+r),]      
Em <- Em.str <-  Y-Z%*%beta_draw%*%t(alpha) -X%*%A_draw
Sigma_draw <- crossprod(Em)/(nrow(Y)-K) #Var-Cov-matrix

Ups_draw <- t(chol(Sigma_draw))/sqrt(diag(Sigma_draw))
diag(Ups_draw) <- 1
omega <- diag(KK)*0.01
omega.mat <- matrix(0.01,KK,M)
zeta.mat <- Ups_draw

scale.start <- 1
theta <- matrix(scale.start,KK,M)
theta.cov  <- lower.tri(Sigma_draw)*1
theta.omega <- matrix(scale.start,KK,M)
theta.zeta <- lower.tri(Sigma_draw)*1

At.full <- array(0,c(T,KK,M))
PI.t_draw <- array(0,c(T,M.beta,M))
PHI.t_draw <- array(0,c(T,K.var,M))
fit_draw <- fit.level_draw <- matrix(0, T, M)

XAt.full1 <- matrix(NA,T,KK*M)
At.full.tilde <- array(NA,c(T,KK,M))
Sig.t_draw <- Ups.t_draw <- Ups.tilde <- array(0,c(T,M,M))
EmUps <- matrix(NA,T,M*(M-1)/2)
for (t in 1:T) Ups.t_draw[t,,] <- diag(M); Ups.tilde[t,,] <- diag(M)

# Storage matrices
A_store <- array(NA,c(nsave, K,M))
A.T_store <- array(NA,c(nsave, K,M))
alpha.T_store <- array(NA, c(nsave, M, r))
beta_store <- array(NA, c(nsave, M.beta, r))
PI_store <- array(NA,c(nsave, M.beta,M))
PI.T_store <- array(NA, c(nsave, M.beta, M))
PHI.T_store <-  array(NA,c(nsave, K.var,M))
sig.T_store <- array(NA,c(nsave , M,M))
Ups.T_store <- array(0, c(nsave, M,M))
Hvar_store <- array(0, c(nsave, T, M))
Em_store <- array(0, c(nsave, T, M))
omega.A_store <- array(0, c(nsave, K, M))
omega.alpha_store <- array(0, c(nsave, M, r))

###--------------------------------------------------------------------------###
###------ START: Markov chain Monte Carlo (MCMC) sampling algorithm ---------###
###--------------------------------------------------------------------------###

to.plot <- seq(1, ntot, 50)
pb <- txtProgressBar(min = 0, max = ntot, style = 3) #start progress bar
start <- Sys.time()
irep <- 1
for (irep in 1:ntot){
###--- Step I: Sample autoregressive parameters
Y_ <- Y
ind.em <- 1

for (mm in 1:M){
# Sample autoregressive coefficients equation by equation
if (mm==1){
  Y.i <- Y_[,mm]-(X.full%*%A.full)[,mm]
  X.i <- X.full%*%diag(omega.mat[,mm])
  a.prior <- A.prior[,mm]
# Draw for the time-varying part At.full
  At.full1 <- KF_fast(t(as.matrix(Y.i)), X.i,as.matrix(tau.fat[,mm]*exp(Hvar[,mm])),t(matrix(1,KK,T)),KK, 1, T, matrix(0,KK,1), diag(KK)*1e-15)
  G <- X.full*t(At.full1)
# Stack X.full and X.full*beta.tilde and normalize with time-varying variances
  normalizer <- as.numeric(exp(-Hvar[,mm]/2)/sqrt(tau.fat[,mm]))
  Xnew <- cbind(X.full,G)*normalizer
  Ynew <- Y[,mm]*normalizer
  
  V.prior <- c(theta[,mm],theta.omega[,mm])
  A.full.i <- get.A(Ynew,Xnew,a.prior,V.prior)
      
  A.full[,mm] <- A.full.i[1:KK]
  diag(omega) <- A.full.i[(KK+1):(2*KK)]
      
# Introduce non-identified element; introduce parameter-expansion, since +/- sqrt(sigma_2) 
  if (runif(1,0,1) > 0.5){
        omega <- -omega
        At.full1 <- -At.full1
  }
  omega.mat[,mm] <- diag(omega)
  for (t in 1:T) Em[t,mm] <- Em.str[t,mm] <- Y[t,mm]-X.full[t,]%*%(as.numeric(A.full[,mm])+diag(omega)*At.full1)[,t]
  At.full[,,mm] <- t((as.numeric(A.full[,mm])+diag(omega)*At.full1))
  XAt.full1[,(KK*(mm-1)+1):(KK*mm)] <-  X.full*t(At.full1)
  At.full.tilde[,,mm] <- t(At.full1)
}else{
  Y.i <- Y_[,mm]-(X.full%*%A.full)[,mm]
  Em.i <- Em[,1:(mm-1),drop = F]
  X.tilde <- cbind(X.full,Em.i)
  Ktilde <- ncol(X.tilde)
  X.i <- X.tilde%*%diag(c(omega.mat[,mm],zeta.mat[mm,1:(mm-1)]))
  a.prior <- matrix(c(A.prior[1:KK,mm],rep(0, Ktilde- KK), A.prior[(KK+1):(2*KK),mm],rep(0, Ktilde- KK)), 2*Ktilde,1)
# Draw for the time-varying part At.full
  At.full1 <- KF_fast(t(as.matrix(Y.i)), X.i,as.matrix(tau.fat[,mm]*exp(Hvar[,mm])),t(matrix(1,Ktilde,T)),Ktilde, 1, T, matrix(0,Ktilde,1), diag(Ktilde)*1/10^15)
  G <- X.tilde*t(At.full1)
# Stack X.full and X.full*beta.tilde and normalize with time-varying variances
  normalizer <- as.numeric(exp(-Hvar[,mm]/2)/sqrt(tau.fat[,mm]))
  Xnew <- cbind(X.tilde,G)*normalizer
  Ynew <- Y_[,mm]*normalizer
  V.prior <- c(theta[,mm],theta.cov[mm,1:(mm-1)],theta.omega[,mm],theta.zeta[mm,1:(mm-1)])
  A.full.i <- get.A(Ynew,Xnew,a.prior,V.prior)
      
# Indicator used to select the appropriate elements of A_draw
  ind_sl_A <- 1:KK
  ind_sl_cov <- (KK+1):Ktilde
  ind_sl_omega <- (Ktilde+1):(2*KK+Ktilde-KK)
  ind_sl_zeta <- (2*KK+Ktilde-KK+1):ncol(Xnew)
      
  A.full[,mm] <- A.full.i[ind_sl_A]
  diag(omega) <- A.full.i[ind_sl_omega] 
  cov.i <- (A.full.i[ind_sl_cov])
  zeta <- A.full.i[ind_sl_zeta] 
  
  if(runif(1,0,1)>0.5){
    omega <- -omega
    zeta <- -zeta
    At.full1 <- -At.full1
  }
      
  omega.mat[,mm] <- diag(omega)
  zeta.mat[mm,1:(mm-1)] <- zeta
      
  for (t in 1:T){
        Em[t,mm] <- Y_[t,mm]-X.full[t,]%*%(as.numeric(A.full[,mm])+diag(omega)*At.full1[ind_sl_A,t])
        Em.str[t,mm] <- Y_[t,mm]-X.full[t,]%*%(as.numeric(A.full[,mm])+diag(omega)*At.full1[ind_sl_A,t])-Em[t,1:(mm-1),drop=FALSE]%*%(as.numeric(cov.i)+zeta*(At.full1[ind_sl_cov,t]))
  }
      
  Ups_draw[mm,1:(mm-1)] <- cov.i
  At.full[,,mm] <- t(as.numeric(A.full[,mm])+diag(omega)*At.full1[ind_sl_A,])
  Ups.t_draw[,mm,1:(mm-1)] <- t(as.numeric(cov.i)+zeta*At.full1[ind_sl_cov,])
  XAt.full1[,(KK*(mm-1)+1):(KK*mm)] <-  X.full*t(At.full1[ind_sl_A,])
  At.full.tilde[,,mm] <- t(At.full1[ind_sl_A,])
  EmUps[,ind.em:(ind.em+ncol(Em.i)-1)] <- Em.i*t(At.full1[ind_sl_cov,,drop = F])
  Ups.tilde[,mm,1:(mm-1)] <- t(At.full1[ind_sl_cov,])
  ind.em <- ind.em+ncol(Em.i)
  }
}
  
# alpha.tilde and alpha.t.tilde denote a* and a*_t (see Koop and Jochmann (2015) appendix)
alpha.tilde <- t(A.full[1:r,])
alpha.t.tilde <- aperm(At.full[,1:r,], c(1,3,2))
A_draw <- A.full[(r+1):KK,]
A.t_draw <- At.full[,(r+1):KK,]

omega.alpha.mat <- t(omega.mat[1:r,])
omega.A.mat <- omega.mat[(r+1):KK,]
  
###--- Step II: Sample Variances (for factors idiosyncratic ones)
for (jj in 1:M){
  Em.str[abs(Em.str[,jj])<1e-10,jj] <- 1e-10
  sv.idio.i <- sv.idio[[jj]]
  sv.idio.i <- svsample_fast_cpp(Em.str[,jj], startpara = sv.idio.i, startlatent = ht[,jj], priorspec = sv_priors)
  sv.idio.i[c("mu", "phi", "sigma", "nu", "rho")] <- as.list(sv.idio.i$para[, c("mu", "phi", "sigma", "nu", "rho")])
  tau.fat.i <- as.numeric(sv.idio.i$tau) 
  ht.i <- t(sv.idio.i$latent)
  ht.i[ht.i < -12] <- -12 
  ht[,jj]      <- ht.i
  tau.fat[,jj] <- tau.fat.i
  Hvar[,jj]    <- ht.i + log(tau.fat.i)
  pars_var[,jj] <- as.numeric(sv.idio.i$para[,c("mu", "phi", "sigma", "nu")])
  sv.idio[[jj]] <- sv.idio.i
}

for (t in 1:T) Sig.t_draw[t,,] <- (Ups.t_draw[t,,])%*%diag(exp(Hvar[t,]))%*%t(Ups.t_draw[t,,])
  
###--- Step III: Sample HS prior
# Constant part of VAR coefficients
hs_draw <- get.hs(bdraw=as.numeric(A.full),lambda.hs=lambda.A,nu.hs=nu.A,tau.hs=tau.A,zeta.hs=zeta.A)
theta <- matrix(hs_draw$psi,KK,M)
lambda.A <- hs_draw$lambda
nu.A <- hs_draw$nu
tau.A <- hs_draw$tau
zeta.A <- hs_draw$zeta

# State innovations of VAR coefficients  
hs_draw <- get.hs(bdraw=as.numeric(omega.mat),lambda.hs=lambda.o,nu.hs=nu.o,tau.hs=tau.o,zeta.hs=zeta.o)
theta.omega <- matrix(hs_draw$psi,KK,M)
lambda.o <- hs_draw$lambda
nu.o <- hs_draw$nu
tau.o <- hs_draw$tau
zeta.o <- hs_draw$zeta
   
# Constant part of covariances 
hs_draw <- get.hs(bdraw=as.numeric(Ups_draw[id_free_cov]),lambda.hs=lambda.cov,nu.hs=nu.cov,tau.hs=tau.cov,zeta.hs=zeta.cov)
theta.cov[id_free_cov] <- hs_draw$psi
lambda.cov <- hs_draw$lambda
nu.cov <- hs_draw$nu
tau.cov <- hs_draw$tau
zeta.cov <- hs_draw$zeta
    
# State innovations of covariances   
hs_draw <- get.hs(bdraw=as.numeric(zeta.mat[id_free_cov]),lambda.hs=lambda.zeta,nu.hs=nu.zeta,tau.hs=tau.zeta,zeta.hs=zeta.zeta)
theta.zeta[id_free_cov] <- hs_draw$psi
lambda.zeta <- hs_draw$lambda
nu.zeta <- hs_draw$nu
tau.zeta <- hs_draw$tau
zeta.zeta <- hs_draw$zeta

###--- Step IV: Cointegration
Y.tilde <- matrix(0, T,r)
alpha.t_draw <- A.Sigt <- array(0, c(T, M, r))
for(tt in 1:T){
  A.Sigt[tt,,] <- chol(solve(Sig.t_draw[tt,,]))%*%alpha.t.tilde[tt,,]
  Y.tilde[tt,] <- (Y[tt,] - X[tt,]%*%A.t_draw[tt,,])%*%solve(Sig.t_draw[tt,,])%*%alpha.t.tilde[tt,,]
}

alphaZ <- matrix(0, M.beta*T, M.beta*r)
for(ii in 1:M){
  for(jj in 1:r){  
    alphaZ[((ii-1)*T+1):(ii*T),((jj-1)*M.beta+1):(jj*M.beta)] <- Z*A.Sigt[,ii,jj]
  }
}
# Draw b* and update cointegration part 
Vp.beta <- solve(crossprod(alphaZ) + solve(kronecker(diag(1,r), Pmat)))
Ap.beta <- Vp.beta%*%as.vector(crossprod(Z,Y.tilde))
beta.tilde <- try(Ap.beta+t(chol(Vp.beta))%*%rnorm(ncol(Vp.beta)),silent=TRUE)
if (is(beta.tilde,"try-error")) beta.tilde <- mvrnorm(1,Ap.beta,Vp.beta)
beta.tilde <- matrix(beta.tilde, M.beta, r)
bb.sqrt <- eigen(crossprod(beta.tilde))
bb.sqrt <- bb.sqrt$vectors%*%diag(sqrt(abs(bb.sqrt$values)))%*%t(bb.sqrt$vectors)
    
beta_draw <- beta.tilde%*%ginv(bb.sqrt)
alpha <- alpha.tilde%*%bb.sqrt 
D_draw <- Z%*%beta.tilde
X.full <- cbind(D_draw,X)

PI_draw <- beta.tilde%*%t(alpha.tilde)
  
for(ttt in 1:T){
  alpha.t_draw[ttt,,] <- alpha.t.tilde[ttt,,]%*%bb.sqrt
  PI.t_draw[ttt,,] <- beta.tilde%*%t(alpha.t.tilde[ttt,,])
  fit_draw[ttt,] <- Z[ttt,]%*%beta.tilde%*%t(alpha.t.tilde[ttt,,]) + X[ttt,]%*%A.t_draw[ttt,,]
  PHI.t_draw[ttt,,] <- VECMtoVAR(PI = PI.t_draw[ttt,,], A = A.t_draw[ttt,,], p = p, M = M, K = K)
  fit.level_draw[ttt,] <- X.levels[ttt,]%*%PHI.t_draw[ttt,,] 
}  
  
###--- Step V: Storage
if(irep %in% save.set){
  jrep <- jrep +1
  Hvar_store[jrep,,] <- Hvar + log(tau.fat)
  Em_store[jrep,,] <- Em
    
  A_store[jrep,,] <- A_draw
  PI_store[jrep,,] <- PI_draw
  A.T_store[jrep,,] <- A.t_draw[T,,]
  PI.T_store[jrep,,] <- PI.t_draw[T,,]
  PHI.T_store[jrep,,] <- PHI.t_draw[T,,]
  alpha.T_store[jrep,,] <- alpha.t.tilde[T,,]
  beta_store[jrep,,] <- beta.tilde
  sig.T_store[jrep,,] <- Sig.t_draw[T,,] 
  Ups.T_store[jrep,,] <- Ups.t_draw[T,,]
  omega.A_store[jrep,,] <- omega.A.mat
  omega.alpha_store[jrep,,] <- omega.alpha.mat
  
}

if (irep %in% to.plot){
  par(mfrow = c(3,4))
  for(pp in 1:M) ts.plot(cbind(fit.level_draw[,pp], Yraw[,pp]), col = c(1,2), lty = c(1,1), main = irep)
}
setTxtProgressBar(pb, irep)

}

end <- Sys.time()
time.min <- (ts(end)-ts(start))/60

message(paste0("Posterior sampling of ", ntot, " draws completed in: ", round(time.min, 3), " mins."))

###--------------------------------------------------------------------------###
###-------- END: Markov chain Monte Carlo (MCMC) sampling algorithm ---------###
###--------------------------------------------------------------------------###

###--------------------------------------------------------------------------###
###------ START: Ex-post sparsification based on obtained MCMC draws --------###
###--------------------------------------------------------------------------###
# Storage
pred.sps_store <- array(NA,c(nsave,M,1))

irep <- 1
for (irep in 1:nsave){
  
###--- Step I: Forecast
A_draw <- A_store[irep,,]
A.T_draw <- A.T_store[irep,,]
PI_draw <- PI_store[irep,,]
PI.T_draw <- PI.T_store[irep,,]
sig.T_draw <- sig.T_store[irep,,]
  
# Sparse coefficients
A.T_sps  <- sparsify.A(X = X, Alpha.draw = A.T_draw, M = M, K = K, kappa = 2, lambda = 1)
PI.T_sps <- sparsify.pi(X = Z, Alpha.draw = PI.T_draw, M = M, K = M.beta, kappa = 2, lambda = 1)
PHI.T_sps  <- VECMtoVAR(PI = PI.T_sps, A = A.T_sps, p = p, M = M, K = K)
 
# Sparse variance covariances
sig.pen <- 0.01/abs(sig.T_draw)^(0.5)
sig.T_sps <- glasso::glasso(sig.T_draw, rho = sig.pen, maxit = 10)$w
  
# Companion of sparse coefficients
get.comp <- get.companion(Beta_= PHI.T_sps[1:(M*p),],varndxv=c(M,0,p))
MM <- get.comp$MM
Jm <- get.comp$Jm

Sig.0 <- Jm%*%sig.T_sps%*%t(Jm)
cholSig.0 <- try(t(chol(Sig.0[1:M,1:M])),silent=TRUE)

# x_t+1
z.tilde <- c(Yraw[T,], X.levels[T,1:(M*(p-1))]) 
X0 <- as.matrix(z.tilde)
X0 <- MM%*%X0 
if(is(cholSig.0,"try-error")){
  pred_sps <- rmvnorm(1,mean=X0[1:M],sigma=Sig.0[1:M,1:M]) 
}else{
  pred_sps <- X0[1:M]+cholSig.0%*%rnorm(M,0,1)
}    

pred.sps_store[irep,,] <- pred_sps

}

###--------------------------------------------------------------------------###
###------ START: Ex-post sparsification based on obtained MCMC draws --------###
###--------------------------------------------------------------------------###

# Save forecast
save(file = paste0(resdir, model, "_pred-", end.in, ".rda"), 
      list = c("pred.sps_store",
               "time.min"))