mlag <- function(X,lag)
{
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(0,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
  }
  return(Xlag)
}

get.A <- function(y,x,A_pr,V_pr){
  V_post <- try(chol2inv(chol(crossprod(x)+diag(1/V_pr))),silent=TRUE)
  if (is(V_post,"try-error")) V_post <- ginv(crossprod(x)+diag(1/V_pr))
  A_post <- V_post%*%(crossprod(x,y)+diag(1/V_pr)%*%A_pr)
  
  A_draw <- try(A_post+t(chol(V_post))%*%rnorm(ncol(x)),silent=TRUE)
  if (is(A_draw,"try-error")) A_draw <- mvrnorm(1,A_post,V_post)
  return(A_draw)
}

get.hs <- function(bdraw,lambda.hs,nu.hs,tau.hs,zeta.hs){
  k <- length(bdraw)
  
  lambda.hs <- invgamma::rinvgamma(k,shape=1,rate=1/nu.hs+bdraw^2/(2*tau.hs))
  tau.hs <- invgamma::rinvgamma(1,shape=(k+1)/2,rate=1/zeta.hs+sum(bdraw^2/lambda.hs)/2)
  nu.hs <- invgamma::rinvgamma(k,shape=1,rate=1+1/lambda.hs)
  zeta.hs <- invgamma::rinvgamma(1,shape=1,rate=1+1/tau.hs)
  
  ret <- list("psi"=(lambda.hs*tau.hs),"lambda"=lambda.hs,"tau"=tau.hs,"nu"=nu.hs,"zeta"=zeta.hs)
  return(ret)
}

VECMtoVAR <- function(PI = PI, A = A, p = p, M = M, K = K){
  p.tilda <- p + 1
  PHI.t <- matrix(NA,M*p.tilda,M)
  for (ll in rev(seq_len(p.tilda))){
  if (ll==p.tilda){
      PHI.00 <- - A[((p-1)*M+1):(p*M),]
    }else if (ll>1 && ll <= p){
      PHI.00 <- A[((ll-1)*M+1):((ll)*M),]- A[((ll-2)*M+1):((ll-1)*M),]
    }else if(ll == 1){
      PHI.00 <- PI + diag(M) + A[((ll-1)*M+1):((ll)*M),]
    }
    PHI.t[((ll-1)*M+1):(ll*M),] <- PHI.00
  }
  
  return(PHI.t)
}

  
norm_vec <- function(x) sqrt(sum(x^2))

sparsify.A <- function(X,Alpha.draw, M,K,kappa, lambda){
  A.draw <- matrix(0, K,M)

  norm.i <- rep(apply(X,2,norm_vec)^2, M)
    
  A.draw.i <- Alpha.draw.i <-  as.vector(Alpha.draw) 
  mu.i <- log(lambda) - kappa*(log(abs(Alpha.draw.i)))
  mu.i <- exp(mu.i)
    
  ind.change <- (abs(Alpha.draw.i)*norm.i)>mu.i 
  A.draw.i[ind.change] <- (sign(Alpha.draw.i)*1/norm.i * ((abs(Alpha.draw.i)*norm.i)-mu.i))[ind.change]
  A.draw.i[!ind.change] <- 0
  A.draw <- matrix(A.draw.i, K, M)
  return(A.draw)
}

sparsify.pi <- function(X,Alpha.draw, M, K, kappa, lambda){
  A.draw <- matrix(0,K,M)

  norm2.Z <- apply(X,2,norm_vec)^2/nrow(X)
  norm.A  <- apply(Alpha.draw,1,norm_vec)
    
  A.draw.i <- Alpha.draw.i <-  as.matrix(Alpha.draw)
  mu.i <- log(lambda) - (kappa)*log(norm.A)
  mu.i <- exp(mu.i)/2
  ind.change <- (norm.A*norm2.Z) > mu.i
  A.draw[!ind.change,] <- 0
  A.draw[ind.change,] <- (Alpha.draw*(1- mu.i/(norm2.Z*norm.A)))[ind.change,]
    
  return(A.draw)
}

get.companion <- function(Beta_,varndxv){
  nn <- varndxv[[1]]
  nd <- varndxv[[2]]
  nl <- varndxv[[3]]
  
  nkk <- nn*nl+nd
  
  Jm <- matrix(0,nkk,nn)
  Jm[1:nn,1:nn] <- diag(nn)
  
  if(nd < 1){
    if(nl == 1){
      MM <- t(Beta_)
    }else{
      MM <- rbind(t(Beta_),cbind(diag((nl-1)*nn), matrix(0,(nl-1)*nn,nn)))
    }
  }else{
   if(nl == 1){
      MM <- rbind(t(Beta_), cbind(matrix(0,nd,nn),diag(1,nd)))
    }else{
      MM <- rbind(t(Beta_),cbind(diag((nl-1)*nn), matrix(0,(nl-1)*nn,nn+nd)),cbind(matrix(0,nd,(nn*nl)),diag(1,nd)))
    }
  }
  return(list(MM=MM,Jm=Jm))
}