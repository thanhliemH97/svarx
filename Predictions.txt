library(MTS)
library(mvtnorm)

# === Commandes prévisions ===
matApower.n <- function(A, n) {
    output <- diag( nrow(A) )
    if(n > 0) {
      for(i in 1:n) {
        output <- A %*% output 
      }
    }
    output
}

mat.subset <- function(MaMat,ind) {
# La matrice MaMat est de dimension n x (mp).
# Autrement dit, p matrices n x m sont les blocs.
# On note que l'on n'a pas besoin de connaître n.
# Le vecteur ind doit être de dimension p, et est constitué de 1 et de 0
# avec 1 pour inclus, et 0 pour non-inclus
  p <- length(ind)
  m <- ncol(MaMat)/p
  output <- NULL
  for(i in 1:p) {
     if(ind[i] == 1) {
       output <- cbind(output, MaMat[,((i-1)*m +1):(i*m)])
     }
  }
  output
}


# --- Fonctions inutiles ---

matC.13 <- function(p,k) {
           rbind( 
                t(mat.Ei(1,p,k)),
                t(mat.Ei(3,p,k))
                )
}

matD.02 <- function(s,kx) {
           rbind( 
                t(mat.Fj(1,s,kx)),
                t(mat.Fj(3,s,kx))
                )
}

# ----------------------------


SigmaSpherique <- function(k, rho) {
 vec.one <- rep(1,k)
 (1-rho)*diag(k) + rho* vec.one %*% t(vec.one)
}

moyenneVARX <- function(p, matPhi, theta0) {
  k <- nrow(matPhi)
  muz <- rep(0,k)
  matSum <- matrix(0, nrow=k, ncol=k)
  for(i in 1:p) {
     matSum <- matSum + matPhi[,((i-1)*k+1):(i*k)]
  }
  muz <- solve( diag(k) - matSum, theta0)
  muz
}

initINDneg <- function(ind, tableau, mu) {
  if(ind > 0) output <- tableau[ind,]
  else output <- mu
  output
}

sumLastLags <- function(ind, p, matA, matX, muX) {
  k <- nrow(matA)
  kx <- ncol(matX)
  sum <- rep(0,k)
  for(i in 1:p) {
    sum <- sum + matA[,((i-1)*kx+1):(i*kx)] %*% initINDneg(ind-i, matX, muX)
  }
  sum
}

simulVARX <- function(n, theta0, 
                      p, matPhi,
                      B0, m, matB,
                      xt, SigmaEps) {
  kx <- ncol(xt)
  mux <- rep(0, kx)
  muz <- moyenneVARX(p, matPhi, theta0)
  k <- nrow(matPhi)
  zt <- matrix(0, nrow=n, ncol=k)
  mesEps <- rmvnorm(n, mean= rep(0, k), sigma = SigmaEps)
  for(i in 1:n) {
    Ui <- sumLastLags(i, p, matPhi, zt, muz)
    Vi <- sumLastLags(i, m, matB, xt, mux)
    zi <- theta0 + B0 %*% xt[i,] + Ui + Vi + mesEps[i,]
    zt[i,] <- zi
  }
  zt
}



matR.VARXp3m2 <- function(k=3, kx=3, ind=c(1,1), indx=c(1,1) ) {
#
# Le modèle postule p=3 maximum pour la partie AR et m=2
# ce qui veut dire que c'est un VARX avec une partie instantanée
# et deux retards. La simulation aura une partie endogène de dimension k
# et exogène de dimension kx.
#
# Il est supposé que dans ind on trouve des 0 (non-inclus) et des 1 (inclus).
#
# Il est supposé un VARX(3,2), ce qui veut dire que le délai 3 est inclus tout le temps,
# et le retard exogène de délai deux est également inclus tout le temps. 
# Les indices intermédiaires peuvent être cependant retirés.
# 
# L'intérêt de la simulation est de ne pas se questionner sur l'identification des indices maximums.
# On suppose que l'on connait (p,m), mais on cherche à voir si on identifie les bons sous-indices.
# On note que la matrice R commence par mettre des 1 partout (et pour p=3, m=2, ils resteront là).
#
#
  p <- 3
  pm1 <- p-1
  m <- 2
  mon.nrow <- 1+p*k + (1+m)*kx
  mon.ncol <- k
  matR <- matrix(1, nrow = mon.nrow, ncol=mon.ncol)
  for(i in 1:pm1) {
    matR[(k*(i-1)+2):(k*i+1),] <- ind[i]
  }
  for(i in 1:m) {
    matR[(p*k+(i-1)*kx+2):(p*k+i*kx+1),] <- indx[i]
  }
  matR
}


matPhiI <- function(matPhi,ind) {
# La matrice matPhi est de dimension k x (kp), avec k le nb de variables
# et p l'ordre autorégressif.
# Le vecteur ind doit être de dimension p, et est constitué de 1 et de 0
# avec 1 pour inclus, et 0 pour non-inclus
  output <- NULL
  k <- nrow(matPhi)
  p <- as.integer(  ncol(matPhi)/k )
  for(i in 1:p) {
     if(ind[i] == 1) {
       output <- cbind(output, matPhi[,((i-1)*k +1):(i*k)])
     }
  }
  output
}

matPhiIstar <- function(matPhi,ind) {
  k <- nrow(matPhi)
  p <- as.integer(  ncol(matPhi)/k )
  a <- sum(ind)
  ma.matPhiI <- matPhiI(matPhi, ind)
  rbind( ma.matPhiI, matrix(0, nrow=k*(p-1), ncol=k*a) )
}

matPhistar <- function(matPhi, ind) {
  k <- nrow(matPhi)
  p <- length(ind)
  a <- sum(ind)
  ma.matPhiIstar <- matPhiIstar(matPhi, ind)
  temp1 <- matrix(0, nrow=k, ncol=k*(p-1) )
  temp2 <- rbind( temp1, diag(k*(p-1)) )
  temp3 <- matrix(0, nrow=k*p, ncol=k)
  temp4 <- cbind( temp2, temp3)
  ma.matPhiIstar %*% matC.I(ind,p,k) + temp4
}


matBJ <- function(matB,ind) {
# La matrice matB est de dimension k x (kx * (s+1)), 
# avec k le nb de variables dans Y,
# avec kx le nb de variables dans X,
# et p l'ordre autorégressif.
# Le vecteur ind doit être de dimension p, et est constitué de 1 et de 0
# avec 1 pour inclus, et 0 pour non-inclus
  output <- NULL
  k <- nrow(matB)
  splus1 <- length(ind)
  s <- splus1 - 1
  kx <- as.integer( ncol(matB)/splus1 )
  for(i in 1:splus1) {
     if(ind[i] == 1) {
       output <- cbind(output, matB[,((i-1)*kx +1):(i*kx)])
     }
  }
  output
}

matBJstar <- function(matB, ind=c(1,0,1), p=3) {
  k <- nrow(matB)
  splus1 <- length(ind)
  s <- splus1 - 1
  b <- sum(ind)
  kx <- as.integer( ncol(matB)/splus1 )
  ma.matBJ <- matBJ(matB, ind)
  rbind( ma.matBJ, matrix(0, nrow=k*(p-1), ncol=kx*b) ) %*% matD.J(ind,s,kx)
}

Ystar <- function(matY, p) {
   n <- nrow(matY)
   k <- ncol(matY)
   output <- matrix(0, nrow=n, ncol=k*p)
   output[1:(p-1),] <- 0
   for(i in p:n) {
     for(j in 1:p) {
          output[i, ((j-1)*k +1):(j*k)] <- matY[i-j+1,]
     }
   }
   output
}

Xstar <- function(matX, s) {
   n <- nrow(matX)
   kx <- ncol(matX)
   output <- matrix(0, nrow=n, ncol=kx*(s+1) )
   output[1:s,] <- 0
   for(i in (s+1):n) {
     for(j in 1:(s+1)) {
          output[i, ((j-1)*kx +1):(j*kx)] <- matX[i-j+1,]
     }
   }
   output
}

predict.subsetVARX <- function(Ynstar, p, k, Xstar.futur, Phistar, Bstar, maxl=12) {
  output <- matrix(0, nrow=maxl, ncol=k)
  for(l in 1:maxl) {
    temp <- rep(0,k)
    for(j in 0:(l-1)) {
        temp <- temp + t(mat.Ei(1,p,k)) %*% matApower.n(Phistar,j) %*% Bstar %*% Xstar.futur[l-j,]
    }
    output[l, ] <- t(mat.Ei(1,p,k)) %*% matApower.n(Phistar,l) %*% Ynstar + temp
  }
  output 
}


predict3.subsetVARX <- function(Ynstar, p, k, Xstar.futur, Phistar, Bstar, maxl=12, 
                                contrast1 = c(1, rep(0,k-1)), 
                                contrast2 = c(rep(0,k-1),1),
                                contrast3 = rep(1/k,k)
                                ) {
  output <- matrix(0, nrow=maxl, ncol=k)
  for(l in 1:maxl) {
    temp <- rep(0,k)
    for(j in 0:(l-1)) {
        temp <- temp + t(mat.Ei(1,p,k)) %*% matApower.n(Phistar,j) %*% Bstar %*% Xstar.futur[l-j,]
    }
    output[l, ] <- t(mat.Ei(1,p,k)) %*% matApower.n(Phistar,l) %*% Ynstar + temp
  }
  output1 <- as.vector( output %*% contrast1 )
  output2 <- as.vector( output %*% contrast2 )
  output3 <- as.vector( output %*% contrast3 )
  list(output.punctual = cbind(output1,output2,output3), output=output)
}

predictIP.subsetVARX <- function(p,k, Phistar, Sigmahat, maxl=12, 
                                contrast1 = c(1, rep(0,k-1)), 
                                contrast2 = c(rep(0,k-1),1),
                                contrast3 = rep(1/k,k)
                                ) {
 outputMSE <- matrix(0, nrow=maxl, ncol=3)
 for(l in 1:maxl) {
     temp <- matrix(0, nrow=k, ncol=k)
     lm1 <- l-1
     for(j in 0:lm1) {
       temp <- temp + t(mat.Ei(1,p,k)) %*% matApower.n(Phistar,j) %*% mat.Ei(1,p,k) %*% 
               Sigmahat %*% t(mat.Ei(1,p,k)) %*% t( matApower.n(Phistar,j) ) %*% mat.Ei(1,p,k)
     }
     outputMSE[l, 1] <- t(contrast1) %*% temp %*% contrast1
     outputMSE[l, 2] <- t(contrast2) %*% temp %*% contrast2
     outputMSE[l, 3] <- t(contrast3) %*% temp %*% contrast3
 }
  outputMSE
}

predictIPcorr.subsetVARX <- function(p,k, 
                                Ynstar, Xstar.futur,
                                Phistar, BJstar, mon.Vhat, Sigmahat, maxl=12, 
                                contrast1 = c(1, rep(0,k-1)), 
                                contrast2 = c(rep(0,k-1),1),
                                contrast3 = rep(1/k,k)
                                ) {
 outputMSE <- matrix(0, nrow=maxl, ncol=3)
 for(l in 1:maxl) {
     temp <- matrix(0, nrow=k, ncol=k)
     lm1 <- l-1
     for(j in 0:lm1) {
       temp <- temp + t(mat.Ei(1,p,k)) %*% matApower.n(Phistar,j) %*% mat.Ei(1,p,k) %*% 
               Sigmahat %*% t(mat.Ei(1,p,k)) %*% t( matApower.n(Phistar,j) ) %*% mat.Ei(1,p,k)
     }
     MSE1 <- temp
     mon.Wnl.out <- mon.Wnl(Ynstar, Xstar.futur, l, k)
     matHIJ.out <- matHIJ.l(l, Phistar, BJstar, k=2,kx=2,p=3,s=2, ind=c(1,0,1), indx=c(1,0,1))
     TermeCorr.l <- mon.Wnl.out %*% matHIJ.out %*% mon.Vhat %*% t(matHIJ.out) %*% t(mon.Wnl.out)
     MSE <- MSE1 + TermeCorr.l
     outputMSE[l, 1] <- t(contrast1) %*% MSE %*% contrast1
     outputMSE[l, 2] <- t(contrast2) %*% MSE %*% contrast2
     outputMSE[l, 3] <- t(contrast3) %*% MSE %*% contrast3
 }
  outputMSE
}



vec.ei <- function(ind,p) {
   output <- rep(0,p)
   output[ind] <- 1
   output
}

vec.fj <- function(ind,s) {
   output <- rep(0,s+1)
   output[ind] <- 1
   output
}

mat.Ei <- function(ind,p,k) {
   kronecker( vec.ei(ind,p), diag(k) )
}

mat.Fj <- function(ind,s,kx) {
   kronecker( vec.fj(ind,s), diag(kx) )
}

matC.I <- function(ind, p, k) {
           output <- NULL
           for(i in 1:p) {
             if(ind[i]==1) {
                 output <- rbind(output, t(mat.Ei(i,p,k)) )
             }
           }
          output
}

matD.J <- function(ind, s, kx) {
           output <- NULL
           splus1 <- s + 1
           for(i in 1:splus1) {
             if(ind[i]==1) {
                 output <- rbind(output, t(mat.Fj(i,s,kx)) )
             }
           }
          output
}

mon.Wnl <- function(Ynstar, MonXstar.futur, ll, k) {
  output <- Ynstar
  for(i in 1:ll) {
    output <- c(output, MonXstar.futur[ll-i+1,])
  }
  kronecker(t(as.vector(output)), diag(k) )
}


poids.Psii <- function(puis, p, k, ma.matPhistar) {
  t(mat.Ei(1,p,k)) %*% matApower.n(ma.matPhistar,puis) %*% mat.Ei(1,p,k)
}

bloc.C <- function(ind,k,s,kx) {
   kronecker( t(matD.J(ind,s,kx)), diag(k) )
}

bloc.B2l <- function(puis,ind,k,p,s,kx,ma.matPhistar) {
   kronecker( t(matD.J(ind,s,kx)), poids.Psii(puis, p, k, ma.matPhistar) )
}

bloc.Al <- function(ll, ma.matPhistar,ind=c(1,0,1),p=3,k=2) {
   aa <- sum(ind)
   output <- matrix(0, nrow=p*k*k, ncol=k*k*aa)
   llmoinsun <- ll- 1
   for(j in 0:llmoinsun) {
       temp <- matApower.n( t(ma.matPhistar), ll-1-j) %*% t(matC.I(ind,p,k))
       Poidsstarj <- poids.Psii(j, p, k, ma.matPhistar)
       output <- output + kronecker(temp, Poidsstarj)
   }
   output
}

bloc.B1l <- function(ll, ma.matPhistar, ma.matBstar, ind=c(1,0,1),p=3,k=2,s=2,kx=2) {
   aa <- sum(ind)
   splus1 <- s+1
   output <- matrix(0, nrow=splus1*k*kx, ncol=k*k*aa)
   llmoinsun <- ll- 1
   for(j in 0:llmoinsun) {
       temp <- t(ma.matBstar) %*% matApower.n( t(ma.matPhistar), ll-1-j) %*% t(matC.I(ind,p,k))
       Poidsstarj <- poids.Psii(j, p, k, ma.matPhistar)
       output <- output + kronecker(temp, Poidsstarj)
   }
   output
}

matHIJ.l1 <- function(ma.matPhistar, k=2,kx=2,p=3,s=2, ind=c(1,0,1), indx=c(1,0,1)) {
 ll <- 1
 sp1 <- s + 1
 a <- sum(ind)
 b <- sum(indx)
 output <- matrix(0, nrow=k*k*p+ll*k*kx*sp1, ncol=k*k*a+k*kx*b)
 output[1:(k*k*p),1:(a*k*k)] <- bloc.Al(1,ma.matPhistar,ind,p,k)
 output[(k*k*p+1):(k*k*p+k*kx*sp1), (a*k*k+1):(a*k*k+b*k*kx)] <- bloc.C(indx,k,s,kx)
 output
}

matHIJ.l6 <- function(ma.matPhistar, ma.matBJstar, k=2,kx=2,p=3,s=2, ind=c(1,0,1), indx=c(1,0,1)) {
 ll <- 6
 llm1 <- ll - 1
 sp1 <- s + 1
 a <- sum(ind)
 b <- sum(indx)
 output <- matrix(0, nrow=k*k*p, ncol=k*k*a+k*kx*b)
 output[1:(k*k*p),1:(a*k*k)] <- bloc.Al(ll,ma.matPhistar,ind,p,k)
 for(j in 1:llm1) {
   temp1 <- bloc.B1l(ll-j, ma.matPhistar, ma.matBJstar, ind=c(1,0,1),p=3,k=2,s=2,kx=2)
   temp2 <- bloc.B2l(ll-j, ind,k,p,s,kx,ma.matPhistar)
   temp <- cbind(temp1, temp2)
   output <- rbind(output,temp)
 }
 temp <- matrix(0, nrow=k*kx*sp1, ncol=k*k*a+k*kx*b)
 temp[,(k*k*a+1):(a*k*k+b*k*kx)] <- bloc.C(indx,k,s,kx)
 output <- rbind(output, temp)
 output
}


matHIJ.l <- function(ll, ma.matPhistar, ma.matBJstar, k=2,kx=2,p=3,s=2, ind=c(1,0,1), indx=c(1,0,1)) {
 llm1 <- ll - 1
 sp1 <- s + 1
 a <- sum(ind)
 b <- sum(indx)
 if(ll==1) {
  output <- matrix(0, nrow=k*k*p+ll*k*kx*sp1, ncol=k*k*a+k*kx*b)
  output[1:(k*k*p),1:(a*k*k)] <- bloc.Al(1,ma.matPhistar,ind,p,k)
  output[(k*k*p+1):(k*k*p+k*kx*sp1), (a*k*k+1):(a*k*k+b*k*kx)] <- bloc.C(indx,k,s,kx)
 }
 if(ll > 1) {
   output <- matrix(0, nrow=k*k*p, ncol=k*k*a+k*kx*b)
   output[1:(k*k*p),1:(a*k*k)] <- bloc.Al(ll,ma.matPhistar,ind,p,k)
   for(j in 1:llm1) {
     temp1 <- bloc.B1l(ll-j, ma.matPhistar, ma.matBJstar, ind=c(1,0,1),p=3,k=2,s=2,kx=2)
     temp2 <- bloc.B2l(ll-j, ind,k,p,s,kx,ma.matPhistar)
     temp <- cbind(temp1, temp2)
     output <- rbind(output,temp)
   }
   temp <- matrix(0, nrow=k*kx*sp1, ncol=k*k*a+k*kx*b)
   temp[,(k*k*a+1):(a*k*k+b*k*kx)] <- bloc.C(indx,k,s,kx)
   output <- rbind(output, temp)
 } 
 output
}

Vhat.VARX <- function(zt, xt, n, Sigmahat,s=2,k=2,ind=c(1,0,1),indx=c(1,0,1)) {
     mon.Ystar <- Ystar(zt, p=3)
     matU.I <- matC.I( ind, p=3, k) %*% t(mon.Ystar[1:(n-1),])
     matV.J <- matD.J( indx, s=2, 2) %*% t(Xstar(xt[2:n,],s=2))
     matZ <- rbind( matU.I, matV.J)
     Vhat <- kronecker( solve( matZ %*% t(matZ) ), Sigmahat)
     Vhat
}



SimulationVARXp3m2.predict <- function(n=128, p=3, m=2, 
                               PhiX, SigmaX, theta0=theta0, matPhi=matPhi, B0=B0, matB=matB, SigmaEps=SigmaEps, alpha=0.05,
                               Nsim=10) {
# Pour la partie exogène on présume un VAR(1). Les paramètres sont PhiX et SigmaX.
# Le processus est un VARX(3,2), avec Phi2 = 0, et B1=0. Ainsi les indices du vrai modèle sont c(1,0,1,1,0,1)
# Le processus exogène est généré à chaque fois.
# 
# Pour les prévisions, on suppose que le bon modèle a été trouvé, ce qui va illustrer la théorie.
# 
set.seed(1)
N <- 2*n + 1 + 12
mon.k <- nrow(matPhi)
mon.kx <- nrow(PhiX)
maxl <- 12
mon.p <- p
mon.m <- m
q.nivconf <- qnorm(1-alpha/2)
# 
# == Matrices pour les EQM de prevision, pas de correction ==
# 
outputIP.l1.lower  <- matrix(0, nrow=Nsim, ncol=3)
outputIP.l1.upper  <- matrix(0, nrow=Nsim, ncol=3)
outputnivIP.l1     <- matrix(0, nrow=Nsim, ncol=3)
outputIP.l6.lower  <- matrix(0, nrow=Nsim, ncol=3)
outputIP.l6.upper  <- matrix(0, nrow=Nsim, ncol=3)
outputnivIP.l6     <- matrix(0, nrow=Nsim, ncol=3)
outputIP.l12.lower <- matrix(0, nrow=Nsim, ncol=3)
outputIP.l12.upper <- matrix(0, nrow=Nsim, ncol=3)
outputnivIP.l12    <- matrix(0, nrow=Nsim, ncol=3)
# 
# == Matrices pour les EQM de prevision, correction a l'ordre 1/n ==
# 
outputIPcorr.l1.lower  <- matrix(0, nrow=Nsim, ncol=3)
outputIPcorr.l1.upper  <- matrix(0, nrow=Nsim, ncol=3)
outputnivIPcorr.l1     <- matrix(0, nrow=Nsim, ncol=3)
outputIPcorr.l6.lower  <- matrix(0, nrow=Nsim, ncol=3)
outputIPcorr.l6.upper  <- matrix(0, nrow=Nsim, ncol=3)
outputnivIPcorr.l6     <- matrix(0, nrow=Nsim, ncol=3)
outputIPcorr.l12.lower <- matrix(0, nrow=Nsim, ncol=3)
outputIPcorr.l12.upper <- matrix(0, nrow=Nsim, ncol=3)
outputnivIPcorr.l12    <- matrix(0, nrow=Nsim, ncol=3)
#
# ==========================================================
#
mat.ztfutur1 <- matrix(0, nrow=Nsim, ncol=3)
mat.pred.ztfutur1 <- matrix(0, nrow=Nsim, ncol=3)
#
# == Contrastes ==
  contrast1 <- c(1, rep(0,mon.k-1))
  contrast2 <- c(rep(0,mon.k-1),1)
  contrast3 <- rep(1/mon.k,mon.k)
# == Définitions des matrices d'output ==
#
# == Simulation du processus ==
for(i in 1:Nsim) {
  ExoX <- VARMAsim(N, arlags = 1, phi=PhiX, sigma=SigmaX)
  xt.out <- ExoX$series
  VARX.out <- simulVARX(N, theta0, mon.p, matPhi, B0, mon.m, matB, xt.out, SigmaEps)
  xt <- xt.out[(n+2):(N-12),]
  zt <- VARX.out[(n+2):(N-12), ]
  zt.futur <- VARX.out[(N-12+1):N, ]
  zt.futur1 <- as.vector(VARX.out[(N-12+1):N, ] %*% contrast1)
  zt.futur2 <- as.vector(VARX.out[(N-12+1):N, ] %*% contrast2)
  zt.futur3 <- as.vector(VARX.out[(N-12+1):N, ] %*% contrast3)
  mat.ztfutur1[i,] <- c( zt.futur1[1], zt.futur2[1], zt.futur3[1])
#
# -- Estimation avec les bonnes matrices --
# 
     OutputY <- zt
     InputX <- xt
     mat.contraintes <- matR.VARXp3m2(k=ncol(OutputY), kx=ncol(InputX), ind=c(1,0),indx=c(1,0))
     model.out       <- VARX(OutputY, p=3, InputX, m=2, fixed=mat.contraintes, include.mean=T, output=F)
     mon.Sigmahat    <- model.out$Sigma
     ma.matPhistar   <- matPhistar(model.out$Phi, c(1,0,1))
     ma.matBJstar    <- matBJstar(model.out$beta, c(1,0,1), p=3)
     mon.Ystar       <- Ystar(zt, 3)
     mon.Ystarn      <- mon.Ystar[n,]
     mon.Xstar       <- Xstar(xt.out,s=2)
     mon.Xstar.futur <- mon.Xstar[(N-12+1):N,]
     mon.Vhat        <- Vhat.VARX(zt, xt, n, mon.Sigmahat)
     predict3.out    <- predict3.subsetVARX(mon.Ystarn, 3,2, mon.Xstar.futur, ma.matPhistar, ma.matBJstar, maxl)
     predictIP.out   <- predictIP.subsetVARX(3, 2, ma.matPhistar, mon.Sigmahat)
     predictIPcorr.out <- predictIPcorr.subsetVARX(3, 2, mon.Ystarn, mon.Xstar.futur, ma.matPhistar, ma.matBJstar, mon.Vhat, mon.Sigmahat)
     mat.pred.ztfutur1[i, ] <- c(predict3.out$output.punctual[1,1], predict3.out$output.punctual[1,2], predict3.out$output.punctual[1,3])
#
# ====== Absence de correction ======
# -- Horizon 1 --
# ---- constraste 1 ----
     outputIP.l1.lower[i,1] <- predict3.out$output.punctual[1,1] - q.nivconf * sqrt(predictIP.out[1,1])
     outputIP.l1.upper[i,1] <- predict3.out$output.punctual[1,1] + q.nivconf * sqrt(predictIP.out[1,1])
     outputnivIP.l1[i,1] <- (zt.futur1[1] > outputIP.l1.lower[i,1]) & (zt.futur1[1] < outputIP.l1.upper[i,1])
# ---- constraste 2 ----
     outputIP.l1.lower[i,2] <- predict3.out$output.punctual[1,2] - q.nivconf * sqrt(predictIP.out[1,2])
     outputIP.l1.upper[i,2] <- predict3.out$output.punctual[1,2] + q.nivconf * sqrt(predictIP.out[1,2])
     outputnivIP.l1[i,2] <- (zt.futur2[1] > outputIP.l1.lower[i,2]) & (zt.futur2[1] < outputIP.l1.upper[i,2])
# ---- constraste 3 ----
     outputIP.l1.lower[i,3] <- predict3.out$output.punctual[1,3] - q.nivconf * sqrt(predictIP.out[1,3])
     outputIP.l1.upper[i,3] <- predict3.out$output.punctual[1,3] + q.nivconf * sqrt(predictIP.out[1,3])
     outputnivIP.l1[i,3] <- (zt.futur3[1] > outputIP.l1.lower[i,3]) & (zt.futur3[1] < outputIP.l1.upper[i,3])
#
# -- Horizon 6 --
# ---- constraste 1 ----
     outputIP.l6.lower[i,1] <- predict3.out$output.punctual[6,1] - q.nivconf * sqrt(predictIP.out[6,1])
     outputIP.l6.upper[i,1] <- predict3.out$output.punctual[6,1] + q.nivconf * sqrt(predictIP.out[6,1])
     outputnivIP.l6[i,1] <- (zt.futur1[6] > outputIP.l6.lower[i,1]) & (zt.futur1[6] < outputIP.l6.upper[i,1])
# ---- constraste 2 ----
     outputIP.l6.lower[i,2] <- predict3.out$output.punctual[6,2] - q.nivconf * sqrt(predictIP.out[6,2])
     outputIP.l6.upper[i,2] <- predict3.out$output.punctual[6,2] + q.nivconf * sqrt(predictIP.out[6,2])
     outputnivIP.l6[i,2] <- (zt.futur2[6] > outputIP.l6.lower[i,2]) & (zt.futur2[6] < outputIP.l6.upper[i,2])
# ---- constraste 3 ----
     outputIP.l6.lower[i,3] <- predict3.out$output.punctual[6,3] - q.nivconf * sqrt(predictIP.out[6,3])
     outputIP.l6.upper[i,3] <- predict3.out$output.punctual[6,3] + q.nivconf * sqrt(predictIP.out[6,3])
     outputnivIP.l6[i,3] <- (zt.futur3[6] > outputIP.l6.lower[i,3]) & (zt.futur3[6] < outputIP.l6.upper[i,3])
#
# -- Horizon 12 --
# ---- constraste 1 ----
     outputIP.l12.lower[i,1] <- predict3.out$output.punctual[12,1] - q.nivconf * sqrt(predictIP.out[12,1])
     outputIP.l12.upper[i,1] <- predict3.out$output.punctual[12,1] + q.nivconf * sqrt(predictIP.out[12,1])
     outputnivIP.l12[i,1] <- (zt.futur1[12] > outputIP.l12.lower[i,1]) & (zt.futur1[12] < outputIP.l12.upper[i,1])
# ---- constraste 2 ----
     outputIP.l12.lower[i,2] <- predict3.out$output.punctual[12,2] - q.nivconf * sqrt(predictIP.out[12,2])
     outputIP.l12.upper[i,2] <- predict3.out$output.punctual[12,2] + q.nivconf * sqrt(predictIP.out[12,2])
     outputnivIP.l12[i,2] <- (zt.futur2[12] > outputIP.l12.lower[i,2]) & (zt.futur2[12] < outputIP.l12.upper[i,2])
# ---- constraste 3 ----
     outputIP.l12.lower[i,3] <- predict3.out$output.punctual[12,3] - q.nivconf * sqrt(predictIP.out[12,3])
     outputIP.l12.upper[i,3] <- predict3.out$output.punctual[12,3] + q.nivconf * sqrt(predictIP.out[12,3])
     outputnivIP.l12[i,3] <- (zt.futur3[12] > outputIP.l12.lower[i,3]) & (zt.futur3[12] < outputIP.l12.upper[i,3])
#
# ====== Avec correction ordre 1/n ======
# -- Horizon 1 --
# ---- constraste 1 ----
     outputIPcorr.l1.lower[i,1] <- predict3.out$output.punctual[1,1] - q.nivconf * sqrt(predictIPcorr.out[1,1])
     outputIPcorr.l1.upper[i,1] <- predict3.out$output.punctual[1,1] + q.nivconf * sqrt(predictIPcorr.out[1,1])
     outputnivIPcorr.l1[i,1] <- (zt.futur1[1] > outputIPcorr.l1.lower[i,1]) & (zt.futur1[1] < outputIPcorr.l1.upper[i,1])
# ---- constraste 2 ----
     outputIPcorr.l1.lower[i,2] <- predict3.out$output.punctual[1,2] - q.nivconf * sqrt(predictIPcorr.out[1,2])
     outputIPcorr.l1.upper[i,2] <- predict3.out$output.punctual[1,2] + q.nivconf * sqrt(predictIPcorr.out[1,2])
     outputnivIPcorr.l1[i,2] <- (zt.futur2[1] > outputIPcorr.l1.lower[i,2]) & (zt.futur2[1] < outputIPcorr.l1.upper[i,2])
# ---- constraste 3 ----
     outputIPcorr.l1.lower[i,3] <- predict3.out$output.punctual[1,3] - q.nivconf * sqrt(predictIPcorr.out[1,3])
     outputIPcorr.l1.upper[i,3] <- predict3.out$output.punctual[1,3] + q.nivconf * sqrt(predictIPcorr.out[1,3])
     outputnivIPcorr.l1[i,3] <- (zt.futur3[1] > outputIPcorr.l1.lower[i,3]) & (zt.futur3[1] < outputIPcorr.l1.upper[i,3])
#
# -- Horizon 6 --
# ---- constraste 1 ----
     outputIPcorr.l6.lower[i,1] <- predict3.out$output.punctual[6,1] - q.nivconf * sqrt(predictIPcorr.out[6,1])
     outputIPcorr.l6.upper[i,1] <- predict3.out$output.punctual[6,1] + q.nivconf * sqrt(predictIPcorr.out[6,1])
     outputnivIPcorr.l6[i,1] <- (zt.futur1[6] > outputIPcorr.l6.lower[i,1]) & (zt.futur1[6] < outputIPcorr.l6.upper[i,1])
# ---- constraste 2 ----
     outputIPcorr.l6.lower[i,2] <- predict3.out$output.punctual[6,2] - q.nivconf * sqrt(predictIPcorr.out[6,2])
     outputIPcorr.l6.upper[i,2] <- predict3.out$output.punctual[6,2] + q.nivconf * sqrt(predictIPcorr.out[6,2])
     outputnivIPcorr.l6[i,2] <- (zt.futur2[6] > outputIPcorr.l6.lower[i,2]) & (zt.futur2[6] < outputIPcorr.l6.upper[i,2])
# ---- constraste 3 ----
     outputIPcorr.l6.lower[i,3] <- predict3.out$output.punctual[6,3] - q.nivconf * sqrt(predictIPcorr.out[6,3])
     outputIPcorr.l6.upper[i,3] <- predict3.out$output.punctual[6,3] + q.nivconf * sqrt(predictIPcorr.out[6,3])
     outputnivIPcorr.l6[i,3] <- (zt.futur3[6] > outputIPcorr.l6.lower[i,3]) & (zt.futur3[6] < outputIPcorr.l6.upper[i,3])
#
# -- Horizon 12 --
# ---- constraste 1 ----
     outputIPcorr.l12.lower[i,1] <- predict3.out$output.punctual[12,1] - q.nivconf * sqrt(predictIPcorr.out[12,1])
     outputIPcorr.l12.upper[i,1] <- predict3.out$output.punctual[12,1] + q.nivconf * sqrt(predictIPcorr.out[12,1])
     outputnivIPcorr.l12[i,1] <- (zt.futur1[12] > outputIPcorr.l12.lower[i,1]) & (zt.futur1[12] < outputIPcorr.l12.upper[i,1])
# ---- constraste 2 ----
     outputIPcorr.l12.lower[i,2] <- predict3.out$output.punctual[12,2] - q.nivconf * sqrt(predictIPcorr.out[12,2])
     outputIPcorr.l12.upper[i,2] <- predict3.out$output.punctual[12,2] + q.nivconf * sqrt(predictIPcorr.out[12,2])
     outputnivIPcorr.l12[i,2] <- (zt.futur2[12] > outputIPcorr.l12.lower[i,2]) & (zt.futur2[12] < outputIPcorr.l12.upper[i,2])
# ---- constraste 3 ----
     outputIPcorr.l12.lower[i,3] <- predict3.out$output.punctual[12,3] - q.nivconf * sqrt(predictIPcorr.out[12,3])
     outputIPcorr.l12.upper[i,3] <- predict3.out$output.punctual[12,3] + q.nivconf * sqrt(predictIPcorr.out[12,3])
     outputnivIPcorr.l12[i,3] <- (zt.futur3[12] > outputIPcorr.l12.lower[i,3]) & (zt.futur3[12] < outputIPcorr.l12.upper[i,3])
# =======================================
}
 list(outputIP.l1.lower=outputIP.l1.lower, outputIP.l1.upper=outputIP.l1.upper, outputnivIP.l1=outputnivIP.l1,
      outputIP.l6.lower=outputIP.l6.lower, outputIP.l6.upper=outputIP.l6.upper, outputnivIP.l6=outputnivIP.l6,
      outputIP.l12.lower=outputIP.l12.lower, outputIP.l12.upper=outputIP.l12.upper, outputnivIP.l12=outputnivIP.l12,
      outputIPcorr.l1.lower=outputIPcorr.l1.lower, outputIPcorr.l1.upper=outputIPcorr.l1.upper, outputnivIPcorr.l1=outputnivIPcorr.l1,
      outputIPcorr.l6.lower=outputIPcorr.l6.lower, outputIPcorr.l6.upper=outputIPcorr.l6.upper, outputnivIPcorr.l6=outputnivIPcorr.l6,
      outputIPcorr.l12.lower=outputIPcorr.l12.lower, outputIPcorr.l12.upper=outputIPcorr.l12.upper, outputnivIPcorr.l12=outputnivIPcorr.l12,
      mat.ztfutur1 = mat.ztfutur1, mat.pred.ztfutur1=mat.pred.ztfutur1)
}



SimulationVARXp3m2.predictpunc <- function(n=128, p=3, m=2, 
                               PhiX, SigmaX, theta0=theta0, matPhi=matPhi, B0=B0, matB=matB, SigmaEps=SigmaEps, 
                               Nsim=10) {
# Pour la partie exogène on présume un VAR(1). Les paramètres sont PhiX et SigmaX.
# Le processus est un VARX(3,2), avec Phi2 = 0, et B1=0. Ainsi les indices du vrai modèle sont c(1,0,1,1,0,1)
# Le processus exogène est généré à chaque fois.
# 
# Pour les prévisions, on suppose que le bon modèle a été trouvé, ce qui va illustrer la théorie.
# 
set.seed(1)
N <- 2*n + 1 + 12
mon.k <- nrow(matPhi)
mon.kx <- nrow(PhiX)
maxl <- 12
mon.p <- p
mon.m <- m
# 
#
# ==========================================================
#
mat.ztfutur1 <- matrix(0, nrow=Nsim, ncol=3)
mat.ztfutur6 <- matrix(0, nrow=Nsim, ncol=3)
mat.ztfutur12 <- matrix(0, nrow=Nsim, ncol=3)
mat.pred.ztfutur1 <- matrix(0, nrow=Nsim, ncol=3)
mat.pred.ztfutur1.full <- matrix(0, nrow=Nsim, ncol=3)
mat.pred.ztfutur6 <- matrix(0, nrow=Nsim, ncol=3)
mat.pred.ztfutur6.full <- matrix(0, nrow=Nsim, ncol=3)
mat.pred.ztfutur12 <- matrix(0, nrow=Nsim, ncol=3)
mat.pred.ztfutur12.full <- matrix(0, nrow=Nsim, ncol=3)
#
# == Contrastes ==
  contrast1 <- c(1, rep(0,mon.k-1))
  contrast2 <- c(rep(0,mon.k-1),1)
  contrast3 <- rep(1/mon.k,mon.k)
# == Définitions des matrices d'output ==
#
# == Simulation du processus ==
for(i in 1:Nsim) {
  ExoX <- VARMAsim(N, arlags = 1, phi=PhiX, sigma=SigmaX)
  xt.out <- ExoX$series
  VARX.out <- simulVARX(N, theta0, mon.p, matPhi, B0, mon.m, matB, xt.out, SigmaEps)
  xt <- xt.out[(n+2):(N-12),]
  zt <- VARX.out[(n+2):(N-12), ]
  zt.futur <- VARX.out[(N-12+1):N, ]
  zt.futur1 <- as.vector(VARX.out[(N-12+1):N, ] %*% contrast1)
  zt.futur2 <- as.vector(VARX.out[(N-12+1):N, ] %*% contrast2)
  zt.futur3 <- as.vector(VARX.out[(N-12+1):N, ] %*% contrast3)
  mat.ztfutur1[i,] <- c( zt.futur1[1], zt.futur2[1], zt.futur3[1])
  mat.ztfutur6[i,] <- c( zt.futur1[6], zt.futur2[6], zt.futur3[6])
  mat.ztfutur12[i,] <- c( zt.futur1[12], zt.futur2[12], zt.futur3[12])

#
# -- Estimation avec les bonnes matrices --
# 
     OutputY <- zt
     InputX <- xt
     mat.contraintes <- matR.VARXp3m2(k=ncol(OutputY), kx=ncol(InputX), ind=c(1,0),indx=c(1,0))
     model.out       <- VARX(OutputY, p=3, InputX, m=2, fixed=mat.contraintes, include.mean=T, output=F)
     mon.Sigmahat    <- model.out$Sigma
     ma.matPhistar   <- matPhistar(model.out$Phi, c(1,0,1))
     ma.matBJstar    <- matBJstar(model.out$beta, c(1,0,1), p=3)
     mon.Ystar       <- Ystar(zt, 3)
     mon.Ystarn      <- mon.Ystar[n,]
     mon.Xstar       <- Xstar(xt.out,s=2)
     mon.Xstar.futur <- mon.Xstar[(N-12+1):N,]
     mon.Vhat        <- Vhat.VARX(zt, xt, n, mon.Sigmahat)
     predict3.out    <- predict3.subsetVARX(mon.Ystarn, 3,2, mon.Xstar.futur, ma.matPhistar, ma.matBJstar, maxl)
     mat.pred.ztfutur1[i, ] <- c(predict3.out$output.punctual[1,1], predict3.out$output.punctual[1,2], predict3.out$output.punctual[1,3])
     mat.pred.ztfutur6[i, ] <- c(predict3.out$output.punctual[6,1], predict3.out$output.punctual[6,2], predict3.out$output.punctual[6,3])
     mat.pred.ztfutur12[i, ] <- c(predict3.out$output.punctual[12,1], predict3.out$output.punctual[12,2], predict3.out$output.punctual[12,3])
#
# -- Estimation du modèle complet --
# 
     OutputY <- zt
     InputX <- xt
     model.full.out       <- VARX(OutputY, p=3, InputX, m=2, include.mean=T, output=F)
     mon.Sigmahat.full    <- model.full.out$Sigma
     ma.matPhistar.full   <- matPhistar(model.full.out$Phi, c(1,1,1))
     ma.matBJstar.full    <- matBJstar(model.full.out$beta, c(1,1,1), p=3)
     mon.Ystar       <- Ystar(zt, 3)
     mon.Ystarn      <- mon.Ystar[n,]
     mon.Xstar       <- Xstar(xt.out,s=2)
     mon.Xstar.futur <- mon.Xstar[(N-12+1):N,]
     mon.Vhat.full        <- Vhat.VARX(zt, xt, n, mon.Sigmahat.full, ind=c(1,1,1), indx=c(1,1,1))
     predict3.full.out    <- predict3.subsetVARX(mon.Ystarn, 3,2, mon.Xstar.futur, ma.matPhistar.full, ma.matBJstar.full, maxl)
     mat.pred.ztfutur1.full[i, ] <- c(predict3.full.out$output.punctual[1,1], predict3.full.out$output.punctual[1,2], predict3.full.out$output.punctual[1,3])
     mat.pred.ztfutur6.full[i, ] <- c(predict3.full.out$output.punctual[6,1], predict3.full.out$output.punctual[6,2], predict3.full.out$output.punctual[6,3])
     mat.pred.ztfutur12.full[i, ] <- c(predict3.full.out$output.punctual[12,1], predict3.full.out$output.punctual[12,2], predict3.full.out$output.punctual[12,3])
}
 list(mat.ztfutur1 = mat.ztfutur1, mat.ztfutur6 = mat.ztfutur6, mat.ztfutur12 = mat.ztfutur12,
      mat.pred.ztfutur1=mat.pred.ztfutur1, 
      mat.pred.ztfutur1.full=mat.pred.ztfutur1.full,
       mat.pred.ztfutur6=mat.pred.ztfutur6,
       mat.pred.ztfutur6.full = mat.pred.ztfutur6.full, 
       mat.pred.ztfutur12=mat.pred.ztfutur12,
       mat.pred.ztfutur12.full = mat.pred.ztfutur12.full
      )
}



# === Expériences numériques ===
# -------------------------------------------------
# --- Modèle 1: Paramètres de la partie exogène ---
# -------------------------------------------------
# --------------
library(MTS)
library(mvtnorm)
# --------------
set.seed(1)
mon.k    <- 2
mon.kx   <- 2
mon.p <- 3
mon.m <- 2
PhiX     <- matrix( c(0.2, -0.6, 0.3, 1), 2,2)
SigmaX   <- matrix( c(4.0, 0.8, 0.8, 1), 2, 2)
Phi1     <- matrix( c(0.5, -0.4, 0.15, 0.5), nrow=2)
Phi2     <- matrix( 0, nrow=2,ncol=2)
Phi3     <- matrix( c(0.2, -0.2, 0.1, 0.25), nrow=2)
matPhi   <- cbind(Phi1, Phi2, Phi3)
SigmaEps <- SigmaSpherique(2, 0.8)
B0       <- matrix( c(0.2, 0.6, -0.3, 1), 2,2)
B1       <- matrix( 0, nrow=2,ncol=2)
B2       <- matrix( c(2.4, -0.9, 0.6, 3.1), 2,2)
matB     <- cbind(B1, B2)
theta0   <- c(0.0, 0.0)

n <- 128
N <- 2*n + 1 + 12
ExoX <- VARMAsim(N, arlags = 1, phi=PhiX, sigma=SigmaX)
xt.out <- ExoX$series
VARX.out <- simulVARX(N, theta0, mon.p, matPhi, B0, mon.m, matB, xt.out, SigmaEps)
xt <- xt.out[(n+2):(N-12),]
zt <- VARX.out[(n+2):(N-12), ]
zt.futur <- VARX.out[(N-12+1):N, ]
mat.contraintes <- matR.VARXp3m2(k=mon.k, kx=mon.kx, ind=c(1,0),indx=c(1,0) )
model.out <- VARX(zt, p= mon.p, xt, m=mon.m, include.mean=F, fixed=mat.contraintes[-1,], output=F)

# --- Estimateur de SigmaEps ---
mon.Sigmahat <- model.out$Sigma

ma.matPhistar <- matPhistar(model.out$Phi, c(1,0,1))
ma.matBJstar <- matBJstar(model.out$beta, c(1,0,1), 3)
mon.Ystar <- Ystar(zt, 3)
mon.Ystarn <- mon.Ystar[n,]
mon.Xstar <- Xstar(xt.out,s=2)
mon.Xstar.futur <- mon.Xstar[(N-12+1):N,]


# --- Estimateur de variance ---
matU.I <- matC.I( c(1,0,1), 3, 2) %*% t(mon.Ystar[1:(n-1),])
matV.J <- matD.J( c(1,0,1), s=2, 2) %*% t(Xstar(xt[2:n,],s=2))
matZ <- rbind( matU.I, matV.J)
Vhat <- kronecker( solve( matZ %*% t(matZ) ), mon.Sigmahat)
sqrt( diag(Vhat) )
model.out$se.coef
Vhat.VARX(zt,xt,n,mon.Sigmahat)


# ---Calcul de la correction ---
ma.matPhiI <- matPhiI(model.out$Phi, c(1,0,1))
ma.matBJ <- matBJ(model.out$beta, c(1,0,1)) 
ma.lambdaIJ <- c( as.vector(ma.matPhiI), as.vector(ma.matBJ) )
mon.Wnl(mon.Ystarn, mon.Xstar.futur, 1, 2)
bloc.Al(1,ma.matPhistar)
bloc.B1l(1,ma.matPhistar,ma.matBJstar)
matHIJ.l1(ma.matPhistar) 
matHIJ.l6(ma.matPhistar, ma.matBJstar) 
matHIJ.l(1, ma.matPhistar, ma.matBJstar) 

# -- Calcul des prévisions ---
predict.subsetVARX <- function(Ynstar, p, k, Xstar.futur, Phistar, Bstar, maxl=12)
predict.subsetVARX(mon.Ystarn, 3,2, mon.Xstar.futur, ma.matPhistar, ma.matBJstar, 12)
predict3.subsetVARX(mon.Ystarn, 3,2, mon.Xstar.futur, ma.matPhistar, ma.matBJstar, 12)
predictIP.subsetVARX(3, 2, ma.matPhistar, mon.Sigmahat)

# -- Comparaisons des prévisions ponctuelles ---
Nsim <- 5000

# ===============
# === n = 120 ===
# ===============
simul120punc.out <- SimulationVARXp3m2.predictpunc(n=120, p=3, m=2, PhiX, SigmaX,
                               theta0=theta0, matPhi=matPhi, B0=B0, matB=matB, SigmaEps=SigmaEps, Nsim)
biais120.subset1 <- apply(simul120punc.out$mat.ztfutur1 - simul120punc.out$mat.pred.ztfutur1, 2, mean)
biais120.full1 <- apply(simul120punc.out$mat.ztfutur1 - simul120punc.out$mat.pred.ztfutur1.full, 2, mean)
EQM120.subset1 <- apply(simul120punc.out$mat.ztfutur1 - simul120punc.out$mat.pred.ztfutur1, 2, var)
EQM120.full1 <- apply(simul120punc.out$mat.ztfutur1 - simul120punc.out$mat.pred.ztfutur1.full, 2, var)

biais120.subset6 <- apply(simul120punc.out$mat.ztfutur6 - simul120punc.out$mat.pred.ztfutur6, 2, mean)
biais120.full6 <- apply(simul120punc.out$mat.ztfutur6 - simul120punc.out$mat.pred.ztfutur6.full, 2, mean)
EQM120.subset6 <- apply(simul120punc.out$mat.ztfutur6 - simul120punc.out$mat.pred.ztfutur6, 2, var)
EQM120.full6 <- apply(simul120punc.out$mat.ztfutur6 - simul120punc.out$mat.pred.ztfutur6.full, 2, var)

biais120.subset12 <- apply(simul120punc.out$mat.ztfutur12 - simul120punc.out$mat.pred.ztfutur12, 2, mean)
biais120.full12 <- apply(simul120punc.out$mat.ztfutur12 - simul120punc.out$mat.pred.ztfutur12.full, 2, mean)
EQM120.subset12 <- apply(simul120punc.out$mat.ztfutur12 - simul120punc.out$mat.pred.ztfutur12, 2, var)
EQM120.full12 <- apply(simul120punc.out$mat.ztfutur12 - simul120punc.out$mat.pred.ztfutur12.full, 2, var)


# ===============
# === n = 240 ===
# ===============
simul240punc.out <- SimulationVARXp3m2.predictpunc(n=240, p=3, m=2, PhiX, SigmaX,
                               theta0=theta0, matPhi=matPhi, B0=B0, matB=matB, SigmaEps=SigmaEps, Nsim)
biais240.subset1 <- apply(simul240punc.out$mat.ztfutur1 - simul240punc.out$mat.pred.ztfutur1, 2, mean)
biais240.full1 <- apply(simul240punc.out$mat.ztfutur1 - simul240punc.out$mat.pred.ztfutur1.full, 2, mean)
EQM240.subset1 <- apply(simul240punc.out$mat.ztfutur1 - simul240punc.out$mat.pred.ztfutur1, 2, var)
EQM240.full1 <- apply(simul240punc.out$mat.ztfutur1 - simul240punc.out$mat.pred.ztfutur1.full, 2, var)

biais240.subset6 <- apply(simul240punc.out$mat.ztfutur6 - simul240punc.out$mat.pred.ztfutur6, 2, mean)
biais240.full6 <- apply(simul240punc.out$mat.ztfutur6 - simul240punc.out$mat.pred.ztfutur6.full, 2, mean)
EQM240.subset6 <- apply(simul240punc.out$mat.ztfutur6 - simul240punc.out$mat.pred.ztfutur6, 2, var)
EQM240.full6 <- apply(simul240punc.out$mat.ztfutur6 - simul240punc.out$mat.pred.ztfutur6.full, 2, var)

biais240.subset12 <- apply(simul240punc.out$mat.ztfutur12 - simul240punc.out$mat.pred.ztfutur12, 2, mean)
biais240.full12 <- apply(simul240punc.out$mat.ztfutur12 - simul240punc.out$mat.pred.ztfutur12.full, 2, mean)
EQM240.subset12 <- apply(simul240punc.out$mat.ztfutur12 - simul240punc.out$mat.pred.ztfutur12, 2, var)
EQM240.full12 <- apply(simul240punc.out$mat.ztfutur12 - simul240punc.out$mat.pred.ztfutur12.full, 2, var)



# ===============
# === n = 360 ===
# ===============
simul360punc.out <- SimulationVARXp3m2.predictpunc(n=360, p=3, m=2, PhiX, SigmaX,
                               theta0=theta0, matPhi=matPhi, B0=B0, matB=matB, SigmaEps=SigmaEps, Nsim)
biais360.subset1 <- apply(simul360punc.out$mat.ztfutur1 - simul360punc.out$mat.pred.ztfutur1, 2, mean)
biais360.full1 <- apply(simul360punc.out$mat.ztfutur1 - simul360punc.out$mat.pred.ztfutur1.full, 2, mean)
EQM360.subset1 <- apply(simul360punc.out$mat.ztfutur1 - simul360punc.out$mat.pred.ztfutur1, 2, var)
EQM360.full1 <- apply(simul360punc.out$mat.ztfutur1 - simul360punc.out$mat.pred.ztfutur1.full, 2, var)

biais360.subset6 <- apply(simul360punc.out$mat.ztfutur6 - simul360punc.out$mat.pred.ztfutur6, 2, mean)
biais360.full6 <- apply(simul360punc.out$mat.ztfutur6 - simul360punc.out$mat.pred.ztfutur6.full, 2, mean)
EQM360.subset6 <- apply(simul360punc.out$mat.ztfutur6 - simul360punc.out$mat.pred.ztfutur6, 2, var)
EQM360.full6 <- apply(simul360punc.out$mat.ztfutur6 - simul360punc.out$mat.pred.ztfutur6.full, 2, var)

biais360.subset12 <- apply(simul360punc.out$mat.ztfutur12 - simul360punc.out$mat.pred.ztfutur12, 2, mean)
biais360.full12 <- apply(simul360punc.out$mat.ztfutur12 - simul360punc.out$mat.pred.ztfutur12.full, 2, mean)
EQM360.subset12 <- apply(simul360punc.out$mat.ztfutur12 - simul360punc.out$mat.pred.ztfutur12, 2, var)
EQM360.full12 <- apply(simul360punc.out$mat.ztfutur12 - simul360punc.out$mat.pred.ztfutur12.full, 2, var)



# ===============
# === n = 480 ===
# ===============
simul480punc.out <- SimulationVARXp3m2.predictpunc(n=480, p=3, m=2, PhiX, SigmaX,
                               theta0=theta0, matPhi=matPhi, B0=B0, matB=matB, SigmaEps=SigmaEps, Nsim)
biais480.subset1 <- apply(simul480punc.out$mat.ztfutur1 - simul480punc.out$mat.pred.ztfutur1, 2, mean)
biais480.full1 <- apply(simul480punc.out$mat.ztfutur1 - simul480punc.out$mat.pred.ztfutur1.full, 2, mean)
EQM480.subset1 <- apply(simul480punc.out$mat.ztfutur1 - simul480punc.out$mat.pred.ztfutur1, 2, var)
EQM480.full1 <- apply(simul480punc.out$mat.ztfutur1 - simul480punc.out$mat.pred.ztfutur1.full, 2, var)

biais480.subset6 <- apply(simul480punc.out$mat.ztfutur6 - simul480punc.out$mat.pred.ztfutur6, 2, mean)
biais480.full6 <- apply(simul480punc.out$mat.ztfutur6 - simul480punc.out$mat.pred.ztfutur6.full, 2, mean)
EQM480.subset6 <- apply(simul480punc.out$mat.ztfutur6 - simul480punc.out$mat.pred.ztfutur6, 2, var)
EQM480.full6 <- apply(simul480punc.out$mat.ztfutur6 - simul480punc.out$mat.pred.ztfutur6.full, 2, var)

biais480.subset12 <- apply(simul480punc.out$mat.ztfutur12 - simul480punc.out$mat.pred.ztfutur12, 2, mean)
biais480.full12 <- apply(simul480punc.out$mat.ztfutur12 - simul480punc.out$mat.pred.ztfutur12.full, 2, mean)
EQM480.subset12 <- apply(simul480punc.out$mat.ztfutur12 - simul480punc.out$mat.pred.ztfutur12, 2, var)
EQM480.full12 <- apply(simul480punc.out$mat.ztfutur12 - simul480punc.out$mat.pred.ztfutur12.full, 2, var)


# -- Présentation des résultats --
# --- n = 120 ---
round(
      c(biais120.subset1[1], biais120.subset6[1], biais120.subset12[1], 
      biais120.subset1[2], biais120.subset6[2], biais120.subset12[2],
      biais120.subset1[3], biais120.subset6[3], biais120.subset12[3]),
     4)
round(
      c(biais120.full1[1], biais120.full6[1], biais120.full12[1], 
      biais120.full1[2], biais120.full6[2], biais120.full12[2],
      biais120.full1[3], biais120.full6[3], biais120.full12[3]),
     4)
round(
     c(EQM120.subset1[1], EQM120.subset6[1], EQM120.subset12[1], 
     EQM120.subset1[2], EQM120.subset6[2], EQM120.subset12[2],
     EQM120.subset1[3], EQM120.subset6[3], EQM120.subset12[3]), 
     4)
round(
    c(EQM120.full1[1], EQM120.full6[1], EQM120.full12[1], 
      EQM120.full1[2], EQM120.full6[2], EQM120.full12[2],
      EQM120.full1[3], EQM120.full6[3], EQM120.full12[3]),
     4)

# --- n = 240 ---
round(
c(biais240.subset1[1], biais240.subset6[1], biais240.subset12[1], 
  biais240.subset1[2], biais240.subset6[2], biais240.subset12[2],
  biais240.subset1[3], biais240.subset6[3], biais240.subset12[3]),
4)
round(
c(biais240.full1[1], biais240.full6[1], biais240.full12[1], 
  biais240.full1[2], biais240.full6[2], biais240.full12[2],
  biais240.full1[3], biais240.full6[3], biais240.full12[3]),
4)
round(
c(EQM240.subset1[1], EQM240.subset6[1], EQM240.subset12[1], 
  EQM240.subset1[2], EQM240.subset6[2], EQM240.subset12[2],
  EQM240.subset1[3], EQM240.subset6[3], EQM240.subset12[3]),
4)
round(
c(EQM240.full1[1], EQM240.full6[1], EQM240.full12[1], 
  EQM240.full1[2], EQM240.full6[2], EQM240.full12[2],
  EQM240.full1[3], EQM240.full6[3], EQM240.full12[3]),
4)

# --- n = 360 ---
round(
c(biais360.subset1[1], biais360.subset6[1], biais360.subset12[1], 
  biais360.subset1[2], biais360.subset6[2], biais360.subset12[2],
  biais360.subset1[3], biais360.subset6[3], biais360.subset12[3]),
4)
round(
  c(biais360.full1[1], biais360.full6[1], biais360.full12[1], 
  biais360.full1[2], biais360.full6[2], biais360.full12[2],
  biais360.full1[3], biais360.full6[3], biais360.full12[3]),
4)
round(
  c(EQM360.subset1[1], EQM360.subset6[1], EQM360.subset12[1], 
  EQM360.subset1[2], EQM360.subset6[2], EQM360.subset12[2],
  EQM360.subset1[3], EQM360.subset6[3], EQM360.subset12[3]),
4)
round(
  c(EQM360.full1[1], EQM360.full6[1], EQM360.full12[1], 
  EQM360.full1[2], EQM360.full6[2], EQM360.full12[2],
  EQM360.full1[3], EQM360.full6[3], EQM360.full12[3]),
4)

# --- n = 480 ---
round(
  c(biais480.subset1[1], biais480.subset6[1], biais480.subset12[1], 
  biais480.subset1[2], biais480.subset6[2], biais480.subset12[2],
  biais480.subset1[3], biais480.subset6[3], biais480.subset12[3]),
4)
round(
  c(biais480.full1[1], biais480.full6[1], biais480.full12[1], 
  biais480.full1[2], biais480.full6[2], biais480.full12[2],
  biais480.full1[3], biais480.full6[3], biais480.full12[3]),
4)
round(
  c(EQM480.subset1[1], EQM480.subset6[1], EQM480.subset12[1], 
  EQM480.subset1[2], EQM480.subset6[2], EQM480.subset12[2],
  EQM480.subset1[3], EQM480.subset6[3], EQM480.subset12[3]),
4)
round(
  c(EQM480.full1[1], EQM480.full6[1], EQM480.full12[1], 
  EQM480.full1[2], EQM480.full6[2], EQM480.full12[2],
  EQM480.full1[3], EQM480.full6[3], EQM480.full12[3]),
4)



> # -- Présentation des résultats --
[1] 0.0126 0.0180 0.0094 0.0202 0.0135 0.0012 0.0164 0.0157 0.0053
[1] 0.0106 0.0158 0.0069 0.0199 0.0149 0.0005 0.0153 0.0153 0.0037
[1] 1.0724 2.0205 2.1369 1.0761 1.8494 2.4923 0.9652 1.2442 1.4639
[1] 1.1100 2.0761 2.1985 1.1218 1.8855 2.5752 1.0036 1.2732 1.5170
> 
[1] -0.0097 -0.0076  0.0022 -0.0134 -0.0254  0.0091 -0.0116 -0.0165  0.0057
[1] -0.0088 -0.0065  0.0041 -0.0114 -0.0248  0.0110 -0.0101 -0.0157  0.0076
[1] 1.0417 1.9174 2.0183 1.0443 1.6504 2.2701 0.9421 1.1620 1.3898
[1] 1.0674 1.9393 2.0472 1.0669 1.6698 2.2931 0.9641 1.1777 1.4082
> 
[1] 0.0001 0.0088 0.0405 0.0002 0.0058 0.0203 0.0001 0.0073 0.0304
[1] -0.0010  0.0077  0.0419 -0.0002  0.0063  0.0174 -0.0006  0.0070  0.0297
[1] 1.0135 1.8934 1.8941 1.0369 1.6498 2.1493 0.9215 1.1563 1.2909
[1] 1.0301 1.9043 1.9087 1.0522 1.6671 2.1599 0.9354 1.1650 1.2993

[1]  0.0052  0.0076  0.0345 -0.0044  0.0103 -0.0209  0.0004  0.0089  0.0068
[1]  0.0042  0.0058  0.0329 -0.0059  0.0101 -0.0219 -0.0008  0.0080  0.0055
[1] 1.0013 1.8475 1.9791 1.0081 1.6544 2.1482 0.9048 1.1270 1.3240
[1] 1.0093 1.8607 1.9860 1.0128 1.6626 2.1592 0.9102 1.1333 1.3287




# --- Simulation ---
Nsim <- 5000

# ===============
# === n = 120 ===
# ===============
simul120niv90.out <- SimulationVARXp3m2.predict(n=120, p=3, m=2, PhiX, SigmaX,
               theta0=theta0, matPhi=matPhi, B0=B0, matB=matB, SigmaEps=SigmaEps, alpha=0.10, Nsim)

simul120niv95.out <- SimulationVARXp3m2.predict(n=120, p=3, m=2, PhiX, SigmaX,
               theta0=theta0, matPhi=matPhi, B0=B0, matB=matB, SigmaEps=SigmaEps, alpha=0.05, Nsim)

# ===============
# === n = 240 ===
# ===============

simul240niv90.out <- SimulationVARXp3m2.predict(n=240, p=3, m=2, PhiX, SigmaX,
               theta0=theta0, matPhi=matPhi, B0=B0, matB=matB, SigmaEps=SigmaEps, alpha=0.10, Nsim)

simul240niv95.out <- SimulationVARXp3m2.predict(n=240, p=3, m=2, PhiX, SigmaX,
               theta0=theta0, matPhi=matPhi, B0=B0, matB=matB, SigmaEps=SigmaEps, alpha=0.05, Nsim)

# ===============
# === n = 360 ===
# ===============

simul360niv90.out <- SimulationVARXp3m2.predict(n=360, p=3, m=2, PhiX, SigmaX,
               theta0=theta0, matPhi=matPhi, B0=B0, matB=matB, SigmaEps=SigmaEps, alpha=0.10, Nsim)

simul360niv95.out <- SimulationVARXp3m2.predict(n=360, p=3, m=2, PhiX, SigmaX,
               theta0=theta0, matPhi=matPhi, B0=B0, matB=matB, SigmaEps=SigmaEps, alpha=0.05, Nsim)

# ===============
# === n = 480 ===
# ===============

simul480niv90.out <- SimulationVARXp3m2.predict(n=480, p=3, m=2, PhiX, SigmaX,
               theta0=theta0, matPhi=matPhi, B0=B0, matB=matB, SigmaEps=SigmaEps, alpha=0.10, Nsim)

simul480niv95.out <- SimulationVARXp3m2.predict(n=480, p=3, m=2, PhiX, SigmaX,
               theta0=theta0, matPhi=matPhi, B0=B0, matB=matB, SigmaEps=SigmaEps, alpha=0.05, Nsim)


=== Présentation des résultats ===
# ===============
# === n = 120 ===
# ===============

n120.l1.niv90      <- apply(simul120niv90.out$outputnivIP.l1, 2, mean)
n120.l1.niv90.corr <- apply(simul120niv90.out$outputnivIPcorr.l1, 2, mean)

n120.l6.niv90      <- apply(simul120niv90.out$outputnivIP.l6, 2, mean)
n120.l6.niv90.corr <- apply(simul120niv90.out$outputnivIPcorr.l6, 2, mean)

n120.l12.niv90      <- apply(simul120niv90.out$outputnivIP.l12, 2, mean)
n120.l12.niv90.corr <- apply(simul120niv90.out$outputnivIPcorr.l12, 2, mean)

n120.l1.niv95      <- apply(simul120niv95.out$outputnivIP.l1, 2, mean)
n120.l1.niv95.corr <- apply(simul120niv95.out$outputnivIPcorr.l1, 2, mean)

n120.l6.niv95      <- apply(simul120niv95.out$outputnivIP.l6, 2, mean)
n120.l6.niv95.corr <- apply(simul120niv95.out$outputnivIPcorr.l6, 2, mean)

n120.l12.niv95      <- apply(simul120niv95.out$outputnivIP.l12, 2, mean)
n120.l12.niv95.corr <- apply(simul120niv95.out$outputnivIPcorr.l12, 2, mean)

# ===============
# === n = 240 ===
# ===============

n240.l1.niv90      <- apply(simul240niv90.out$outputnivIP.l1, 2, mean)
n240.l1.niv90.corr <- apply(simul240niv90.out$outputnivIPcorr.l1, 2, mean)

n240.l6.niv90      <- apply(simul240niv90.out$outputnivIP.l6, 2, mean)
n240.l6.niv90.corr <- apply(simul240niv90.out$outputnivIPcorr.l6, 2, mean)

n240.l12.niv90      <- apply(simul240niv90.out$outputnivIP.l12, 2, mean)
n240.l12.niv90.corr <- apply(simul240niv90.out$outputnivIPcorr.l12, 2, mean)

n240.l1.niv95      <- apply(simul240niv95.out$outputnivIP.l1, 2, mean)
n240.l1.niv95.corr <- apply(simul240niv95.out$outputnivIPcorr.l1, 2, mean)

n240.l6.niv95      <- apply(simul240niv95.out$outputnivIP.l6, 2, mean)
n240.l6.niv95.corr <- apply(simul240niv95.out$outputnivIPcorr.l6, 2, mean)

n240.l12.niv95      <- apply(simul240niv95.out$outputnivIP.l12, 2, mean)
n240.l12.niv95.corr <- apply(simul240niv95.out$outputnivIPcorr.l12, 2, mean)



# ===============
# === n = 360 ===
# ===============


n360.l1.niv90      <- apply(simul360niv90.out$outputnivIP.l1, 2, mean)
n360.l1.niv90.corr <- apply(simul360niv90.out$outputnivIPcorr.l1, 2, mean)

n360.l6.niv90      <- apply(simul360niv90.out$outputnivIP.l6, 2, mean)
n360.l6.niv90.corr <- apply(simul360niv90.out$outputnivIPcorr.l6, 2, mean)

n360.l12.niv90      <- apply(simul360niv90.out$outputnivIP.l12, 2, mean)
n360.l12.niv90.corr <- apply(simul360niv90.out$outputnivIPcorr.l12, 2, mean)

n360.l1.niv95      <- apply(simul360niv95.out$outputnivIP.l1, 2, mean)
n360.l1.niv95.corr <- apply(simul360niv95.out$outputnivIPcorr.l1, 2, mean)

n360.l6.niv95      <- apply(simul360niv95.out$outputnivIP.l6, 2, mean)
n360.l6.niv95.corr <- apply(simul360niv95.out$outputnivIPcorr.l6, 2, mean)

n360.l12.niv95      <- apply(simul360niv95.out$outputnivIP.l12, 2, mean)
n360.l12.niv95.corr <- apply(simul360niv95.out$outputnivIPcorr.l12, 2, mean)


# ===============
# === n = 480 ===
# ===============


n480.l1.niv90      <- apply(simul480niv90.out$outputnivIP.l1, 2, mean)
n480.l1.niv90.corr <- apply(simul480niv90.out$outputnivIPcorr.l1, 2, mean)

n480.l6.niv90      <- apply(simul480niv90.out$outputnivIP.l6, 2, mean)
n480.l6.niv90.corr <- apply(simul480niv90.out$outputnivIPcorr.l6, 2, mean)

n480.l12.niv90      <- apply(simul480niv90.out$outputnivIP.l12, 2, mean)
n480.l12.niv90.corr <- apply(simul480niv90.out$outputnivIPcorr.l12, 2, mean)

n480.l1.niv95      <- apply(simul480niv95.out$outputnivIP.l1, 2, mean)
n480.l1.niv95.corr <- apply(simul480niv95.out$outputnivIPcorr.l1, 2, mean)

n480.l6.niv95      <- apply(simul480niv95.out$outputnivIP.l6, 2, mean)
n480.l6.niv95.corr <- apply(simul480niv95.out$outputnivIPcorr.l6, 2, mean)

n480.l12.niv95      <- apply(simul480niv95.out$outputnivIP.l12, 2, mean)
n480.l12.niv95.corr <- apply(simul480niv95.out$outputnivIPcorr.l12, 2, mean)

# =================== Tableaux =================
# --- n = 120 ---
c(n120.l1.niv90[1], n120.l6.niv90[1], n120.l12.niv90[1], 
  n120.l1.niv90[2], n120.l6.niv90[2], n120.l12.niv90[2],
  n120.l1.niv90[3], n120.l6.niv90[3], n120.l12.niv90[3])
c(n120.l1.niv90.corr[1], n120.l6.niv90.corr[1], n120.l12.niv90.corr[1], 
  n120.l1.niv90.corr[2], n120.l6.niv90.corr[2], n120.l12.niv90.corr[2],
  n120.l1.niv90.corr[3], n120.l6.niv90.corr[3], n120.l12.niv90.corr[3])

c(n120.l1.niv95[1], n120.l6.niv95[1], n120.l12.niv95[1], 
  n120.l1.niv95[2], n120.l6.niv95[2], n120.l12.niv95[2],
  n120.l1.niv95[3], n120.l6.niv95[3], n120.l12.niv95[3])
c(n120.l1.niv95.corr[1], n120.l6.niv95.corr[1], n120.l12.niv95.corr[1], 
  n120.l1.niv95.corr[2], n120.l6.niv95.corr[2], n120.l12.niv95.corr[2],
  n120.l1.niv95.corr[3], n120.l6.niv95.corr[3], n120.l12.niv95.corr[3])

# --- n = 240 ---
c(n240.l1.niv90[1], n240.l6.niv90[1], n240.l12.niv90[1], 
  n240.l1.niv90[2], n240.l6.niv90[2], n240.l12.niv90[2],
  n240.l1.niv90[3], n240.l6.niv90[3], n240.l12.niv90[3])
c(n240.l1.niv90.corr[1], n240.l6.niv90.corr[1], n240.l12.niv90.corr[1], 
  n240.l1.niv90.corr[2], n240.l6.niv90.corr[2], n240.l12.niv90.corr[2],
  n240.l1.niv90.corr[3], n240.l6.niv90.corr[3], n240.l12.niv90.corr[3])

c(n240.l1.niv95[1], n240.l6.niv95[1], n240.l12.niv95[1], 
  n240.l1.niv95[2], n240.l6.niv95[2], n240.l12.niv95[2],
  n240.l1.niv95[3], n240.l6.niv95[3], n240.l12.niv95[3])
c(n240.l1.niv95.corr[1], n240.l6.niv95.corr[1], n240.l12.niv95.corr[1], 
  n240.l1.niv95.corr[2], n240.l6.niv95.corr[2], n240.l12.niv95.corr[2],
  n240.l1.niv95.corr[3], n240.l6.niv95.corr[3], n240.l12.niv95.corr[3])


# --- n = 360 ---
c(n360.l1.niv90[1], n360.l6.niv90[1], n360.l12.niv90[1], 
  n360.l1.niv90[2], n360.l6.niv90[2], n360.l12.niv90[2],
  n360.l1.niv90[3], n360.l6.niv90[3], n360.l12.niv90[3])
c(n360.l1.niv90.corr[1], n360.l6.niv90.corr[1], n360.l12.niv90.corr[1], 
  n360.l1.niv90.corr[2], n360.l6.niv90.corr[2], n360.l12.niv90.corr[2],
  n360.l1.niv90.corr[3], n360.l6.niv90.corr[3], n360.l12.niv90.corr[3])

c(n360.l1.niv95[1], n360.l6.niv95[1], n360.l12.niv95[1], 
  n360.l1.niv95[2], n360.l6.niv95[2], n360.l12.niv95[2],
  n360.l1.niv95[3], n360.l6.niv95[3], n360.l12.niv95[3])
c(n360.l1.niv95.corr[1], n360.l6.niv95.corr[1], n360.l12.niv95.corr[1], 
  n360.l1.niv95.corr[2], n360.l6.niv95.corr[2], n360.l12.niv95.corr[2],
  n360.l1.niv95.corr[3], n360.l6.niv95.corr[3], n360.l12.niv95.corr[3])



# --- n = 480 ---
c(n480.l1.niv90[1], n480.l6.niv90[1], n480.l12.niv90[1], 
  n480.l1.niv90[2], n480.l6.niv90[2], n480.l12.niv90[2],
  n480.l1.niv90[3], n480.l6.niv90[3], n480.l12.niv90[3])
c(n480.l1.niv90.corr[1], n480.l6.niv90.corr[1], n480.l12.niv90.corr[1], 
  n480.l1.niv90.corr[2], n480.l6.niv90.corr[2], n480.l12.niv90.corr[2],
  n480.l1.niv90.corr[3], n480.l6.niv90.corr[3], n480.l12.niv90.corr[3])

c(n480.l1.niv95[1], n480.l6.niv95[1], n480.l12.niv95[1], 
  n480.l1.niv95[2], n480.l6.niv95[2], n480.l12.niv95[2],
  n480.l1.niv95[3], n480.l6.niv95[3], n480.l12.niv95[3])
c(n480.l1.niv95.corr[1], n480.l6.niv95.corr[1], n480.l12.niv95.corr[1], 
  n480.l1.niv95.corr[2], n480.l6.niv95.corr[2], n480.l12.niv95.corr[2],
  n480.l1.niv95.corr[3], n480.l6.niv95.corr[3], n480.l12.niv95.corr[3])



Nsim = 5000
> # =================== Tableaux =================
> # --- n = 120 ---
[1] 0.8672 0.8660 0.8660 0.8704 0.8662 0.8566 0.8678 0.8638 0.8624
[1] 0.8822 0.8934 0.8996 0.8832 0.8948 0.8928 0.8796 0.8938 0.8946
[1] 0.9266 0.9216 0.9246 0.9318 0.9218 0.9186 0.9288 0.9248 0.9224
[1] 0.9390 0.9440 0.9464 0.9404 0.9396 0.9450 0.9388 0.9438 0.9470
> 
> # --- n = 240 ---
[1] 0.8862 0.8858 0.8814 0.8824 0.8898 0.8876 0.8848 0.8880 0.8814
[1] 0.8930 0.8990 0.8984 0.8884 0.9048 0.9032 0.8902 0.9026 0.8982
[1] 0.9398 0.9384 0.9338 0.9420 0.9462 0.9382 0.9394 0.9452 0.9354
[1] 0.9450 0.9496 0.9474 0.9454 0.9566 0.9496 0.9434 0.9550 0.9456
> 
> # --- n = 360 ---
[1] 0.8892 0.8864 0.8932 0.8860 0.8932 0.8966 0.8894 0.8870 0.8942
[1] 0.8936 0.8958 0.9040 0.8894 0.9030 0.9074 0.8936 0.8936 0.9056
[1] 0.9480 0.9388 0.9438 0.9388 0.9428 0.9496 0.9404 0.9402 0.9468
[1] 0.9512 0.9452 0.9518 0.9414 0.9498 0.9540 0.9446 0.9484 0.9532
> 
> # --- n = 480 ---
[1] 0.8992 0.8924 0.8884 0.8964 0.8912 0.8962 0.8970 0.8962 0.8968
[1] 0.9038 0.8984 0.8960 0.8996 0.8978 0.9038 0.8994 0.9036 0.9030
[1] 0.9440 0.9430 0.9376 0.9486 0.9450 0.9462 0.9472 0.9446 0.9452
[1] 0.9462 0.9480 0.9448 0.9514 0.9510 0.9512 0.9500 0.9496 0.9506





== Tableau pour LaTeX ==
$n=120$                                &  \multicolumn{9}{c}{Nominal confidence level: 90\%} \\
P.~I. Coverage                         & 0.8672 & 0.8660 & 0.8660 & 0.8704 & 0.8662 & 0.8566 & 0.8678 & 0.8638 & 0.8624 \\
P.~I. Coverage, $\bfO_P(n^{-1})$ C.~F. & 0.8822 & 0.8934 & 0.8996 & 0.8832 & 0.8948 & 0.8928 & 0.8796 & 0.8938 & 0.8946 \\
                                       &  \multicolumn{9}{c}{Nominal confidence level: 95\%} \\
P.~I. Coverage                         & 0.9266 & 0.9216 & 0.9246 & 0.9318 & 0.9218 & 0.9186 & 0.9288 & 0.9248 & 0.9224 \\
P.~I. Coverage, $\bfO_P(n^{-1})$ C.~F. & 0.9390 & 0.9440 & 0.9464 & 0.9404 & 0.9396 & 0.9450 & 0.9388 & 0.9438 & 0.9470 \\
$n=240$                                &  \multicolumn{9}{c}{Nominal confidence level: 90\%} \\
P.~I. Coverage                         & 0.8862 & 0.8858 & 0.8814 & 0.8824 & 0.8898 & 0.8876 & 0.8848 & 0.8880 & 0.8814 \\
P.~I. Coverage, $\bfO_P(n^{-1})$ C.~F. & 0.8930 & 0.8990 & 0.8984 & 0.8884 & 0.9048 & 0.9032 & 0.8902 & 0.9026 & 0.8982 \\
                                       &  \multicolumn{9}{c}{Nominal confidence level: 95\%} \\
P.~I. Coverage                         & 0.9398 & 0.9384 & 0.9338 & 0.9420 & 0.9462 & 0.9382 & 0.9394 & 0.9452 & 0.9354 \\
P.~I. Coverage, $\bfO_P(n^{-1})$ C.~F. & 0.9450 & 0.9496 & 0.9474 & 0.9454 & 0.9566 & 0.9496 & 0.9434 & 0.9550 & 0.9456 \\
$n=360$                                &  \multicolumn{9}{c}{Nominal confidence level: 90\%} \\
P.~I. Coverage                         & 0.8892 & 0.8864 & 0.8932 & 0.8860 & 0.8932 & 0.8966 & 0.8894 & 0.8870 & 0.8942 \\
P.~I. Coverage, $\bfO_P(n^{-1})$ C.~F. & 0.8936 & 0.8958 & 0.9040 & 0.8894 & 0.9030 & 0.9074 & 0.8936 & 0.8936 & 0.9056 \\
                                       &  \multicolumn{9}{c}{Nominal confidence level: 95\%} \\
P.~I. Coverage                         & 0.9480 & 0.9388 & 0.9438 & 0.9388 & 0.9428 & 0.9496 & 0.9404 & 0.9402 & 0.9468 \\
P.~I. Coverage, $\bfO_P(n^{-1})$ C.~F. & 0.9512 & 0.9452 & 0.9518 & 0.9414 0& .9498 & 0.9540 & 0.9446 & 0.9484 & 0.9532 \\
$n=480$                                &  \multicolumn{9}{c}{Nominal confidence level: 90\%} \\
P.~I. Coverage                         & 0.8992 & 0.8924 & 0.8884 & 0.8964 & 0.8912 & 0.8962 & 0.8970 & 0.8962 & 0.8968 \\
P.~I. Coverage, $\bfO_P(n^{-1})$ C.~F. & 0.9038 & 0.8984 & 0.8960 & 0.8996 & 0.8978 & 0.9038 & 0.8994 & 0.9036 & 0.9030 \\
                                       &  \multicolumn{9}{c}{Nominal confidence level: 95\%} \\
P.~I. Coverage                         & 0.9440 & 0.9430 & 0.9376 & 0.9486 & 0.9450 & 0.9462 & 0.9472 & 0.9446 & 0.9452 \\
P.~I. Coverage, $\bfO_P(n^{-1})$ C.~F. & 0.9462 & 0.9480 & 0.9448 & 0.9514 & 0.9510 & 0.9512 & 0.9500 & 0.9496 & 0.9506 \\





