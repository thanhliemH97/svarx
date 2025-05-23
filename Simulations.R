# -----------------------------------------------------------
# -- Ne pas oublier de charger les packages MTS et mvtnorm --
# -----------------------------------------------------------

VARpStable <- function(k,p, matPhi) {
  kp1 <- k+1
  pm1 <- p-1
  kfoisp <- k*p
  AR <- matrix(0, nrow= kfoisp, ncol= kfoisp)
  AR[1:k, ] <- matPhi
  AR[kp1:kfoisp, 1:(pm1*k)] <- diag( pm1*k)
  eigen.out <- eigen(AR)$values
  Mod.out <- Mod(eigen.out)
  list(AR=AR, eigen=eigen.out, Mod=Mod.out)
}

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

TousLesVARXp3m2 <- function() {
#
# Cette fonction liste tous les 2^2*2^2 modèles VARX(3,2), avec le délai p=3 toujours inclus et m=2 toujours inclus.
#
  mon.row <- 2^2
  output.var <- matrix(0,nrow=mon.row, ncol=2)
  no.ligne <- 1
  for(i1 in 0:1) {
   for(i2 in 0:1) {
     output.var[no.ligne,] <- c(i1,i2)
     no.ligne <- no.ligne+1
    }
  }
  mon.row.exo <- 2^2
  output.exo <- matrix(0,nrow=mon.row.exo, ncol=2)
  no.ligne <- 1
  for(i1 in 0:1) {
   for(i2 in 0:1) {
      output.exo[no.ligne,] <- c(i1,i2)
      no.ligne <- no.ligne+1
    }
  }
  list(output.var=output.var, output.exo=output.exo)
}

AjusterTous <- function(OutputY, mon.p= 3, InputX, mon.m=2) {
#
# Dans l'ajustement de tous les modèles, avec p=3, m=2, en conservant
# toujours les indices maximum, alors il y a 2^2 * 2^2 = 16 modèles.
# L'output va inclure les indices du modèle, et les trois colonnes
# suivantes les trois critères.
  output <- matrix(0, nrow=2^2*2^2, ncol=3+3+3)
  nbobs <- nrow(OutputY)
  TousLesCas <- TousLesVARXp3m2()
  no.mod <- 1 
  for(i in 1:2^2) {
   for(j in 1:2^2) {
     indi <- TousLesCas$output.var[i,]
     indix <- TousLesCas$output.exo[j,]
     mat.contraintes <- matR.VARXp3m2(k=ncol(OutputY), kx=ncol(InputX), ind=indi,indx=indix)
     model.out <- VARX(OutputY, p= mon.p, InputX, m=mon.m, fixed=mat.contraintes, output=F)
     npar.model <- length(model.out$coef) - sum(model.out$coef== 0)
     nobe.model <- nbobs - max(mon.p,mon.m)
     AICi <- log(det(model.out$Sigma))+  2*npar.model/nobe.model
     HQi  <-  log(det(model.out$Sigma))+ 2*log(log(nobe.model))*npar.model/nobe.model
     BICi <- log(det(model.out$Sigma))+ log(nobe.model)*npar.model/nobe.model
     output[no.mod,] <- c(indi, 1, indix, 1, AICi, HQi, BICi)
     no.mod <- no.mod+1
   }
  }
  output
}

ListerTous <- function() {
  output <- matrix(0, nrow=16,ncol=6)
  TousLesCas <- TousLesVARXp3m2()
  no.mod <- 1 
  for(i in 1:2^2) {
   for(j in 1:2^2) {
     indi <- TousLesCas$output.var[i,]
     indix <- TousLesCas$output.exo[j,]
     output[no.mod,] <- c(indi, 1, indix, 1)
     no.mod <- no.mod+1
   }
  }
  output
}

SimulationVARXp3m2 <- function(n=128, p=3, m=2, 
                               PhiX, SigmaX, theta0=theta0, matPhi=matPhi, B0=B0, matB=matB, 
                               Nsim=10) {
# Pour la partie exogène on présume un VAR(1). Les paramètres sont PhiX et SigmaX.
# Le processus est un VARX(3,2), avec Phi2 = 0, et B1=0. Ainsi les indices du vrai modèle sont c(1,0,1,1,0,1)
# Le processus exogène est généré à chaque fois.
set.seed(1)
N <- 2*n + 1
mon.k <- nrow(matPhi)
mon.kx <- nrow(PhiX)
outputOrdreAIC <- matrix(0, nrow=Nsim, ncol=2)
outputOrdreHQ <- matrix(0, nrow=Nsim, ncol=2)
outputOrdreBIC <- matrix(0, nrow=Nsim, ncol=2)
BonOrdreMaxAIC <- rep(0,Nsim)
BonOrdreMaxHQ <- rep(0,Nsim)
BonOrdreMaxBIC <- rep(0,Nsim)
outputIndicesAIC     <- matrix(0, nrow=Nsim, ncol=6)
outputModeleChoixAIC <- rep(0, Nsim)
outputBonChoixAIC    <- matrix(0, nrow=Nsim, ncol=6)
outputIndicesHQ      <- matrix(0, nrow=Nsim, ncol=6)
outputModeleChoixHQ  <- rep(0, Nsim)
outputBonChoixHQ     <- matrix(0, nrow=Nsim, ncol=6)
outputIndicesBIC     <- matrix(0, nrow=Nsim, ncol=6)
outputModeleChoixBIC <- rep(0, Nsim)
outputBonChoixBIC    <- matrix(0, nrow=Nsim, ncol=6)
for(i in 1:Nsim) {
  ExoX <- VARMAsim(N, arlags = 1, phi=PhiX, sigma=SigmaX)
  xt <- ExoX$series
  VARX.out <- simulVARX(N, theta0, p, matPhi, B0, m, matB, xt, SigmaEps)
  xt <- xt[(n+2):N,]
  zt <- VARX.out[(n+2):N,]
#
# --- Choix des ordres maximum ---
#
  VARXorder.out <- VARXorder(zt,xt,maxp=6, maxm=3,output=F)
  outputOrdreAIC[i,] <- VARXorder.out$aicor
  outputOrdreHQ[i,]  <- VARXorder.out$hqor
  outputOrdreBIC[i,] <- VARXorder.out$bicor
  BonOrdreMaxAIC[i] <- all(outputOrdreAIC[i,] == c(3,2))
  BonOrdreMaxHQ[i] <-  all(outputOrdreHQ[i,] == c(3,2))
  BonOrdreMaxBIC[i] <- all(outputOrdreBIC[i,] == c(3,2))
#
# --- Ajuster tous les sous-modèles ---
#
  AjusterTous.out <- AjusterTous(zt, p, xt, m)
#
# --- Selection AIC ---
#
  outputModeleChoixAIC[i] <- which.min(AjusterTous.out[,7])
  outputIndicesAIC[i,] <- AjusterTous.out[outputModeleChoixAIC[i],1:6]
  outputBonChoixAIC[i,] <- (outputIndicesAIC[i,] == c(1,0,1,1,0,1))
#
# --- Selection HQ ---
#
  outputModeleChoixHQ[i] <- which.min(AjusterTous.out[,8])
  outputIndicesHQ[i,] <- AjusterTous.out[outputModeleChoixHQ[i],1:6]
  outputBonChoixHQ[i,] <- (outputIndicesHQ[i,] == c(1,0,1,1,0,1))
#
# --- Selection BIC ---
#
  outputModeleChoixBIC[i] <- which.min(AjusterTous.out[,9])
  outputIndicesBIC[i,] <- AjusterTous.out[outputModeleChoixBIC[i],1:6]
  outputBonChoixBIC[i,] <- (outputIndicesBIC[i,] == c(1,0,1,1,0,1))
 }
 list(outputOrdreAIC = outputOrdreAIC,
      outputOrdreHQ = outputOrdreHQ,
      outputOrdreBIC = outputOrdreBIC,
      BonOrdreMaxAIC=BonOrdreMaxAIC,
      BonOrdreMaxHQ=BonOrdreMaxHQ,
      BonOrdreMaxBIC=BonOrdreMaxBIC,
      outputIndicesAIC=outputIndicesAIC, 
      outputModeleChoixAIC=outputModeleChoixAIC, 
      outputBonChoixAIC=outputBonChoixAIC,
      outputIndicesHQ=outputIndicesHQ, 
      outputModeleChoixHQ=outputModeleChoixHQ, 
      outputBonChoixHQ=outputBonChoixHQ,
      outputIndicesBIC=outputIndicesBIC, 
      outputModeleChoixBIC=outputModeleChoixBIC, 
      outputBonChoixBIC=outputBonChoixBIC)
}


# =====================
# ==== Simulations ====
# =====================
# Nsim  <- 1000
Nsim  <- 5000
mon.p <- 3
mon.m <- 2

# -------------------------------------------------
# --- Modèle 1: Paramètres de la partie exogène ---
# -------------------------------------------------
mon.k    <- 2
mon.kx   <- 2
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
theta0   <- c(1.0, 2.0)

n <- 12*10
Mod1n120.out <- SimulationVARXp3m2(n, mon.p, mon.m, PhiX, SigmaX, theta0, matPhi, B0, matB, Nsim)
n <- 12*20
Mod1n240.out <- SimulationVARXp3m2(n, mon.p, mon.m, PhiX, SigmaX, theta0, matPhi, B0, matB, Nsim)
n <- 12*30
Mod1n360.out <- SimulationVARXp3m2(n, mon.p, mon.m, PhiX, SigmaX, theta0, matPhi, B0, matB, Nsim)
n <- 12*40
Mod1n480.out <- SimulationVARXp3m2(n, mon.p, mon.m, PhiX, SigmaX, theta0, matPhi, B0, matB, Nsim)


# -------------------------------------------------
# --- Modèle 2: Paramètres de la partie exogène ---
# -------------------------------------------------
mon.k    <- 3
mon.kx   <- 3
PhiX     <- matrix( c(0.2, -0.6, 0.1, 0.4, 0.3, 0.7, -0.4, 0.3, 0.2), 3,3)
SigmaX   <- matrix( c(4.0, 0.8, 0.6, 0.8, 2.0, 0.1, 0.6, 0.1, 1), 3, 3)
Phi1     <- matrix( c(0.32, 0.04, -0.20,
                    0.15, -0.15, 0.23,
                    1.00, 0.30, 0.25), nrow=3, ncol=3)  
Phi2     <- matrix( 0, nrow=3,ncol=3)
Phi3     <- matrix( c(0.12, 0.35, -0.41,
                  -0.1, -0.07, 0.12,
                  -0.1, -0.10, 0.09), nrow=3, ncol=3)  
matPhi   <- cbind(Phi1, Phi2, Phi3)
SigmaEps <- SigmaSpherique(3, 0.8)
B0       <- matrix( c(1.0, 0.8, 0.6,
                  1.2, -0.9, 1.4,
                  1.0, 1.4, 1.8), nrow=3, ncol=3) 
B1       <- matrix( 0, nrow=3,ncol=3)
B2       <- matrix( c(0.8, 0.4, 0.3,
                  1.0, -0.9, 1.4,
                  0.4, 0.4, 0.8), nrow=3, ncol=3) 
matB <- cbind(B1, B2)
theta0 <- c(2.0, 1.0, 0.5)

n <- 12*10
Mod2n120.out <- SimulationVARXp3m2(n, mon.p, mon.m, PhiX, SigmaX, theta0, matPhi, B0, matB, Nsim)
n <- 12*20
Mod2n240.out <- SimulationVARXp3m2(n, mon.p, mon.m, PhiX, SigmaX, theta0, matPhi, B0, matB, Nsim)
n <- 12*30
Mod2n360.out <- SimulationVARXp3m2(n, mon.p, mon.m, PhiX, SigmaX, theta0, matPhi, B0, matB, Nsim)
n <- 12*40
Mod2n480.out <- SimulationVARXp3m2(n, mon.p, mon.m, PhiX, SigmaX, theta0, matPhi, B0, matB, Nsim)


# -------------------------------------------------
# --- Modèle 3: Paramètres de la partie exogène ---
# -------------------------------------------------
mon.k <- 5
mon.kx <- 2
PhiX <- matrix( c(0.2, -0.6, 0.3, 1), 2,2)
SigmaX <- matrix( c(4.0, 0.8, 0.8, 1), 2, 2)
Phi1 <- matrix( c(0.42, 0.14,   -0.2,  0.22,  0.18,
                  0.15, -0.15,  0.23,  0.11,  -0.1,
                  1.0,  0.30,   0.25,  0.21, 0.19,
                  0.3,  -0.18,    -0.15,  0.17,  -0.17,
                  0.2,  0.15,   0.11,  0.18,  0.11), nrow=5, ncol=5)  
Phi2 <- matrix( 0, nrow=5,ncol=5)
Phi3 <- matrix( c(0.12,  -0.15, -0.14, 0.11, 0.12,
                  0.35, -0.07,  0.10,  0.14,  -0.13,
                  -0.41, 0.12,  0.09, 0.11,  0.19,
                  0.19, 0.17,  -0.14,  0.11,  0.12,
                  0.11,  0.11,  0.12, 0.12,  0.15), nrow=5, ncol=5)  
matPhi <- cbind(Phi1, Phi2, Phi3)
round(VARpStable(k=5,p=3, matPhi)$Mod, 4)
SigmaEps <- SigmaSpherique(5, 0.8)
B0 <- matrix( c(0.2, 0.6, -0.3, 1.0, 1.2, 
                1.1, 1.8, 2.4, 5.6, 2.2), nrow=5, ncol=2)
B1 <- matrix( 0, nrow=5,ncol=2)
B2 <- matrix( c(2.4, -0.9, 0.6, 3.1, 2.4, 
                1.8,  2.0, 1.0, 2.2, 3.1), nrow=5,ncol=2)
matB <- cbind(B1, B2)
theta0 <- c(1.0, 2.0, 1.4, 1.8, 2.2)

n <- 12*10
Mod3n120.out <- SimulationVARXp3m2(n, mon.p, mon.m, PhiX, SigmaX, theta0, matPhi, B0, matB, Nsim)
n <- 12*20
Mod3n240.out <- SimulationVARXp3m2(n, mon.p, mon.m, PhiX, SigmaX, theta0, matPhi, B0, matB, Nsim)
n <- 12*30
Mod3n360.out <- SimulationVARXp3m2(n, mon.p, mon.m, PhiX, SigmaX, theta0, matPhi, B0, matB, Nsim)
n <- 12*40
Mod3n480.out <- SimulationVARXp3m2(n, mon.p, mon.m, PhiX, SigmaX, theta0, matPhi, B0, matB, Nsim)


# -------------------------------------------------
# --- Modèle 4: Paramètres de la partie exogène ---
# -------------------------------------------------
mon.k  <- 5
mon.kx <- 3
PhiX   <- matrix( c(0.2, -0.6, 0.1, 0.4, 0.3, 0.7, -0.4, 0.3, 0.2), 3,3)
SigmaX <- matrix( c(4.0, 0.8, 0.6, 0.8, 2.0, 0.1, 0.6, 0.1, 1), 3, 3)
Phi1   <- matrix( c(0.42, 0.14,   -0.2,  0.22,  0.18,
                  0.15, -0.15,  0.23,  0.11,  -0.1,
                  1.0,  0.30,   0.25,  0.21, 0.19,
                  0.3,  -0.18,    -0.15,  0.17,  -0.17,
                  0.2,  0.15,   0.11,  0.18,  0.11), nrow=5, ncol=5)  
Phi2   <- matrix( 0, nrow=5,ncol=5)
Phi3   <- matrix( c(0.12,  -0.15, -0.14, 0.11, 0.12,
                  0.35, -0.07,  0.10,  0.14,  -0.13,
                  -0.41, 0.12,  0.09, 0.11,  0.19,
                  0.19, 0.17,  -0.14,  0.11,  0.12,
                  0.11,  0.11,  0.12, 0.12,  0.15), nrow=5, ncol=5)  
matPhi <- cbind(Phi1, Phi2, Phi3)
round(VARpStable(k=5,p=3, matPhi)$Mod, 4)
SigmaEps <- SigmaSpherique(5, 0.8)
B0       <- matrix( c(1.2, -0.3, -0.6, 1.8, 1.2, 
                      1.1, 1.8, 2.4, 5.6, 2.2,
                      0.9, 0.7, 0.3, -0.2, -0.9), nrow=5, ncol=3)
B1       <- matrix( 0, nrow=5,ncol=3)
B2       <- matrix( c(1.4, 0.9, 0.4, 1.8, 2.4, 
                     -1.6, -2.0, 1.0, 2.2, 3.1,
                      3.0, 1.8, 1.2, 0.6, 2.2), nrow=5,ncol=3)
matB <- cbind(B1, B2)
theta0 <- c(1.0, 2.0, 1.4, 1.8, 2.2)

n <- 12*10
Mod4n120.out <- SimulationVARXp3m2(n, mon.p, mon.m, PhiX, SigmaX, theta0, matPhi, B0, matB, Nsim)
n <- 12*20
Mod4n240.out <- SimulationVARXp3m2(n, mon.p, mon.m, PhiX, SigmaX, theta0, matPhi, B0, matB, Nsim)
n <- 12*30
Mod4n360.out <- SimulationVARXp3m2(n, mon.p, mon.m, PhiX, SigmaX, theta0, matPhi, B0, matB, Nsim)
n <- 12*40
Mod4n480.out <- SimulationVARXp3m2(n, mon.p, mon.m, PhiX, SigmaX, theta0, matPhi, B0, matB, Nsim)


# === Traitement des résultats ===
#
# =================================
# Sélection des bons ordres maximum
# =================================
# 
# === Modèle 1 ===
PropBonOrdreMaxDGP1.AIC <- c(mean(Mod1n120.out$BonOrdreMaxAIC),
                            mean(Mod1n240.out$BonOrdreMaxAIC),
                            mean(Mod1n360.out$BonOrdreMaxAIC),
                            mean(Mod1n480.out$BonOrdreMaxAIC)
                            )
PropBonOrdreMaxDGP1.HQ <- c(mean(Mod1n120.out$BonOrdreMaxHQ),
                            mean(Mod1n240.out$BonOrdreMaxHQ),
                            mean(Mod1n360.out$BonOrdreMaxHQ),
                            mean(Mod1n480.out$BonOrdreMaxHQ)
                            )
PropBonOrdreMaxDGP1.BIC <- c(mean(Mod1n120.out$BonOrdreMaxBIC),
                            mean(Mod1n240.out$BonOrdreMaxBIC),
                            mean(Mod1n360.out$BonOrdreMaxBIC),
                            mean(Mod1n480.out$BonOrdreMaxBIC)
                            )
PropBonOrdreMaxDGP1 <- rbind(PropBonOrdreMaxDGP1.AIC,
                             PropBonOrdreMaxDGP1.HQ,
                             PropBonOrdreMaxDGP1.BIC)

# === Modèle 2 ===
PropBonOrdreMaxDGP2.AIC <- c(mean(Mod2n120.out$BonOrdreMaxAIC),
                            mean(Mod2n240.out$BonOrdreMaxAIC),
                            mean(Mod2n360.out$BonOrdreMaxAIC),
                            mean(Mod2n480.out$BonOrdreMaxAIC)
                            )
PropBonOrdreMaxDGP2.HQ <- c(mean(Mod2n120.out$BonOrdreMaxHQ),
                            mean(Mod2n240.out$BonOrdreMaxHQ),
                            mean(Mod2n360.out$BonOrdreMaxHQ),
                            mean(Mod2n480.out$BonOrdreMaxHQ)
                            )
PropBonOrdreMaxDGP2.BIC <- c(mean(Mod2n120.out$BonOrdreMaxBIC),
                            mean(Mod2n240.out$BonOrdreMaxBIC),
                            mean(Mod2n360.out$BonOrdreMaxBIC),
                            mean(Mod2n480.out$BonOrdreMaxBIC)
                            )
PropBonOrdreMaxDGP2 <- rbind(PropBonOrdreMaxDGP2.AIC,
                             PropBonOrdreMaxDGP2.HQ,
                             PropBonOrdreMaxDGP2.BIC)

# === Modèle 3 ===
PropBonOrdreMaxDGP3.AIC <- c(mean(Mod3n120.out$BonOrdreMaxAIC),
                            mean(Mod3n240.out$BonOrdreMaxAIC),
                            mean(Mod3n360.out$BonOrdreMaxAIC),
                            mean(Mod3n480.out$BonOrdreMaxAIC)
                            )
PropBonOrdreMaxDGP3.HQ <- c(mean(Mod3n120.out$BonOrdreMaxHQ),
                            mean(Mod3n240.out$BonOrdreMaxHQ),
                            mean(Mod3n360.out$BonOrdreMaxHQ),
                            mean(Mod3n480.out$BonOrdreMaxHQ)
                            )
PropBonOrdreMaxDGP3.BIC <- c(mean(Mod3n120.out$BonOrdreMaxBIC),
                            mean(Mod3n240.out$BonOrdreMaxBIC),
                            mean(Mod3n360.out$BonOrdreMaxBIC),
                            mean(Mod3n480.out$BonOrdreMaxBIC)
                            )
PropBonOrdreMaxDGP3 <- rbind(PropBonOrdreMaxDGP3.AIC,
                             PropBonOrdreMaxDGP3.HQ,
                             PropBonOrdreMaxDGP3.BIC)

# === Modèle 4 ===
PropBonOrdreMaxDGP4.AIC <- c(mean(Mod4n120.out$BonOrdreMaxAIC),
                            mean(Mod4n240.out$BonOrdreMaxAIC),
                            mean(Mod4n360.out$BonOrdreMaxAIC),
                            mean(Mod4n480.out$BonOrdreMaxAIC)
                            )
PropBonOrdreMaxDGP4.HQ <- c(mean(Mod4n120.out$BonOrdreMaxHQ),
                            mean(Mod4n240.out$BonOrdreMaxHQ),
                            mean(Mod4n360.out$BonOrdreMaxHQ),
                            mean(Mod4n480.out$BonOrdreMaxHQ)
                            )
PropBonOrdreMaxDGP4.BIC <- c(mean(Mod4n120.out$BonOrdreMaxBIC),
                            mean(Mod4n240.out$BonOrdreMaxBIC),
                            mean(Mod4n360.out$BonOrdreMaxBIC),
                            mean(Mod4n480.out$BonOrdreMaxBIC)
                            )
PropBonOrdreMaxDGP4 <- rbind(PropBonOrdreMaxDGP4.AIC,
                             PropBonOrdreMaxDGP4.HQ,
                             PropBonOrdreMaxDGP4.BIC)

PropBonOrdreMax <- rbind(PropBonOrdreMaxDGP1,
                         PropBonOrdreMaxDGP2,
                         PropBonOrdreMaxDGP3,
                         PropBonOrdreMaxDGP4)
PropBonOrdreMax
write.table(PropBonOrdreMax, "PropBonOrdreMax.txt")


# ================================================================
# Tableau contient ici les pourcentages de sélection du bon modèle.
# ================================================================
#
# === Modèle 1 ===
#
PropBonModeleDGP1.AIC <- c( mean(Mod1n120.out$outputModeleChoixAIC == 11),
                            mean(Mod1n240.out$outputModeleChoixAIC == 11),
                             mean(Mod1n360.out$outputModeleChoixAIC == 11),
                             mean(Mod1n480.out$outputModeleChoixAIC == 11)
                            )
PropBonModeleDGP1.HQ <- c( mean(Mod1n120.out$outputModeleChoixHQ == 11),
                             mean(Mod1n240.out$outputModeleChoixHQ == 11),
                             mean(Mod1n360.out$outputModeleChoixHQ == 11),
                             mean(Mod1n480.out$outputModeleChoixHQ == 11)
                            )
PropBonModeleDGP1.BIC <- c( mean(Mod1n120.out$outputModeleChoixBIC == 11),
                             mean(Mod1n240.out$outputModeleChoixBIC == 11),
                             mean(Mod1n360.out$outputModeleChoixBIC == 11),                             
                             mean(Mod1n480.out$outputModeleChoixBIC == 11)
                            )
PropBonModeleDGP1 <- rbind(PropBonModeleDGP1.AIC,PropBonModeleDGP1.HQ,PropBonModeleDGP1.BIC)
#
# === Modèle 2 ===
#
PropBonModeleDGP2.AIC <- c( mean(Mod2n120.out$outputModeleChoixAIC == 11),
                           mean(Mod2n240.out$outputModeleChoixAIC == 11),
                           mean(Mod2n360.out$outputModeleChoixAIC == 11),
                           mean(Mod2n480.out$outputModeleChoixAIC == 11)
                          )
PropBonModeleDGP2.HQ <-  c(mean(Mod2n120.out$outputModeleChoixHQ == 11),
                           mean(Mod2n240.out$outputModeleChoixHQ == 11),
                           mean(Mod2n360.out$outputModeleChoixHQ == 11),
                           mean(Mod2n480.out$outputModeleChoixHQ == 11)
                          )
PropBonModeleDGP2.BIC <- c(mean(Mod2n120.out$outputModeleChoixBIC == 11),
                           mean(Mod2n240.out$outputModeleChoixBIC == 11),
                           mean(Mod2n360.out$outputModeleChoixBIC == 11),
                           mean(Mod2n480.out$outputModeleChoixBIC == 11)
                          )
PropBonModeleDGP2 <- rbind(PropBonModeleDGP2.AIC,PropBonModeleDGP2.HQ,PropBonModeleDGP2.BIC)
#
# === Modèle 3 ===
#
PropBonModeleDGP3.AIC <- c(mean(Mod3n120.out$outputModeleChoixAIC == 11),
                           mean(Mod3n240.out$outputModeleChoixAIC == 11),
                           mean(Mod3n360.out$outputModeleChoixAIC == 11),
                           mean(Mod3n480.out$outputModeleChoixAIC == 11)
                          )
PropBonModeleDGP3.HQ <-  c(mean(Mod3n120.out$outputModeleChoixHQ == 11),
                           mean(Mod3n240.out$outputModeleChoixHQ == 11),
                           mean(Mod3n360.out$outputModeleChoixHQ == 11),
                           mean(Mod3n480.out$outputModeleChoixHQ == 11)
                          )
PropBonModeleDGP3.BIC <- c(mean(Mod3n120.out$outputModeleChoixBIC == 11),
                           mean(Mod3n240.out$outputModeleChoixBIC == 11),
                           mean(Mod3n360.out$outputModeleChoixBIC == 11),
                           mean(Mod3n480.out$outputModeleChoixBIC == 11)
                          )
PropBonModeleDGP3 <- rbind(PropBonModeleDGP3.AIC,PropBonModeleDGP3.HQ,PropBonModeleDGP3.BIC)
#
# === Modèle 4 ===
#
PropBonModeleDGP4.AIC <- c(mean(Mod4n120.out$outputModeleChoixAIC == 11),
                           mean(Mod4n240.out$outputModeleChoixAIC == 11),
                           mean(Mod4n360.out$outputModeleChoixAIC == 11),
                           mean(Mod4n480.out$outputModeleChoixAIC == 11)
                          )
PropBonModeleDGP4.HQ <-  c(mean(Mod4n120.out$outputModeleChoixHQ == 11),
                           mean(Mod4n240.out$outputModeleChoixHQ == 11),
                           mean(Mod4n360.out$outputModeleChoixHQ == 11),
                           mean(Mod4n480.out$outputModeleChoixHQ == 11)
                          )
PropBonModeleDGP4.BIC <- c(mean(Mod4n120.out$outputModeleChoixBIC == 11),
                           mean(Mod4n240.out$outputModeleChoixBIC == 11),
                           mean(Mod4n360.out$outputModeleChoixBIC == 11),
                           mean(Mod4n480.out$outputModeleChoixBIC == 11)
                          )
PropBonModeleDGP4 <- rbind(PropBonModeleDGP4.AIC,PropBonModeleDGP4.HQ,PropBonModeleDGP4.BIC)

PropBonModele <- rbind(PropBonModeleDGP1,
                       PropBonModeleDGP2,
                       PropBonModeleDGP3,
                       PropBonModeleDGP4)
PropBonModele
write.table(PropBonModele,"PropBonModele.txt")

#
# ================================================================
# Tableau contient ici les pourcentages de sélection par matrice.
# ================================================================
#
#
# === Modèle 1 ===
#
PropBonChoixMatDGP1.AIC <- c(apply(Mod1n120.out$outputBonChoixAIC,2,mean),
                             apply(Mod1n240.out$outputBonChoixAIC,2,mean),
                             apply(Mod1n360.out$outputBonChoixAIC,2,mean),
                             apply(Mod1n480.out$outputBonChoixAIC,2,mean)
                            )
PropBonChoixMatDGP1.HQ <- c( apply(Mod1n120.out$outputBonChoixHQ,2,mean),
                             apply(Mod1n240.out$outputBonChoixHQ,2,mean),
                             apply(Mod1n360.out$outputBonChoixHQ,2,mean),
                             apply(Mod1n480.out$outputBonChoixHQ,2,mean)
                            )
PropBonChoixMatDGP1.BIC <- c(apply(Mod1n120.out$outputBonChoixBIC,2,mean),
                             apply(Mod1n240.out$outputBonChoixBIC,2,mean),
                             apply(Mod1n360.out$outputBonChoixBIC,2,mean),
                             apply(Mod1n480.out$outputBonChoixBIC,2,mean)
                            )
PropBonChoixMatDGP1 <- rbind(PropBonChoixMatDGP1.AIC,PropBonChoixMatDGP1.HQ,PropBonChoixMatDGP1.BIC)
#
# === Modèle 2 ===
#
PropBonChoixMatDGP2.AIC <- c(apply(Mod2n120.out$outputBonChoixAIC,2,mean),
                             apply(Mod2n240.out$outputBonChoixAIC,2,mean),
                             apply(Mod2n360.out$outputBonChoixAIC,2,mean),
                             apply(Mod2n480.out$outputBonChoixAIC,2,mean)
                            )
PropBonChoixMatDGP2.HQ <- c( apply(Mod2n120.out$outputBonChoixHQ,2,mean),
                             apply(Mod2n240.out$outputBonChoixHQ,2,mean),
                             apply(Mod2n360.out$outputBonChoixHQ,2,mean),
                             apply(Mod2n480.out$outputBonChoixHQ,2,mean)
                            )
PropBonChoixMatDGP2.BIC <- c( apply(Mod2n120.out$outputBonChoixBIC,2,mean),
                             apply(Mod2n240.out$outputBonChoixBIC,2,mean),
                             apply(Mod2n360.out$outputBonChoixBIC,2,mean),
                             apply(Mod2n480.out$outputBonChoixBIC,2,mean)
                            )
PropBonChoixMatDGP2 <- rbind(PropBonChoixMatDGP2.AIC,PropBonChoixMatDGP2.HQ,PropBonChoixMatDGP2.BIC)
#
# === Modèle 3 ===
#
PropBonChoixMatDGP3.AIC <- c( apply(Mod3n120.out$outputBonChoixAIC,2,mean),
                             apply(Mod3n240.out$outputBonChoixAIC,2,mean),
                             apply(Mod3n360.out$outputBonChoixAIC,2,mean),
                             apply(Mod3n480.out$outputBonChoixAIC,2,mean)
                            )
PropBonChoixMatDGP3.HQ <- c( apply(Mod3n120.out$outputBonChoixHQ,2,mean),
                             apply(Mod3n240.out$outputBonChoixHQ,2,mean),
                             apply(Mod3n360.out$outputBonChoixHQ,2,mean),
                             apply(Mod3n480.out$outputBonChoixHQ,2,mean)
                            )
PropBonChoixMatDGP3.BIC <- c(apply(Mod3n120.out$outputBonChoixBIC,2,mean),
                             apply(Mod3n240.out$outputBonChoixBIC,2,mean),
                             apply(Mod3n360.out$outputBonChoixBIC,2,mean),
                             apply(Mod3n480.out$outputBonChoixBIC,2,mean)
                            )
PropBonChoixMatDGP3 <- rbind(PropBonChoixMatDGP3.AIC,PropBonChoixMatDGP3.HQ,PropBonChoixMatDGP3.BIC)
#
# === Modèle 4 ===
#
PropBonChoixMatDGP4.AIC <- c( apply(Mod4n120.out$outputBonChoixAIC,2,mean),
                             apply(Mod4n240.out$outputBonChoixAIC,2,mean),
                             apply(Mod4n360.out$outputBonChoixAIC,2,mean),
                             apply(Mod4n480.out$outputBonChoixAIC,2,mean)
                            )
PropBonChoixMatDGP4.HQ <- c(apply(Mod4n120.out$outputBonChoixHQ,2,mean),
                             apply(Mod4n240.out$outputBonChoixHQ,2,mean),
                             apply(Mod4n360.out$outputBonChoixHQ,2,mean),
                             apply(Mod4n480.out$outputBonChoixHQ,2,mean)
                            )
PropBonChoixMatDGP4.BIC <- c(apply(Mod4n120.out$outputBonChoixBIC,2,mean),
                             apply(Mod4n240.out$outputBonChoixBIC,2,mean),
                             apply(Mod4n360.out$outputBonChoixBIC,2,mean),
                             apply(Mod4n480.out$outputBonChoixBIC,2,mean)
                            )
PropBonChoixMatDGP4 <- rbind(PropBonChoixMatDGP4.AIC,PropBonChoixMatDGP4.HQ,PropBonChoixMatDGP4.BIC)


PropBonChoixMat <- rbind(PropBonChoixMatDGP1,
                       PropBonChoixMatDGP2,
                       PropBonChoixMatDGP3,
                       PropBonChoixMatDGP4)
PropBonChoixMat
write.table(PropBonChoixMat,"PropBonChoixMat.txt")

# ========= Analyses marginales parmi les modèles correctement identifiés ========
#
#
# === Modèle 1 ===
#
IterBonMod.Mod3n120.AIC <- (1:5000)[Mod3n120.out$BonOrdreMaxAIC == 1]
mean(Mod3n120.out$outputModeleChoixAIC[IterBonMod.Mod3n120.AIC] == 11)


# =================== Output =====================
> PropBonOrdreMax
                          [,1]   [,2]   [,3]   [,4]
PropBonOrdreMaxDGP1.AIC 0.4074 0.4622 0.4764 0.4942
PropBonOrdreMaxDGP1.HQ  0.7214 0.8114 0.8362 0.8586
PropBonOrdreMaxDGP1.BIC 0.9136 0.9642 0.9742 0.9808
PropBonOrdreMaxDGP2.AIC 0.1390 0.1864 0.1960 0.2072
PropBonOrdreMaxDGP2.HQ  0.4202 0.5870 0.6480 0.6878
PropBonOrdreMaxDGP2.BIC 0.7576 0.9104 0.9414 0.9598
PropBonOrdreMaxDGP3.AIC 0.3456 0.4516 0.4868 0.5064
PropBonOrdreMaxDGP3.HQ  0.7212 0.8692 0.9084 0.9306
PropBonOrdreMaxDGP3.BIC 0.9600 0.9926 0.9964 0.9986
PropBonOrdreMaxDGP4.AIC 0.0498 0.1090 0.1362 0.1408
PropBonOrdreMaxDGP4.HQ  0.2960 0.5490 0.6386 0.6744
PropBonOrdreMaxDGP4.BIC 0.7174 0.9322 0.9712 0.9820


> PropBonModele
                        [,1]   [,2]   [,3]   [,4]
PropBonModeleDGP1.AIC 0.7558 0.8030 0.8128 0.8200
PropBonModeleDGP1.HQ  0.9492 0.9770 0.9822 0.9848
PropBonModeleDGP1.BIC 0.9936 0.9990 0.9996 0.9996
PropBonModeleDGP2.AIC 0.8484 0.9010 0.9090 0.9164
PropBonModeleDGP2.HQ  0.9926 0.9984 0.9992 1.0000
PropBonModeleDGP2.BIC 1.0000 1.0000 1.0000 1.0000
PropBonModeleDGP3.AIC 0.8978 0.9442 0.9504 0.9636
PropBonModeleDGP3.HQ  0.9958 0.9994 0.9996 1.0000
PropBonModeleDGP3.BIC 1.0000 1.0000 1.0000 1.0000
PropBonModeleDGP4.AIC 0.9202 0.9618 0.9716 0.9802
PropBonModeleDGP4.HQ  0.9988 1.0000 1.0000 1.0000
PropBonModeleDGP4.BIC 1.0000 1.0000 1.0000 1.0000

> PropBonChoixMat
                        [,1]   [,2] [,3] [,4]   [,5] [,6] [,7]   [,8] [,9] [,10]  [,11] [,12] [,13]  [,14] [,15] [,16]
PropBonChoixMatDGP1.AIC    1 0.8674    1    1 0.8674    1    1 0.8970    1     1 0.8912     1     1 0.9018     1     1
PropBonChoixMatDGP1.HQ     1 0.9770    1    1 0.9714    1    1 0.9898    1     1 0.9872     1     1 0.9930     1     1
PropBonChoixMatDGP1.BIC    1 0.9958    1    1 0.9978    1    1 0.9994    1     1 0.9996     1     1 0.9996     1     1
PropBonChoixMatDGP2.AIC    1 0.9198    1    1 0.9190    1    1 0.9484    1     1 0.9494     1     1 0.9468     1     1
PropBonChoixMatDGP2.HQ     1 0.9964    1    1 0.9962    1    1 0.9992    1     1 0.9992     1     1 0.9994     1     1
PropBonChoixMatDGP2.BIC    1 1.0000    1    1 1.0000    1    1 1.0000    1     1 1.0000     1     1 1.0000     1     1
PropBonChoixMatDGP3.AIC    1 0.9766    1    1 0.9174    1    1 0.9926    1     1 0.9506     1     1 0.9932     1     1
PropBonChoixMatDGP3.HQ     1 1.0000    1    1 0.9958    1    1 1.0000    1     1 0.9994     1     1 1.0000     1     1
PropBonChoixMatDGP3.BIC    1 1.0000    1    1 1.0000    1    1 1.0000    1     1 1.0000     1     1 1.0000     1     1
PropBonChoixMatDGP4.AIC    1 0.9736    1    1 0.9412    1    1 0.9916    1     1 0.9690     1     1 0.9948     1     1
PropBonChoixMatDGP4.HQ     1 1.0000    1    1 0.9988    1    1 1.0000    1     1 1.0000     1     1 1.0000     1     1
PropBonChoixMatDGP4.BIC    1 1.0000    1    1 1.0000    1    1 1.0000    1     1 1.0000     1     1 1.0000     1     1
                         [,17] [,18] [,19]  [,20] [,21] [,22]  [,23] [,24]
PropBonChoixMatDGP1.AIC 0.8962     1     1 0.9000     1     1 0.9082     1
PropBonChoixMatDGP1.HQ  0.9892     1     1 0.9908     1     1 0.9938     1
PropBonChoixMatDGP1.BIC 1.0000     1     1 0.9998     1     1 0.9998     1
PropBonChoixMatDGP2.AIC 0.9582     1     1 0.9562     1     1 0.9578     1
PropBonChoixMatDGP2.HQ  0.9998     1     1 1.0000     1     1 1.0000     1
PropBonChoixMatDGP2.BIC 1.0000     1     1 1.0000     1     1 1.0000     1
PropBonChoixMatDGP3.AIC 0.9570     1     1 0.9970     1     1 0.9664     1
PropBonChoixMatDGP3.HQ  0.9996     1     1 1.0000     1     1 1.0000     1
PropBonChoixMatDGP3.BIC 1.0000     1     1 1.0000     1     1 1.0000     1
PropBonChoixMatDGP4.AIC 0.9764     1     1 0.9966     1     1 0.9834     1
PropBonChoixMatDGP4.HQ  1.0000     1     1 1.0000     1     1 1.0000     1
PropBonChoixMatDGP4.BIC 1.0000     1     1 1.0000     1     1 1.0000     1
> 



# =============== Menage ===============
rm(Mod1n120.out, Mod1n240.out, Mod1n360.out, Mod1n480.out,
Mod2n120.out,Mod2n240.out,Mod2n360.out,Mod2n480.out,
Mod3n120.out,Mod3n240.out,Mod3n360.out,Mod3n480.out,
Mod4n120.out,Mod4n240.out,Mod4n360.out,Mod4n480.out,
PropBonOrdreMaxDGP1.AIC,PropBonOrdreMaxDGP2.AIC,PropBonOrdreMaxDGP3.AIC,
PropBonOrdreMaxDGP4.AIC,
PropBonOrdreMaxDGP1.HQ,PropBonOrdreMaxDGP2.HQ,PropBonOrdreMaxDGP3.HQ,
PropBonOrdreMaxDGP4.HQ,
PropBonOrdreMaxDGP1.BIC,PropBonOrdreMaxDGP2.BIC,PropBonOrdreMaxDGP3.BIC,
PropBonOrdreMaxDGP4.BIC,
PropBonOrdreMaxDGP1,PropBonOrdreMaxDGP2,PropBonOrdreMaxDGP3,PropBonOrdreMaxDGP4,
PropBonOrdreMax,
PropBonModeleDGP1.AIC,PropBonModeleDGP2.AIC,PropBonModeleDGP3.AIC,
PropBonModeleDGP4.AIC,
PropBonModeleDGP1.HQ,PropBonModeleDGP2.HQ,PropBonModeleDGP3.HQ,PropBonModeleDGP4.HQ,
PropBonModeleDGP1.BIC,PropBonModeleDGP2.BIC,PropBonModeleDGP3.BIC,PropBonModeleDGP4.BIC,
PropBonModeleDGP1,PropBonModeleDGP2,PropBonModeleDGP3,PropBonModeleDGP4,
PropBonModele,
PropBonChoixMatDGP1.AIC,PropBonChoixMatDGP2.AIC,PropBonChoixMatDGP3.AIC,
PropBonChoixMatDGP4.AIC,
PropBonChoixMatDGP1.HQ,PropBonChoixMatDGP2.HQ,PropBonChoixMatDGP3.HQ,
PropBonChoixMatDGP4.HQ,
PropBonChoixMatDGP1.BIC,PropBonChoixMatDGP2.BIC,PropBonChoixMatDGP3.BIC,
PropBonChoixMatDGP4.BIC,
PropBonChoixMatDGP1,PropBonChoixMatDGP2,PropBonChoixMatDGP3,PropBonChoixMatDGP4,
PropBonChoixMat)
