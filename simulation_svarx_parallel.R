source("~/Simulations/fctSVARX.R")
library(parallel)

#SVARX({1,3},{0,2})
#DGP1 (voir table 2)
#Les séries y et x sont bidimensionnelles.
#Il y a le délai 1 et 3 pour y et délai 0 et 2 pour x.

Phi1 = matrix(c(.5, .15, -.4, .5), nrow=2, byrow=T)

Phi2 = matrix(c(.2, .1, -.2, .25), 2, byrow=T)

###############
#Confirmer la causalité du processus

bigPhi = rbind(cbind(Phi1,0,0,Phi2),
               cbind(diag(1,4),0,0))
round(Mod(eigen(bigPhi)$values),4) # tous inférieurs à 1
###############


phi1 = matrix(c(1,2),2,1)

Beta1 = matrix(c(.2, -.3, .6, 1), 2, byrow=T) 
Beta2 = matrix(c(2.4, .6, -.9, 3.1), 2, byrow=T) 

bigBeta = rbind(cbind(Beta1,0,0, Beta2), 0,0,0,0)

rho = .8; 

Sigma1 = diag((1-rho),2) + rho*matrix(1,2,2)

Phi1.X = matrix(c(.2, .3, -.6, 1), nrow=2, byrow=T)

Sigma1.X = matrix(c(4,.8,.8,1), 2)

Phi3 = matrix(c(.32, .15, 1, 
                .4, -.15, .3,
                -.2, .23, .25), 3, byrow=T)
Phi4 = matrix(c(.12, -.1, -.1, 
                .35, -.07, -.1,
                -.41, .12, .09), 3, byrow=T)
##############
#Confirmer la causalité du processus
bigPhi = rbind(cbind(Phi3,0,0,0,Phi4),
               cbind(diag(1,6),0,0,0))
round(Mod(eigen(bigPhi)$values),4) # tous inférieurs à 1
##############
phi2 = matrix(c(2,1,.5),3,1)

Beta3 = matrix(c(1, 1.2, 1,
                 .8, -.9, 1.4,
                 .6, 1.4, 1.8), 3, byrow=T)
Beta4 = matrix(c(.8, 1, .4,
                 .4, -.09, .4,
                 .3, 1.4, .8), 3, byrow=T)

rho = .8; 
Sigma2 = diag((1-rho),3) + rho*matrix(1,3,3)

Phi2.X = matrix(c(.2, .4, -.4,
                  -.6, .3, .3,
                  .1, .7, .2), nrow=3, byrow=T)

Sigma2.X = matrix(c(4,.8,.6,
                    .8,2,.1,
                    .6,.1,1),3)

#DGP3 (voir table 2)
#Les séries y et x sont 5 dimensions et 2 dimensions.
#Il y a le délai 1 et 3 pour y et délai 0 et 2 pour x.

Phi5 = matrix(c(.42, .15, 1, .3, .2, 
                .14, -.15, .3, -.18, .15,
                -.2, .23, .25, -.15, .11,
                .22, .11, .21, .17, .18,
                .18, -.1, .19, -.17, .11), 5, byrow=T)
Phi6 = matrix(c(.12, .35, -.41, .19, .11, 
                -.15, -.17, .12, .17, .11,
                -.11, .1, .09, -.14, .12,
                .11, .14, .11, .11, .12,
                .12, -.13, .19, .12, .15), 5, byrow=T)

phi3 = matrix(c(1,2,1.4, 1.8, 2.2),5,1)

Beta5 = matrix(c(2, 1.1,
                 .6, 1.8,
                 -.3, 2.4,
                 1, 5.6,
                 1.2, 2.2), 5,2, byrow=T)
Beta6 = matrix(c(2.4, 1.8,
                 -.9, 2,
                 -.6, 1,
                 3.1, 2.2,
                 2.4, 3.1), 5,2, byrow=T)

rho = .8

#Confirmer la causalité du processus
bigPhi = rbind(cbind(Phi5,0,0,0,0,0,Phi6),
               cbind(diag(1,10),0,0,0,0,0))
round(Mod(eigen(bigPhi)$values),4) # tous inférieurs à 1
##############

Sigma3 = diag((1-rho),5) + rho*matrix(1,5,5)

Phi1.X = matrix(c(.2, .3, -.6, 1), nrow=2, byrow=T)

Sigma1.X = matrix(c(4,.8,.8,1), 2)

#DPG4
#Les séries y et x sont 5 dimensions et 3 dimensions.
#Il y a le délai 1 et 3 pour y et délai 0 et 2 pour x.

Beta7 = matrix(c(1.2, 1.1, .9,
                 -.3, 1.8, .7,
                 .6, 2.4, .3,
                 1.8, 5.6, -.2,
                 1.2, 2.2, -.9), 5,3, byrow=T)
Beta8 = matrix(c(1.4, -1.6, 3,
                 .9, -2, 1.8,
                 .4, 1, 1.2,
                 1.8, 2.2, .6,
                 2.4, 3.1, 2.2), 5,3, byrow=T)

#Expérience 1 : Trouver p et q minimisant le critère de sélection
#Trouver p et q à l'aide de VARXorder version corrigée

phi = list(Phi1,Phi2)
beta = list(Beta1,Beta2)
sigma = Sigma1
phiX = Phi1.X
sigmaX = Sigma1.X
phi0 = phi1

sim1a = function(sim,n) {
  out = genSVARX(n, phi, beta, sigma, c(1, 3), c(0, 2), phiX, sigmaX, phi0)
  obj = VARXorder(out$y, out$x, 6, 3, output=F)
  result = c(paste(obj$aicor,collapse = "|"), paste(obj$hqcor,collapse = "|"), 
             paste(obj$bicor,collapse = "|"))
  return(result)
}
sim1b = function(sim,n) {
  out = genSVARX(n, phi, beta, sigma, c(1, 3), c(0, 2), phiX, sigmaX, phi0)
  
  # Calcul des critères pour chaque matrice
  vecAIC = aic.svarx(y=out$y,x=out$x, p=3, s=2)$IJ
  vecHQC = aic.svarx(y=out$y,x=out$x, p=3, s=2, crit="hqc")$IJ
  vecBIC = aic.svarx(y=out$y,x=out$x, p=3, s=2, crit="bic")$IJ
  
  # Sélection du meilleur modèle pour chaque critère
  res = c(vecAIC,vecHQC,vecBIC)
  return(res)
}

# Configuration des clusters
cl = makeCluster(detectCores()-1) #on crée des copies de R qui exécutent
#en même temps. Le -1 pour laisser un coeur libre pour l'ordinateur ne ralentit pas

# Charger les packages nécessaires sur chaque noeud
clusterEvalQ(cl, library(MASS))

nsim = 5000

n = 120 #changer n pour 120, 240, 360, 480
# Exporter toutes les variables et fonctions nécessaires
clusterExport(cl, 
              varlist = c("genSVARX", "n", "phi", "beta", "sigma", "phiX",
                          "sigmaX", "phi0","VARXorder","subsets.lags",
                          "subsets.lags.exo","aic.svarx","sim1a","sim1b"))

clusterSetRNGStream(cl, iseed = 1)
tic = Sys.time()
whatselect_a = t(parSapply(cl, 1:nsim, function(i) {sim1a(i, n)}))
tac = Sys.time() - tic

aicSelection_a = sort(table(whatselect_a[,1]),decreasing = T)
hqcSelection_a = sort(table(whatselect_a[,2]),decreasing = T)
bicSelection_a = sort(table(whatselect_a[,3]),decreasing = T)

#sink(file = "SVARX-DPG1.txt", append = T)
cat("Expérience 1 :","\n")
cat("Nombre de simulations :", nsim, "\n")
cat("Nombre d'observations :", n, "\n")
tac
print("AIC")
#aicSelection_a
aicSelection_a/nsim
print("HQC")
#hqcSelection_a
hqcSelection_a/nsim
print("BIC")
#bicSelection_a
bicSelection_a/nsim
tac
#cat("\n")
#sink()
#tac

tic = Sys.time()
clusterSetRNGStream(cl, iseed = 1)
whatselect_b = t(parSapply(cl, 1:nsim, function(i) {sim1b(i, n)}))
tac = Sys.time() - tic

aicSelection_b = sort(table(whatselect_b[,1]),decreasing = T)
hqcSelection_b = sort(table(whatselect_b[,2]),decreasing = T)
bicSelection_b = sort(table(whatselect_b[,3]),decreasing = T)

#sink(file = "SVARX-DPG1.txt", append = T)
cat("Expérience 2 :","\n")
cat("Nombre de simulations :", nsim, "\n")
cat("Nombre d'observations :", n, "\n")
tac
print("AIC")
#aicSelection_b
aicSelection_b/nsim
print("HQC")
#hqcSelection_b
hqcSelection_b/nsim
print("BIC")
#bicSelection_b
bicSelection_b/nsim
cat("\n")
#sink()
tac

stopCluster(cl)
