source("~/these_doctorat/Simulations/fctVARXcoint.R")

phi1 = matrix(c(0.6,.12,1,.7),2,2)
v0 = matrix(c(2,4),2,1)
sigma.e = matrix(c(25,5.4,5.4,9.0),2,2)
phix.1 = matrix(0.8)
sigma.b = matrix(4)

#vérifier les valeurs propres de la matrice companion
eigen(phi1)$values #valeur propre 1 <=> racine 1 pour det{\Phi(z)} = 0
                    #valeur propre 0.3 <=> racine  10/3 pour det{\Phi(z)} = 0

sim1 = function(i,n) {
  out = genVARX(n, phi1, v0, sigma.e, p=1, s=0, phix.1, sigma.b, pX=1)
  Y = out$y
  X = out$x
  
  W=diff(Y)
  Ylag1 = Y[1:(n-1),]
  Z = cbind(Ylag1, X[2:n,])
  
  #Estimateurs des moindres carrés de rang plein
  
  para.ls = t(solve(t(Z)%*%Z, t(Z)%*%W)) 
  resi = W - Z %*% t(para.ls)
  cov.ls = (t(resi) %*% resi) / nrow(W)
  
  #Estimation de rang réduit
  out2 = reducedRank(para.ls, cov.ls)

  #Estimateurs de vraisemblance maximale
  
  f = function(par) {
    
    A = matrix(par[1:2],2,1)
    B = matrix(c(1,par[3]),1,2)
    V0 = matrix(par[4:5],2,1)
    
    Epsilon =  W - Ylag1 %*% t(A%*%B) - X[2:n,] %*% t(V0)
    
    cov = (t(Epsilon) %*% Epsilon) / nrow(W)
    
    loglike = log(det(cov))
    return(loglike)
  }
  
  #départ avec les estimés des moindres carrés
  start.ls = c(out2$A,out2$B0, para.ls[,3])
  
  # estimés du maximum de vraisemblance
  out3 = nlminb(start = start.ls, objective = f)
  par.mle = out3$par
  
  # covariance par maximum de vraisemblance
  A.est = matrix(par.mle[1:2],2,1)
  B.est = matrix(c(1,par.mle[3]),1,2)
  V0.est = matrix(par.mle[4:5],2,1)

  #mêmes donnees, seuls les paramètres changent  
  Epsilon =  W - Ylag1 %*% t(A.est%*%B.est) - X[2:n,] %*% t(V0.est) 
  
  cov.mle = (t(Epsilon) %*% Epsilon) / nrow(W)
  
  res = c(start.ls, cov.ls, par.mle, cov.mle)
  return(unname(c(start.ls, cov.ls, par.mle, cov.mle)))
}

library(parallel)
cl = makeCluster(detectCores()-1)

n = 50 #changer pour 50, 100, 200, 300, 400
clusterExport(cl, 
              varlist = c("sim1", "n", "phi1", "v0", "sigma.e", "phix.1",
                          "sigma.b","genVARX","reducedRank"))

clusterSetRNGStream(cl, iseed = 1)
tic = Sys.time()
nsim = 10000
res = t(parSapply(cl, 1:nsim, function(i) {sim1(i, n)}))
(tac = Sys.time() - tic)
round(colMeans(res),5)
round(colSD(res),5)

#brouillon

phi1 = matrix(c(0.6,.12,1,.7),2,2)
v0 = matrix(c(2,4),2,1)
sigma.e = matrix(c(25,5.4,5.4,9.0),2,2)
phix.1 = matrix(0.8)
sigma.b = matrix(4)

n = 100
set.seed(1)
out = genVARX(n, phi1, v0, sigma.e, p=1, s=0, phix.1, sigma.b, pX=1)
Y = out$y
X = out$x

W=diff(Y)
Ylag1 = Y[1:(n-1), ]
Xtrunc = X[2:n, ]
Z = cbind(Ylag1, Xtrunc)

#Estimateurs des moindres carrés à rang plein

(C.ls = t(ols(Z,W)$parameters) )
(cov.ls = ols(Z,W)$covariance)

#Le bon C est :
C = matrix(c(-.4,.12,1,-.3),2,2)

#Estimation de rang réduit

out2 = reducedRank(C.ls, cov.ls)

#Les vrais A et B sont :
(A = matrix(c(-.4,.12),2,1))
(B = matrix(c(1,-2.5),1,2))
(C = A%*%B)

#Estimateurs de vraisemblance maximale

f = function(par) {
  
  A = matrix(par[1:2],2,1)
  B = matrix(c(1,par[3]),1,2)
  V0 = matrix(par[4:5],2,1)
  Epsilon =  W - Ylag1 %*% t(A%*%B) - X[2:n,] %*% t(V0)
  
  cov = (t(Epsilon) %*% Epsilon) / nrow(W)
  
  loglike = log(det(cov))
  return(loglike)
}

#le départ avec les estimées par moindres carrés
(start.ls = c(out2$A,out2$B0, para.ls[,3]))
out3 = nlminb(start = start.ls, objective = f)
(par.mle=out3$par)

#les vrais paramètres :
start = c (-.4, .12, -2.5, 2, 4)
A = matrix(par.mle[1:2],2,1)
B = matrix(c(1,par.mle[3]),1,2)
V0 = matrix(par.mle[4:5],2,1)

#nouveau sigma
Epsilon =  W - Ylag1 %*% t(A%*%B) - X[2:n,] %*% t(V0)

(cov.mle = (t(Epsilon) %*% Epsilon) / nrow(W))

#Comparaison avec la fonction personnelle qui exploite les équations de
#Newton-Raphson pour estimation de rang reduit selon Yap et Reinsel 

newtonRaphson1(start.ls)

#Vérifier la convergence du vecteur aléatoire du théorème 4 de Ahn et Reinsel

sim2 = function(i,n) {
  out = genVARX(n, phi1, v0, sigma.e, p=1, s=0, phix.1, sigma.b, pX=1)
  Y = out$y
  X = out$x
  
  W=diff(Y)
  Ylag1 = Y[1:(n-1),]
  Z = cbind(Ylag1, X[2:n,])
  
  #Estimateurs des moindres carrés de rang plein
  
  C.ls = t(solve(t(Z)%*%Z, t(Z)%*%W)) 
  resi = W - Z %*% t(C.ls)
  cov.ls = (t(resi) %*% resi) / nrow(W)
  
  #Estimation de rang réduit
  
  out2 = reducedRank(C.ls, cov.ls)
  
  #Estimateurs de vraisemblance maximale
  
  f = function(par) {
    
    A = matrix(par[1:2],2,1)
    B = matrix(c(1,par[3]),1,2)
    V0 = matrix(par[4:5],2,1)
    
    Epsilon =  W - Ylag1 %*% t(A%*%B) - X[2:n,] %*% t(V0)
    
    cov = (t(Epsilon) %*% Epsilon) / nrow(W)
    
    loglike = log(det(cov))
    return(loglike)
  }
  
  #départ avec les estimés des moindres carrés
  start.ls = c(out2$A,out2$B0, C.ls[,3])
  
  # estimés du maximum de vraisemblance
  out3 = nlminb(start = start.ls, objective = f)
  par.mle = out3$par
  
  b0 = par.mle[3]
  stat.lse = sqrt(sum(Ylag1[,2]^2))*(out2$B0+2.5)
  stat.mle = sqrt(sum(Ylag1[,2]^2))*(par.mle[3]+2.5)
  return(c(stat.lse,stat.mle))
}

library(parallel)
cl = makeCluster(detectCores()-1)

n = 300 #changer pour 50, 100, 200, 300, 400
clusterExport(cl, 
              varlist = c("sim2", "n", "phi1", "v0", "sigma.e", "phix.1",
                          "sigma.b","genVARX","reducedRank"))

clusterSetRNGStream(cl, iseed = 1)
tic = Sys.time()
nsim = 10000

res = t(parSapply(cl, 1:nsim, function(i) {sim2(i, n)}))
(tac = Sys.time() - tic)
colVar(res)

#lse      mle
109.9659; 110.2971 #n = 50
97.38322; 97.40790 #n = 100
91.51731; 91.53219 #n = 200
88.86611; 88.87049 #n = 300
84.98013; 84.98368 #n = 400
85.72940; 85.72739 #n = 500
85.24994; 85.24906 #n = 600
83.94052; 83.93837 #n = 700
86.37483; 86.37478 #n = 800
86.41822; 86.41990 #n = 900
84.18582; 84.18950 #n = 1000

#vraie variance : (t(A)%*%\sigma.e^-1%*%A))^-1

A = matrix(c(-.4,.12), 2,1)
B = matrix(c(1,-2.5), 1,2)
sigma.e = matrix(c(25,5.4,5.4,9.0), 2,2)
solve(t(A)%*%solve(sigma.e)%*%A) 
84.47205
