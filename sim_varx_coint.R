source("~/these_doctorat/Simulations/fctVARXcoint.R")

phi1 = matrix(c(0.6,.12,1,.7),2,2)
v0 = matrix(c(2,4),2,1)
sigma.e = matrix(c(25,5.4,5.4,9.0),2,2)
phix.1 = matrix(0.8)
sigma.b = matrix(4)

#vérifier les valeurs propres de la matrice compagnon
eigen(phi1)$values #valeur propre 1 <=> racine 1 pour det{\Phi(z)} = 0
                    #valeur propre 0.3 <=> racine  10/3 pour det{\Phi(z)} = 0
qr(phi1)$rank #rang 2
qr(diag(1,2)-phi1)$rank #rang 1

sim1 = function(i,n) {
  out = genVARX(n, phi1, v0, sigma.e, p=1, s=0, phix.1, sigma.b, pX=1)
  Y = out$y
  X = out$x

  W = diff(Y)
  Ylag1 = Y[1:(n-1),]
  Xlag0 = X[2:n,]
  Z = cbind(Ylag1, Xlag0)
  
  #Estimateurs des moindres carrés de rang plein
  
  para.ls = t(solve(t(Z)%*%Z, t(Z)%*%W)) 
  resi = W - Z %*% t(para.ls)
  cov.ls = (t(resi) %*% resi) / nrow(resi)
  
  #Estimation de rang réduit
  
  C.ls = para.ls[,1:2]
  out2 = reducedRank(C.ls, cov.ls)

  #Estimateurs de vraisemblance maximale
  
  f = function(par) {
    
    A = matrix(par[1:2],2,1)
    B = matrix(c(1,par[3]),1,2)
    V0 = matrix(par[4:5],2,1)
    
    Epsilon =  W - Ylag1 %*% t(A%*%B) - Xlag0 %*% t(V0)
    
    cov = (t(Epsilon) %*% Epsilon) / nrow(Epsilon)
    
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
  Epsilon =  W - Ylag1 %*% t(A.est%*%B.est) - Xlag0 %*% t(V0.est) 
  cov.mle = (t(Epsilon) %*% Epsilon) / nrow(Epsilon)

  return(unname(c(start.ls, cov.ls, par.mle, cov.mle)))
}

library(parallel)
cl = makeCluster(detectCores()-1)

n = 400 #changer pour 50, 100, 200, 300, 400
clusterExport(cl, 
              varlist = c("n", "phi1", "v0", "sigma.e", "phix.1","sigma.b",
                          "genVARX","reducedRank","sim1"))

clusterSetRNGStream(cl, iseed = 1)
tic = Sys.time()
nsim = 10000
res = t(parSapply(cl, 1:nsim, function(i) {sim1(i, n)}))
(tac = Sys.time() - tic)
round(colMeans(res),3)
round(colSD(res),3)

#Vérifier la convergence du vecteur aléatoire du théorème 4 de Ahn et Reinsel
#généralisé par l'inclusion de la variable exogène

sim2 = function(i,n) {
  out = genVARX(n, phi1, v0, sigma.e, p=1, s=0, phix.1, sigma.b, pX=1)
  Y = out$y
  X = out$x
  
  W=diff(Y)
  Ylag1 = Y[1:(n-1),]
  Xlag0 = X[2:n,]
  Z = cbind(Ylag1, Xlag0)
  
  #Estimateurs des moindres carrés de rang plein
  
  para.ls = t(solve(t(Z)%*%Z, t(Z)%*%W)) 
  resi = W - Z %*% t(para.ls)
  cov.ls = (t(resi) %*% resi) / nrow(resi)
  
  #Estimation de rang réduit
  
  C.ls = para.ls[,1:2]
  out2 = reducedRank(C.ls, cov.ls)
  
  #Estimateurs de vraisemblance maximale
  
  f = function(par) {
    
    A = matrix(par[1:2],2,1)
    B = matrix(c(1,par[3]),1,2)
    V0 = matrix(par[4:5],2,1)
    
    Epsilon =  W - Ylag1 %*% t(A%*%B) - Xlag0 %*% t(V0)
    
    cov = (t(Epsilon) %*% Epsilon) / nrow(Epsilon)
    
    loglike = log(det(cov))
    return(loglike)
  }
  
  #départ avec les estimés des moindres carrés
  
  start.ls = c(out2$A,out2$B0, para.ls[,3])
  
  #estimés du maximum de vraisemblance
  
  out3 = nlminb(start = start.ls, objective = f)
  par.mle = out3$par
  
  #vecteur aléatoire du théorème 4
  b0 = par.mle[3]
  stat.lse = sqrt(sum(Ylag1[,2]^2))*(out2$B0+2.5) #car b0 = -2.5
  stat.mle = sqrt(sum(Ylag1[,2]^2))*(b0+2.5)
  
  return(c(stat.lse,stat.mle))
}

library(parallel)
cl = makeCluster(detectCores()-1)

n = 3000 #changer pour 50, 100, 200, 300, 400
clusterExport(cl, 
              varlist = c("n", "phi1", "v0", "sigma.e", "phix.1","sigma.b",
                          "genVARX","reducedRank","sim2"))

clusterSetRNGStream(cl, iseed = 1)
tic = Sys.time()
nsim = 10000

res = t(parSapply(cl, 1:nsim, function(i) {sim2(i, n)}))
(tac = Sys.time() - tic)
round(colMeans(res),3)
round(colVar(res),3)

#means
#lse      mle
 0.363;  0.212 #n=50
 0.129;  0.064 #n=100
 0.182;  0.152 #n=200
-0.046; -0.066 #n=300
-0.106; -0.121 #n=400
 0.096;  0.084 #n=500
 0.036;  0.026 #n=600
-0.011; -0.019 #n=700
-0.034; -0.041 #n=800
-0.040; -0.046 #n=900
-0.042; -0.049 #n=1000
 0.079;  0.074 #n=1500
 0.015;  0.008 #n=2000

#variances
#lse      mle
109.966;110.297 #n=50
 97.383; 97.408 #n=100
 91.517; 91.532 #n=200
 88.866; 88.870 #n=300
 84.980; 84.984 #n=400
 85.729; 85.727 #n=500
 85.250; 85.249 #n=600
 83.941; 83.938 #n=700
 86.375; 86.375 #n=800
 86.418; 86.420 #n=900
 84.186; 84.190 #n=1000
 84.074; 84.072 #n=1500
 83.660; 83.658 #n=2000

#vraie variance : (t(A)%*%\sigma.e^-1%*%A))^-1

A = matrix(c(-.4,.12), 2,1)
sigma.e = matrix(c(25,5.4,5.4,9.0), 2,2)
solve(t(A)%*%solve(sigma.e)%*%A) 
84.472