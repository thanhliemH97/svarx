source("~/these_doctorat/Simulations/fctVARXcoint.R")

D = diag(c(1, .5, .25))

P = matrix(c(4, 0, 1,
             0, 1, 1,
             1, 1, 0), 3, 3, byrow=T)
det(P) # = 0

(sumofPhi = P %*% D %*% solve(P))
sumofPhi = round(sumofPhi,2)
eigen(sumofPhi)$values #une valeur propre 1, les autres inférieurs à 1 

(Phi1 = matrix(c(0.4, -0.7, -0.2,
                -0.1, 0.5, -0.1,
                 0.1, -0.1, 0.3),3,3, byrow=T))
(Phi2 = round(sumofPhi - Phi1,2))

(C = round(sumofPhi - diag(1,3),2))
round(eigen(C)$values,4) # une valeur propre 0
qr(C)$rank #rang 2 = 3-1

#les vrais A et B
A = C[,1:2]
B = cbind(diag(1,2), matrix(c(-4,0),2,1))
A%*%B == C

bigPhi = rbind(cbind(Phi1,Phi2),
               cbind(diag(1,3),0,0,0))

round(Mod(eigen(bigPhi)$values),4) #une valeur propre 1 et 5 inférieurs à 1. 


V02 = matrix(c(0.4,0.6,0,0.3,0.6,-0.2,0,-0.1,0.7),3,3)
sigma.e2 = matrix(c(25,5,3,5,9,4,3,4,16),3,3)
phix.2 = matrix(c( 0.2,-.6,.1,.4,.3,.7,-.4,.3,.2),3,3) 
sigma.b2 = matrix(c(4,-.8,-.6,-.8, 2,-.1,-.6,-.1,1),3,3)

Mod(eigen(phix.2)$values) #valeurs propres inférieurs à 1

cbind(C,-Phi2,V02)

sim1 = function(i,n) {
  out = genVARX(n, list(Phi1, Phi2), V02, sigma.e2, p=2, s=0, phix.2, sigma.b2,
                pX=1)
  Y = out$y
  X = out$x
  W = diff(Y)
  Wlag0 = W[2:(n-1),]; Wlag1 = W[1:(n-2),]
  Ylag1 = Y[2:(n-1),]
  Z = cbind(Ylag1, Wlag1, X[3:n,])
  
  #Estimateurs des moindres carrés de rang plein
  
  para.ls = t(solve(t(Z)%*%Z, t(Z)%*%Wlag0)) 
  resi = Wlag0 - Z %*% t(para.ls)
  cov.ls = (t(resi) %*% resi) / nrow(Wlag0)
  
  #Estimation de rang réduit
  C.ls = para.ls[,1:3]
  out2 = reducedRank(C.ls, cov.ls, r=2)
  
  #Estimateurs de vraisemblance maximale
  
  f = function(par) {
    
    A = matrix(par[1:6],3,2)
    B = matrix(cbind(diag(1,2),par[7:8]),2,3)
    Phi1star = matrix(par[9:17],3,3)
    V0 = matrix(par[18:26],3,3)
    
    Epsilon =  Wlag0 - Ylag1 %*% t(A%*%B) - Wlag1 %*% t(Phi1star) - 
      X[3:n,] %*% t(V0)
    
    cov = (t(Epsilon) %*% Epsilon) / nrow(Wlag0)
    
    loglike = log(det(cov))
    return(loglike)
  }
  
  #départ avec les estimés des moindres carrés
  start.ls = c(out2$A,out2$B0, para.ls[,-(1:3)])
  
  # estimés du maximum de vraisemblance
  out3 = nlminb(start = start.ls, objective = f)
  par.mle = out3$par
  
  # covariance par maximum de vraisemblance
  A.est = matrix(par.mle[1:6],3,2)
  B.est = matrix(cbind(diag(1,2),par.mle[7:8]),2,3)
  Phi1star.est = matrix(par.mle[9:17],3,3)
  V0.est = matrix(par.mle[18:26],3,3)
  
  #mêmes donnees, seuls les paramètres changent  
  Epsilon =  Wlag0 - Ylag1 %*% t(A.est %*% B.est) - Wlag1 %*% t(Phi1star.est) - 
    X[3:n,] %*% t(V0.est)
  
  cov.mle = (t(Epsilon) %*% Epsilon) / nrow(Wlag0)
  
  res = c(start.ls, cov.ls, par.mle, cov.mle)
  return(unname(c(start.ls, cov.ls, par.mle, cov.mle)))
}

library(parallel)
cl = makeCluster(detectCores()-1)

n = 400 #changer pour 50, 100, 200, 300, 400
clusterExport(cl, 
              varlist = c("sim1", "n", "Phi1", "Phi2", "V02", "sigma.e2", 
                          "phix.2", "sigma.b2", "genVARX", "reducedRank"))

clusterSetRNGStream(cl, iseed = 1)
tic = Sys.time()
nsim = 10000
res = t(parSapply(cl, 1:nsim, function(i) {sim1(i, n)}))
(tac = Sys.time() - tic)
round(colMeans(res),5)
round(colSD(res),5)

#Vérifier la convergence du vecteur aléatoire du théorème 4 de Ahn et Reinsel

sim2 = function(i,n) {
  out = genVARX(n, list(Phi1, Phi2), V02, sigma.e2, p=2, s=0, phix.2, sigma.b2,
                pX=1)
  Y = out$y
  X = out$x
  rm(out)
  W = diff(Y)
  Wlag0 = W[2:(n-1),]; Wlag1 = W[1:(n-2),]
  Ylag1 = Y[2:(n-1),]
  Z = cbind(Ylag1, Wlag1, X[3:n,])
  
  #Estimateurs des moindres carrés de rang plein
  
  para.ls = t(solve(t(Z)%*%Z, t(Z)%*%Wlag0)) 
  resi = Wlag0 - Z %*% t(para.ls)
  cov.ls = (t(resi) %*% resi) / nrow(Wlag0)
  
  #Estimation de rang réduit
  C.ls = para.ls[,1:3]
  out2 = reducedRank(C.ls, cov.ls, r=2)
  
  #Estimateurs de vraisemblance maximale
  
  f = function(par) {
    
    A = matrix(par[1:6],3,2)
    B = matrix(cbind(diag(1,2),par[7:8]),2,3)
    Phi1star = matrix(par[9:17],3,3)
    V0 = matrix(par[18:26],3,3)
    
    Epsilon =  Wlag0 - Ylag1 %*% t(A%*%B) - Wlag1 %*% t(Phi1star) - 
      X[3:n,] %*% t(V0)
    
    cov = (t(Epsilon) %*% Epsilon) / nrow(Wlag0)
    
    loglike = log(det(cov))
    return(loglike)
  }
  
  #départ avec les estimés des moindres carrés
  start.ls = c(out2$A,out2$B0, para.ls[,-(1:3)])
  
  # estimés du maximum de vraisemblance
  out3 = nlminb(start = start.ls, objective = f)
  par.mle = out3$par
  # vecteur aléatoire du théorème 4
  B0.lse = c(out2$B0)
  B0.mle = par.mle[7:8]
  B0 = c(-4,0)
  stat.lse = sqrt(sum(Ylag1[,3]^2))*(B0.lse - B0)
  stat.mle = sqrt(sum(Ylag1[,3]^2))*(B0.mle - B0)
  rm(X,Y,Z,W,Wlag0,Wlag1,Ylag1,Xlag0, out,out2,out3)
  return(c(stat.lse,stat.mle))
}

library(parallel)
cl = makeCluster(detectCores()-1)

n = 2500 #changer pour 50, 100, 200, 300, 400
clusterExport(cl, 
              varlist = c("sim2", "n", "Phi1", "Phi2", "V02", "sigma.e2",
                          "phix.2", "sigma.b2","genVARX","reducedRank"))

clusterSetRNGStream(cl, iseed = 1)
tic = Sys.time()
nsim = 10000
res = t(parSapply(cl, 1:nsim, function(i) {sim2(i, n)}))
(tac = Sys.time() - tic)
round(colMeans(res),3)
round(colVar(res),3)
round(cov(res[,1],res[,2]),3) 
round(cov(res[,3],res[,4]),3)

#moyennes
#lse                mle
12.158; -1.368; -3.372; -0.003 #n=50
 5.783; -0.721; -3.958;  0.429 #n=100
 2.996; -0.228; -0.579;  0.225 #n=200
 1.950; -0.274; -0.250; -0.001 #n=300
 1.526; -0.213;  0.001; -0.021 #n=400
 0.161; -0.014; -0.178;  0.029 #n=1500 
-0.333;  0.033; -0.579;  0.064 #n=2000
 0.245;  0.033;  0.048;  0.057 #n=2500 
#variances et covariances
#lse                             mle
2903.167; 87.387; -339.161;110376.473; 2870.452;-12912.275 #n=50
1661.526; 54.838; -198.408;  9623.476;  216.442; -1277.851 #n=100
1187.552; 39.787; -135.574;  1234.447;   42.238; -144.020 #n=200
1044.153; 35.359; -120.542;  1051.496;   35.752; -121.148 #n=300
 968.219; 33.608; -111.762;   965.996;   33.843; -111.859 #n=400
 951.257; 32.613; -112.377;   955.702;   32.855; -113.255 #n=500
 919.043; 31.739; -105.858;   923.603;   31.869; -106.393 #n=600
 886.289; 31.275; -105.326;   883.408;   31.240; -104.854 #n=700
 898.671; 30.587; -104.735;   899.737;   30.643; -104.936 #n=800
 865.032; 29.537;  -99.101;   865.988;   29.621;  -99.389 #n=900
 881.389; 31.147; -107.248;   881.633;   31.195; -107.399 #n=1000
 841.241; 29.053;  -97.740;   841.698;   29.061;  -97.770 #n=1500
 808.989; 28.835;  -95.401;   809.211;   28.849;  -95.453 #n=2000
 835.448; 28.925;  -97.497;   835.740;   28.915;  -97.467 #n=2500
 #vraie matrice de covariance : (t(A)%*%\sigma.e^-1%*%A))^-1

solve(t(A) %*% solve(sigma.e2,A)) 
794.433; -92.936
-92.936;  27.591