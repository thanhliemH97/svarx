library(MTS)
genVARX = function(n, phi, beta, sigma, p, s, phiX, sigmaX, pX){
  
  #Le modèle est dans le sens Y = XB + e plutôt que Y = BX + e
  
  #n : nombre d'observations
  #phi : liste de matrices de coefficients autorégressifs
  #beta : liste de matrices de coefficients régressifs exogènes
  #sigma : matrice de covariance de bruits blancs
  #phiX : array de VAR(pX) pour la variable exogène
  #p : ordre autorégressif maximal
  #s : ordre exogène maximal
  #pX : ordre autorégressif maximal de la variable exogène
  
  #Si avoir une seule matrice phi, changer sa classe pour une liste.
  if (class(phi)[1] == "matrix") phi = list(phi) 
  
  if (class(beta)[1] == "matrix") beta = list(beta)
  
  if (class(phiX)[1] == "numeric") phiX = as.matrix(phiX)
  
  k = ncol(phi[[1]]) # dimension de la série chronologique endogène
  l = ncol(phiX[[1]]) # dimension de la série chronologique exogène
  
  #Matrice covariance identité par défaut
  if (missing(sigma)) sigma = diag(1,k)
  
  if (missing(sigmaX)) sigmaX = diag(1,l)
  
  chauffe = n+1 #temps de chauffe de "n+1" 
  nn = n + chauffe 
  d = max(p,s) #max pour éviter des erreurs de récursion dans la simulation
  w = MASS::mvrnorm(nn, mu=rep(0,k), Sigma=sigma) #bruit blanc vectoriel
  x = MTS::VARMAsim(nn, arlags = 1:pX, phi = phiX, sigma = sigmaX)$series
  x = as.matrix(x)
  y = matrix(0, ncol=k, nrow=nn) #variable endogène initialisée à 0
  
  #Récursion 
  for (t in (d+1):nn) {
    sum = matrix(0, ncol=k, nrow=1)
    #Yt ~ SVAR(I,J)
    for (i in 1:p) sum = sum + y[t-i,] %*% t(phi[[i]])
    for (j in 0:s) sum = sum + x[t-j,] %*% t(beta[[j+1]]) #car ne peut pas
    # avoir beta[[0]]
    
    y[t,] =  sum + w[t,]
  }
  #enlever les "chauffe" premiers nombres aléatoires
  
  y = y[-(1:chauffe), ]
  x = x[-(1:chauffe), ]
  return(list(y = as.matrix(y), x = as.matrix(x)))
}

ols = function(X,Y){ 
  # X: variables explicatives, Y: variables expliquées

  if (class(X)[1] == "numeric") X = matrix(X)
  if (class(Y)[1] == "numeric") Y = matrix(Y)
  param = solve(a=t(X)%*%X, b=t(X)%*%Y)
  resi = W - Z %*% param
  cov = (t(resi) %*% resi) / nrow(Y)
  return(list(parameters=param, covariance=cov))
}

reducedRank = function(C, sigma) {
  A = matrix(C[,1], 2,  1)
  C2 = matrix(C[,2], 2, 1)
  
  part1 = t(A) %*% solve(sigma, A)  
  part2 = t(A) %*% solve(sigma, C2)
  
  B0 = solve(part1, part2)
  
  return(list(A=A, B0=B0))
}

colSD = function(matrix){ #écart-type de chaque colonne d'une matrice
  apply(matrix,2,function(y) sd(y)) #"2" indique que la fonction s'applique aux colonnes
}

newtonRaphson1 = function(par){
  
  condition = Inf; i = 1
  
  while(condition>1e-10){
    
    A = matrix(par[1:2],2,1) 
    B = matrix(c(1,par[3]),1,2)
    V0 = matrix(par[4:5],2,1)
    
    Epsilon =  W - Ylag1 %*% t(A%*%B) - Xtrunc %*% t(V0)
    cov = (t(Epsilon) %*% Epsilon) / nrow(W)
    
    sum1  = matrix(0, 5, 5); sum2 = matrix(0, 5, 1)
    
    for (t in 2:n) {
      
      bloc1 = Y[t-1,2]  %x% A
      Ut = cbind(t(B%*%Y[t-1,]), t(X[t,]))
      bloc2 = Ut %x% diag(1,2)
      Zstar_t = cbind(bloc1,bloc2) 
      sum1 = sum1 + t(Zstar_t) %*% solve(cov, Zstar_t)
      sum2 = sum2 + t(Zstar_t) %*% solve(cov, Epsilon[t-1,])
    }
    remain = solve(sum1, sum2)
    newpar = par + remain[c(2,3,1,4,5),]
    par = newpar
    condition = sqrt(sum(remain^2))
    i = i+1
  }
  print(i)
  return(par)
}