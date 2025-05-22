library(MASS) #plus rapide que library(mvtnorm)
library(MTS)
library(Matrix)
########################################################
###### Générer une série chronologique SVARX     #######
########################################################

genSVARX = function(n,phi,beta,sigma,lagY,lagX,phiX,sigmaX,phi0){
  
  #Le modèle est dans le sens Y = XB + e plutôt que Y = BX + e
  
  #n : nombre d'observations
  #phi : liste de matrices de coefficients autorégressifs
  #beta : liste de matrices de coefficients régressifs exogènes
  #sigma : matrice de covariance de bruits blancs
  #phi0 : vecteur colonne constant 
  #phiX : matrice de VAR(1) pour la variable exogène
  #lagY : vecteur de délais endogènes
  #lagX : vecteur de délais exogènes
  
  #Si avoir une seule matrice phi, changer sa classe pour une liste.
  if (class(phi)[1] == "matrix") phi = list(phi) 
  
  if (class(beta)[1] == "matrix") beta = list(beta)

  k = ncol(phi[[1]]) # dimension de la série chronologique endogène
  #l = ncol(beta[[1]]) # dimension de la série chronologique exogène
  
  #Matrice covariance identité par défaut
  if (missing(sigma)) sigma = diag(1,k)
  
  if (missing(phi0)) phi0 = matrix(0, nrow = k, ncol = 1)
  
  chauffe = n+1 #temps de chauffe de "n+1" 
  nn = n + chauffe 
  p = max(lagY) #ordre autorégressif 
  q = max(lagX) #ordre autorégressif exogène (il y a un délai 0)
  d = max(p,q) #max pour éviter des erreurs de récursion dans la simulation
  w = MASS::mvrnorm(nn, mu = rep(0,k), Sigma = sigma) #Wt ~ bruit blanc vectoriel
  x = MTS::VARMAsim(nn, arlags = 1, phi = phiX, sigma = sigmaX)$series # Xt ~ VAR(1)
  y = matrix(0, ncol = k, nrow = nn) #variable endogène initialisée à 0
  
  #Récursion Yt ~ SVAR(I,J)
  for (t in (d+1):nn) {
    sum = matrix(0, ncol = k, nrow = 1)
    
    for (i in 1:length(lagY)) sum = sum + y[t-lagY[i], ] %*% t(phi[[i]])
    for (j in 1:length(lagX)) sum = sum + x[t-lagX[j], ] %*% t(beta[[j]])

    y[t,] =  sum + w[t,] + t(phi0)
  }
  
  #Enlever les "chauffe" premiers nombres aléatoires
  y = y[-(1:chauffe), ]
  x = x[-(1:chauffe), ]
  return(list(y = as.matrix(y), x = as.matrix(x)))
}

#Sélection de p et s

VARXorder2 = function(y, x, max_p = 13, max_s = 3, output = T) {
  
  y = as.matrix(y)
  x = as.matrix(x)
  n = dim(y)[1]
  k = dim(y)[2]
  n_x = dim(x)[1]
  l = dim(x)[2]
  if (max_p < 1) max_p = 1
  
  if (n_x != n) {
    cat("Adjustment made for different nobs:", c(n, n_x), "\n")
    n = min(n, n_x)
  }

  aic = matrix(0, max_p + 1, max_s + 1)
  rownames(aic) = paste0("p=", 0:max_p)
  colnames(aic) = paste0("s=", 0:max_s)
  bic = aic
  hqc = aic
  
  for (s in 0:max_s) {
    
    n_prime_x = n - s
    y_trunc = y[(s+1):n, ]
    z = rep(1, n_prime_x)

    for (j in 0:s) z = cbind(z, x[(s+1 - j):(n - j), ])
      
    ztz = t(z) %*% z
    zty = t(z) %*% y_trunc
    beta = solve(ztz, zty)
    resi = y_trunc - z %*% beta
    sigma = t(resi) %*% resi / n_prime_x
    ln_ds = log(det(sigma))
      
    npar_x = k*l*(s+1)
      
    aic[1, s + 1] = ln_ds + 2 * npar_x / n_prime_x
    bic[1, s + 1] = ln_ds + log(n_prime_x) * npar_x / n_prime_x
    hqc[1, s + 1] = ln_ds + 2 * log(log(n_prime_x)) * npar_x / n_prime_x
    
    for (p in 1:max_p) {
      
      d = max(p, s)
      n_prime  = n - d
      y_trunc = y[(d+1):n, ]
      z = rep(1, n_prime)
      
      for (i in 1:p) z = cbind(z, y[(d+1 - i):(n - i), ])
      
      for (j in 0:s) z = cbind(z, x[(d+1 - j):(n - j), ])
      
      ztz = t(z) %*% z
      zty = t(z) %*% y_trunc
      beta = solve(ztz, zty)
      resi = y_trunc - z %*% beta
      sigma = t(resi) %*% resi / n_prime
      ln_ds = log(det(sigma))
      
      npar = k*k*p + k*l*(s+1)
      
      aic[p + 1, s + 1] = ln_ds + 2 * npar / n_prime
      bic[p + 1, s + 1] = ln_ds + log(n_prime) * npar / n_prime
      hqc[p + 1, s + 1] = ln_ds + 2 * log(log(n_prime)) * npar / n_prime
    }
  }
  ind.min = function(A) {
    index = which(A == min(A), arr.ind = TRUE)
    p = index[1,1]
    s = index[1,2]
    return(c(p,s))
  }
  
  aic_order = ind.min(aic) - 1
  bic_order = ind.min(bic) - 1
  hqc_order = ind.min(hqc) - 1
  
  if (output) {
    cat("selected order(p,s): aic = ", aic_order, "\n")
    cat("selected order(p,s): hqc = ", hqc_order, "\n")
    cat("selected order(p,s): bic = ", bic_order, "\n")
  }
  VARXorder2 = list(aic = aic, aicor = aic_order,  
                    hqc = hqc, hqcor = hqc_order,
                    bic = bic, bicor = bic_order)
}

#Sélection d'indices

subsets.lags = function(p) {
  set = list()
  for (i in 0:(p-1)) {
    set = c(set, combn(p-1, i, function(x) c(x,p), simplify = F))
  }
  return(set)
}

subsets.lags.exo = function(s) {
  set = list()
  for (j in 0:s) {
    set = c(set,combn(0:(s-1), j, function(x) c(x,s), simplify = F))
  }
  return(set)
}

aic.svarx = function(y, x, p, s, crit, intercept=T){
  
  if (class(x)[1] == "numeric") x = matrix(x)
  if (class(y)[1] == "numeric") y = matrix(y)
  n = nrow(y)
  k = ncol(y)
  l = ncol(x)
  d = max(p,s)
  n_prime = n-d
  y_trunc = y[(d+1):n, ]
  
  if (missing(crit)) crit = 2   #par défaut : AIC
  if (crit == "bic") crit = log(n) #BIC
  if (crit == "hqc") crit = 2*log(log(n)) #HQC : Hannan-Quinn
  
  subsetI = subsets.lags(p)
  subsetJ = subsets.lags.exo(s)
  nSubsetI = length(subsetI)
  nSubsetJ = length(subsetJ)
  minValue = Inf
  bestI = NULL
  bestJ = NULL
  
  for (rj in 1:nSubsetJ){
    for (ri in 1:nSubsetI) {
      
      if (intercept) {
        z = rep(1, n_prime)
      } else {
        z = NULL
      }
      for (i in subsetI[[ri]]) z = cbind(z, y[(d+1 - i):(n - i), ])
      
      for (j in subsetJ[[rj]]) z = cbind(z, x[(d+1 - j):(n - j), ])
      
      ztz = t(z) %*% z
      zty = t(z) %*% y_trunc
      beta = solve(ztz, zty)
      resi = y_trunc - z %*% beta
      sigma = (t(resi) %*% resi) / n_prime
      ln_ds = log(det(sigma))
      npar = k * (k*length(subsetI[[ri]]) + l*length(subsetJ[[rj]]))
      
      value = ln_ds + (crit * npar) / n_prime
      if (value < minValue) {
        minValue = value
        bestI = ri
        bestJ = rj
      }
    }
  }
  I = subsetI[[bestI]]; J = subsetJ[[bestJ]]
  IJ = paste0("{", paste(I, collapse = ","), "}",",{",paste(J,collapse = ","),"}")
  subsetIJ = list(I = I, J = J, min = minValue, IJ = IJ)
  return(subsetIJ)
}

SVARX = function(y, x, lagY, lagX, intercept=T, print=T){
  if (class(y)[1] == "numeric") y = matrix(y)
  if (class(x)[1] == "numeric") x = matrix(x)
  
  p = max(lagY)
  q = max(lagX)
  n = nrow(y)
  d = max(p,q)
  k = ncol(y)
  l = ncol(x)
  a = length(lagY)
  b = length(lagX)
  n_prime  = n-d
  Y = y[(d+1):n, ]
  
  if (intercept) {Z = rep(1, n_prime)}
  else {Z = NULL}
  
  for (i in lagY) Z = cbind(Z, y[(d+1 - i):(n - i), ])
  
  for (j in lagX) Z = cbind(Z, x[(d+1 - j):(n - j), ])
  
  ztz = t(Z) %*% Z
  zty = t(Z) %*% Y
  beta = solve(ztz, zty)
  resi = Y - Z %*% beta
  npar = k * (k*a + l*b)
  
  sigma = t(resi) %*% resi / n_prime
  ln_ds = log(det(sigma))

  
  loglik = ln_ds
  aic = ln_ds + (2 * npar) / n_prime
  bic = ln_ds + (log(n_prime) * npar) / n_prime
  hqc = ln_ds + (2 * log(log(n_prime)) * npar) / n_prime
  
  covBeta = solve(ztz) %x% sigma
  se.Beta = matrix(sqrt(diag(covBeta)),nrow=k,byrow=F)
  
  if (intercept) {
    intercepts = matrix(beta[1,],nrow=k)
    remainder = t(beta[-1,])
  } else {
    intercepts = matrix(0,nrow=k)
    remainder = t(beta)
  }
  
  colnames(intercepts) = "Intercept"
  rownames(intercepts) = paste0("y", 1:k)

  ar = remainder[, 1:(k*a)]
  colnames(ar) = paste0("AR", rep(lagY,each=k))
  rownames(ar) = paste0("y", 1:k)
  exo = remainder[, -(1:(k*a))]
  colnames(exo) = paste0("Exo", rep(lagX,each=l))
  rownames(exo) = paste0("y", 1:k)
  coef = list(intercept = intercepts, ar=ar, exo=exo)
  lag = list(y=lagY, x=lagX)
  
  if (print) print(cbind(intercepts,ar,exo))
  
  SVARX = list(coefficients=coef,covariance=sigma,loglik=loglik,aic=aic,
               bic=bic,hqc=hqc,covCoef=covBeta,stdCoef=se.Beta,resid = resi,
               y=y,x=x,lag=lag)
}

#Prévision

SVARX.pred = function(model, new.x, h, level=0.95, cf = T, print=T){
  y = model$y
  x = rbind(model$x, new.x)
  n = nrow(y)
  p = max(model$lag$y)
  s = max(model$lag$x)
  k = ncol(y)
  l = ncol(x)
  lagY = model$lag$y
  lagX = model$lag$x
  qn = qnorm(1-(1-level)/2)
  nh = n+h
  
  #data
  yStar =  function(y,p){
    n = nrow(y)
    return(matrix( c(t(y[n:(n-p+1),])), ncol = 1 ))
  }
  xStarLagged = function(x,j,s){
    nx = nrow(x)
    return(matrix(c(t(x[(nx-j):(nx-s-j), ])), ncol = 1))
  }
  bigW = function(h, outputY, inputX, p, s, intercept=T) {
    n = nrow(outputY)
    k = ncol(outputY)
    truncX = inputX[1:(n+h),]
    if (intercept) {
      tempo = cbind( matrix(1, ncol = h), t(yStar(outputY, p)) )
    } else {
      tempo = t(yStar(outputY, p))
    }
    
    for (j in 0:(h-1)) tempo = cbind( tempo, t(xStarLagged(truncX, j, s)) )
    tempo %x% diag(k)
  }
  C.I = function(lagY,k) {
    p = max(lagY)
    output = NULL
    for (i in lagY){
      vecteur = rep(0, p)
      vecteur[i] = 1
      Ei = kronecker(vecteur, diag(k))
      output = rbind(output, t(Ei))
    }
    rownames(output) = paste0("AR", rep(lagY, each = k))
    colnames(output) = paste0("AR", rep(1:p, each = k))
    output
  }
  D.J = function(lagX,l) {
    s = max(lagX)
    output = NULL
    
    for (j in lagX){
      vecteur = rep(0, s+1)
      vecteur[j+1] = 1
      Fj = kronecker(vecteur, diag(l))
      output = rbind(output, t(Fj))
    }
    rownames(output) = paste0("Exo", rep(lagX, each = l))
    colnames(output) = paste0("Exo", rep(0:s, each = l))
    output
  }
  pow = function(A, n) {
    if (!is.matrix(A)) return(A**n)
    else {
      output = diag(1, nrow(A))
      
      if(n == 0) return(output)
      if(n == 1) return(A)
      
      while (n > 0) {
        if (n %% 2 == 1) output = A %*% output
        A = A %*% A
        n = n %/% 2
      }
      return(output)
    }
  }
  
  #point forecast 
  
  phi0 = model$coefficients$intercept
  phi0Star = rbind(phi0, matrix(0, k*(p-1), 1))
  
  PhiStar = rbind(model$coefficients$ar %*% C.I(lagY, k),
                  cbind(diag(1, k*(p-1)), matrix(0, k*(p-1), k)))
  
  Beta2Star = model$coefficients$exo %*% D.J(lagX, l)  
  BetaStar = rbind(Beta2Star, matrix(0, k*(p-1), l*(s+1)))
  
  PsiStar = function(A, n, k) (pow(A, n)) [1:k, 1:k]
  
  if (any(phi0)==0){
    add_0 = 0
  } else {
  add_0 = Reduce("+", lapply(0:(h-1), function(j) PsiStar(PhiStar, j, k))) %*% phi0
  }
  add_x = Reduce("+", lapply(0:(h-1), function(j) 
    PsiStar(PhiStar, j, k) %*% Beta2Star %*% xStarLagged( x[1:nh, ], j,s) ))
  
  add_y = matrix( (pow(PhiStar, h) %*% yStar(y, p))[1:k, ], ncol = 1 ) 
  
  ypred = add_0 + add_y + add_x
  colnames(ypred) = paste0("hor",h)
  
  #conditionnal mean square forecast error
  
  bigH = function(h, phi0Star, PhiStar, BetaStar, lagY, lagX, k, l, p, s){
    CI = C.I(lagY, k)
    DJ = D.J(lagX, k)
    
    #i = 1,...., h-1
    blocB = function(phi0Star, PhiStar, i, k) {
      result = Reduce("+", lapply(0:(i-1), function(j) 
        t(CI %*% pow(PhiStar, i-j-1) %*% phi0Star) %x% PsiStar(PhiStar, j, k) ))
      return(result)
    }
    
    #i = h
    blocC = function(PhiStar, i, k) {
      result = Reduce("+", lapply(0:(i-1), function(j) 
        t(CI %*% pow(PhiStar, i-j-1)) %x% PsiStar(PhiStar, j, k) ))
      return(result)
    }
    
    #i = 1,...,h-1
    blocD = function(PhiStar, BetaStar, i, k){
      result = Reduce("+", lapply(0:(i-1), function(j) 
        ((t(CI %*% pow(PhiStar, i-j-1) %*% BetaStar)) %x% PsiStar(PhiStar, j, k)) ))
      return(result)
    }
    
    #i = 0,...,h-1
    blocE = function(PhiStar, i, k) t(DJ) %x% PsiStar(PhiStar, i, k)
    
    
    if (h==1){
      tempoA = diag(1, k)
      tempoC = t(CI) %x% diag(1, k)
      tempoE = t(DJ) %x% diag(1, l)
      
      M = ifelse(any(phi0Star)==0,bdiag(tempoC, tempoE),
                 bdiag(tempoA, tempoC, tempoE))
      return(as.matrix(M))
    } else {
      tempoC = blocC(PhiStar, h, k)
      
      tempoD = NULL
      for (i in 1:(h-1)){
        tempoD = rbind(tempoD, blocD(PhiStar, BetaStar, i, k))
      }
      
      tempoE = NULL
      for (i in 0:(h-1)){
        tempoE = rbind(tempoE, blocE(PhiStar, i, k))
      }
      
      if (any(phi0Star)==0){

        nrow1 = nrow(blocE(PhiStar, 0, k))
        nrow2 = nrow(tempoC)
        
        ncol1 = ncol(tempoC)
        ncol2 = ncol(tempoE)
        
        H = cbind(rbind(tempoC, matrix(0, nrow1, ncol1), tempoD),
                  rbind(matrix(0, nrow2, ncol2), tempoE))
      } else {
        tempoA = diag(1,k)
        for (i in 1:(h-1)){
          tempoA = rbind(tempoA, PsiStar(PhiStar, i, k))
        }
        
        tempoB = NULL
        for (i in 1:(h-1)){
          tempoB = rbind(tempoB, blocB(phi0Star, PhiStar, i, k))
        }
        nrow1 = nrow(tempoC) + nrow(tempoE)
        nrow2 = length(phi0)
        nrow3 = nrow(blocE(PhiStar, 0, k))
        nrow4 = nrow(tempoA) + nrow(tempoC)
        
        ncol1 = ncol(tempoA)
        ncol2 = ncol(tempoB)
        ncol3 = ncol(tempoE)
        
        H = cbind(rbind(tempoA, matrix(0, nrow1, ncol1)),
                  rbind(matrix(0, nrow2, ncol2), tempoB, tempoC,
                        matrix(0, nrow3, ncol2), tempoD),
                  rbind(matrix(0, nrow4, ncol3), tempoE))
      }
      return(H)
    }          
  }
  
  sigma = model$covariance
  varBeta = model$covCoef
  
  eqmp = Reduce("+", lapply(0:(h-1), function(j) 
    PsiStar(PhiStar, j, k) %*% sigma %*% t(PsiStar(PhiStar, j, k))))
  
  
  if (any(phi0)==0){
    W_nh = bigW(h, y, x, p, s, F)
  } else {
    W_nh = bigW(h, y, x, p, s)
  }
  
  H_IJ = bigH(h, phi0Star, PhiStar, BetaStar, lagY, lagX, k, l, p, s)
  
  correcteur = W_nh %*% H_IJ %*% varBeta %*% t(H_IJ) %*% t(W_nh)
  
  if (cf==T) eqmp = eqmp + correcteur
  
  IP = matrix(0,nrow=k,ncol=2)
  
  for (i in 1:k) IP[i,] = ypred[i,] + c(-1,1) * qn * sqrt(eqmp[i,i])
  
  colnames(IP) = c("L95 PI", "U95 PI")
  
  if (print==T) print(cbind(ypred,IP))
  
  SVARX.pred = list(pred=ypred, ip=IP, eqmp=eqmp)
}


#SVAR

VARorder2 = function(y, max_p = 13, output = T) {
  y = as.matrix(y)
  n = dim(y)[1]
  k = dim(y)[2]
  d = max_p
  
  n_prime = n - d
  y_trunc = y[(d+1):n, ]
  
  aic = rep(0, max_p)
  names(aic) = paste0("p=", 1:max_p)
  bic = aic
  hqc = aic
  
  for (p in 1:max_p) {
    
    z = rep(1, n_prime)
    
    for (i in 1:p) z = cbind(z, y[(d+1 - i):(n - i), ])
    
    ztz = t(z) %*% z
    zty = t(z) %*% y_trunc
    beta = solve(ztz, zty)
    resi = y_trunc - z %*% beta
    sigma = (t(resi) %*% resi) / n_prime
    ln_ds = log(det(sigma))
    npar = k * (p*k + 1)
    
    aic[p] = ln_ds + (2 * npar) / n_prime
    bic[p] = ln_ds + (log(n_prime) * npar) / n_prime
    hqc[p] = ln_ds + (2 * log(log(n_prime)) * npar) / n_prime
  }
  
  aic_order = which.min(aic)
  bic_order = which.min(bic)
  hqc_order = which.min(hqc)
  
  if (output) {
    cat("selected order(p): aic = ", aic_order, "\n")
    cat("selected order(p): hqc = ", hqc_order, "\n")
    cat("selected order(p): bic = ", bic_order, "\n")
  }
  VARorder2 = list(aic = aic, aicor = aic_order,  
                   hqc = hqc, hqcor = hqc_order,
                   bic = bic, bicor = bic_order)
}

aic.svar = function(y, p, crit, intercept=T){
  
  if (class(x)[1] == "numeric") x = matrix(x)
  n = nrow(y)
  k = ncol(y)
  d = p
  n_prime = n-d
  y_trunc = y[(d+1):n, ]
  
  if (missing(crit)) crit = 2   #par défaut : AIC
  if (crit == "bic") crit = log(n) #BIC
  if (crit == "hqc") crit = 2*log(log(n)) #HQC : Hannan-Quinn
  
  subsetI = subsets.lags(p)
  nSubsetI = length(subsetI)
  minValue = Inf
  bestI = NULL
  
  for (ri in 1:nSubsetI) {

    z = ifelse(intercept, rep(1, n_prime), NULL)
      
    for (i in subsetI[[ri]]) z = cbind(z, y[(d+1 - i):(n - i), ])
      
    ztz = t(z) %*% z
    zty = t(z) %*% y_trunc
    beta = solve(ztz, zty)
    resi = y_trunc - z %*% beta
    sigma = (t(resi) %*% resi) / n_prime
    ln_ds = log(det(sigma))
    npar = k*k*length(subsetI[[ri]])
      
    value = ln_ds + (crit * npar) / n_prime
    
    if (value < minValue) {
      minValue = value
      bestI = ri
    }
  }
  subsetI = list(I = subsetI[[bestI]], min = minValue)
  return(subsetI)
}


#test diagnostic Duchesne-Roy
qs <- function(x) {
  cte <- sqrt(5/3)
  ind0 <- (x == 0)
  ind1 <- (x != 0)
  out <- rep(1, length(x))
  out[ind0] <- 1
  out[ind1] <- 9/(5 * x[ind1]^2 * pi^2) * (sin(cte * pi * x[ind1])/(cte * pi * x[ind1]) - cos(cte * pi * x[ind1]))
  out
}
daniell <- function(x) {
  ind0 <- (x == 0)
  ind1 <- (x != 0)
  out <- x
  out[ind0] <- 1
  out[ind1] <- sin(pi * x[ind1])/(pi * x[ind1])
  out
}

parzen <- function(x) {
  out <- rep(0, length(x))
  A <- (abs(x) > 6/pi)
  B <- (abs(x) > 3/pi) & (abs(x) <= 6/pi)
  C <- (abs(x) <= 3/pi)
  out[A] <- 0
  out[B] <- 2 * (1 - abs((pi * x[B])/6))^3
  out[C] <- 1 - 6 * ((pi * x[C])/6)^2 + 6 * abs((pi * x[C])/6)^3
  out
}


bartlett <- function(x) {
  out <- rep(0, length(x) )
  A <- (abs(x) <= 1)
  B <- (abs(x) > 1)
  out[A] <- 1 - abs(x[A])
  out[B] <- 0
  out
}


tronque <- function(x) {
  out <- rep(0, length(x))
  A <- (abs(x) <= 1)
  B <- (abs(x) > 1)
  out[A] <- 1
  out[B] <- 0
  out
}
##fonction de calcul de l'autocovariance d'ordre k###
acf.lagk.mult <- function(k,X){
  # la matrice des donnees X est de dimension n x d
  n <- nrow(X)
  d <- ncol(X)
  meanX <- apply(X,2,mean)
  output <- matrix(0, nrow=d, ncol=d)
  for (t in (k+1): n){
    mat <- (X[t,] - meanX) %*% t(X[t-k, ] - meanX)
    output <- output+mat/n      
  }
  output
}

##########la fonction de calcul de la statistique#########
trmat <- function(A) {
  sum(diag(A))
}
statTn <- function(res, noyau, Pn){
  n <- nrow(res)
  nm1 <- n-1
  d <- ncol(res)
  Mn <- 0
  Vn <- 0
  Qn <- 0
  sigmam1 <- solve(acf.lagk.mult(0,res))
  for (j in 1:nm1){
    acf.lagj <- acf.lagk.mult(j, res)
    mat <- t( acf.lagj ) %*%  sigmam1  %*%  acf.lagj  %*%  sigmam1
    Mn <- Mn + (1-j/n)*( noyau(j/Pn) )^2
    Vn <- Vn + (1- j/n)*(1-(j+1)/n)*( noyau(j/Pn) )^4
    Qn <- Qn+trmat(mat)*( noyau(j/Pn) )^2
  }
  Tn <- ( n*Qn-(d^2)*Mn )/sqrt(2*d*d*Vn)
  return(Tn)
}

########################## Menage ########################
# rm(list=ls())

