VARXorder = function(y, x, maxp = 13, maxs = 3, output = T){
  
  y = as.matrix(y)
  x = as.matrix(x)
  n = dim(y)[1]
  k = dim(y)[2]
  n_x = dim(x)[1]
  m = dim(x)[2]
  if (maxp < 1) maxp = 1
  
  if (n_x != n) {
    cat("Adjustment made for different nobs:", c(n, n_x), "\n")
    n = min(n, n_x)
  }
  
  aic = matrix(0, maxp + 1, maxs + 1)
  rownames(aic) = paste0("p=", 0:maxp)
  colnames(aic) = paste0("s=", 0:maxs)
  bic = aic
  hqc = aic
  
  for (s in 0:maxs) {
    
    y_trunc = y[(s+1):n, ]
    n_prime = n - s
    z = rep(1, n_prime)
    
    for (j in 0:s) {
      z = cbind(z, x[(s+1 - j):(n - j), ])
    }
    
    ztz = t(z) %*% z
    zty = t(z) %*% y_trunc
    beta = solve(ztz, zty)
    resi = y_trunc - z %*% beta
    
    mse = t(resi) %*% resi / n_prime
    ln_ds = log(det(mse))
    npar = k * m * (s+1)
    
    aic[1, s + 1] = ln_ds + 2 * npar / n_prime
    bic[1, s + 1] = ln_ds + log(n_prime) * npar / n_prime
    hqc[1, s + 1] = ln_ds + 2 * log(log(n_prime)) * npar / n_prime
    
    for (p in 1:maxp) {
      
      d = max(p,s)
      n_prime = n - d
      y_trunc = y[(d+1):n, ]
      z = rep(1, n_prime)
      
      for (i in 1:p) {
        z = cbind(z, y[(d+1 - i):(n - i), ])
      }
      for (j in 0:s) {
        z = cbind(z, x[(d+1 - j):(n - j), ])
      }
      
      ztz = t(z) %*% z
      zty = t(z) %*% y_trunc
      beta = solve(ztz, zty)
      resi = y_trunc - z %*% beta
      
      mse = t(resi) %*% resi / n_prime
      ln_ds = log(det(mse))
      npar2 = k*k * p + k*m * (s+1)
      
      aic[p + 1, s + 1] = ln_ds + 2 * npar2 / n_prime
      bic[p + 1, s + 1] = ln_ds + log(n_prime) * npar2 / n_prime
      hqc[p + 1, s + 1] = ln_ds + 2 * log(log(n_prime)) * npar2 / n_prime
    }
  }
  ind.min = function(A) {
    index = which(A == min(A), arr.ind = TRUE)
    p = index[1,1]
    s = index[1,2]
    return(c(p,s))
  }
  aicor = ind.min(aic) - 1
  bicor = ind.min(bic) - 1
  hqcor = ind.min(hqc) - 1
  
  if (output) {
    cat("selected order(p,s): aic = ", aicor, "\n")
    cat("selected order(p,s): bic = ", bicor, "\n")
    cat("selected order(p,s): hqc = ", hqcor, "\n")
  }
  VARXorder = list(aic = aic, aicor = aicor, bic = bic, bicor = bicor, 
                   hqc = hqc, hqcor = hqcor)
}