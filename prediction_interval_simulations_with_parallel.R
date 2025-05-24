Phi1 = matrix(c(.5, .15, -.4, .5), nrow=2, byrow=T)
Phi2 = matrix(c(.2, .1, -.2, .25), 2, byrow=T)
phi1 = matrix(c(1,2),2,1)
Beta1 = matrix(c(.2, -.3, .6, 1), 2, byrow=T) 
Beta2 = matrix(c(2.4, .6, -.9, 3.1), 2, byrow=T) 
bigBeta = rbind(cbind(Beta1,0,0, Beta2), 0,0,0,0)
rho = .8
Sigma1 = diag((1-rho),2) + rho*matrix(1,2,2)
Phi1.X = matrix(c(.2, .3, -.6, 1), nrow=2, byrow=T)
Sigma1.X = matrix(c(4,.8,.8,1), 2)
listPHI = list(Phi1,Phi2)
listBETA = list(Beta1,Beta2)

source("~/these_doctorat/Simulations/fctSVARX.R")

sim1c = function(sim,n,hmax,h) {
  out = genSVARX(n+hmax, listPHI, listBETA, Sigma1, c(1, 3), c(0, 2), Phi1.X, 
                 Sigma1.X, matrix(c(0,0),2,1))
  
  trainY = out$y[1:n,]
  trainX = out$x[1:n,]
  testY = out$y[(n+1):(n+hmax),]
  testX = out$x[(n+1):(n+hmax),]
  
  regsub = SVARX(trainY,trainX,c(1,3),c(0,2), intercept=F, print=F)

  pred1 = SVARX.pred(regsub, testX, h, cf = F, print = F)
  pred2 = SVARX.pred(regsub, testX, h, cf = T, print = F)
  
  ypred = pred1$pred
  
  #longueur avec vs sans CF
  eqmp1 = pred1$mspe
  eqmp2 = pred2$mspe
  
  qn = qnorm(0.975)
  
  se_est_1 = sqrt(eqmp1[1,1])
  se_est_2 = sqrt(eqmp1[2,2])
  se_est_moy = sqrt(c(1/2,1/2) %*% eqmp1 %*% c(1/2,1/2))
  
  true_1 = abs(testY[h,1]-ypred[1]) <= qn*se_est_1
  true_2 = abs(testY[h,2]-ypred[2]) <= qn*se_est_2  
  true_mean = abs(mean(testY[h,])-mean(ypred)) <= qn*se_est_moy
  
  se_est_1_c = sqrt(eqmp2[1,1])
  se_est_2_c = sqrt(eqmp2[2,2])
  se_est_moy_c = sqrt(c(1/2,1/2) %*% eqmp2 %*% c(1/2,1/2))
  
  true_1_c = abs(testY[h,1]-ypred[1]) <= qn*se_est_1_c
  true_2_c = abs(testY[h,2]-ypred[2]) <= qn*se_est_2_c  
  true_mean_c = abs(mean(testY[h,])-mean(ypred)) <= qn*se_est_moy_c
  
  res = c(true_1, true_2, true_mean, true_1_c, true_2_c, true_mean_c)
  
  return(res)
}

library(parallel)
n = 480; hmax = 12; h = 1; nsim = 5000
cl = makeCluster(detectCores()-1)
clusterEvalQ(cl, library(MASS, Matrix))
clusterExport(cl, 
              varlist = c("genSVARX", "n", "hmax", "h",
                          "listPHI", "listBETA", "Sigma1", "Phi1.X",
                          "Sigma1.X", "phi1","subsets.lags",
                          "subsets.lags.exo","SVARX","SVARX.pred",
                          "aic.svarx","sim1c"))

clusterSetRNGStream(cl, iseed = 1)
tic = Sys.time()
couverture = t(parSapply(cl, 1:nsim, function(i) {sim1c(i, n, hmax, h)}))
(tac = Sys.time() - tic)
colMeans(couverture)
#sink(file = "SVARX-DPG1-pred.txt", append = T)
#cat("Expérience 3 : SVARX({1,3},{0,2})","\n")
#cat("Nombre de simulations :", nsim, "\n")
#cat("Nombre d'observations :", n, "\n")
#cat("Horizon :", h, "\n")
#tac
#cat("Taux nominal : 0,95","\n")
#cat("Taux de couverture :","n")
#cat("Yn(1), Yn(2), 1/2*(Yn(1)+Yn(2)) :","\n")
#colMeans(couverture)
#cat("\n")
#sink()
