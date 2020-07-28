if (!require("expm")) install.packages("expm")
library(expm)
if (!require("metafor")) install.packages("metafor")
library(metafor)
if (!require("tsDyn")) install.packages("tsDyn")
library(tsDyn)
if (!require("vars")) install.packages("vars")
library(vars)
#
set.seed(123)

q = 2 # Number of variables
Phi_pop <- matrix(c(0.50, 0.15,
                    0.25, 0.40), nrow=q, byrow=TRUE) # population Phi matrix (i.e., lagged effects matrix)
A_pop <- logm(Phi_pop) # underlying population drift matrix
Phi_pop <- expm(A_pop)
vecPhi_pop <- as.vector(t(Phi_pop)) # vectorized population Phi
Gamma_pop <- matrix(c(1,    0.3,
                      0.3,    1), nrow=q, byrow=TRUE) # population stationary covariance matrix
# Since Gamma is now a correlation matrix, Phi is already standardized
SigmaVAR_pop <- Gamma_pop - Phi_pop %*% Gamma_pop %*% t(Phi_pop) # population residual covariance matrix

#######################################

# Used from Work engagement & burnout paper
N <- matrix(c(
  643,
  651,
  473,
  387,
  187,
  209,
  2897,
  160,
  1964,
  848,
  926,
  274,
  433,
  256,
  409,
  926,
  162,
  262,
  247,
  102,
  171,
  201,
  309,
  77,
  67), byrow=TRUE)

S <- length(N)

var1 <- matrix(NA, ncol=(q^2), nrow=S)
for(s in 1:S){
  var1[s,] <- rep(1/(N[s] - 3), q^2)
}

# Time intervals ('time lags'), using Work engagement & burnout paper
TI <- matrix(c(
  12,
  12,
  12,
  4,
  9,
  12,
  12,
  2,
  36,
  12,
  12,
  8,
  24,
  12*(1/365),
  24,
  12,
  2,
  14,
  12,
  10,
  48,
  12,
  12,
  1,
  8), byrow=TRUE)
TI = TI/12 # now in years


# Needed in plot later on
Step = 0.1
Min = 0 #min(TI) 
Max = max(TI)
TIs<-seq(Min,Max,by=Step)
PhiTIs<-array(data=NA,dim=c(q,q,length(TIs)))
for(i in 1:length(TIs)){
  PhiTIs[,,i]<-expm(A*TIs[i])
}


## Make dummy variables ##
repTI <- matrix(rep(TI, each = q^2), nrow = (q^2 * dim(TI)[1]), ncol=1)
if (!require("fastDummies")) install.packages("fastDummies")
library(fastDummies)
G <- length(unique(TI))
D <- dummy_cols(TI)
repD <- dummy_cols(repTI)
# First column contain data itself, so take that out: 
if(dim(D)[2] > G){
  D <- D[,-1]
  repD <- repD[,-1]
}


## Determine population values for the G = 12 unique time intervals ##
vecPhi_G_pop <- matrix(NA, G, q^2)
for(g in 1:G){
  vecPhi_G_pop[g,] <- as.vector(t(expm(A_pop * unique(TI)[g])))
}


#######################################################################################################################################

## Example ##

# Sample - for each study s in the meta-analysis - data based on population values, N_s, and DeltaT_s
S <- length(N)
vecPhi <- array(NA, dim = c(S*q*q))
CovMx <- array(0, dim = c(S*q*q, S*q*q))
vecPhi1 <- array(NA, dim = c(S*q*q))
CovMx1 <- array(0, dim = c(S*q*q, S*q*q))
#
G <- length(unique(TI))
vecPhi_G <- array(NA, dim = c(S*q*q, G))
CovMx_G <- array(0, dim = c(S*q*q, S*q*q, G))
#
s <- 1 # index for number of studies, with in total S studies.
while(s <= S){
  Y_mu <- VAR.sim(Phi_pop, n = N[s], lag = 1, include = c("none"), starting = NULL, varcov = SigmaVAR_pop)
  Y <- scale(Y_mu, scale = FALSE) # substracts the means (so , it centers, does not standardize now because of scale = F)
  outcomeVAR <- VAR(Y, p = 1)
  Phi_VARest <- Acoef(outcomeVAR)[[1]]
  if(any(Re(eigen(Phi_VARest)$values) < 0)){
    s <- s # No CTM-equivalent, so do not proceed
  }else{  
    SigmaVAR1_VARest <- cov(resid(outcomeVAR))
    Gamma_VARest <- cov(Y)
    # Standardize parameters! 
    Sxy <- sqrt(diag(diag(Gamma_VARest)))
    Gamma_VARest <- solve(Sxy) %*% Gamma_VARest %*% solve(Sxy) 
    Phi_VARest <- solve(Sxy) %*% Phi_VARest %*% Sxy
    SigmaVAR1_VARest <- solve(Sxy) %*% SigmaVAR1_VARest %*% solve(Sxy) 
    #
    invGamma <- solve(Gamma_VARest)
    B_VARest <- -logm(Phi_VARest)/1
    #
    Phi_VARest <- expm(-B_VARest*TI[s])
    vecPhi[((s-1)*(q*q)+1):(s*q*q)] <- as.vector(t(Phi_VARest))
    SigmaVAR_VARest <- Gamma_VARest - Phi_VARest %*% Gamma_VARest %*% t(Phi_VARest)
    CovMx[((s-1)*q^2+1):(s*q^2),((s-1)*q^2+1):(s*q^2)] <- kronecker(SigmaVAR_VARest, invGamma) / (N[s]-q)
    if(any( eigen( CovMx[((s-1)*q^2+1):(s*q^2),((s-1)*q^2+1):(s*q^2)] )$values < 0 )){
      s <- s # Cov mx should be pos def
    }else{
      #
      Phi_VARest <- expm(-B_VARest*1)
      vecPhi1[((s-1)*(q*q)+1):(s*q*q)] <- as.vector(t(Phi_VARest))
      SigmaVAR_VARest <- Gamma_VARest - Phi_VARest %*% Gamma_VARest %*% t(Phi_VARest)
      CovMx1[((s-1)*q^2+1):(s*q^2),((s-1)*q^2+1):(s*q^2)] <- kronecker(SigmaVAR_VARest, invGamma) / (N[s]-q)
      if(any( eigen( CovMx1[((s-1)*q^2+1):(s*q^2),((s-1)*q^2+1):(s*q^2)] )$values < 0 )){
        s <- s # Cov mx should be pos def
      }else{
        for(g in 1:G){
          Phi_VARest <- expm(-B_VARest*unique(TI)[g])
          vecPhi_G[((s-1)*(q*q)+1):(s*q*q),g] <- as.vector(t(Phi_VARest))
          SigmaVAR_VARest <- Gamma_VARest - Phi_VARest %*% Gamma_VARest %*% t(Phi_VARest)
          CovMx_G[((s-1)*q^2+1):(s*q^2),((s-1)*q^2+1):(s*q^2),g] = kronecker(SigmaVAR_VARest, invGamma) / (N[s]-q)
        }
        if(any( apply(x <- CovMx_G[((s-1)*q^2+1):(s*q^2),((s-1)*q^2+1):(s*q^2),], 3, function(x){eigen(x)$values}) < 0 )){
          s <- s # Cov mx should be pos def
        }else{
          s <- s+1
        }
      }
    }
  }
}

any(eigen(CovMx1)$val < 0)
# FALSE = pos def; like we want! & TRUE = not pos def!
#
any(eigen(CovMx)$val < 0)


###########################################################################################################


# Multivariate Meta-analyses

# Needed in the GLS/multivariate meta-an
sub = NULL
for(i in 1:q){
  sub = c(sub, paste(i, 1:q, sep=""))
}
outcome <- rep(sub, S) 




# Ignore time interval
metaan <- rma.mv(yi=vecPhi, V=CovMx, mods = ~ outcome - 1, method = "FE")
Phi_N <- coef(metaan)
sePhi_N <- metaan$se 
CovMxPhi_N <- metaan$vb
# elliptical 95%CI 
# Determine points on 95% LL contour
vecPhi_metaan <- Phi_N
CovMx_metaan <- CovMxPhi_N
eigenCovMx <- eigen(CovMx_metaan) 
lambda <- eigenCovMx$val
E <- eigenCovMx$vec
Chi2 <- qchisq(p=0.05, df=(q*q), lower.tail=FALSE)
LB_vecPhi <- matrix(NA, nrow=q*q, ncol =q*q)
UB_vecPhi <- matrix(NA, nrow=q*q, ncol =q*q)
LL <- matrix(NA, nrow=q*q, ncol=2)
teller = 0
for(row in 1:q){
  for(column in 1:q){
    teller = teller + 1
    LB_vecPhi[teller,] <- vecPhi_metaan - sqrt(Chi2 * lambda[teller]) * E[,teller]
    UB_vecPhi[teller,] <- vecPhi_metaan + sqrt(Chi2 * lambda[teller]) * E[,teller]  
    LL[teller,1] <- t(LB_vecPhi[teller,]-vecPhi_metaan) %*% solve(CovMx_metaan) %*% (LB_vecPhi[teller,]-vecPhi_metaan)
    LL[teller,2] <- t(UB_vecPhi[teller,]-vecPhi_metaan) %*% solve(CovMx_metaan) %*% (UB_vecPhi[teller,]-vecPhi_metaan)
    Phi_LB_t <- matrix(LB_vecPhi[teller,], ncol=q, nrow=q)
    Phi_UB_t <- matrix(UB_vecPhi[teller,], ncol=q, nrow=q)
  }
}
# Note that UB can be smaller than LB and LB larger than UB!
minPhi_N <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, min)
maxPhi_N <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, max)
#
# It is actually the estimate for all time intervals, so compare to all!
(minPhi_N < vecPhi_G_pop) + (maxPhi_N > vecPhi_G_pop) - 1 # coverage
# unique(TI) # order of the matrices


# Linear moderator
G <- length(unique(TI))
minPhi_Mod <- matrix(NA, ncol=(q^2), nrow=G)
maxPhi_Mod <- matrix(NA, ncol=(q^2), nrow=G)
Phi_Mod <- matrix(NA, ncol=(q^2), nrow=G)
sePhi_Mod <- matrix(NA, ncol=(q^2), nrow=G)
CovMxPhi_Mod <- array(NA, dim = c(G, q^2, q^2))
for(g in 1:G){
  metaan <- rma.mv(yi=vecPhi, V=CovMx, mods = ~ outcome + outcome:I(repTI-unique(repTI)[g]) - 1, method = "FE")
  Phi_Mod[g,] <- coef(metaan)[1:q^2]
  sePhi_Mod[g,] <- metaan$se[1:q^2]
  CovMxPhi_Mod[g,,] <- metaan$vb[1:q^2,1:q^2]
  #
  #
  # Determine points on 95% LL contour
  vecPhi_metaan <- Phi_Mod[g,]
  CovMx_metaan <- CovMxPhi_Mod[g,,]
  eigenCovMx <- eigen(CovMx_metaan) 
  lambda <- eigenCovMx$val
  E <- eigenCovMx$vec
  Chi2 <- qchisq(p=0.05, df=(q*q), lower.tail=FALSE)
  LB_vecPhi <- matrix(NA, nrow=q*q, ncol =q*q)
  UB_vecPhi <- matrix(NA, nrow=q*q, ncol =q*q)
  LL <- matrix(NA, nrow=q*q, ncol=2)
  teller = 0
  for(row in 1:q){
    for(column in 1:q){
      teller = teller + 1
      LB_vecPhi[teller,] <- vecPhi_metaan - sqrt(Chi2 * lambda[teller]) * E[,teller]
      UB_vecPhi[teller,] <- vecPhi_metaan + sqrt(Chi2 * lambda[teller]) * E[,teller]  
      LL[teller,1] <- t(LB_vecPhi[teller,]-vecPhi_metaan) %*% solve(CovMx_metaan) %*% (LB_vecPhi[teller,]-vecPhi_metaan)
      LL[teller,2] <- t(UB_vecPhi[teller,]-vecPhi_metaan) %*% solve(CovMx_metaan) %*% (UB_vecPhi[teller,]-vecPhi_metaan)
      Phi_LB_t <- matrix(LB_vecPhi[teller,], ncol=q, nrow=q)
      Phi_UB_t <- matrix(UB_vecPhi[teller,], ncol=q, nrow=q)
    }
  }
  # Note that UB can be smaller than LB and LB larger than UB!
  minPhi_Mod[g,] <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, min)
  maxPhi_Mod[g,] <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, max)
}
(minPhi_Mod < vecPhi_G_pop) + (maxPhi_Mod > vecPhi_G_pop) - 1
#unique(TI)


# Linear & quadratic moderator 
#
# metaan <- rma.mv(yi=vecPhi, V=CovMx, mods = ~ outcome + outcome:I(repTI-mean(repTI)) + outcome:I((repTI-mean(repTI))^2) - 1, method = "FE")
# To obtain estimates with cov mx / se's per time interval, we should substract TI[s]-mean(TI). Then mean(TI) cancels out!
#
G <- length(unique(TI))
minPhi_Quadr <- matrix(NA, ncol=(q^2), nrow=G)
maxPhi_Quadr <- matrix(NA, ncol=(q^2), nrow=G)
Phi_Quadr <- matrix(NA, ncol=(q^2), nrow=G)
sePhi_Quadr <- matrix(NA, ncol=(q^2), nrow=G)
CovMxPhi_Quadr <- array(NA, dim = c(G, q^2, q^2))
for(g in 1:G){
  metaan <- rma.mv(yi=vecPhi, V=CovMx, mods = ~ outcome + outcome:I(repTI-unique(repTI)[g]) + outcome:I((repTI-unique(repTI)[g])^2) - 1, method = "FE")
  Phi_Quadr[g,] <- coef(metaan)[1:q^2]
  sePhi_Quadr[g,] <- metaan$se[1:q^2]
  CovMxPhi_Quadr[g,,] <- metaan$vb[1:q^2,1:q^2]
  #
  #
  # Determine points on 95% LL contour
  vecPhi_metaan <- Phi_Quadr[g,]
  CovMx_metaan <- CovMxPhi_Quadr[g,,]
  eigenCovMx <- eigen(CovMx_metaan) 
  lambda <- eigenCovMx$val
  E <- eigenCovMx$vec
  Chi2 <- qchisq(p=0.05, df=(q*q), lower.tail=FALSE)
  LB_vecPhi <- matrix(NA, nrow=q*q, ncol =q*q)
  UB_vecPhi <- matrix(NA, nrow=q*q, ncol =q*q)
  LL <- matrix(NA, nrow=q*q, ncol=2)
  teller = 0
  for(row in 1:q){
    for(column in 1:q){
      teller = teller + 1
      LB_vecPhi[teller,] <- vecPhi_metaan - sqrt(Chi2 * lambda[teller]) * E[,teller]
      UB_vecPhi[teller,] <- vecPhi_metaan + sqrt(Chi2 * lambda[teller]) * E[,teller]  
      LL[teller,1] <- t(LB_vecPhi[teller,]-vecPhi_metaan) %*% solve(CovMx_metaan) %*% (LB_vecPhi[teller,]-vecPhi_metaan)
      LL[teller,2] <- t(UB_vecPhi[teller,]-vecPhi_metaan) %*% solve(CovMx_metaan) %*% (UB_vecPhi[teller,]-vecPhi_metaan)
      Phi_LB_t <- matrix(LB_vecPhi[teller,], ncol=q, nrow=q)
      Phi_UB_t <- matrix(UB_vecPhi[teller,], ncol=q, nrow=q)
    }
  }
  # Note that UB can be smaller than LB and LB larger than UB!
  minPhi_Quadr[g,] <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, min)
  maxPhi_Quadr[g,] <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, max)
}
(minPhi_Quadr < vecPhi_G_pop) + (maxPhi_Quadr > vecPhi_G_pop) - 1
#unique(TI)


# Dummies
repD <- as.matrix(repD)
metaan <- rma.mv(yi=vecPhi, V=CovMx, mods = ~ outcome:repD - 1, method = 'FE') 
#summary(metaan)
Phi_D <- matrix(coef(metaan), length(unique(TI)), q^2, byrow=T) 
sePhi_D <- matrix(metaan$se, length(unique(TI)), q^2, byrow=T) 
CovMxPhi_D <- metaan$vb
#
# elliptical 95%CI
# Determine points on 95% LL contour - Has to be done per unique time interval
minPhi_D <- matrix(NA, length(unique(TI)), q^2)
maxPhi_D <- matrix(NA, length(unique(TI)), q^2)
tellerCovMx <- 1
for(i in 1:length(unique(TI))){
  vecPhi_metaan <- Phi_D[i,]
  CovMx_metaan <- CovMxPhi_D[tellerCovMx:(tellerCovMx+3),tellerCovMx:(tellerCovMx+3)]
  tellerCovMx <- tellerCovMx + 4
  eigenCovMx <- eigen(CovMx_metaan) 
  lambda <- eigenCovMx$val
  E <- eigenCovMx$vec
  Chi2 <- qchisq(p=0.05, df=(q*q), lower.tail=FALSE)
  LB_vecPhi <- matrix(NA, nrow=q*q, ncol =q*q)
  UB_vecPhi <- matrix(NA, nrow=q*q, ncol =q*q)
  LL <- matrix(NA, nrow=q*q, ncol=2)
  teller = 0
  for(row in 1:q){
    for(column in 1:q){
      teller = teller + 1
      LB_vecPhi[teller,] <- vecPhi_metaan - sqrt(Chi2 * lambda[teller]) * E[,teller]
      UB_vecPhi[teller,] <- vecPhi_metaan + sqrt(Chi2 * lambda[teller]) * E[,teller]  
      LL[teller,1] <- t(LB_vecPhi[teller,]-vecPhi_metaan) %*% solve(CovMx_metaan) %*% (LB_vecPhi[teller,]-vecPhi_metaan)
      LL[teller,2] <- t(UB_vecPhi[teller,]-vecPhi_metaan) %*% solve(CovMx_metaan) %*% (UB_vecPhi[teller,]-vecPhi_metaan)
      Phi_LB_t <- matrix(LB_vecPhi[teller,], ncol=q, nrow=q)
      Phi_UB_t <- matrix(UB_vecPhi[teller,], ncol=q, nrow=q)
    }
  }
  #
  # Note that UB can be smaller than LB and LB larger than UB!
  minPhi_D[i,] <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, min)
  maxPhi_D[i,] <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, max)
}
(minPhi_D < vecPhi_G_pop) + (maxPhi_D > vecPhi_G_pop) - 1
# unique(TI)


# CTmeta
G <- length(unique(TI))
Phi_Trans <- matrix(NA, ncol=(q^2), nrow=G)
sePhi_Trans <- matrix(NA, ncol=(q^2), nrow=G)
CovMxPhi_Trans <- array(NA, dim = c(G, q^2, q^2))
for(g in 1:G){
  metaan <- rma.mv(yi=vecPhi_G[,g], V=CovMx_G[,,g], mods = ~ outcome - 1, method = "FE") 
  Phi_Trans[g,] <- coef(metaan)
  sePhi_Trans[g,] <- metaan$se
  CovMxPhi_Trans[g,,] <- metaan$vb
}
# Determine points on 95% LL contour - Has to be done per unique time interval
minPhi_Trans <- matrix(NA, length(unique(TI)), q^2)
maxPhi_Trans <- matrix(NA, length(unique(TI)), q^2)
tellerCovMx <- 1
for(i in 1:length(unique(TI))){
  vecPhi_metaan <- Phi_Trans[i,]
  CovMx_metaan <- CovMxPhi_Trans[i,,]
  tellerCovMx <- tellerCovMx + 4
  eigenCovMx <- eigen(CovMx_metaan) 
  lambda <- eigenCovMx$val
  E <- eigenCovMx$vec
  Chi2 <- qchisq(p=0.05, df=(q*q), lower.tail=FALSE)
  LB_vecPhi <- matrix(NA, nrow=q*q, ncol =q*q)
  UB_vecPhi <- matrix(NA, nrow=q*q, ncol =q*q)
  LL <- matrix(NA, nrow=q*q, ncol=2)
  teller = 0
  for(row in 1:q){
    for(column in 1:q){
      teller = teller + 1
      LB_vecPhi[teller,] <- vecPhi_metaan - sqrt(Chi2 * lambda[teller]) * E[,teller]
      UB_vecPhi[teller,] <- vecPhi_metaan + sqrt(Chi2 * lambda[teller]) * E[,teller]  
      LL[teller,1] <- t(LB_vecPhi[teller,]-vecPhi_metaan) %*% solve(CovMx_metaan) %*% (LB_vecPhi[teller,]-vecPhi_metaan)
      LL[teller,2] <- t(UB_vecPhi[teller,]-vecPhi_metaan) %*% solve(CovMx_metaan) %*% (UB_vecPhi[teller,]-vecPhi_metaan)
      Phi_LB_t <- matrix(LB_vecPhi[teller,], ncol=q, nrow=q)
      Phi_UB_t <- matrix(UB_vecPhi[teller,], ncol=q, nrow=q)
    }
  }
  #
  # Note that UB can be smaller than LB and LB larger than UB!
  minPhi_Trans[i,] <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, min)
  maxPhi_Trans[i,] <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, max)
}
(minPhi_Trans < vecPhi_G_pop) + (maxPhi_Trans > vecPhi_G_pop) - 1
(minPhi_Trans < 0) + (maxPhi_Trans > 0) - 1 # you want this to be 0: then, 0 is not contained in 95CI
# unique(TI)



#--- Save and Load (for new data)
save(list = ls(), file = "Example_CTmeta.RData")


################

# Plot
Col <- c(3, "red", 8, 1, 4, "orange", "purple")

## Plot Neglecting and Lin & Quadr
Col <- c(3, "red", 8, 1, 4, "orange", "purple")
for(teller in 1:q^2){
  op <- par(mfrow=c(1,1)) 
  #title <- as.list(expression(paste("How does ", phi["jk"](Delta[t]), " and its overall estimate vary ", "as a function of the time-interval"))) 
  j = (teller+q-1)%/%q
  k = teller-(j-1)*q
  plot(y=PhiTIs[j,k,], x=TIs, type="l", 
       ylim=
         c(min(minPhi_N[teller], minPhi_Mod[,teller], minPhi_Quadr[,teller],
               maxPhi_N[teller], maxPhi_Mod[,teller], maxPhi_Quadr[,teller]), 
           max(minPhi_N[teller], minPhi_Mod[,teller], minPhi_Quadr[,teller],
               maxPhi_N[teller], maxPhi_Mod[,teller], maxPhi_Quadr[,teller])), 
       ylab = expression(paste(phi["jk"](Delta[t]), " values")), xlab = expression(paste("Time-interval (", Delta[t], ")", sep="")), 
       col=Col[1], lwd=2, lty=1
       #,
       #main=mtext(do.call(expression, title), side=3) #, line = c(2,1,0), cex = 1 )
  )
  #
  #
  lines(y=Phi_Mod[,teller], x=unique(TI), col=Col[3], lwd=0.5, lty=1, type = "p", pch = 19)
  LB <- minPhi_Mod[,teller]
  UB <- maxPhi_Mod[,teller]
  lines(y=LB,x=unique(TI), col=Col[3], lwd=1, lty=1, type = "p", pch = 24) #, lwd=10, lty=1, type = "p", "-")
  lines(y=UB,x=unique(TI), col=Col[3], lwd=1, lty=1, type = "p", pch = 25) #, lwd=10, lty=1, type = "p", "+")
  #
  lines(y=LB[order(unique(TI))],x=unique(TI)[order(unique(TI))], col=Col[3], lwd=1, lty=1, type = "l")
  lines(y=UB[order(unique(TI))],x=unique(TI)[order(unique(TI))], col=Col[3], lwd=1, lty=1, type = "l")
  #
  #
  lines(y=Phi_Quadr[,teller], x=unique(TI), col=Col[6], lwd=0.5, lty=1, type = "p", pch = 19)
  LB <- minPhi_Quadr[,teller]
  UB <- maxPhi_Quadr[,teller]
  lines(y=LB,x=unique(TI), col=Col[6], lwd=1, lty=1, type = "p", pch = 24) 
  lines(y=UB,x=unique(TI), col=Col[6], lwd=1, lty=1, type = "p", pch = 25)
  #
  lines(y=LB[order(unique(TI))],x=unique(TI)[order(unique(TI))], col=Col[6], lwd=1, lty=1, type = "l")
  lines(y=UB[order(unique(TI))],x=unique(TI)[order(unique(TI))], col=Col[6], lwd=1, lty=1, type = "l")
  #
  #
  lines(y=rep(Phi_N[teller],12),x=unique(TI)[order(unique(TI))], col=Col[2], lwd=0.5, lty=1, type = "p", pch = 19)
  lines(y=rep(minPhi_N[teller],12),x=unique(TI)[order(unique(TI))], col=Col[2], lwd=1, lty=1, type = "p", pch = 24)
  lines(y=rep(maxPhi_N[teller],12),x=unique(TI)[order(unique(TI))], col=Col[2], lwd=1, lty=1, type = "p", pch = 25)
  lines(y=rep(minPhi_N[teller],12),x=unique(TI)[order(unique(TI))], col=Col[2], lwd=0.5, lty=1, type = "l")
  lines(y=rep(maxPhi_N[teller],12),x=unique(TI)[order(unique(TI))], col=Col[2], lwd=0.5, lty=1, type = "l")
  #
  #
  e1 <- "True values"
  e5 <- "Overall estimates ignoring the time interval"
  e3 <- "Overall estimates using a linear relationship"
  e6 <- "Overall estimates using a quadratic relationship"
  legendT2 = c(e1,e5,e3,e6)
  #if(j==3 & k==1){
  #  pos <- "bottomright"
  #}else{
  #  pos <- "topright"
  #}
  pos <- "bottomleft"
  legend(pos,
         legend = legendT2, cex=1,
         bty = "n",
         lty=rep(1,4), # gives the legend appropriate symbols (lines)
         lwd=rep(0,4), 
         pch = rep(19,4),
         col=Col[c(1,2,3,6)] # gives the legend lines the correct color and width
  )
  dev.copy(pdf,paste0('PhiPlot_Example_Phi', j, k, '_95CI_N_L_Q.pdf'))
  dev.off()
  par(op)
}
# NOTE: scale differs per plot!



## Plot Dummy and trans 
Col <- c(3, "red", 8, 1, 4, "orange", "purple")
for(teller in 1:q^2){
  op <- par(mfrow=c(1,1)) 
  #title <- as.list(expression(paste("How does ", phi["jk"](Delta[t]), " and its overall estimate vary ", "as a function of the time-interval"))) 
  j = (teller+q-1)%/%q
  k = teller-(j-1)*q
  plot(y=PhiTIs[j,k,], x=TIs, type="l", 
       ylim=
         c(min(minPhi_D[,teller], minPhi_Trans[,teller],
               maxPhi_D[,teller], maxPhi_Trans[,teller]), 
           max(minPhi_D[,teller], minPhi_Trans[,teller],
               maxPhi_D[,teller], maxPhi_Trans[,teller])), 
       ylab = expression(paste(phi["jk"](Delta[t]), " values")), xlab = expression(paste("Time-interval (", Delta[t], ")", sep="")), 
       col=Col[1], lwd=2, lty=1
       #,
       #main=mtext(do.call(expression, title), side=3) #, line = c(2,1,0), cex = 1 )
  )
  #
  #
  lines(y=Phi_D[,teller], x=unique(TI), col=Col[4], lwd=0.5, lty=1, type = "p", pch = 19)
  LB <- minPhi_D[,teller]
  UB <- maxPhi_D[,teller]
  lines(y=LB,x=unique(TI), col=Col[4], lwd=1, lty=1, type = "p", pch = 24)
  lines(y=UB,x=unique(TI), col=Col[4], lwd=1, lty=1, type = "p", pch = 25) 
  #
  lines(y=LB[order(unique(TI))],x=unique(TI)[order(unique(TI))], col=Col[4], lwd=0.5, lty=1, type = "l")
  lines(y=UB[order(unique(TI))],x=unique(TI)[order(unique(TI))], col=Col[4], lwd=0.5, lty=1, type = "l")
  #
  #
  lines(y=Phi_Trans[,teller], x=unique(TI), col=Col[5], lwd=0.5, lty=1, type = "p", pch = 19)
  LB <- minPhi_Trans[,teller]
  UB <- maxPhi_Trans[,teller]
  lines(y=LB,x=unique(TI), col=Col[5], lwd=1, lty=1, type = "p", pch = 24) 
  lines(y=UB,x=unique(TI), col=Col[5], lwd=1, lty=1, type = "p", pch = 25) 
  #
  ###lines(y=Phi_Trans[order(unique(TI)),teller],x=unique(TI)[order(unique(TI))], col=Col[5], lwd=2, lty=1, type = "l")
  ##lines(y=LB[order(unique(TI))],x=unique(TI)[order(unique(TI))], col=Col[5], lwd=1, lty=1, type = "l")
  ##lines(y=UB[order(unique(TI))],x=unique(TI)[order(unique(TI))], col=Col[5], lwd=1, lty=1, type = "l")
  ## Get smooth lines:
  end1 <- 1 
  end <- G 
  x <- unique(TI)[order(unique(TI))][end1:end]
  y <- LB[order(unique(TI))][end1:end]
  lo <- loess(y~x)
  xl <- seq(min(x),max(x), (max(x) - min(x))/1000)
  lines(xl, predict(lo,xl), col=Col[5], lwd=1, lty=1, type = "l")
  x <- unique(TI)[order(unique(TI))][end1:end]
  y <- UB[order(unique(TI))][end1:end]
  lo <- loess(y~x)
  xl <- seq(min(x),max(x), (max(x) - min(x))/1000)
  lines(xl, predict(lo,xl), col=Col[5], lwd=1, lty=1, type = "l")
  #
  #
  e1 <- "True values"
  e4 <- "Overall estimates using dummies"
  e2 <- "Overall estimates using CTmeta"
  legendT2 = c(e1,e4,e2)
  if(j==3 & k==1){
    pos <- "bottomright"
  }else{
    pos <- "topright"
  }
  legend(pos,
         legend = legendT2, cex=1,
         bty = "n",
         lty=rep(1,3), # gives the legend appropriate symbols (lines)
         lwd=rep(0,3), 
         pch = rep(19,3),
         col=Col[c(1,4,5)] # gives the legend lines the correct color and width
  )
  dev.copy(pdf,paste0('PhiPlot_Example_Phi', j, k, '_95CI_D_T.pdf'))
  dev.off()
  par(op)
}
# NOTE: scale differs per plot!


