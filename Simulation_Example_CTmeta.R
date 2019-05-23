
# Simulation based on Example_CTmeta.R
NrSim <- 1000
set.seed(123)
#set.seed(124)

# Prelimanary - based on Example_CTmeta.R
if (!require("expm")) install.packages("expm")
library(expm)
if (!require("metafor")) install.packages("metafor")
library(metafor)
if (!require("tsDyn")) install.packages("tsDyn")
library(tsDyn)
if (!require("vars")) install.packages("vars")
library(vars)
#
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


# Used from Work engagement paper
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
#N <- matrix(100,S) # CHECK this once

var1 <- matrix(NA, ncol=(q^2), nrow=S)
for(s in 1:S){
  var1[s,] <- rep(1/(N[s] - 3), q^2)
}

# Time intervals ('time lags'), using Work engagement paper
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
DeltaT <- matrix(c(0.75, 1, 2), byrow=TRUE)


vecPhi_G_pop <- matrix(NA, G, q^2)
for(g in 1:G){
  vecPhi_G_pop[g,] <- as.vector(t(expm(A_pop * unique(TI)[g])))
}


## Dummies ##
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

#######################################################################################################################################

## Simulation based on Example ##
Phi_true <- array(NA, dim=c(q,q,G, NrSim))
#
minPhi_N <- array(NA, dim=c(NrSim, q^2))
maxPhi_N <- array(NA, dim=c(NrSim, q^2))
minPhi_Mod <- array(NA, dim=c(NrSim, G, q^2))
maxPhi_Mod <- array(NA, dim=c(NrSim, G, q^2))
minPhi_Quadr <- array(NA, dim=c(NrSim, G, q^2))
maxPhi_Quadr <- array(NA, dim=c(NrSim, G, q^2))
minPhi_D <- array(NA, dim=c(NrSim, G, q^2))
maxPhi_D <- array(NA, dim=c(NrSim, G, q^2))
minPhi_Trans <- array(NA, dim=c(NrSim, G, q^2))
maxPhi_Trans <- array(NA, dim=c(NrSim, G, q^2))
#
Coverage_LB_Mod <- array(NA, dim=c(q,q,G, NrSim))
Coverage_UB_Mod <- array(NA, dim=c(q,q,G, NrSim))
Coverage_Mod <- array(NA, dim=c(q,q,G, NrSim))
Coverage_LB_Quadr <- array(NA, dim=c(q,q,G, NrSim))
Coverage_UB_Quadr <- array(NA, dim=c(q,q,G, NrSim))
Coverage_Quadr <- array(NA, dim=c(q,q,G, NrSim))
Coverage_LB_D <- array(NA, dim=c(q,q,G, NrSim))
Coverage_UB_D <- array(NA, dim=c(q,q,G, NrSim))
Coverage_D <- array(NA, dim=c(q,q,G, NrSim))
Coverage_LB_D_var1 <- array(NA, dim=c(q,q,G, NrSim))
Coverage_UB_D_var1 <- array(NA, dim=c(q,q,G, NrSim))
Coverage_D_var1 <- array(NA, dim=c(q,q,G, NrSim))
Coverage_LB_Trans <- array(NA, dim=c(q,q,G, NrSim))
Coverage_UB_Trans <- array(NA, dim=c(q,q,G, NrSim))
Coverage_Trans <- array(NA, dim=c(q,q,G, NrSim))
#
Coverage_0_Mod <- array(NA, dim=c(q,q,G, NrSim))
Coverage_0_Quadr <- array(NA, dim=c(q,q,G, NrSim))
Coverage_0_D <- array(NA, dim=c(q,q,G, NrSim))
Coverage_0_Trans <- array(NA, dim=c(q,q,G, NrSim))
#
MeanPhi_Mod <- array(0, dim=c(G,q^2))
MeanPhi_Quadr <- array(0, dim=c(G,q^2))
MeanPhi_D <- array(0, dim=c(G,q^2))
MeanPhi_Trans <- array(0, dim=c(G,q^2))
MeanSqDiffPhi_Mod <- array(0, dim=c(G,q^2))
MeanSqDiffPhi_Quadr <- array(0, dim=c(G,q^2))
MeanSqDiffPhi_D <- array(0, dim=c(G,q^2))
MeanSqDiffPhi_Trans <- array(0, dim=c(G,q^2))
MeanAbsDiffPhi_Mod <- array(0, dim=c(G,q^2))
MeanAbsDiffPhi_Quadr <- array(0, dim=c(G,q^2))
MeanAbsDiffPhi_D <- array(0, dim=c(G,q^2))
MeanAbsDiffPhi_Trans <- array(0, dim=c(G,q^2))
#
SeTransLower <- array(NA, dim=c(G,NrSim))


# Needed in the meta-an
sub = NULL
for(i in 1:q){
  sub = c(sub, paste(i, 1:q, sep=""))
}
outcome <- rep(sub, S) 


# Simulate data and do the meta-analyses
NotPosDefCovMx <- array(NA, dim=c(NrSim,S))
for(sim in 1:NrSim){

  # Example
  
  # Sample for each study s in the meta-analysis data based on population values and N_s and DeltaT_s
  S <- length(N)
  vecA <- array(NA, dim = c(S*q*q))
  vecPhi <- array(NA, dim = c(S*q*q))
  CovMx <- array(0, dim = c(S*q*q, S*q*q))
  vecPhi1 <- array(NA, dim = c(S*q*q))
  CovMx1 <- array(0, dim = c(S*q*q, S*q*q))
  #
  G <- length(unique(TI))
  vecPhi_G <- array(NA, dim = c(S*q*q, G))
  CovMx_G <- array(0, dim = c(S*q*q, S*q*q, G))
  #
  s <- 1
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
      vecA[((s-1)*(q*q)+1):(s*q*q)] <- as.vector(t(-B_VARest))
      #
      Phi_VARest <- expm(-B_VARest*TI[s])
      vecPhi[((s-1)*(q*q)+1):(s*q*q)] <- as.vector(t(Phi_VARest))
      SigmaVAR_VARest <- Gamma_VARest - Phi_VARest %*% Gamma_VARest %*% t(Phi_VARest)
      CovMx[((s-1)*q^2+1):(s*q^2),((s-1)*q^2+1):(s*q^2)] <- kronecker(SigmaVAR_VARest, invGamma) / (N[s]-q)
      if(any( eigen( CovMx[((s-1)*q^2+1):(s*q^2),((s-1)*q^2+1):(s*q^2)] )$values < 0 )){
        s <- s # Cov mx should be pos def
      }else{
        #
        Phi_VARest <- expm(-B_VARest)
        vecPhi1[((s-1)*(q*q)+1):(s*q*q)] <- as.vector(t(Phi_VARest))
        SigmaVAR1_VARest <- Gamma_VARest - Phi_VARest %*% Gamma_VARest %*% t(Phi_VARest)
        CovMx1[((s-1)*q^2+1):(s*q^2),((s-1)*q^2+1):(s*q^2)] <- kronecker(SigmaVAR1_VARest, invGamma) / (N[s]-q)
        # this equates
        #CovMx[((s-1)*q^2+1):(s*q^2),((s-1)*q^2+1):(s*q^2)] * kronecker(SigmaVAR1_VARest / SigmaVAR_VARest, matrix(1,q,q))
        if(any( eigen( CovMx1[((s-1)*q^2+1):(s*q^2),((s-1)*q^2+1):(s*q^2)] )$values < 0 )){
          s <- s # Cov mx should be pos def
        }else{
          for(g in 1:G){
            Phi_VARest <- expm(-B_VARest*unique(TI)[g])
            vecPhi_G[((s-1)*(q*q)+1):(s*q*q),g] <- as.vector(t(Phi_VARest))
            SigmaVAR_VARest <- Gamma_VARest - Phi_VARest %*% Gamma_VARest %*% t(Phi_VARest)
            CovMx_G[((s-1)*q^2+1):(s*q^2),((s-1)*q^2+1):(s*q^2),g] = kronecker(SigmaVAR_VARest, invGamma) / (N[s]-q)
            # this equates
            #CovMx1[((s-1)*q^2+1):(s*q^2),((s-1)*q^2+1):(s*q^2)] * kronecker(SigmaVAR_VARest / SigmaVAR1_VARest, matrix(1,q,q))
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
  
  ###########################################################################################################
  
  # Multivariate Meta-analyses
  if( any( abs(CovMx[lower.tri(diag(S*q^2))] - CovMx[upper.tri(diag(S*q^2))]) < 0.1^10 ) ){
    CovMx <- (CovMx + t(CovMx))/2
  }
  
  # Ignore time interval
  metaan <- rma.mv(yi=vecPhi, V=CovMx, mods = ~ outcome - 1, method = "FE")
  #summary(metaan)
  Phi_N <- coef(metaan)
  sePhi_N <- metaan$se 
  CovMxPhi_N <- metaan$vb
  #
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
  #
  # Note that UB can be smaller than LB and LB larger than UB!
  minPhi_N[sim,] <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, min)
  maxPhi_N[sim,] <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, max)
  
  
  # Linear moderator
  G <- length(unique(TI))
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
    minPhi_Mod[sim,g,] <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, min)
    maxPhi_Mod[sim,g,] <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, max)
  }
  
  # Linear & quadratic moderator 
  # metaan <- rma.mv(yi=vecPhi, V=CovMx, mods = ~ outcome + outcome:I(repTI-mean(repTI)) + outcome:I((repTI-mean(repTI))^2) - 1, method = "FE")
  # To obtain estimates with cov mx / se's per time interval, we should substract TI[s]-mean(TI). Then mean(TI) cancels out!
  G <- length(unique(TI))
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
    minPhi_Quadr[sim,g,] <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, min)
    maxPhi_Quadr[sim,g,] <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, max)
  }

  
  # Dummies
  repD <- as.matrix(repD)
  metaan <- rma.mv(yi=vecPhi, V=CovMx, mods = ~ outcome:repD - 1, method = 'FE') 
  
  #summary(metaan)
  Phi_D <- matrix(coef(metaan), length(unique(TI)), q^2, byrow=T) #coef(metaan)
  sePhi_D <- matrix(metaan$se, length(unique(TI)), q^2, byrow=T) #metaan$se 
  CovMxPhi_D <- metaan$vb
  #
  #
  # Determine points on 95% LL contour - Has to be done per unique time interval
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
    minPhi_D[sim,i,] <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, min)
    maxPhi_D[sim,i,] <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, max)
  }
  
  
  # CTmeta
  G <- length(unique(TI))
  Phi_Trans <- matrix(NA, ncol=(q^2), nrow=G)
  sePhi_Trans <- matrix(NA, ncol=(q^2), nrow=G)
  CovMxPhi_Trans <- array(NA, dim = c(G, q^2, q^2))
  for(g in 1:G){
    if( any( abs(CovMx_G[,,g][lower.tri(diag(S*q^2))] - CovMx_G[,,g][upper.tri(diag(S*q^2))]) < 0.1^10 ) ){
      CovMx_G[,,g] <- (CovMx_G[,,g] + t(CovMx_G[,,g]))/2
    }
    metaan <- rma.mv(yi=vecPhi_G[,g], V=CovMx_G[,,g], mods = ~ outcome - 1, method = "FE") 
    Phi_Trans[g,] <- coef(metaan)
    sePhi_Trans[g,] <- metaan$se
    CovMxPhi_Trans[g,,] <- metaan$vb
  }
  #
  # Determine points on 95% LL contour - Has to be done per unique time interval
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
    minPhi_Trans[sim,i,] <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, min)
    maxPhi_Trans[sim,i,] <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, max)
  }
  # Coverage
  
  #> unique(TI)
  #[,1]
  #[1,] 1.000000000
  #[2,] 0.333333333
  #[3,] 0.750000000
  #[4,] 0.166666667
  #[5,] 3.000000000
  #[6,] 0.666666667
  #[7,] 2.000000000
  #[8,] 0.002739726
  #[9,] 1.166666667
  #[10,] 0.833333333
  #[11,] 4.000000000
  #[12,] 0.083333333
  
  #LB_Phi_metaan_N[j,k]
  #UB_Phi_metaan_N[j,k]
  
  for(g in 1:G){
    Phi_true[,,g,sim] <- expm(A_pop * unique(TI)[g]) # is same per sim of course!
    #
    #LB < pop & UB > pop 
    Coverage_Mod[,,g,sim] <- (minPhi_Mod[sim,g,] < t(Phi_true[,,g,sim])) + (maxPhi_Mod[sim,g,] > t(Phi_true[,,g,sim])) - 1
    #LB < 0 & UB > 0 
    Coverage_0_Mod[,,g,sim] <- (minPhi_Mod[sim,g,] < 0) + (maxPhi_Mod[sim,g,] > 0) - 1
    # Sum of parameters - at the end divide by NrSim
    MeanPhi_Mod[g,] <- MeanPhi_Mod[g,] + Phi_Mod[g,]
    MeanSqDiffPhi_Mod[g,] <- MeanSqDiffPhi_Mod[g,] + t((Phi_Mod[g,] - Phi_true[,,g,1])^2)
    MeanAbsDiffPhi_Mod[g,] <- MeanAbsDiffPhi_Mod[g,] + t(abs(Phi_Mod[g,] - Phi_true[,,g,1]))
    #
    #
    #LB < pop & UB > pop 
    Coverage_Quadr[,,g,sim] <- (minPhi_Quadr[sim,g,] < t(Phi_true[,,g,sim])) + (maxPhi_Quadr[sim,g,] > t(Phi_true[,,g,sim])) - 1
    #LB < 0 & UB > 0 
    Coverage_0_Quadr[,,g,sim] <- (minPhi_Quadr[sim,g,] < 0) + (maxPhi_Quadr[sim,g,] > 0) - 1
    # Sum of parameters - at the end divide by NrSim
    MeanPhi_Quadr[g,] <- MeanPhi_Quadr[g,] + Phi_Quadr[g,]
    MeanSqDiffPhi_Quadr[g,] <- MeanSqDiffPhi_Quadr[g,] + t((Phi_Quadr[g,] - Phi_true[,,g,1])^2)
    MeanAbsDiffPhi_Quadr[g,] <- MeanAbsDiffPhi_Quadr[g,] + t(abs(Phi_Quadr[g,] - Phi_true[,,g,1]))
    #
    #
    #LB < pop & UB > pop 
    Coverage_D[,,g,sim] <- (minPhi_D[sim,g,] < t(Phi_true[,,g,sim])) + (maxPhi_D[sim,g,] > t(Phi_true[,,g,sim])) - 1
    #LB < 0 & UB > 0 
    Coverage_0_D[,,g,sim] <- (minPhi_D[sim,g,] < 0) + (maxPhi_D[sim,g,] > 0) - 1
    # Sum of parameters - at the end divide by NrSim
    MeanPhi_D[g,] <- MeanPhi_D[g,] + Phi_D[g,]
    MeanSqDiffPhi_D[g,] <- MeanSqDiffPhi_D[g,] + t((Phi_D[g,] - Phi_true[,,g,1])^2)
    MeanAbsDiffPhi_D[g,] <- MeanAbsDiffPhi_D[g,] + t(abs(Phi_D[g,] - Phi_true[,,g,1]))
    #
    #
    #LB < pop & UB > pop 
    Coverage_Trans[,,g,sim] <- (minPhi_Trans[sim,g,] < t(Phi_true[,,g,sim])) + (maxPhi_Trans[sim,g,] > t(Phi_true[,,g,sim])) - 1
    #LB < 0 & UB > 0 
    Coverage_0_Trans[,,g,sim] <- (minPhi_Trans[sim,g,] < 0) + (maxPhi_Trans[sim,g,] > 0) - 1
    # Sum of parameters - at the end divide by NrSim
    MeanPhi_Trans[g,] <- MeanPhi_Trans[g,] + Phi_Trans[g,]
    MeanSqDiffPhi_Trans[g,] <- MeanSqDiffPhi_Trans[g,] + t((Phi_Trans[g,] - Phi_true[,,g,1])^2)
    MeanAbsDiffPhi_Trans[g,] <- MeanAbsDiffPhi_Trans[g,] + t(abs(Phi_Trans[g,] - Phi_true[,,g,1]))
    #
    #
    SeTransLower[g,sim] <- all(sePhi_Trans[g,] < sePhi_D[g,])
    #
    #
    # Transform overallA to Phis
    Phi_overallA <- expm(matrix(overallA, byrow=T, ncol=q) * unique(TI)[g])
    #LB < pop & UB > pop 
    Coverage_overallA[,,g,sim] <- (minPhi_overallA[sim,g,] < t(Phi_true[,,g,sim])) + (maxPhi_overallA[sim,g,] > t(Phi_true[,,g,sim])) - 1
    #LB < 0 & UB > 0 
    Coverage_0_overallA[,,g,sim] <- (minPhi_overallA[sim,g,] < 0) + (maxPhi_overallA[sim,g,] > 0) - 1
    # As a guess for LB&UB:
    minPhi_overallA_2[sim,g,] <- as.vector(t(expm(matrix(min_overallA[sim,], byrow=T, ncol=q) * unique(TI)[g])))
    maxPhi_overallA_2[sim,g,] <- as.vector(t(expm(matrix(max_overallA[sim,], byrow=T, ncol=q) * unique(TI)[g])))
    #LB < pop & UB > pop 
    Coverage_overallA_2[,,g,sim] <- (minPhi_overallA_2[sim,g,] < t(Phi_true[,,g,sim])) + (maxPhi_overallA_2[sim,g,] > t(Phi_true[,,g,sim])) - 1
    #LB < 0 & UB > 0 
    Coverage_0_overallA_2[,,g,sim] <- (minPhi_overallA_2[sim,g,] < 0) + (maxPhi_overallA_2[sim,g,] > 0) - 1
    # Sum of parameters - at the end divide by NrSim
    MeanPhi_overallA[g,] <- MeanPhi_overallA[g,] + Phi_Trans[g,]
    MeanSqDiffPhi_overallA[g,] <- MeanSqDiffPhi_overallA[g,] + t((Phi_overallA - Phi_true[,,g,1])^2)
    MeanAbsDiffPhi_overallA[g,] <- MeanAbsDiffPhi_overallA[g,] + t(abs(Phi_overallA - Phi_true[,,g,1]))
  }

} # End simulation

# Coverage
Cov_Mod_ <- apply(Coverage_Mod, c(1,2,3), mean, na.rm=TRUE)
Cov_Quadr_ <- apply(Coverage_Quadr, c(1,2,3), mean, na.rm=TRUE)
Cov_D_ <- apply(Coverage_D, c(1,2,3), mean, na.rm=TRUE)
Cov_Trans_ <- apply(Coverage_Trans, c(1,2,3), mean, na.rm=TRUE)
#
Cov_Mod <- array(NA, c(G, q^2))
Cov_Quadr <- array(NA, c(G, q^2))
Cov_D <- array(NA, c(G, q^2))
Cov_Trans <- array(NA, c(G, q^2))
Phi_true_G <- array(NA, c(G, q^2))
for(g in 1:G){
  Cov_Mod[g,] <- t(Cov_Mod_[,,g])
  Cov_Quadr[g,] <- t(Cov_Quadr_[,,g])
  Cov_D[g,] <- t(Cov_D_[,,g])
  Cov_Trans[g,] <- t(Cov_Trans_[,,g])
  Phi_true_G[g,] <- t(Phi_true[,,g,1])
}
rownames(Cov_Mod) <- unique(TI)
rownames(Cov_Quadr) <- unique(TI)
rownames(Cov_D) <- unique(TI)
rownames(Cov_Trans) <- unique(TI)
#
Cov_Mod[order(unique(TI)),]
Cov_Quadr[order(unique(TI)),]
Cov_D[order(unique(TI)),]
Cov_Trans[order(unique(TI)),]
apply((abs(Cov_Trans - 0.95) < abs(Cov_D - 0.95))[order(unique(TI)),], 2 ,mean)
apply((abs(Cov_Trans - 0.95) == abs(Cov_D - 0.95))[order(unique(TI)),], 2 ,mean)
apply((abs(Cov_Trans - 0.95) > abs(Cov_D - 0.95))[order(unique(TI)),], 2 ,mean)


# Power: 0 (not) contained in 95%
Cov_0_Mod_ <- apply(Coverage_0_Mod, c(1,2,3), mean, na.rm=TRUE)
Cov_0_Quadr_ <- apply(Coverage_0_Quadr, c(1,2,3), mean, na.rm=TRUE)
Cov_0_D_ <- apply(Coverage_0_D, c(1,2,3), mean, na.rm=TRUE)
Cov_0_Trans_ <- apply(Coverage_0_Trans, c(1,2,3), mean, na.rm=TRUE)
#
Cov_0_Mod <- array(NA, c(G, q^2))
Cov_0_Quadr <- array(NA, c(G, q^2))
Cov_0_D <- array(NA, c(G, q^2))
Cov_0_Trans <- array(NA, c(G, q^2))
Phi_true_G <- array(NA, c(G, q^2))
for(g in 1:G){
  Cov_0_Mod[g,] <- t(Cov_0_Mod_[,,g])
  Cov_0_Quadr[g,] <- t(Cov_0_Quadr_[,,g])
  Cov_0_D[g,] <- t(Cov_0_D_[,,g])
  Cov_0_Trans[g,] <- t(Cov_0_Trans_[,,g])
  Phi_true_G[g,] <- t(Phi_true[,,g,1])
}
rownames(Cov_0_Mod) <- unique(TI)
rownames(Cov_0_Quadr) <- unique(TI)
rownames(Cov_0_D) <- unique(TI)
rownames(Cov_0_Trans) <- unique(TI)
#
apply(Cov_0_Trans, 2, mean)
apply(Cov_0_D, 2, mean)
#
Cov_0_Mod[order(unique(TI)),] # one wants this to be zero, at least for specific time intervals
Cov_0_Quadr[order(unique(TI)),] # one wants this to be zero, at least for specific time intervals
Cov_0_D[order(unique(TI)),]
Cov_0_Trans[order(unique(TI)),]
#
apply((Cov_0_Trans < Cov_0_D)[order(unique(TI)),], 2, mean)
apply((Cov_0_Trans == Cov_0_D)[order(unique(TI)),], 2, mean)
apply((Cov_0_Trans > Cov_0_D)[order(unique(TI)),], 2, mean)


# Compare standard errors of parameters
apply(SeTransLower, c(1), mean, na.rm=TRUE)


apply(maxPhi_Mod - minPhi_Mod, c(2,3), mean)[order(unique(TI)),]
apply(maxPhi_Quadr - minPhi_Quadr, c(2,3), mean)[order(unique(TI)),]
apply(maxPhi_D - minPhi_D, c(2,3), mean)[order(unique(TI)),]
apply(maxPhi_Trans - minPhi_Trans, c(2,3), mean)[order(unique(TI)),]
apply(maxPhi_Trans - minPhi_Trans, c(2,3), mean)[order(unique(TI)),] < apply(maxPhi_D - minPhi_D, c(2,3), mean)[order(unique(TI)),]
apply(maxPhi_D - minPhi_D, c(2,3), mean)[order(unique(TI)),] / apply(maxPhi_Trans - minPhi_Trans, c(2,3), mean)[order(unique(TI)),]


# Sum of parameters - at the end divide by NrSim
##MeanPhi_Mod / NrSim
##MeanPhi_D / NrSim
##MeanPhi_Trans / NrSim
##Phi_true_G
#MeanPhi_Mod / NrSim - Phi_true_G
#MeanPhi_D / NrSim - Phi_true_G
#MeanPhi_Trans / NrSim - Phi_true_G


# MSE
(MeanSqDiffPhi_Mod / NrSim)[order(unique(TI)),]
(MeanSqDiffPhi_Quadr / NrSim)[order(unique(TI)),]
(MeanSqDiffPhi_D / NrSim)[order(unique(TI)),]
(MeanSqDiffPhi_Trans / NrSim)[order(unique(TI)),]


apply(MeanSqDiffPhi_D / NrSim, 2, mean) / vecPhi_pop^2
apply(MeanSqDiffPhi_Trans / NrSim, 2, mean) / vecPhi_pop^2
# The MSE can be written as the sum of the variance of the estimator and the squared bias of the estimator


# RMSE
# Taking the square root of the average squared errors has some interesting implications for RMSE. 
# Since the errors are squared before they are averaged, 
# the RMSE gives a relatively high weight to large errors. 
# This means the RMSE should be more useful when large errors are particularly undesirable.
sqrt(MeanSqDiffPhi_Mod / NrSim)[order(unique(TI)),]
sqrt(MeanSqDiffPhi_Quadr / NrSim)[order(unique(TI)),]
sqrt(MeanSqDiffPhi_D / NrSim)[order(unique(TI)),]
sqrt(MeanSqDiffPhi_Trans / NrSim)[order(unique(TI)),]
apply(sqrt(MeanSqDiffPhi_Trans / NrSim)[order(unique(TI)),] < sqrt(MeanSqDiffPhi_D / NrSim)[order(unique(TI)),], 2, mean)
apply(sqrt(MeanSqDiffPhi_Trans / NrSim)[order(unique(TI)),] > sqrt(MeanSqDiffPhi_D / NrSim)[order(unique(TI)),], 2, mean)


# MAE
#From an interpretation standpoint, MAE is clearly the winner. 
#RMSE does not describe average error alone and has other implications that are more difficult to tease out and understand.
(MeanAbsDiffPhi_Mod / NrSim)[order(unique(TI)),]
(MeanAbsDiffPhi_Quadr / NrSim)[order(unique(TI)),]
(MeanAbsDiffPhi_D / NrSim)[order(unique(TI)),]
(MeanAbsDiffPhi_Trans / NrSim)[order(unique(TI)),]
#
apply(MeanSqDiffPhi_Trans < MeanSqDiffPhi_D, 2, mean)
apply(MeanSqDiffPhi_Trans > MeanSqDiffPhi_D, 2, mean)


MeanSqDiffPhi_D[order(unique(TI)),] / MeanSqDiffPhi_Trans[order(unique(TI)),] 
MeanAbsDiffPhi_D[order(unique(TI)),] / MeanAbsDiffPhi_Trans[order(unique(TI)),] 




#--- Save and Load (for new data)
save(list = ls(), file = "Sim_CTmeta.RData")  # saves all objects including data in a compressed format; use extensions like .rda, .Rda, .RData etc.

