#Purpose: 
#MCMC for SVR model: 
#p=2, q = 1, with covariates, in the individual form, involving two mean functions M1 and M2 and using dlm package
#sigma_xi (i.e. sigma_U in the paper) following a lognorm prior with mean mu and variance sigma2
#

#Input:
#Y: n*m observation matrix, where n is number of obervations and m is number of subjects
#X: m*k design matrix of k covariates, including intercept
#obs.idx.m: n*m observed index matrix,  1 yes; 0 no 
#group1.idx.m: m*1 if subject belongs to group 1, 1 yes;0 no
#times: (n+1)*1 all possible times, including the intitial time
#sigma2.eps: 1*1
#sigma2.zeta: 2*1   (i.e. sigma_M in the paper)
#sigma2.xi: m*1     (i.e. sigma_U in the paper)
#sigma2.U0: 1*1      variance of initial values U0
#betas: 1*k
#sigma2: 1*1 
#n.iter: number of iteration of MCMC
#burn.in: number of burn-in stage
#thin: thin


#Output:
#sims:  MCMC output 
#for M1(t_j), N1(t_j),
#    M2(t_j), N2(t_j),
#    U_i(t_j), i = 1...m, j=1...n
#    M1(t_0), N1(t_0), M2(t_0), N2(t_0),U_i(t_0).
#    sigma2.eps, sigma2.zeta, sigma2.xi, betas, sigma2

#######################################################
#Auxiliary functions for MCMC steps 
#######################################################
LogPosteriorPDFofOneSigma2Xi <- function(sigma2.xi.i, sigma2.xi.i.tmp.sum, 
                                         mu, sigma2, n){
  return(-(n/2 + 1) * log(sigma2.xi.i) - 
         sigma2.xi.i.tmp.sum / (2 * sigma2.xi.i) - 
         (log(sigma2.xi.i) - mu)^2 / (2 * sigma2))
}


LogProposedPDFofOneSigma2Xi <- function(sigma2.xi.i, alpha, beta.){
  return(-(alpha + 1) * log(sigma2.xi.i) - beta. / sigma2.xi.i)
}

SampleSigma2Xi <- function(theta.sample, Kobs.idx.m, Kdeltas.i, Kns,  sigma2.xi, times,  mus, sigma2, n, m){
  for(i in 1:m){
    theta.sample.i <- theta.sample[c(1, Kobs.idx.m[[i]] + 1), i + 4]
    sigma2.xi.i.tmp.sum = sum((theta.sample.i[-1]-theta.sample.i[-(Kns[i]+1)])^2 / Kdeltas.i[[i]])
    
    prop <- rgamma(1, shape = Kns[i] / 2 + 0.3, 
                      rate  = sigma2.xi.i.tmp.sum / 2 + 0.3)^-1
    old  <- sigma2.xi[i]
    
    r <- exp(LogPosteriorPDFofOneSigma2Xi(prop, sigma2.xi.i.tmp.sum, mus[i], sigma2, Kns[i])+
             LogProposedPDFofOneSigma2Xi(old, Kns[i] / 2 + 0.3, sigma2.xi.i.tmp.sum / 2 + 0.3)-
             LogPosteriorPDFofOneSigma2Xi(old, sigma2.xi.i.tmp.sum, mus[i], sigma2, Kns[i])-
             LogProposedPDFofOneSigma2Xi(prop, Kns[i] / 2 + 0.3, sigma2.xi.i.tmp.sum / 2 + 0.3))
    sigma2.xi[i] <- ifelse(runif(1) < r, prop, old)
  }
  return(sigma2.xi)
}

SampleBetaSigma2 <- function(sigma2.xi, KVBeta, X, m, k){
  log.sigma2.xi <- log(sigma2.xi)
  
  beta.hat <- KVBeta %*% t(X) %*% log.sigma2.xi
  sigma2.hat <- sum((log.sigma2.xi - X %*% beta.hat)^2) / (m - k)
  
  nu <- m - k
  sigma2 <- nu * sigma2.hat * rchisq(1, nu, ncp = 0)^-1
  betas <- rmvnorm(n = 1, mean = beta.hat, sigma = sigma2 * KVBeta)
  
  return(list(betas = betas,sigma2 = sigma2))
}

#######################################################
#The main function
#######################################################
MCMC_SVR_q0_sub4 <- function(Y, X, obs.idx.m, group1.idx.m,  times,  
                             sigma2.eps, sigma2.zeta, sigma2.xi, sigma2.U0,
                             betas, sigma2,
                             n.iter, burn.in, thin){
#Constants
  deltas <- as.matrix(diff(times))  #interval between adjancnet times 
  n <- dim(Y)[1]
  m <- dim(Y)[2]
  m1 <- sum(group1.idx.m)    # number of group 1
  m2 <- m - m1               # number of group 2
  k <- dim(X)[2]             # number of covariates
  deltas.m <- cbind(deltas^3 / 3, deltas^2 / 2, deltas) # matrix of deltas
  
  KGG <- array(0, dim=c(2, 2, n))  # G matrix in the state equation of SVM.WN of function M
  Kinv.W <- array(0, dim=c(2, 2, n)) # inverse of the W matrix in the state equation of SVM.WN
  for(j in 1 : n){
    KGG[, , j] <- matrix(c(1, deltas[j], 0, 1), 2, 2, byrow = T)
    W.i <- matrix(c(deltas[j]^3 / 3, deltas[j]^2 / 2, 
                    deltas[j]^2 / 2, deltas[j]), 2, 2, byrow = T)
    Kinv.W[, , j] <- solve(W.i) 
  }
  
  KVBeta <- solve(t(X) %*% X)
  
  Kobs.idx.m <- vector("list", m) # idx of observations for each subject
  for(i in 1:m){
    Kobs.idx.m[[i]] <- which(obs.idx.m[, i] == 1)
  }
  Kns <- unlist(lapply(Kobs.idx.m, length)) # number of observations for each subject 
  
  Kgroup1.idx.m <- which(group1.idx.m == 1) # idx of subject from gruoup 1
  Kgroup2.idx.m <- which(group1.idx.m == 0) # idx of subject from gruoup 2
  
  Kdeltas.i <- vector("list",m) # time intervals for each subject.
  for(i in 1:m){
   Kdeltas.i[[i]] <- diff(times[c(1,Kobs.idx.m[[i]]+1)]) 
  } 
  
  Kn.obs <- sum(obs.idx.m)   

#The MCMC results storage
  sims.M <- matrix(NA, (n.iter - burn.in) / thin, 2 * n)
  sims.N <- matrix(NA, (n.iter - burn.in) / thin, 2 * n)
  sims.init <- matrix(NA, (n.iter - burn.in) / thin, 2 * 2 + m)
  sims.para <- matrix(NA, (n.iter - burn.in) / thin, 1 + 2 + m + 1 + k + 1)

  dimnames(sims.para) <- list(NULL, c("sigma2.eps", "sigma2.zeta1", "sigma2.zeta2",
                                      paste("sigma2.xi_", 1 : m, sep = ""),
                                      "sigma2.U0", paste("beta_", 1 : k, sep = ""), "sigma2"))

  sims.U <-vector("list", m)
  for(i in 1:m){
    sims.U[[i]] <- matrix(NA, (n.iter - burn.in) / thin, Kns[[i]])
    #dimnames(sims.U[[i]]) <- list(NULL, paste("U_", 1 : Kns[[i]], sep = ""))
  }  
  
  sigma2.zeta.tmp <- matrix(-9999, n, 2)  #specify some temporary matrix
  theta.sample <- matrix(-9999, n + 1, m + 2 + 2)
  sigma2.eps.tmp <- rep(0, m)
  
  #specify SSMs
  #for U_i(t): SPM.WN
  SPM.WN<-dlm(FF  = 1,
              V   = sigma2.eps, 
              GG  = 1,
              W   = -9,
              m0  = 0,
              C0  = 10^4, 
              JW  = 1,
              X   = sigma2.xi[1] * deltas 
             )
            
  #for M1(t): SVM.RepeatedIW1
  SVM.RepeatedIW1 <- dlmRandom(m1, 2)

  FF(SVM.RepeatedIW1) <- cbind(rep(-9, m1), 0)
  JFF(SVM.RepeatedIW1) <- cbind(1 : m1 + 4, 0)
  V(SVM.RepeatedIW1) <- diag(99, m1)
  JV(SVM.RepeatedIW1) <- diag(1 : m1 + m1 + 4)

  GG(SVM.RepeatedIW1) <- matrix(c(1, -9, 0, 1), 2, 2, byrow = T)
  JGG(SVM.RepeatedIW1) <- matrix(c(0, 4, 0, 0), 2, 2, byrow = T)
  W(SVM.RepeatedIW1) <- matrix(c(99, 99, 99, 99), 2, 2, byrow = T)
  JW(SVM.RepeatedIW1) <- matrix(c(1, 2, 2, 3), 2, 2, byrow = T)

  m0(SVM.RepeatedIW1) <- c(0, 0)
  C0(SVM.RepeatedIW1) <- diag(c(10^4, 10^4))
  
  #for M2(t): SVM.RepeatedIW2
  SVM.RepeatedIW2 <- dlmRandom(m2, 2)

  FF(SVM.RepeatedIW2) <- cbind(rep(-9, m2), 0)
  JFF(SVM.RepeatedIW2) <- cbind(1 : m2 + 4, 0)
  V(SVM.RepeatedIW2) <- diag(99, m2)
  JV(SVM.RepeatedIW2) <- diag(1 : m2 + m2 + 4)

  GG(SVM.RepeatedIW2) <- matrix(c(1, -9, 0, 1), 2, 2, byrow = T)
  JGG(SVM.RepeatedIW2) <- matrix(c(0, 4, 0, 0), 2, 2, byrow = T)
  W(SVM.RepeatedIW2) <- matrix(c(99, 99, 99, 99), 2, 2, byrow = T)
  JW(SVM.RepeatedIW2) <- matrix(c(1, 2, 2, 3), 2, 2, byrow = T)

  m0(SVM.RepeatedIW2) <- c(0, 0)
  C0(SVM.RepeatedIW2) <- diag(c(10^4, 10^4))
  
  theta.sample[-1, 1] = rowMeans(Y[, Kgroup1.idx.m])   # group 1 mean M1
  theta.sample[-1, 3] = rowMeans(Y[, Kgroup2.idx.m])   # group 2 mean M2
  
  for(iter in 1 : n.iter){
    
    #print(iter)
    
    #Sample U_i
    for(i in 1 : m){
      if(group1.idx.m[i] == 1){ # for the group 1
        theta.sample[c(1, Kobs.idx.m[[i]] + 1), i + 4] = 
        SampleW(Y.i = Y[Kobs.idx.m[[i]], i] - theta.sample[Kobs.idx.m[[i]]+1, 1],
                n.i = Kns[[i]],
                deltas.i = Kdeltas.i[[i]],
                SPM.WN, sigma2.eps, sigma2.xi[i], sigma2.U0)  
      }
      else{
        theta.sample[c(1, Kobs.idx.m[[i]] + 1), i + 4] = 
        SampleW(Y.i = Y[Kobs.idx.m[[i]], i] - theta.sample[Kobs.idx.m[[i]]+1, 3],
                 n.i = Kns[[i]],
                deltas.i = Kdeltas.i[[i]],
                SPM.WN, sigma2.eps, sigma2.xi[i], sigma2.U0)
      } 
    }
    
    # Sample M1, M2, N1 and N2
    theta.sample[, 1 : 2] =
    SampleRepeatedIW(Y[, Kgroup1.idx.m] - theta.sample[-1, Kgroup1.idx.m + 4], n,
    obs.idx.m[, Kgroup1.idx.m], deltas.m, SVM.RepeatedIW1, sigma2.eps, sigma2.zeta[1])
    
    theta.sample[, 3 : 4] =
    SampleRepeatedIW(Y[, Kgroup2.idx.m] - theta.sample[-1, Kgroup2.idx.m + 4], n, 
    obs.idx.m[, Kgroup2.idx.m], deltas.m, SVM.RepeatedIW2, sigma2.eps, sigma2.zeta[2])
  
    
    
    for(i in 1:m){
      if(group1.idx.m[i] == 1){
        sigma2.eps.tmp[i] <- sum((Y[Kobs.idx.m[[i]], i]-
                                  theta.sample[Kobs.idx.m[[i]]+1, 1]- 
                                  theta.sample[Kobs.idx.m[[i]]+1, i + 4])^2)
      }
      else{
        sigma2.eps.tmp[i] <- sum((Y[Kobs.idx.m[[i]], i]-
                                  theta.sample[Kobs.idx.m[[i]]+1, 3]- 
                                  theta.sample[Kobs.idx.m[[i]]+1, i + 4])^2)
      }
    }
        
    #sample sigma2.eps
    sigma2.eps <- rgamma(1, shape = Kn.obs / 2 + 0.3, 
                            rate  = sum(sigma2.eps.tmp) / 2 + 0.3)^-1
    
    #sample sigma2.zeta
    for(j in 1 : n){                     
      xi.MN1 <- theta.sample[j + 1, 1 : 2] -  KGG[, , j] %*% theta.sample[j, 1 : 2] 
      xi.MN2 <- theta.sample[j + 1, 3 : 4] -  KGG[, , j] %*% theta.sample[j, 3 : 4] 
      
      sigma2.zeta.tmp[j, 1] <- t(xi.MN1) %*% Kinv.W[, , j] %*% xi.MN1
      sigma2.zeta.tmp[j, 2] <- t(xi.MN2) %*% Kinv.W[, , j] %*% xi.MN2
    }

    sigma2.zeta[1] <- rgamma(1, shape = n + 0.3, 
                                rate  = sum(sigma2.zeta.tmp[, 1]) / 2 + 0.3)^-1
    sigma2.zeta[2] <- rgamma(1, shape = n + 0.3, 
                                rate  = sum(sigma2.zeta.tmp[, 2]) / 2 + 0.3)^-1   
    
    #sample sigma2.xi, betas and sigma2
    mus <- X %*% t(betas)                                                                          
    sigma2.xi <- SampleSigma2Xi(theta.sample, Kobs.idx.m, Kdeltas.i, Kns, sigma2.xi, times,  mus, sigma2, n, m)
    
    betas.sigma2 <- SampleBetaSigma2(sigma2.xi, KVBeta, X, m, k)
    betas <- betas.sigma2$betas
    sigma2 <- betas.sigma2$sigma2                                 
    
    #sample sigma2.U0 and sigma2.V0                                                                           
    sigma2.U0 <- rgamma(1, shape = m / 2 + 0.3, 
                           rate = sum(theta.sample[1, 1 : m + 4 ]^2) / 2 + 0.3)^-1
    #save results:
    if(iter > burn.in & iter %% thin == 0){
      sims.M[(iter - burn.in) / thin, ] <- c(theta.sample[-1, 1], theta.sample[-1, 3])
      sims.N[(iter - burn.in) / thin, ] <- c(theta.sample[-1, 2], theta.sample[-1, 4])
      #sims.U[(iter - burn.in) / thin, ] <- as.vector(theta.sample[-1, 1 : m + 4])
      for(i in 1:m){
        sims.U[[i]][(iter - burn.in) / thin, ] <- theta.sample[Kobs.idx.m[[i]]+1, 4 + i]
      }
      sims.init[(iter - burn.in) / thin, ] <- theta.sample[1, ]
      sims.para[(iter - burn.in) / thin, ] <- c(sigma2.eps, sigma2.zeta, 
                                                sigma2.xi, sigma2.U0,
                                                betas, sigma2)
      }
    
  }
  return(list(sims.M = sims.M, sims.N = sims.N, sims.U = sims.U, 
              sims.init = sims.init, sims.para = sims.para))

}

print("MCMC_SVR_q0_sub4 and auxiliary functions loaded")

