SimuW <- function(n, deltas, sigma2, init){
#Simulate Wiener process only
  theta <- matrix(NA, n+1, 1)
  eps <- rnorm(n, 0, sqrt(sigma2 * deltas))
  theta[1] <- init
  theta[-1] <- cumsum(eps) + init
 
  return(theta)
  }

SampleW <- function(Y.i, n.i, deltas.i, SPM.WN, sigma2.eps, sigma2.xi.i, sigma2.U0){
#Sample latent process W from the SPM.WN  
  V(SPM.WN) <- sigma2.eps
  X(SPM.WN) <- sigma2.xi.i * deltas.i
  C0(SPM.WN) <- sigma2.U0
  
  #slower:              
  #theta = dlmBSample(dlmFilter(Y.i, SPM.WN)[-1])
  
  #faster:
  theta.all <- SimuW(n.i, deltas.i, sigma2.xi.i, rnorm(1, 0, sqrt(sigma2.U0)))

  theta.hat <- dlmSmooth(Y.i, SPM.WN)$s
  theta.star <- dlmSmooth(theta.all[-1] + rnorm(n.i, 0, sqrt(sigma2.eps)), SPM.WN)$s
  
  theta = theta.hat + theta.all - theta.star
  
  return(theta)    
}