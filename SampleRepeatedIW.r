SampleRepeatedIW <- function(Y.i, n,  obs.idx.m, deltas.m, SVM.RepeatedIW, sigma2.eps, sigma2.zeta){
  
  X(SVM.RepeatedIW) <- cbind(sigma2.zeta * deltas.m, deltas.m[, 3], obs.idx.m,  sigma2.eps * obs.idx.m)
               
  #may be improved by the simulation smoother
  theta = dlmBSample(dlmFilter(Y.i, SVM.RepeatedIW)[-1])
  
  return(theta)    
}