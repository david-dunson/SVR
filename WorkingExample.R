library(dlm)
source("MCMC_SVR_q0_sub4.R")  
source("SampleW.R")
source("SampleRepeatedIW.R")


#step1: simulate the case II at paper
set.seed(95545)  

m <- 100
n <-20
deltas <- rep(0.2,n)
times <- cumsum(deltas)


sigma2.eps <- 1

M_U1 <- 10*(times+sin(times))
M_U2 <- 10*(times+cos(times))

eigen1.1 <- 0.6*cos(times*pi/10)
eigne1.2<- 0.2*sin(times*pi/10)

eigen2.1 <- 0.5*cos(times*pi/10)
eigne2.2<- 0.3*sin(times*pi/10)

#simulate the Y and theta
Y <- matrix(NA, nrow = n, ncol = m)
M_U <- matrix(NA, nrow = n, ncol = m)

group1.idx.m <- rbinom(m, 1, 0.5)

for(i in 1:m)
{
  if(group1.idx.m[i] == 1){
    M_U[, i] <- M_U1 + 2*rnorm(1)*eigen1.1 + rnorm(1)*eigne1.2
    Y[, i] <- M_U[, i] + rnorm(n,1) 
  }else{
    M_U[, i] <- M_U2 + 2*rnorm(1)*eigen2.1 + rnorm(1)*eigne2.2
    Y[, i] <- M_U[, i] + rnorm(n,1)
    
  }
}

missing.idx.m <- matrix(rbinom(n * m, 1, 0.2), n, m) # create 20% missing
obs.idx.m <- 1 - missing.idx.m

#plot
# matplot(times,Y, type="l",lty=1, xlab="Time")
# lines(times, M_U1, type="l",lty=1, col="green", lwd=4)
# lines(times, M_U2, type="l",lty=1, col="yellow", lwd=4)

simu <- list(Y = Y, M_U = M_U,  times = times, obs.idx.m = obs.idx.m, group1.idx.m = group1.idx.m)

#step2:run MCMC
times <- simu$times
sims = MCMC_SVR_q0_sub4 (Y=simu$Y, X = as.matrix(rep(1,m)),  obs.idx.m = simu$obs.idx.m, group1.idx.m = simu$group1.idx.m,
                                                    times = c(0, simu$times),  
                                                    sigma2.eps = 1, sigma2.zeta = c(2,2), 
                                                    sigma2.xi = rep(1,m), sigma2.U0 = 4,
                                                    betas = 0, sigma2 = 1, 
                                                    n.iter =15, burn.in = 5, thin = 1)
                                                    #n.iter =15000, burn.in = 5000, thin = 10)


