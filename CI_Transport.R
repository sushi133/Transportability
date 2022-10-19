library(tictoc)
library(parallel)

sim<-function(i){
  
  #study
  n1<-5000
  a0<-2
  a1<-4
  u_x1<-1
  sd_x1<-2
  x1<-rnorm(n1,u_x1,sd_x1)
  x1_sd<-sd(x1)
  x1_mean<-mean(x1)
  y1<-rnorm(n1,a0+a1*x1,1)
  dat1<-as.data.frame(cbind(y1,x1))
  te1<-a0+a1*u_x1
  ste1<-mean(y1)
  
  #target
  n0<-1000
  b0<-2
  b1<-4
  u_x0<--1
  sd_x0<-4
  x0<-rnorm(n0,u_x0,sd_x0)
  x0_sd<-sd(x0)
  x0_mean<-mean(x0)
  y0<-rnorm(n0,b0+b1*x0,1)
  dat0<-as.data.frame(cbind(y0,x0))
  te0<-b0+b1*u_x0
  ste0<-mean(y0)

  fun <- function(beta){
    f<-sum(exp(x1*beta)-(x0_mean*exp(x1*beta)/x1))
    return(f)
  }
  output <- optim(par = 0, fn = fun, lower=-Inf, upper=Inf, method="L-BFGS-B", hessian = T)
  #standard errors
  sd <- sqrt(diag(solve(output$hessian)))
  #95% CI
  beta_l<-output$par-qnorm(0.975)*sd
  beta_u<-output$par+qnorm(0.975)*sd
  
  funl <- function(beta){
    f<-sum(exp(x1*beta)-((x0_mean-qnorm(0.975)*x0_sd/sqrt(n0))*exp(x1*beta)/x1))
    return(f)
  }
  output <- optim(par = 0, fn = funl, lower=-Inf, upper=Inf, method="L-BFGS-B", hessian = T)
  #standard errors
  sd <- sqrt(diag(solve(output$hessian)))
  #95% CI
  betal_l<-output$par-qnorm(0.95)*sd
  betal_u<-output$par+qnorm(0.95)*sd
  
  funu <- function(beta){
    f<-sum(exp(x1*beta)-((x0_mean+qnorm(0.975)*x0_sd/sqrt(n0))*exp(x1*beta)/x1))
    return(f)
  }
  output <- optim(par = 0, fn = funu, lower=-Inf, upper=Inf, method="L-BFGS-B", hessian = T)
  #standard errors
  sd <- sqrt(diag(solve(output$hessian)))
  #95% CI
  betau_l<-output$par-qnorm(0.95)*sd
  betau_u<-output$par+qnorm(0.95)*sd
  
  te10l_l<-sum(dat1$y1*exp(betal_l*dat1$x1))/sum(exp(betal_l*dat1$x1))
  te10u_u<-sum(dat1$y1*exp(betau_u*dat1$x1))/sum(exp(betau_u*dat1$x1))
  te10_l<-sum(dat1$y1*exp(beta_l*dat1$x1))/sum(exp(beta_l*dat1$x1))
  te10_u<-sum(dat1$y1*exp(beta_u*dat1$x1))/sum(exp(beta_u*dat1$x1))
  
  cr<-c(I(te10_l<=te0 & te0<=te10_u),te10_l,te10_u,I(te10l_l<=te0 & te0<=te10u_u),te10l_l,te10u_u)
}

#time start
tic()
runs<-1000
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
save <- mclapply(1:runs, sim, mc.cores = 8, mc.set.seed = TRUE)
toc()
#time end
results<-round(apply(array(as.numeric(unlist(save)),dim = c(6, runs)),1,mean),3)
results



