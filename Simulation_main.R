##### Simulated Data #####

rm(list = ls())

library(LearnBayes)
library(mvtnorm)
library(coda)  
library(ggplot2)
library(ggpubr)
library(Matrix)
library(MASS)
library(plotly)
library(SoDA)
library(stringr)
library(gmp)
library(dlm)

set.seed(1997)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))





##### Simulate from Multivariate Ornstein Uhlenbeck #####
# n is the number of observations
# mu is the mean of the gaussian
# lambda is the covariance matrix of the gaussian
# d is the number of dimensions

MultivariateOU <- function(n,mu,lambda,gamma,d=3){
  x <- matrix(nrow = n, ncol = d)
  x[1,]=mvrnorm(n = 1, mu, lambda)
  for (i in 2:n) {
    x[i,] = mvrnorm(n = 1, mu+gamma*(x[i-1]-mu), (1-gamma^2)*lambda)
  }
  return(x)
}





##### Data Simulation Parameters for changes in mean #####

N1=50
mu1=c(0,0,0)
lambda1=matrix(c(1,0,0,0,1,0,0,0,1), nrow = 3, ncol = 3)
gamma1=0.4

N2=35
mu2=c(4,4,4)
lambda2=matrix(c(1,0,0,0,1,0,0,0,1), nrow = 3, ncol = 3)
gamma2=0.4

N3=65
mu3=c(4,4,-4)
lambda3=matrix(c(1,0,0,0,1,0,0,0,1), nrow = 3, ncol = 3)
gamma3=0.4






##### Data Simulation Parameters for changes in variance #####

N1=50
mu1=c(0,0,0)
lambda1=matrix(c(1,0,0,0,1,0,0,0,1), nrow = 3, ncol = 3)
gamma1=0.4

N2=35
mu2=c(0,0,0)
lambda2=matrix(c(5,0,0,0,5,0,0,0,5), nrow = 3, ncol = 3)
gamma2=0.4

N3=65
mu3=c(0,0,0)
lambda3=matrix(c(5,4,4,4,5,4,4,4,5), nrow = 3, ncol = 3)
gamma3=0.4






##### Data Simulation Parameters for changes in mean,variance and correlation #####

N1=50
mu1=c(0,0,0)
lambda1=matrix(c(1,0,0,0,1,0,0,0,1), nrow = 3, ncol = 3)
gamma1=0.4

N2=35
mu2=c(3,3,3)
lambda2=matrix(c(3,0,0,0,3,0,0,0,3), nrow = 3, ncol = 3)
gamma2=0.7

N3=65
mu3=c(0,0,0)
lambda3=matrix(c(2,1,1,1,2,1,1,1,2), nrow = 3, ncol = 3)
gamma3=0.9





##### Real Change Points #####
true_change_points = cumsum(c(N1,N2))+1


##### Data Simulation Generation #####
y1=MultivariateOU(N1, mu1, lambda1, gamma1)
y2=MultivariateOU(N2, mu2, lambda2, gamma2)
y3=MultivariateOU(N3, mu3, lambda3, gamma3)
y = rbind(y1,y2,y3)





##### load functions  and hyper parameters ####
source('simulation_parameters.R')
source('functions.R')


##### call for algorithm #####
solution=multi_split_merge_MCMC(start_k,start_n,q,nrep,burnin,start_sigma,start_theta,start_gamma)

post_prob_change_points = posterior_prob_change_points(solution)





##### Estimates of posterior #####
solution_rhon=solution$n
solution_k=solution$k
solution_gamma=solution$gamma
solution_theta=solution$theta
solution_sigma=solution$sigma





##### Gamma diagnostic #####
mcmc_gamma_post <- mcmc(solution_gamma, start = burnin + 1, end = burnin+nrep)
summary(mcmc_gamma_post) 
plot(mcmc_gamma_post)
acfplot(mcmc_gamma_post)
cumuplot(mcmc_gamma_post)
effectiveSize(mcmc_gamma_post)
plot(ts(solution_gamma[17000:20000]))


##### Theta diagnostic #####
mcmc_theta_post <- mcmc(solution_theta, start = burnin + 1, end = burnin+nrep)
summary(mcmc_theta_post) 
plot(mcmc_theta_post)
acfplot(mcmc_theta_post)
cumuplot(mcmc_theta_post)
effectiveSize(mcmc_theta_post)
plot(ts(solution_theta[17000:20000]))


##### Sigma diagnostic #####
mcmc_sigma_post <- mcmc(solution_sigma, start = burnin + 1, end = burnin+nrep)
summary(mcmc_sigma_post) 
plot(mcmc_sigma_post)
acfplot(mcmc_sigma_post)
cumuplot(mcmc_sigma_post)
effectiveSize(mcmc_sigma_post)
plot(ts(solution_sigma[17000:20000]))





##### Estimates for change points #####

# Order the results
adata_ordered = pretty_results(solution_rhon)

# Computation of the mode
rhon_mode=adata_ordered[1,1][[1]]
kcoherent=length(adata_ordered[1,1][[1]])
occurrences_of_mode=adata_ordered[1,2]
prob_rhon=occurrences_of_mode/nrep

# Computation of change points
sum_prec=1
change_points=numeric(length(rhon_mode)-1)
for(i in 1:(length(rhon_mode)-1)){
  sum_prec=sum_prec+rhon_mode[i]
  change_points[i]=sum_prec
}
change_points





##### Estimates for k #####

# K mode 
kmode=Mode(solution_k)
kmode
kmode_occurrences=sort(table(solution_k),decreasing=TRUE)[1]
prob_k=kmode_occurrences/nrep

# K coherent with change points estimate
kcoherent=length(adata_ordered[1,1][[1]])
kcoherent





##### Plot of change points #####
layout(matrix(c(1,1,2,3),2,2,byrow = TRUE))
 
for(i in 1:d){
  plot(y[,i],col=i,type = 'l')
  for(i in 1:length(change_points)){
    abline(v=change_points[i], col="red")
  }
  for(i in 1:length(true_change_points)){
    abline(v=true_change_points[i], col = "blue")
  }
}


change_points
true_change_points
