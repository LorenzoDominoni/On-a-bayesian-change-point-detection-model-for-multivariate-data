##### Finance Data Application #####

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

set.seed(16091997)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

ENEL=read.csv("ENEL.MI_2020-01-01_to_2020-11-29.csv")
ENI=read.csv("ENI.MI_2020-01-01_to_2020-11-29.csv")
FERRARI=read.csv("RACE.MI_2020-01-01_to_2020-11-29.csv")
STMcompleta=read.csv("STM.MI_2020-01-01_to_2020-11-29.csv")
INTESA=read.csv("ISP.MI_2020-01-01_to_2020-11-29.csv")
FTSEMIBcompleta=read.csv("FTSEMIB.MI_2020-01-01_to_2020-11-29.csv")

#FTSEMIB has 2 observations more than the others
FTSEMIBcompleta=FTSEMIBcompleta[-140,]
FTSEMIB=FTSEMIBcompleta[-183,]

#STM has an observation more than the others (20 july)
STM=STMcompleta[-140,]





##### Creation of closing price dataset for companies #####
yraw <- data.frame("dummy"=1:nrow(ENEL))

yraw$ENEL=ENEL$Close
yraw$INTESA=INTESA$Close
yraw$FERRARI=FERRARI$Close
yraw$ENI=ENI$Close
yraw$STM=STM$Close
yraw=yraw[,-1]  


##### Creation of closing price dataset for FTSEMIB #####
yraw_tot <- data.frame("dummy"=1:nrow(FTSEMIB))
yraw_tot$FTSEMIB=FTSEMIB$Close
yraw_tot=yraw_tot[,-1]  





# Plot data
layout(matrix(c(1,2,3,4,5,6),2,3))
plot(yraw_tot,col='red',type = 'l',ylab='Price',xlab='Time', main='FTSEMIB')
plot(yraw[,1],col='blue',type = 'l',ylab='Price',xlab='Time', main='ENEL')
plot(yraw[,2],col='blue',type = 'l',ylab='Price',xlab='Time', main='INTESA')
plot(yraw[,3],col='blue',type = 'l',ylab='Price',xlab='Time', main='FERRARI')
plot(yraw[,4],col='blue',type = 'l',ylab='Price',xlab='Time', main='ENI')
plot(yraw[,5],col='blue',type = 'l',ylab='Price',xlab='Time', main='STM')


# Scaling and preparing for the algorithm
yraw=as.matrix(yraw)
y=scale(yraw)
yraw_tot=as.matrix(yraw_tot)
y_tot=scale(yraw_tot)





##### Algorithm application ####


##### load functions and hyper parameters #####
source('finance_parameters.R')
source('functions.R')


##### call for algorithm #####
solution=multi_split_merge_MCMC(start_k,start_n,q,nrep,burnin,start_sigma,start_theta,start_gamma,debug)

#solution=readRDS("Finance_Multivariate_Solution.rds") #Load the solution when already computed

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

adata_ordered = pretty_results(solution_rhon)


# Mode (1st estimate) with its posterior probability
rhon_mode=adata_ordered[1,1][[1]]
occurrences_of_mode=adata_ordered[1,2]
prob_rhon=occurrences_of_mode/nrep

sum_prec=1
change_points1=numeric(length(rhon_mode)-1)
for(i in 1:(length(rhon_mode)-1)){
  sum_prec=sum_prec+rhon_mode[i]
  change_points1[i]=sum_prec
}

change_points1
# 24  45  63 109 183 214
change_dates1=ENEL[change_points1,"Date"]
change_dates1
# 2020-02-04 2020-03-04 2020-03-30 2020-06-05 2020-09-21 2020-11-03


# Time steps with a posterior probability higher than a treshold 0.2 (2nd estimate)
treshold=0.2
change_points2=find_best(post_prob_change_points, treshold, 10)
change_points2

change_points2
# 23  44  62 108 182 213
change_dates2=ENEL[change_points2,"Date"]
change_dates2
# 2020-02-03 2020-03-03 2020-03-27 2020-06-04 2020-09-17 2020-11-02





##### Estimates for the number of different regimes #####

# K coherent with first estimate
kcoherent1=length(adata_ordered[1,1][[1]])
kcoherent1
# 7


# K coherent with second estimate
kcoherent2=length(change_points2)+1
kcoherent2
# 7


# K mode 
kmode=Mode(solution_k)
kmode
# 7





##### Plot of change points 2 with posterior probability #####

layout(matrix(c(1,2),2,1))
plot(post_prob_change_points,type = 'l',ylab='Probability',xlab='Time', main='Posterior probability of Change Points')
abline(h=treshold, col="orange")
plot(yraw_tot,col='red',type = 'l',ylab='Price',xlab='Time', main='FTSEMIB Change Points')
for(i in 1:length(change_points2)){
  abline(v=change_points2[i], col="black")
}





##### Plot of change points 2 for all companies #####

layout(matrix(c(1,2,3,4,5,6),2,3))

plot(yraw_tot,col='red',type = 'l',ylab='Price',xlab='Time', main='FTSEMIB')
for(i in 1:length(change_points2)){
  abline(v=change_points2[i], col="black")
}
plot(yraw[,1],col='blue',type = 'l',ylab='Price',xlab='Time', main='ENEL')
for(i in 1:length(change_points2)){
  abline(v=change_points2[i], col="black")
}
plot(yraw[,2],col='blue',type = 'l',ylab='Price',xlab='Time', main='INTESA')
for(i in 1:length(change_points2)){
  abline(v=change_points2[i], col="black")
}
plot(yraw[,3],col='blue',type = 'l',ylab='Price',xlab='Time', main='FERRARI')
for(i in 1:length(change_points2)){
  abline(v=change_points2[i], col="black")
}
plot(yraw[,4],col='blue',type = 'l',ylab='Price',xlab='Time', main='ENI')
for(i in 1:length(change_points2)){
  abline(v=change_points2[i], col="black")
}
plot(yraw[,5],col='blue',type = 'l',ylab='Price',xlab='Time', main='STM')
for(i in 1:length(change_points2)){
  abline(v=change_points2[i], col="black")
}





# Focus Intesa Ferrari Enel Eni
layout(matrix(c(1,2,3,4),2,2))

plot(yraw[,1],col='blue',type = 'l',ylab='Price',xlab='Time', main='ENEL')
for(i in 1:length(change_points2)){
  abline(v=change_points2[i], col="black")
}
abline(v=change_points2[2], col="orange")
abline(v=change_points2[4], col="red")
abline(v=change_points2[6], col="green")
plot(yraw[,2],col='blue',type = 'l',ylab='Price',xlab='Time', main='INTESA')
for(i in 1:length(change_points2)){
  abline(v=change_points2[i], col="black")
}
abline(v=change_points2[2], col="orange")
abline(v=change_points2[4], col="red")
abline(v=change_points2[6], col="green")
plot(yraw[,3],col='blue',type = 'l',ylab='Price',xlab='Time', main='FERRARI')
for(i in 1:length(change_points2)){
  abline(v=change_points2[i], col="black")
}
abline(v=change_points2[2], col="orange")
abline(v=change_points2[4], col="red")
abline(v=change_points2[6], col="green")
plot(yraw[,4],col='blue',type = 'l',ylab='Price',xlab='Time', main='ENI')
for(i in 1:length(change_points2)){
  abline(v=change_points2[i], col="black")
}
abline(v=change_points2[2], col="orange")
abline(v=change_points2[4], col="red")
abline(v=change_points2[6], col="green")





##### Comparison of mode and highest probability change points #####
layout(matrix(c(1,2),2,1))
plot(yraw_tot,col='black',type = 'l',ylab='Price',xlab='Time', main='FTSEMIB mode estimate')
for(i in 1:length(change_points1)){
  abline(v=change_points1[i], col="red")
}
plot(yraw_tot,col='black',type = 'l',ylab='Price',xlab='Time', main='FTSEMIB highest probability estimate')
for(i in 1:length(change_points2)){
  abline(v=change_points2[i], col="blue")
}





##### Plot of posterior probabilities of change points #####
layout(matrix(c(1),1,1))
plot(post_prob_change_points,type = 'l',ylab='Probability',xlab='Time', main='Posterior probability of Change Points')
abline(h=0.6, col="blue")
abline(h=0.2, col="red")
legend("topleft", legend=c("High Treshold", "Low Treshold"), col=c("blue", "red"), lty=1, cex=0.8)

# Save Solution
saveRDS(solution, file="Finance_Multivariate_Solution.rds")





##### Comparison with FTSE MIB #####
set.seed(16091997)
y=y_tot


##### load functions and parameters ####
source('finance_parameters.R')
source('functions.R')


##### call for algorithm #####
solution=multi_split_merge_MCMC(start_k,start_n,q,nrep,burnin,start_sigma,start_theta,start_gamma,debug)

#solution=readRDS("Finance_Univariate_Solution.rds") #Load the solution when already computed

post_prob_change_points_uni = posterior_prob_change_points(solution)





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

adata_ordered = pretty_results(solution_rhon)


# Time steps with a posterior probability higher than a treshold 0.2
treshold_uni=0.25
change_points_uni=find_best(post_prob_change_points_uni, treshold_uni, 7)
change_points_uni

change_points_uni
# 47 58 121 182 200 217
change_dates_uni=ENEL[change_points_uni,"Date"]
change_dates_uni
# 2020-03-06 2020-03-23 2020-06-23 2020-09-17 2020-10-14 2020-11-06





##### Plot of posterior probabilities of change points #####
layout(matrix(c(1,2),2,1))

plot(post_prob_change_points_uni,col='black',type = 'l',ylab='Probability',xlab='Time', main='Posterior probability of Univariate Change Points')
abline(h=treshold_uni, col="red")

plot(post_prob_change_points,col='black',type = 'l',ylab='Probability',xlab='Time', main='Posterior probability of Multivariate Change Points')
abline(h=treshold_uni, col="blue")





##### Plot of change points #####
layout(matrix(c(1,2),2,1))

plot(yraw_tot,col='black',type = 'l',ylab='Cases',xlab='Time', main='Univariate Change Points')
for(i in 1:length(change_points_uni)){
  abline(v=change_points_uni[i], col="red")
}
plot(yraw_tot,col='black',type = 'l',ylab='Cases',xlab='Time', main='Multivariate Change Points')
for(i in 1:length(change_points2)){
  abline(v=change_points2[i], col="blue")
}


# Save Solution
saveRDS(solution, file="Finance_Univariate_Solution.rds")

