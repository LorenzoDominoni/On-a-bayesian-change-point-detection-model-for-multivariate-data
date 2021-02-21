##### Covid Data Application #####

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

df=read.csv("dpc-covid19-ita-regioni.csv")





##### Data for the daily new positives in regions #####

# Create dataset for regions
data <- data.frame("dummy"=1:length(unique(df$data)))
for(i in 1:length(unique(df$denominazione_regione))){
  region=unique(df$denominazione_regione)[i]
  data[[as.character(region)]]=df[which(df$denominazione_regione==region),"nuovi_positivi"]
}

# Set north, center, south
north=c("Emilia-Romagna","Friuli Venezia Giulia","Liguria","Lombardia","Piemonte","P.A. Bolzano","P.A. Trento","Valle d'Aosta","Veneto")
center=c("Lazio","Marche","Toscana","Umbria")
south=c("Abruzzo","Basilicata","Calabria","Campania","Molise","Puglia","Sardegna","Sicilia")

# Group data by macroregion
yraw <- data.frame("dummy"=1:length(unique(df$data)))
yraw$north=rowSums(data[north])
yraw$center=rowSums(data[center])
yraw$south=rowSums(data[south])
yraw=yraw[,-1]  

# Goup all data together
y_all=rowSums(yraw)
y_all=as.matrix(y_all)





# Plot data together
layout(matrix(c(1),1,1))
plot(yraw[,1],col='red',type = 'l',ylab='Cases',xlab='Time', main='Covid daily cases')
lines(yraw[,2],col='blue',type = 'l',ylab='Center',xlab='Time')
lines(yraw[,3],col='green',type = 'l',ylab='South',xlab='Time')
legend("topleft", legend=c("North", "Center", "South"), col=c("red", "blue", "green"), lty=1, cex=0.8)

# Plot data divided
layout(matrix(c(1,2,3,4),2,2))
plot(y_all,col='red',type = 'l',ylab='Cases',xlab='Time', main='Italy')
plot(yraw[,1],col='blue',type = 'l',ylab='Cases',xlab='Time', main='North')
plot(yraw[,2],col='blue',type = 'l',ylab='Cases',xlab='Time', main='Center')
plot(yraw[,3],col='blue',type = 'l',ylab='Cases',xlab='Time', main='South')

# Scaling and preparing for the algorithm
yraw=as.matrix(yraw)
y=scale(yraw)





##### Algorithm application ####


##### load functions and parameters ####
source('covid_parameters.R')
source('functions.R')


##### call for algorithm #####
solution=multi_split_merge_MCMC(start_k,start_n,q,nrep,burnin,start_sigma,start_theta,start_gamma,debug)

solution=readRDS("Covid_Multivariate_Solution.rds") #Load the solution when already computed

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
# 15  43  87 196 237
change_dates1=unique(df$data)[change_points1]
change_dates1
# 2020-03-09 2020-04-06 2020-05-20 2020-09-06 2020-10-17


# Time steps with a posterior probability higher than a high treshold 0.35 (2nd estimate)
treshold1=0.35
change_points2=find_best(post_prob_change_points, treshold1, 5)
change_points2

change_points2
# 14  86 196 236
change_dates2=unique(df$data)[change_points2]
change_dates2
# 2020-03-08 2020-05-19 2020-09-06 2020-10-16


# Time steps with a posterior probability higher than a low treshold 0.2 (3rd estimate)
treshold2=0.2
change_points3=find_best(post_prob_change_points, treshold2, 5)
change_points3

change_points3
# 14  42  86 196 236 279
change_dates3=unique(df$data)[change_points3]
change_dates3
# 2020-03-08 2020-04-05 2020-05-19 2020-09-06 2020-10-16 2020-11-28





##### Estimates for the number of different regimes #####

# K coherent with first estimate
kcoherent1=length(adata_ordered[1,1][[1]])
kcoherent1
# 6


# K coherent with second estimate
kcoherent2=length(change_points2)+1
kcoherent2
# 5


# K coherent with third estimate
kcoherent3=length(change_points3)+1
kcoherent3
# 7


# K mode 
kmode=Mode(solution_k)
kmode
# 7





##### Plot of change points 1 #####
layout(matrix(c(1,2,3),3,1))
plot(yraw[,1],col='red',type = 'l',ylab='Cases',xlab='Time', main='North')
for(i in 1:length(change_points1)){
  abline(v=change_points1[i], col="black")
}
plot(yraw[,2],col='blue',type = 'l',ylab='Cases',xlab='Time', main='Center')
for(i in 1:length(change_points1)){
  abline(v=change_points1[i], col="black")
}
plot(yraw[,3],col='green',type = 'l',ylab='Cases',xlab='Time', main='South')
for(i in 1:length(change_points1)){
  abline(v=change_points1[i], col="black")
}





# Plot change points 1 together
layout(matrix(c(1),1,1))
plot(yraw[,1],col='red',type = 'l',ylab='Cases',xlab='Time', main='Covid daily cases')
lines(yraw[,2],col='blue',type = 'l',ylab='Center',xlab='Time')
lines(yraw[,3],col='green',type = 'l',ylab='South',xlab='Time')
for(i in 1:length(change_points1)){
  abline(v=change_points1[i], col="grey")
}
abline(v=change_points1[5], col="black")
legend("topleft", legend=c("North", "Center", "South"), col=c("red", "blue", "green"), lty=1, cex=0.8)





##### Plot of change points 2 #####
layout(matrix(c(1,2,3),3,1))
plot(yraw[,1],col='red',type = 'l',ylab='Cases',xlab='Time', main='North')
for(i in 1:length(change_points2)){
  abline(v=change_points2[i], col="black")
}
plot(yraw[,2],col='blue',type = 'l',ylab='Cases',xlab='Time', main='Center')
for(i in 1:length(change_points2)){
  abline(v=change_points2[i], col="black")
}
plot(yraw[,3],col='green',type = 'l',ylab='Cases',xlab='Time', main='South')
for(i in 1:length(change_points2)){
  abline(v=change_points2[i], col="black")
}





##### Plot of change points 3 #####
layout(matrix(c(1,2,3),3,1))
plot(yraw[,1],col='red',type = 'l',ylab='Cases',xlab='Time', main='North')
for(i in 1:length(change_points3)){
  abline(v=change_points3[i], col="black")
}
plot(yraw[,2],col='blue',type = 'l',ylab='Cases',xlab='Time', main='Center')
for(i in 1:length(change_points3)){
  abline(v=change_points3[i], col="black")
}
plot(yraw[,3],col='green',type = 'l',ylab='Cases',xlab='Time', main='South')
for(i in 1:length(change_points3)){
  abline(v=change_points3[i], col="black")
}




##### Plot of posterior probabilities of change points #####
layout(matrix(c(1),1,1))
plot(post_prob_change_points,type = 'l',ylab='Probability',xlab='Time', main='Posterior probability of Change Points')
abline(h=treshold1, col="blue")
abline(h=treshold2, col="red")
legend("topleft", legend=c("High Treshold", "Low Treshold"), col=c("blue", "red"), lty=1, cex=0.8)

# Save Solution
saveRDS(solution, file="Covid_Multivariate_Solution.rds")





##### Comparison with the univariate sum #####
set.seed(16091997)


##### Scaling and preparing for the algorithm #####
y=scale(y_all)


# Plot data
layout(matrix(c(1),1,1))
plot(y_all,col='blue',type = 'l',ylab='Italy',xlab='Time', main='Covid national daily cases')


##### load functions and parameters ####
source('covid_parameters.R')
source('functions.R')


##### call for algorithm #####
solution=multi_split_merge_MCMC(start_k,start_n,q,nrep,burnin,start_sigma,start_theta,start_gamma,debug)

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


# Time steps with a posterior probability higher than a medium treshold 0.3 
treshold_uni=0.3
change_points_uni=find_best(post_prob_change_points_uni, treshold_uni, 5)
change_points_uni

change_points_uni
# 15  85 184 220 226
change_dates_uni=unique(df$data)[change_points_uni]
change_dates_uni
# 2020-03-09 2020-05-18 2020-08-25 2020-09-30 2020-10-06


##### Estimates for the number of different regimes #####

# K coherent with univariate estimate
kcoherent_uni=length(adata_ordered[1,1][[1]])
kcoherent_uni
# 6





##### Plot of posterior probabilities of change points #####
layout(matrix(c(1,2),2,1))

plot(post_prob_change_points_uni,col='black',type = 'l',ylab='Probability',xlab='Time', main='Posterior probability of Univariate Change Points')
abline(h=treshold_uni, col="red")

plot(post_prob_change_points,col='black',type = 'l',ylab='Probability',xlab='Time', main='Posterior probability of Multivariate Change Points')
abline(h=treshold_uni, col="blue")





##### Plot of change points #####
layout(matrix(c(1,2),2,1))

plot(y_all,col='black',type = 'l',ylab='Italy',xlab='Time', main='Univariate Change Points')
for(i in 1:length(change_points_uni)){
  abline(v=change_points_uni[i], col="red")
}
plot(y_all,col='black',type = 'l',ylab='Italy',xlab='Time', main='Multivariate Change Points')
for(i in 1:length(change_points1)){
  abline(v=change_points1[i], col="blue")
}


# Save Solution
saveRDS(solution, file="Covid_Univariate_Solution.rds")

  