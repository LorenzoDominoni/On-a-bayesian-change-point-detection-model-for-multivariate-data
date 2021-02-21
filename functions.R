##### Functions #####
# These are the functions that make the algorithm work and show its results
# All functions below are in log scale because of numerical issues (too large numbers)





##### Helper Functions #####
# Useful functions not already implemented in R

##### Logpochhammer Function #####
# Computes the log of x pochhammer m

logpochhammer = function(x,m){
  if(m==0){
    result=0
  }
  else{
    result=sum(log(x:(x+m-1)))
  }
  return(result)
}


##### Logmultigamma Function #####
# Computes the logarithm of the multivariate gamma function with parameter a in d dimesnsions

log_Multi_Gamma = function(a,d){
  log(pi)*(d*(d-1)/4)+sum(lgamma(a+(1-(1:d))/2))
}


##### Mode Function #####
# Computes the mode of a vector

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}





##### Posterior Distribution Functions #####
# Functions to evaluate the posterior density of the random parameters sigma, theta, gamma

#### LogSigma Posterior Distribution ####
# Computes the logarithm of the posterior density of sigma parameter of the prior

log_post_sigma = function(k,n,sigma,theta){
  result = log(sigma)*(a1-1)+log(1-sigma)*(b1-1)+log(theta+sigma)*(c1-1)+(-d1*sigma)
  partial=0
  if(k>1){
    for(j in 1:(k-1)){
      partial = partial+log(theta+j*sigma)
    }
  }
  partial1=0
  for(j in 1:k){
    nj=n[j]
    partial1 = partial1+(logpochhammer(1-sigma,nj-1))
  }
  return(result+partial+partial1)
}


#### LogTheta Posterior Distribution ####
# Computes the logarithm of the posterior density of theta parameter of the prior

log_post_theta = function(k,n,sigma,theta){
  partial=0
  if(k>1){
    partial=sum(log(theta+(1:(k-1))*sigma))
  }
  result=log(theta+sigma)*(c1-1)+(-d1*sigma)-logpochhammer(theta+1,size-1)+partial
  return(result)
}


#### LogGamma Posterior Distribution ####
# Computes the logarithm of the posterior density of gamma parameter for correlation, it's the likelihood of the data

log_post_gamma = function(k,n,sigma,theta,gamma){ 
  result=0
  for (j in 1:k){
    nj=n[j]
    
    if(j==1){
      startingdatapoint=1
    }
    else{
      startingdatapoint = sum(n[1:(j-1)])+1
    }
    yj = y[startingdatapoint:(startingdatapoint+nj-1),]
    yj=as.matrix(yj)
    if (nj==1){
      yj=t(as.matrix(yj))
    }
    
    partial=rep(0,d)
    if(nj!=1){ 
      for(i in 2:nj){
        partial=partial+yj[i,]-gamma*yj[i-1,]
      }
    }
    
    partial1=matrix(data=rep(0,d*d),nrow=d,ncol=d)
    if(nj!=1){ 
      for(i in 2:nj){
        temp=yj[i,]-gamma*yj[i-1,]
        partial1=partial1+(temp%*%t(temp)/(1-gamma^2))
      }
    }
    
    nunj=nu0+nj
    
    knj=k0+(1-gamma)^2/(1-gamma^2)*(nj-1)+1
    
    mnj=1/knj*(m0*k0+yj[1,]+(1-gamma)/(1-gamma^2)*partial)
    
    psinj=psi0 + yj[1,]%*%t(yj[1,]) + partial1 + k0*(m0%*%t(m0)) + knj*(mnj%*%t(mnj)) 
    
    current= log(knj)*(-d/2)-log(k0)*(-d/2) + log(pi)*(-nj*d/2)-log(1-gamma^2)*((nj-1)/2) + log_Multi_Gamma(nunj/2,d)- log_Multi_Gamma(nu0/2,d) + log(det(psi0))*(nu0/2)-log(det(psinj))*(nunj/2) 
    
    result  = result + current
  }
  return(result)
}





##### Prior Functions #####
# Functions to evaluate the ratio of the general prior for split, merge or shuffle

##### prior for split ratio #####
prior_split_ratio = function(k,nj,l,sigma,theta,gamma){ 
  if(l>=nj-l){
    result=(theta+k*sigma)/(k+1)*exp(lchoose(nj,l)+logpochhammer(1-sigma,nj-l-1)-logpochhammer(l-sigma,nj-l))
  }
  else{
    result=(theta+k*sigma)/(k+1)*exp(lchoose(nj,nj-l)+logpochhammer(1-sigma,l-1)-logpochhammer(nj-l-sigma,l))
  }
  return(result)
}


##### prior for merge ratio #####
prior_merge_ratio = function(k,nj,nj1,sigma,theta,gamma){ 
  if(nj>=nj1){
    result=k/(theta+(k-1)*sigma)*exp(-lchoose(nj+nj1,nj)+logpochhammer(nj-sigma,nj1)-logpochhammer(1-sigma,nj1-1))
  }
  else{
    result=k/(theta+(k-1)*sigma)*exp(-lchoose(nj+nj1,nj1)+logpochhammer(nj1-sigma,nj)-logpochhammer(1-sigma,nj-1))
  }
  return(result)
}


##### prior for shuffle ratio #####
prior_shuffle_ratio = function(ni,ni1,j,sigma,theta,gamma){ 
  imax=max(ni,ni1)
  imin=min(ni,ni1)
  jmax=max(j,ni+ni1-j)
  jmin=min(j,ni+ni1-j)
  maxtot=max(imax,jmax)
  secondtot=min(imax,jmax)
  thirdtot=max(imin,jmin)
  mintot=min(imin,jmin)
  if(imax==jmax){
    result=1
  }
  else if(jmax>=imax){
    result=exp(sum(log((mintot+1):thirdtot))-sum(log((secondtot+1):maxtot))+logpochhammer(secondtot-sigma,maxtot-secondtot)-logpochhammer(mintot-sigma,maxtot-secondtot))
  }
  else {
    result=1/exp(sum(log((mintot+1):thirdtot))-sum(log((secondtot+1):maxtot))+logpochhammer(secondtot-sigma,maxtot-secondtot)-logpochhammer(mintot-sigma,maxtot-secondtot))
  }
  return(result)
}





##### Likelihood Functions #####
# Functions to evaluate the ratio of the integrated likelihood for split, merge or shuffle

##### multivariate likelihood ratio for split ####
likelihood_split_ratio = function(k,nj,l,startingdatapoint,sigma,theta,gamma){
  
  yj = y[startingdatapoint:(startingdatapoint+nj-1),]
  yj=as.matrix(yj)
  if (nj==1){
    yj=t(as.matrix(yj))
  }
  
  partialj=rep(0,d)
  if(nj!=1){
    for(i in 2:nj){
      partialj=partialj+yj[i,]-gamma*yj[i-1,]
    }
  }
  
  partialj1=matrix(data=rep(0,d*d),nrow=d,ncol=d)
  if(nj!=1){
    for(i in 2:nj){
      temp=yj[i,]-gamma*yj[i-1,]
      partialj1=partialj1+(temp%*%t(temp)/(1-gamma^2))
    }
  }
  
  yl = y[startingdatapoint:(startingdatapoint+l-1),]
  yl=as.matrix(yl)
  if (l==1){
    yl=t(as.matrix(yl))
  }
  
  partiall=rep(0,d)
  if(l!=1){
    for(i in 2:l){
      partiall=partiall+yl[i,]-gamma*yl[i-1,]
    }
  }
  
  partiall1=matrix(data=rep(0,d*d),nrow=d,ncol=d)
  if(l!=1){
    for(i in 2:l){
      temp=yl[i,]-gamma*yl[i-1,]
      partiall1=partiall1+(temp%*%t(temp)/(1-gamma^2))
    }
  }
  
  yjl = y[(startingdatapoint+l):(startingdatapoint+nj-1),]
  yjl=as.matrix(yjl)
  if (nj-l==1){
    yjl=t(as.matrix(yjl))
  }
  
  partialjl=rep(0,d)
  if(nj-l!=1){
    for(i in 2:(nj-l)){
      partialjl=partialjl+yjl[i,]-gamma*yjl[i-1,]
    }
  }
  
  partialjl1=matrix(data=rep(0,d*d),nrow=d,ncol=d)
  if(nj-l!=1){
    for(i in 2:(nj-l)){
      temp=yjl[i,]-gamma*yjl[i-1,]
      partialjl1=partialjl1+(temp%*%t(temp)/(1-gamma^2))
    }
  }
  
  nul=nu0+l
  nunjl=nu0+nj-l
  nunj=nu0+nj
  
  kl=k0+(1-gamma)^2/(1-gamma^2)*(l-1)+1
  knjl=k0+(1-gamma)^2/(1-gamma^2)*(nj-l-1)+1
  knj=k0+(1-gamma)^2/(1-gamma^2)*(nj-1)+1
  
  ml=1/kl*(m0*k0+yl[1,]+(1-gamma)/(1-gamma^2)*partiall)
  mnjl=1/knjl*(m0*k0+yjl[1,]+(1-gamma)/(1-gamma^2)*partialjl)
  mnj=1/knj*(m0*k0+yj[1,]+(1-gamma)/(1-gamma^2)*partialj)
  
  psil=psi0 + yl[1,]%*%t(yl[1,]) + partiall1 + k0*(m0%*%t(m0)) + kl*(ml%*%t(ml))
  psinjl=psi0 + yjl[1,]%*%t(yjl[1,]) + partialjl1 + k0*(m0%*%t(m0)) + knjl*(mnjl%*%t(mnjl))
  psinj=psi0 + yj[1,]%*%t(yj[1,]) + partialj1 + k0*(m0%*%t(m0)) + knj*(mnj%*%t(mnj))
  
  loggammapart=sum(lgamma(nul/2+(1-1:d)/2)+lgamma(nunjl/2+(1-1:d)/2)-lgamma(nunj/2+(1-1:d)/2)-lgamma(nu0/2+(1-1:d)/2))
  logdetpart=(nu0/2)*log(det(psi0))+(nunj/2)*log(det(psinj))-(nul/2)*log(det(psil))-(nunjl/2)*log(det(psinjl))
  result=(kl*knjl/knj/k0)^(-d/2)/(1-gamma^2)^(-1/2)*exp(loggammapart+logdetpart)
  return(result)
}


##### multivariate likelihood ratio for merge ####
likelihood_merge_ratio = function(k,nj,nj1,startingdatapoint,sigma,theta,gamma){
  
  yj = y[startingdatapoint:(startingdatapoint+nj-1),]
  yj=as.matrix(yj)
  if (nj==1){
    yj=t(as.matrix(yj))
  }
  
  partialj=rep(0,d)
  if(nj!=1){
    for(i in 2:nj){
      partialj=partialj+yj[i,]-gamma*yj[i-1,]
    }
  }
  
  partialj1=matrix(data=rep(0,d*d),nrow=d,ncol=d)
  if(nj!=1){
    for(i in 2:nj){
      temp=yj[i,]-gamma*yj[i-1,]
      partialj1=partialj1+(temp%*%t(temp)/(1-gamma^2))
    }
  }
  
  yj1 = y[(startingdatapoint+nj):(startingdatapoint+nj+nj1-1),]
  yj1=as.matrix(yj1)
  if (nj1==1){
    yj1=t(as.matrix(yj1))
  }
  
  partialnj1=rep(0,d)
  if(nj1!=1){
    for(i in 2:nj1){
      partialnj1=partialnj1+yj1[i,]-gamma*yj1[i-1,]
    }
  }
  
  partialnj11=matrix(data=rep(0,d*d),nrow=d,ncol=d)
  if(nj1!=1){
    for(i in 2:nj1){
      temp=yj1[i,]-gamma*yj1[i-1,]
      partialnj11=partialnj11+(temp%*%t(temp)/(1-gamma^2))
    }
  }
  
  ytot = y[startingdatapoint:(startingdatapoint+nj+nj1-1),]
  ytot=as.matrix(ytot)
  if (nj+nj1==1){
    ytot=t(as.matrix(ytot))
  }
  
  partialtot=rep(0,d)
  if(nj+nj1!=1){
    for(i in 2:(nj+nj1)){
      partialtot=partialtot+ytot[i,]-gamma*ytot[i-1,]
    }
  }
  
  partialtot1=matrix(data=rep(0,d*d),nrow=d,ncol=d)
  if(nj+nj1!=1){
    for(i in 2:(nj+nj1)){
      temp=ytot[i,]-gamma*ytot[i-1,]
      partialtot1=partialtot1+(temp%*%t(temp)/(1-gamma^2))
    }
  }
  
  nunj1=nu0+nj1
  nutot=nu0+nj+nj1
  nunj=nu0+nj
  
  knj1=k0+(1-gamma)^2/(1-gamma^2)*(nj1-1)+1
  ktot=k0+(1-gamma)^2/(1-gamma^2)*(nj+nj1-1)+1
  knj=k0+(1-gamma)^2/(1-gamma^2)*(nj-1)+1
  
  mnj1=1/knj1*(m0*k0+yj1[1,]+(1-gamma)/(1-gamma^2)*partialnj1)
  mtot=1/ktot*(m0*k0+ytot[1,]+(1-gamma)/(1-gamma^2)*partialtot)
  mnj=1/knj*(m0*k0+yj[1,]+(1-gamma)/(1-gamma^2)*partialj)
  
  psinj1=psi0 + yj1[1,]%*%t(yj1[1,]) + partialnj11 + k0*(m0%*%t(m0)) + knj1*(mnj1%*%t(mnj1))
  psitot=psi0 + ytot[1,]%*%t(ytot[1,]) + partialtot1 + k0*(m0%*%t(m0)) + ktot*(mtot%*%t(mtot))
  psinj=psi0 + yj[1,]%*%t(yj[1,]) + partialj1 + k0*(m0%*%t(m0)) + knj*(mnj%*%t(mnj))
  
  loggammapart=sum(lgamma(nutot/2+(1-1:d)/2)+lgamma(nu0/2+(1-1:d)/2)-lgamma(nunj/2+(1-1:d)/2)-lgamma(nunj1/2+(1-1:d)/2))
  logdetpart=(nunj/2)*log(det(psinj))+(nunj1/2)*log(det(psinj1))-(nutot/2)*log(det(psitot))-(nu0/2)*log(det(psi0))
  result=(ktot*k0/knj/knj1)^(-d/2)/(1-gamma^2)^(1/2)*exp(loggammapart+logdetpart)
  return(result)
}


##### multivariate likelihood ratio for shuffle ####
likelihood_shuffle_ratio = function(ni,ni1,j,startingdatapoint,sigma,theta,gamma){
  
  yi = y[startingdatapoint:(startingdatapoint+ni-1),]
  yi=as.matrix(yi)
  if (ni==1){
    yi=t(as.matrix(yi))
  }
  
  partiali=rep(0,d)
  if(ni!=1){
    for(i in 2:ni){
      partiali=partiali+yi[i,]-gamma*yi[i-1,]
    }
  }
  
  partiali1=matrix(data=rep(0,d*d),nrow=d,ncol=d)
  if(ni!=1){
    for(i in 2:ni){
      temp=yi[i,]-gamma*yi[i-1,]
      partiali1=partiali1+(temp%*%t(temp)/(1-gamma^2))
    }
  }
  
  yi1 = y[(startingdatapoint+ni):(startingdatapoint+ni+ni1-1),]
  yi1=as.matrix(yi1)
  if (ni1==1){
    yi1=t(as.matrix(yi1))
  }
  
  partialni1=rep(0,d)
  if(ni1!=1){
    for(i in 2:ni1){
      partialni1=partialni1+yi1[i,]-gamma*yi1[i-1,]
    }
  }
  
  partialni11=matrix(data=rep(0,d*d),nrow=d,ncol=d)
  if(ni1!=1){
    for(i in 2:ni1){
      temp=yi1[i,]-gamma*yi1[i-1,]
      partialni11=partialni11+(temp%*%t(temp)/(1-gamma^2))
    }
  }
  
  yj = y[startingdatapoint:(startingdatapoint+j-1),]
  yj=as.matrix(yj)
  if (j==1){
    yj=t(as.matrix(yj))
  }
  
  partialj=rep(0,d)
  if(j!=1){
    for(i in 2:(j)){
      partialj=partialj+yj[i,]-gamma*yj[i-1,]
    }
  }
  
  partialj1=matrix(data=rep(0,d*d),nrow=d,ncol=d)
  if(j!=1){
    for(i in 2:(j)){
      temp=yj[i,]-gamma*yj[i-1,]
      partialj1=partialj1+(temp%*%t(temp)/(1-gamma^2))
    }
  }
  
  yj1 = y[(startingdatapoint+j):(startingdatapoint+ni+ni1-1),]
  yj1=as.matrix(yj1)
  if (ni+ni1-j==1){
    yj1=t(as.matrix(yj1))
  }
  
  partialnj1=rep(0,d)
  if(ni+ni1-j!=1){
    for(i in 2:(ni+ni1-j)){
      partialnj1=partialnj1+yj1[i,]-gamma*yj1[i-1,]
    }
  }
  
  partialnj11=matrix(data=rep(0,d*d),nrow=d,ncol=d)
  if(ni+ni1-j!=1){
    for(i in 2:(ni+ni1-j)){
      temp=yj1[i,]-gamma*yj1[i-1,]
      partialnj11=partialnj11+(temp%*%t(temp)/(1-gamma^2))
    }
  }
  
  nui=nu0+ni
  nui1=nu0+ni1
  nuj=nu0+j
  nuj1=nu0+ni+ni1-j
  
  ki=k0+(1-gamma)^2/(1-gamma^2)*(ni-1)+1
  ki1=k0+(1-gamma)^2/(1-gamma^2)*(ni1-1)+1
  kj=k0+(1-gamma)^2/(1-gamma^2)*(j-1)+1
  kj1=k0+(1-gamma)^2/(1-gamma^2)*(ni+ni1-j-1)+1
  
  mi=1/ki*(m0*k0+yi[1,]+(1-gamma)/(1-gamma^2)*partiali)
  mi1=1/ki1*(m0*k0+yi1[1,]+(1-gamma)/(1-gamma^2)*partialni1)
  mj=1/kj*(m0*k0+yj[1,]+(1-gamma)/(1-gamma^2)*partialj)
  mj1=1/kj1*(m0*k0+yj1[1,]+(1-gamma)/(1-gamma^2)*partialnj1)
  
  psii=psi0 + yi[1,]%*%t(yi[1,]) + partiali1 + k0*(m0%*%t(m0)) + ki*(mi%*%t(mi))
  psii1=psi0 + yi1[1,]%*%t(yi1[1,]) + partialni11 + k0*(m0%*%t(m0)) + ki1*(mi1%*%t(mi1))
  psij=psi0 + yj[1,]%*%t(yj[1,]) + partialj1 + k0*(m0%*%t(m0)) + kj*(mj%*%t(mj))
  psij1=psi0 + yj1[1,]%*%t(yj1[1,]) + partialnj11 + k0*(m0%*%t(m0)) + kj1*(mj1%*%t(mj1))
  
  loggammapart=sum(lgamma(nuj/2+(1-1:d)/2)+lgamma(nuj1/2+(1-1:d)/2)-lgamma(nui/2+(1-1:d)/2)-lgamma(nui1/2+(1-1:d)/2))
  logdetpart=(nui/2)*log(det(psii))+(nui1/2)*log(det(psii1))-(nuj/2)*log(det(psij))-(nuj1/2)*log(det(psij1))
  result=(kj*kj1/ki/ki1)^(-d/2)*exp(loggammapart+logdetpart)
  return(result)
}





##### Split, Merge, Shuffle Functions #####
# Functions to compute the acceptance probabilities and update the partition accordingly

##### multivariate split function #####
# It decides where to split, computes the split probability and updates the partition accordingly

multi_split_partition = function(k,n,sigma,theta,gamma){ 
  
  if(length(which(n>1))==1){
    j=which(n>1)
  }
  else{
    j = sample(which(n>1),1)
  }
  nj=n[j]
  l = sample(1:(nj-1),1)
  treshold = runif(1,0,1)
  ngk = length(which(n>1))
  
  if(j==1){
    startingdatapoint=1
  }
  else{
    startingdatapoint = sum(n[1:(j-1)])+1
  }
  
  if(k>1){
    alpha = min(1, (1-q)/q*prior_split_ratio(k,nj,l,sigma,theta,gamma)*likelihood_split_ratio(k,nj,l,startingdatapoint,sigma,theta,gamma)*ngk*(nj-1)/k)
  }else{
    alpha = min(1, (1-q)*(size-1)*prior_split_ratio(k,nj,l,sigma,theta,gamma)*likelihood_split_ratio(k,nj,l,startingdatapoint,sigma,theta,gamma)) 
  }
  
  if(alpha > treshold){ 
    
    if(k==1){
      temp_n = c(l,nj-l) 
    } 
    else if(j==1){
      temp_n = c(l,nj-l,n[(j+1):k]) 
    } 
    else if(j==k){
      temp_n = c(n[1:(j-1)],l,nj-l) 
    }
    else{
      temp_n = c(n[1:(j-1)],l,nj-l,n[(j+1):k]) 
    }
    
    n = temp_n
    k = k+1
  }
  return(list("n" = n, "k" = k)) 
}


##### multivariate merge function #####
# It decides where to merge, computes the merge probability and updates the partition accordingly

multi_merge_partition = function(k,n,sigma,theta,gamma){ 
  
  j = sample(1:(k-1),1)
  nj=n[j]
  nj1=n[j+1]
  treshold = runif(1,0,1)
  
  if(k==2){
    temp_n = c(nj+nj1)
  }
  else if(j==1){
    temp_n=c(nj+nj1,n[(j+2):k])
  }
  else if(j==(k-1)){
    temp_n = c(n[1:(j-1)],nj+nj1)
  }
  else{
    temp_n = c(n[1:(j-1)],nj+nj1,n[(j+2):k])
  }
  
  ngk_1 = length(which(temp_n>1)) 
  
  if(j==1){
    startingdatapoint=1
  }
  else{
    startingdatapoint = sum(n[1:(j-1)])+1
  }
  
  if(k<size){
    alpha = min(q/(1-q)*prior_merge_ratio(k,nj,nj1,sigma,theta,gamma)*likelihood_merge_ratio(k,nj,nj1,startingdatapoint,sigma,theta,gamma)/ngk_1/(nj+nj1-1),1) 
  }else{
    alpha = min(q*(size-1)*prior_merge_ratio(k,nj,nj1,sigma,theta,gamma)*likelihood_merge_ratio(k,nj,nj1,startingdatapoint,sigma,theta,gamma),1) 
  }
  
  if(alpha > treshold){
    
    n = temp_n
    k = k-1
  }
  return(list("n" = n, "k" = k)) 
}

##### multivariate shuffle function #####
# It decides where to shuffle, computes the shuffle probability and updates the partition accordingly

multi_shuffle_partition = function(k,n,sigma,theta,gamma){
  
  i = sample(1:(k-1),1)
  ni=n[i]
  ni1=n[i+1]
  j = sample(1:(ni+ni1-1),1)
  
  if(i==1){
    startingdatapoint=1
  }
  else{
    startingdatapoint = sum(n[1:(i-1)])+1
  }
  
  alpha = min(prior_shuffle_ratio(ni,ni1,j,sigma,theta,gamma)*likelihood_shuffle_ratio(ni,ni1,j,startingdatapoint,sigma,theta,gamma),1)
  treshold = runif(1,0,1)
  
  if(alpha > treshold){ 
    
    if(k==2){
      temp_n = c(j,ni+ni1-j)
    }
    else if(i==1){
      temp_n = c(j,ni+ni1-j,n[(i+2):k])
    }
    else if(i==(k-1)){
      temp_n = c(n[1:(i-1)],j,ni+ni1-j)
    }
    else{
      temp_n = c(n[1:(i-1)],j,ni+ni1-j,n[(i+2):k])
    }
    
    n = temp_n
  }
  return(list("n" = n, "k" = k)) 
}






##### Algorithm Main Function #####

# This is the function to call in the main, it performs the mcmc 

# k is the starting number of partitions (length(n))
# n is the starting partition 
# q is the split proposal probability of each iteration (1-q is the merge proposal probability)
# nrep is the useful number of iterations of the MCMC
# burnin is the number of burnin iterations of the MCMC
# sigma is the starting value for sigma
# theta is the starting value for theta
# gamma is the starting value for gamma
# if debug is true the algorithm prints the proposed partition at each iteration

multi_split_merge_MCMC = function(k,n,q,nrep,burnin,sigma,theta,gamma,debug=FALSE){
  
  results_rhon= vector("list", nrep)
  results_k = numeric(nrep)
  results_sigma = numeric(nrep)
  results_theta = numeric(nrep)
  results_gamma = numeric(nrep)
  
  pb <- txtProgressBar(min = 1, max = burnin+nrep, initial = 1, style = 3) 
  
  for(iter in 1:(burnin+nrep)){ 
    u = runif(1,0,1)
    
    if((k==1) || (u < q && k<size) ){ 
      splitted=multi_split_partition(k,n,sigma,theta,gamma) 
      n=splitted[["n"]]
      k=splitted[["k"]]
      
    }else{
      merged=multi_merge_partition(k,n,sigma,theta,gamma)
      n=merged[["n"]]
      k=merged[["k"]]
    }
    if(k>1){
      shuffled=multi_shuffle_partition(k,n,sigma,theta,gamma)
      n=shuffled[["n"]]
      k=shuffled[["k"]]
    }
    if(debug){
      print(n)
    }
    sigmanew = arms(sigma, function(sigma) log_post_sigma(k,n,sigma,theta), function(sigma) (sigma>max(-theta,0))*(sigma<1), 1)
    thetanew = arms(theta, function(theta) log_post_theta(k,n,sigmanew,theta), function(theta) (theta>-sigmanew)*(theta<size), 1)
    gammanew = arms(gamma, function(gamma) log_post_gamma(k,n,sigma,theta,gamma), function(gamma) (gamma>0)*(gamma<1), 1)
    theta=thetanew
    sigma=sigmanew
    gamma=gammanew
    
    if(iter>burnin){
      results_k[iter-burnin] = k
      results_rhon[[iter-burnin]] = n
      results_sigma[iter-burnin] = sigma
      results_theta[iter-burnin] = theta
      results_gamma[iter-burnin] = gamma
    }
    setTxtProgressBar(pb, iter)
  }
  close(pb)
  return(list("n" = results_rhon, "k" = results_k, "sigma" = results_sigma ,"theta" = results_theta, "gamma" = results_gamma ))
}





##### Ordered results #####
# Function to have all partitions sampled ordered with their correspondent number of occurrences

pretty_results = function(solution_rhon){
  
  ux = unique(solution_rhon)
  matched = match(solution_rhon, ux)
  a=cbind(partition = ux, occurrences = tabulate(matched))
  adata=as.data.frame(a)
  for(i in 1:nrow(adata)){
    adata[i,"final_occurrences"]=adata[i,"occurrences"][[1]]
  }
  adata=adata[,-2] 
  adata_ordered <- adata[order(-adata$final_occurrences),]
  
  return(adata_ordered)
}





##### posterior probability of the change points #####
# Function to compute the posterior probability for every time step to be a change point

posterior_prob_change_points = function(solution){
  
  prob_change_points = rep(0,size)
  
  for(i in 1:length(solution[['n']])){
    
    n_act = solution[['n']][[i]]
    current_change_point = 0
    for(j in n_act){
      current_change_point = current_change_point + j
      prob_change_points[current_change_point] = prob_change_points[current_change_point] + 1
    }
    
  }
  prob_change_points[length(prob_change_points)] = 0
  prob_change_points  = prob_change_points/length(solution[['n']])
  
  
  return(prob_change_points)
}





##### best partition estimator #####
# Function to compute the partition with change points detected as the ones with a posterior probability above a treshold

find_best = function(post_prob_change_points, treshold, k){
  
  rhon_highest_posterior_raw=which(post_prob_change_points>treshold)
  rhon_highest_posterior <- c()
  
  for (i in 1:length(rhon_highest_posterior_raw)){
    if (i!=1 && rhon_highest_posterior_raw[i]-rhon_highest_posterior_raw[i-1]<k){
      if (post_prob_change_points[rhon_highest_posterior_raw[i-1]]<post_prob_change_points[rhon_highest_posterior_raw[i]]){
        rhon_highest_posterior=rhon_highest_posterior[-length(rhon_highest_posterior)]
        rhon_highest_posterior <- c(rhon_highest_posterior, rhon_highest_posterior_raw[i])
      }
    }
    else{
      rhon_highest_posterior <- c(rhon_highest_posterior, rhon_highest_posterior_raw[i])
    }
  }
  
  return(rhon_highest_posterior)
}


