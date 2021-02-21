##### Parameters for Finance Data #####
# Standard choices since the algorithm is robust, very important k0 parameter





##### General Parameters #####

d=ncol(y) # dimensions
size=nrow(y) # time steps


##### Prior Parameters #####

a1=1 # a parameter in the prior for sigma
b1=1 # b parameter in the prior for sigma
c1=1 # c parameter in the prior for theta
d1=1 # d parameter in the prior for theta
start_theta=0.1897 # starting value for theta
start_sigma=0.1 # starting value for sigma


##### Likelihood Parameters #####

m0=colMeans(y) # m0 parameter in the integrated likelihood
k0=10 #0.3 # k0 parameter in the integrated likelihood
nu0=5 # nu0 parameter in the integrated likelihood
psi0=var(y) # psi0 parameter in the integrated likelihood
start_gamma=0.5 # starting value for gamma


##### Algorithm Parameters #####

start_n = c(58,58,58,58) # starting partition 
start_k = length(start_n) # starting number of partitions (not really important)
q = 0.5 # split proposal probability of each iteration (1-q is the merge proposal probability)
burnin = 5000 # number of burnin iterations
nrep = 20000 # number of useful iterations
debug=TRUE # shows the partition at each iteration


