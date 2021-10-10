library("survival")
library("MASS")
library("zoo")
library("PermAlgo")
library("dplyr")

# parameter settings
n = 100
p = 3
rho = 0.7
d = 4  # total number of dep and indep covariates

## generate time-indept covariate
indep_input1 <- rbinom(n,1,0.55)
# center the data
indep_input1 <- indep_input1-mean(indep_input1)

indep_input2 <- rbinom(n,1,0.55)
# center the data
indep_input2 <- indep_input2-mean(indep_input2)


## generate time-dept covariate

## autoregressive covariance matrix
covmat <- function(rho, p) {
  rho^(abs(outer(seq(p), seq(p), "-")))
}
sigma<-covmat(rho,p)


## We assume that the time-dependent variable remains constant between measurements.
# p = 3 here represents 3 time intervals [0, 10) [10, 20), [20, --]


# U is exponential
U <- rexp(n=n, 0.05) 

# generate a N(0, sigma) (n by p)
dep_input1  <-  mvrnorm(n = n, rep(0,p), sigma) + matrix(0.2*U-4,n,p,byrow=FALSE)
dep_input2  <-  mvrnorm(n = n, rep(0,p), sigma) + matrix(0.2*U-4,n,p,byrow=FALSE)

a1 <- matrix(dep_input1[,1], nrow = n, ncol = 10, byrow = F)
a2 <- matrix(dep_input1[,2], nrow = n, ncol = 10, byrow = F)
a3 <- matrix(dep_input1[,3], nrow = n, ncol = 80, byrow = F)

b1 <- matrix(dep_input2[,1], nrow = n, ncol = 10, byrow = F)
b2 <- matrix(dep_input2[,2], nrow = n, ncol = 10, byrow = F)
b3 <- matrix(dep_input2[,3], nrow = n, ncol = 80, byrow = F)


# organize the time-dependent data into a long format
# each row represents a measurement for a person at a specific time
# first two entries of each row corresponds to time-independent covariates
# last two engries of each row corresponds to time-dependent covariates
# each person has m time measurements, the first m rows of Xmat are data for
# the 1st person, row m+1 to row 2m are data for the second person, etc.

m=100   # maximum time interval considered
# colnames(dep_input)<-paste("t",1:3,sep="")
Xmat = matrix(0, ncol = d, nrow = m * n)
Xmat[,1] = rep(indep_input1,each = m)
Xmat[,2] = rep(indep_input2,each = m)
Xmat[,3] = as.double(t(cbind(a1, a2, a3)))
Xmat[,4] = as.double(t(cbind(b1, b2, b3)))

# generate random event time, PermAlgo will match it with covariates.
eventRandom_V <- rexp(n, 0.01) + 10

# censorRandom is set to m+1, so there is no censoring
censorRandom = rep(m+1, n)  # 



#### generating V

gamma_true_V <- c(0.4,0.4,0.2,0.2)
data_for_V <- permalgorithm(n, m, Xmat, 
                            XmatNames=c("t.indep1", "t.indep2","t.dep1", "t.dep2"), 
                            eventRandom = eventRandom_V, censorRandom=censorRandom, 					  
                            betas = gamma_true_V, groupByD=FALSE)

## extract unique V time
temp=""
data_temp <- lapply(1:n, 
                    function(i) temp[i]=(unique((dplyr::filter(data_for_V, Id==i))$Fup)))

# break the tie
V_temp <- unlist(data_temp) - abs(rnorm(length(data_temp), 0,  1e-5))

# break the tie for V
for(i in 1:n){
	data_for_V[data_for_V$Id==i,]$Fup <- V_temp[i] 	
} 


## generate censoring status C  
# generate the time varying treatment variable
trt  <- unlist(lapply(V_temp, function(x) c(rep(1, round(x)),rep(0,m-round(x)))))
# append treatment to the data matrix
Xmat_for_C <- cbind(Xmat, trt)
# generate true gamma for c
gamma_true_C <- c(0.4,0.4,0.2,0.2,0.5) 

eventRandom_C <- rexp(n, 0.01) + 10
# generate C using the covariates and treatment 
data_for_C <- permalgorithm(n, m, Xmat_for_C, XmatNames=c("t.indep1", "t.indep2","t.dep1", "t.dep2", "trt"), eventRandom = eventRandom_C, censorRandom=censorRandom, betas = gamma_true_C, groupByD=FALSE )


temp=""
data_temp <- lapply(1:n, function(i) temp[i]=(unique((dplyr::filter(data_for_C, Id==i))$Fup)))

# break the tie
C_temp <- unlist(data_temp) - abs(rnorm(length(data_temp), 0,  1e-5))

# break the tie for C
for(i in 1:n){
	data_for_C[data_for_C$Id==i,]$Fup <- C_temp[i] 	
} 



## generate T according to a SFTM
psi_true <- 0.5
T1 <- U * exp(-psi_true)
T <- (T1 < V_temp) * T1 + (T1 >= V_temp) * (U + V_temp - V_temp * exp(psi_true))


# check whether there is tie in V or in min(T,C)
# T_censored is death or censor
# V_temp is discontinuition 
T_censored <- apply(cbind(T,C_temp),1,min) 
check_tie <- diff(sort(c(V_temp, T_censored)))
if(min(check_tie)<1e-10){
  print("Please re-generate the data in order to avoid the same time points for V and U")
}

# gammaV is indicator of discontinuation before time to death or censoring
# V is death, discontinuation or censor
V <- apply(cbind(V_temp,C_temp,T),1,min)
gammaV <- (V_temp < T_censored)



deltaT <- ( T < C_temp ) 

## *************end of simulation ******************

## step 1

fit_V = coxph(Surv(Start, Stop, Event) ~ t.indep1 + t.indep2
              + t.dep1 + t.dep2, data = data_for_V)
gammahat_V <- fit_V$coefficients
data_for_V$exp_gammahat_V_by_g <- exp(as.matrix(data_for_V[,6:9]) %*% gammahat_V)


cum_hazard <- -log(survfit(fit_V)$surv) 
lambda0hat_V <- cum_hazard-c(0,cum_hazard[1:(length(cum_hazard)-1)])  

##  multiply the <same id exp_gammahat_V_by_g> terms with lambda accordingly
temp=list()
data_temp <- lapply(1:n,
       function(i){
         slice = ((dplyr::filter(data_for_V, Id==i))$exp_gammahat_V_by_g)
         temp[[i]]= slice * lambda0hat_V[1:nrow(slice)] 
         })
data_for_V$exp_gammahat_V_g_lambda<- unlist(data_temp)


############### warning -- For real data we should use floor() function
## here instead of ceiling function

ceiling(V)
temp=list()
data_temp_trucated <- lapply(1:n,
       function(i){
         temp[[i]]= dplyr::filter(data_for_V, Id==i)[1:ceiling(V)[i], ]
         })
		
tmp <- matrix(0,m,n)		
for(i in 1:n){
	dat <- dplyr::filter(data_for_V, Id==i)
	tmp[1:nrow(dat), i] <- dat$exp_gammahat_V_g_lambda
	if(nrow(dat)+1 <= m) tmp[(nrow(dat)+1):m, i] <- 0
}	
exp_gammahat_V_g_lambda_YV <- as.vector(tmp)


tmp <- matrix(0,m,n)
for(i in 1:n){
	tmp[floor(V_temp)[i] , i] <- 1
}		 
dN_V <- as.vector(tmp)

dM_V <- dN_V - exp_gammahat_V_g_lambda_YV


## step2 
fit_C = coxph(Surv(Start, Stop, Event) ~ t.indep1 + t.indep2
              + t.dep1 + t.dep2 + trt, data = data_for_C)
gammahat_C <- fit_C$coefficients


## gammahat_C is a list of length5, since we have an additional coefficients for "trt"
## cropped gammahat_C[1:4] to match length of variables 
data_for_C$exp_gammahat_C_by_g <- exp(as.matrix(data_for_C[,6:10]) %*% gammahat_C)

cum_hazard <- -log(survfit(fit_C)$surv) 
lambda0hat_C <- cum_hazard-c(0,cum_hazard[1:(length(cum_hazard)-1)])


# to get exp_gammahat_C_g_lambda

temp=list()
data_temp <- lapply(1:n,
                    function(i){
                      slice = ((dplyr::filter(data_for_C, Id==i))$exp_gammahat_C_by_g)
                      temp[[i]]= slice * lambda0hat_C[1:nrow(slice)]
                    })
data_for_C$exp_gammahat_C_g_lambda<- unlist(data_temp)


######## Using Prod(1+g(u)du) = exp int g(du) du 

K_c <- rep(0,n)
for(i in 1:n){
	tmp <- dplyr::filter(data_for_C, Id==i)[1:ceiling(T_censored)[i], ]
	K_c[i] <- exp(sum(-tmp$exp_gammahat_C_g_lambda))
}

######## here is an alternative way to compute K_c, using the direct computation of the product 

# K_c <- rep(0,n)
# for(i in 1:n){
# 	tmp <- dplyr::filter(data_for_C, Id==i)[1:ceiling(T_censored)[i], ]
# 	K_c[i] <- prod(1-tmp$exp_gammahat_C_g_lambda)
# }



############### warning -- For real data we should use floor() function
## here instead of ceiling function



