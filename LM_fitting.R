library(copula)
library(nleqslv)
library(survival)

####################################################################################
### Purpose: To generate data from landmark model using file 'LM_Data_Generation'
###.         of survival and multivariate longitudianl outcomes
###          special case theta obeys Beta distribution
###          M: number of longitudinal biomarkers
###
###
### Both linear and nonlinear trajectory of longitudinal biomarkers can be generated, 
### by using different 'y.pdf.input'
### 'y.pdf.input.beta' for linear trajectory, 'y.pdf.input.beta.nonlinear' for nonlinear
###
###
####################################################################################
source('LM_Data_Generation.R') 

#difference between two landmark times
repeated.interval = 1
#number of landmark times
n.landmark.time = 6
#landmark times
landmark.time = repeated.interval * seq(0, n.landmark.time - 1) 

##########################################################################
###  Used to generate input data, Marginal distribution of theta(s) obeys
###  Beta distribution to get linear trajectory of longitudinal biomarkers
##########################################################################
b_beta_theta =  2.5 + 0.05 * landmark.time 
a_beta_theta = 2 - 0.02 * landmark.time  

###  Marginal distribution of theta(s) obeys Beta distribution 
###  PDf of Beta(S)
beta.distribution =  function(theta.var, s){   #s: index
  gamma(a_beta_theta[s] + b_beta_theta[s])/(gamma(a_beta_theta[s]) * 
                                              gamma(b_beta_theta[s])) * theta.var^(a_beta_theta[s] - 1) * 
    (1- theta.var)^(b_beta_theta[s] - 1)
}

##generate data
upper.theta <- 1                          # upper bound of theta;
delta.theta <- 0.002                      # grid to discretize theta;
theta.beta <- seq(delta.theta/2, 1, by = delta.theta)
y.pdf.input.beta  = matrix(min(theta.beta), length(landmark.time), length(theta.beta))
for(s in 1: length(landmark.time)){
  y.pdf.input.beta[s,] = beta.distribution(theta.beta, s)
}
##########################################################################
###  Used to generate input data, Marginal distribution of theta(s) obeys
###  Beta distribution to get nonlinear trajectory of longitudinal biomarkers
##########################################################################
b_beta_theta = 5  + 4 * sin(0.5 * landmark.time) 
a_beta_theta = 3 * (2 - 0.2 * landmark.time)

###  Marginal distribution of theta(s) obeys Beta distribution 
###  PDf of Beta(S)
beta.distribution =  function(theta.var, s){   #s: index
  gamma(a_beta_theta[s] + b_beta_theta[s])/(gamma(a_beta_theta[s]) * 
                                              gamma(b_beta_theta[s])) * theta.var^(a_beta_theta[s] - 1) * 
    (1 - theta.var)^(b_beta_theta[s] - 1)
}

##generate data
upper.theta <- 1                          # upper bound of theta;
delta.theta <- 0.002                      # grid to discretize theta;
theta.beta <- seq(delta.theta/2, 1, by = delta.theta)
y.pdf.input.beta.nonlinear  = matrix(min(theta.beta), length(landmark.time), length(theta.beta))
for(s in 1: length(landmark.time)){
  y.pdf.input.beta.nonlinear[s,] = beta.distribution(theta.beta, s)
}

####################################################################################
### 
###  Parameters and Data Generation of
###  landmark model and model fitting with one biomarker and two baseline covariates
###  Can be easily extended to multiple biomarkers by changing the value of M
###
####################################################################################
theta = theta.beta 

### you can choose different y.pdf.input to use linear trajectory or nonlinear trajectory
y.pdf.input = y.pdf.input.beta 
              #y.pdf.input.beta.nonlinear

#number of subjects
n.subjects = 1500

#number of biomarkers and baseline covariates
M = 1; KK = 2
rho1 = 0.6 # for copula

## true parameter
beta = matrix(1.2, length(landmark.time), M)
beta_cor = matrix(0, length(landmark.time), KK)
beta_cor[, 1] = rep(1.3, n.landmark.time) 
beta_cor[, 2] = rep(0.55, n.landmark.time) 

##parameters for baseline hazard
lambda = 0.3; kappa = 1.5
index = 1    ##used for random seed

## censor status = 1 have censor, = 0 without censor  
censor.status = 0

seed <- 20 * (index - 1)
### data used to train the model
data.raw.sim = LandmarkGeneDataGeneral(theta, y.pdf.input, n.subjects, 
                                       M = 1, KK = 2, rho1, landmark.time, beta, 
                                       beta_cor,
                                       repeated.interval,
                                       lambda, kappa, seed, censor.status)

### Fit landmark model using coxph 
beta.est <- matrix(0, length(landmark.time),  M + KK)
baselinehazard <- list()
for (s in 1: 6){ 
  fit <- coxph(Surv(Ts, status) ~ longitudinal1 + covariate1  + covariate2,
               data = data.raw.sim, 
               subset = time == (s - 1) * repeated.interval, ties = "breslow")
  beta.est[s,] <- as.numeric(fit$coef)
  baselinehazard[[s]] = basehaz(fit, centered = FALSE)
}    

## first col for beta, second and third col for beta_cor
beta.est

