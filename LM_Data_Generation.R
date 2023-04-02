
###############################################################################################
##  Purpose: Function used to generate data of landmark model of survival and 
##           multivariate longitudianl outcomes
##           
##
## ---------------------------------------------------------------------------           
##  Date:    04/02/2023
##  Author:  Wenhao Li
##           E-mail: wenhao.li@uth.tmc.edu
###############################################################################################


####################################################################################
## 
##  Function of parameters and generate simulations data for M longitudinal based on 
##  landmark model of special case
##  (marginal distribution of theta: general case)
##
####################################################################################

LandmarkGeneDataGeneral = function(theta, y.pdf.input, n.subjects = 50, M = 1, 
                                   KK = 1, rho1 = 0.2, landmark.time = c(0, 1, 2, 3, 4, 5),
                                   beta = c(0.3, 1.4, 0.5, 1.05, 0.4, 0.1), 
                                   beta_cor = c(0.3, 1.4, 0.5, 1.05, 0.4, 0.1), 
                                   repeated.interval = 1, lambda = 0.3, kappa = 1.5, 
                                   seed =1, censor.status = 1){
  set.seed(seed)
  
  ##########################################################################
  ###  Solve for Baseline Hazard H_0(u, s) ###
  ##########################################################################
  
  ###pre-defined H_0(u, 0) Weibull dist
  H_0_base <- function(t) {               
    (lambda * t)^kappa                   
  }
  ###function to get H_0(u, s)
  #s: index of landmark time, landmark.time[s]: landmark time
  #u: residual life time = T - s, where T is time to event
  H_0_solve = function(u, s){
    H_0_s_u = sum(exp(-H_0_base(landmark.time[s] + u) * theta) * y.pdf.input[1,] * 
                    delta.theta)
    H_0_s = sum(exp(-H_0_base(landmark.time[s]) * theta) * y.pdf.input[1,] * 
                  delta.theta)
    H_0_general = function(H_0){
      sum(exp(-H_0 * theta) * y.pdf.input[s,] * delta.theta) - H_0_s_u/H_0_s
    }
    
    H_0_start = 1  ## initial value
    return(nleqslv(H_0_start, H_0_general)$x)
  }
  
  ##########################################################################
  ###  Marginal distribution of survival time T ###
  ##########################################################################
  
  ##generate survival function S(T)
  S_T = function(u, s){
    S.T =  sum(exp(-H_0_base(u) * theta) * y.pdf.input[1,] * 
                 delta.theta)
    return(S.T)
  }
  
  ##discrete time to event T used to discrete the survival function
  upper.T <- 30                               # upper bound of time (T);
  delta.T <- 0.002
  grid.T <- seq(0, upper.T, by=delta.T)
  
  ##grid S(T)
  S.T.grid = c()
  for(i in 1:length(grid.T)){
    S.T.grid[i] = S_T(grid.T[i], 1)
  }
  S.T.grid[length(grid.T)] = 0
  
  ##generate uniform random variable
  T.uniform <- runif(n.subjects, min = 0, max = 1)
  
  ##generate time to event
  T_surv_sapply = function(x){
    min(grid.T[S.T.grid <= x]) 
  }
  time100 <- proc.time()
  T.surv = sapply(T.uniform, T_surv_sapply)
  sim.time.sapply <-proc.time() - time100
  
  ## censoring: maximum of measurements per subject
  if(censor.status == 1){
    C =  runif(n.subjects , min = 5, max = 30) #
  }else{
    C =  rep(65, n.subjects)
  }
  
  ## generate real survival time
  T.real = rep(0,n.subjects)
  delta = rep(0,n.subjects)
  for (i in 1: n.subjects){
    if (T.surv[i] > C[i]){
      T.real[i] = C[i]  ##censored
      delta[i] = 0
    }else{
      T.real[i] = T.surv[i]
      delta[i] = 1
    }
  }
  
  ## censor rate
  censor.rate = 1 - sum(delta/length(delta))
  
  ##ni responses in the ith subject,i = 1,2...n. measurements every 1 months,
  n.visits = floor(T.real / repeated.interval) + 1 
  
  ##landmark time
  t = list()
  for (i in 1:n.subjects){
    t[[i]] = seq(from = 0, to = (n.visits[i] - 1) * repeated.interval, 
                 by = repeated.interval)
  }
  
  ##residual life time
  T.resi = list()
  for (i in 1: n.subjects){
    T.resi[[i]] = T.real[i] - t[[i]]
  }
  
  ##########################################################################
  ###  Generate conditional distribution of theta(s) given W(s) ###
  ###########################################################################
  
  ##grid PDF of theta(s) given W(S)
  pdf_theta_conditional_W = function(s){
    pdf_w_conditional_theta = exp( -H_0_solve(T.resi[[i]][s], s) * theta) * 
      H_0_solve(T.resi[[i]][s], s) * theta
    pdf_w = sum(exp( -H_0_solve(T.resi[[i]][s], s) * theta) * H_0_solve(T.resi[[i]][s], s) * 
                  theta * y.pdf.input[s,] * delta.theta)
    pdf_theta_conditional_W = pdf_w_conditional_theta * 
      y.pdf.input[s,] / pdf_w
    return(pdf_theta_conditional_W)
  }
  
  ##generate random uniform variable with and without copula
  theta.random.copula = list()
  for (i in 1:n.subjects){
    if (n.visits[i]>1){
      ##random number from copula
      cp <- normalCopula(param = rho1, dim = n.visits[i] , dispstr = "ex")
      theta.random.copula[[i]] = rCopula(n = 1, copula = cp)
      ##random number from uniform distribution
    }else{
      theta.random.copula[[i]] <- runif(n = 1)
    }
  }
  
  ##discrete the PDF of theta(s) given W(S)
  ##generate discrete CDF of theta(s) given W(S)
  ##generate randon variable theta(s) obeys PDF of theta(s) given W(S)
  ##generate M biomarders from theta(s)
  theta.surv = list()
  X = list()  #Y.long longitudianl data: X
  X_cor = list()  ## baseline covariates
  for(i in 1:n.subjects){
    theta.surv[[i]] = rep(0, n.visits[i] )
    X[[i]] = matrix(0, n.visits[i], M) 
    
    if(KK == 0){
      X_cor[[i]] = 0
    }else{
      X_cor[[i]] = rep(0, KK)
      if(KK == 1){
        X_cor[[i]] = rbinom(1, 1, 0.5)
      }else{
        for(kk in 1:(KK-1) ){
          X_cor[[i]][kk] = rbinom(1, 1, 0.5)
          #rnorm(1, 0, 4)
        }
        X_cor[[i]][KK] = rnorm(1, 0, 4)
      }
    }
    for(s in 1: min(n.visits[i], length(landmark.time)) ){ ##s: index of landmark time
      pdf.theta.grid = pdf_theta_conditional_W(s)
      cdf.theta.grid = cumsum(pdf.theta.grid * delta.theta)
      cdf.theta.grid[1] <- 0
      F.theta = runif(1, min = 0, max = 1)      
      #with copula
      theta.surv[[i]][s] =  max(min(theta[cdf.theta.grid >=        
                                            theta.random.copula[[i]][s] ]) - 
                                  delta.theta/2, delta.theta/4) 
      if(M == 1) {
        ###one biomarker
        X[[i]][s, 1] = (log(theta.surv[[i]][s]) - sum(X_cor[[i]] * 
                                                        as.matrix(beta_cor)[s, ]) )/beta[s]
      } else {
        ###multiple biomarkers
        for(m in 1:(M - 1)){
          ##longitudianl biomarker increases
          X[[i]][s, m] = 0.5 * landmark.time[s] + rnorm(1, 0, 0.5)
          
        }
        X[[i]][s, M] = (log(theta.surv[[i]][s]) - sum(X[[i]][s, 1:M - 1] * 
                                                        beta[s, 1:M - 1]) - sum(X_cor[[i]] * 
                                                                                  as.matrix(beta_cor)[s, ]) )/beta[s, M]
      }
      
    }
  }
  
  ##########################################################################
  ### Data frame
  ##########################################################################
  time=c()
  for(ii in 1:n.subjects){
    time=cbind(time, cbind(t(t[[ii]])) )
  }
  start=c()
  for(ii in 1:n.subjects){
    start=cbind(start, cbind(t(t[[ii]])) )
  }
  fuyrs=c()
  for(ii in 1:n.subjects){
    fuyrs=cbind(fuyrs, cbind(t(rep(T.real[[ii]],n.visits[[ii]]) )) )
  }
  Ts=c()
  for(ii in 1:n.subjects){
    Ts=cbind(Ts, cbind(t(T.resi[[ii]])) )
  }
  stop=c()
  for(ii in 1:length(start)){
    stop[ii]=min(start[ii]+1,fuyrs[ii])
  }
  status=c()
  for(ii in 1:n.subjects){
    status=cbind(status, cbind(t(rep(delta[[ii]],n.visits[[ii]]) )) )
  }
  longitudinal.data = vector("list", M)
  for(m in 1: M){   
    for(ii in 1:n.subjects){
      longitudinal.data[[m]] = cbind(longitudinal.data[[m]], cbind(t(X[[ii]][,m]) )) 
    }
  }
  if(KK == 0){
    KK = 0
  }else{
    covariate.data = vector("list", KK)
    for(kk in 1: KK){ 
      for(ii in 1:n.subjects){
        covariate.data[[kk]]=cbind(covariate.data[[kk]], cbind(t(rep(X_cor[[ii]][kk], n.visits[[ii]]) )) )
      }
    }
  }
  theta.data = c()
  for(ii in 1:n.subjects){
    theta.data = cbind(theta.data, cbind(t(theta.surv[[ii]]) )) 
  }
  num=c()
  for(ii in 1:n.subjects){
    num=cbind(num, t(rep(ii, n.visits[[ii]]))  )
  }
  
  ###***###
  data.raw3=data.frame(
    num=c(num), 
    time=c(time),
    fuyrs=as.numeric(c(fuyrs)),
    Ts=c(Ts),
    status=as.numeric(c(status)),
    start=as.numeric(c(start)),
    stop=as.numeric(c(stop)),
    censor_rate = censor.rate,
    theta = as.numeric(c(theta.data))
  )
  for(m in 1:M) {
    data.raw3[paste("longitudinal", m, sep = "")] = c(longitudinal.data[[m]])
  }
  if(KK == 0){
    KK = 0
  }else{
    for(kk in 1:KK) {
      data.raw3[paste("covariate", kk, sep = "")] = c(covariate.data[[kk]])
    }
  }
  
  ## before max landmark time data
  data.raw.sim = data.raw3[data.raw3$time<= max(landmark.time),]
  rownames(data.raw.sim)<-1:nrow(data.raw.sim)
  
  return(data.raw.sim)
}


