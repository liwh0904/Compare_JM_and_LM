library(survival)
library(pROC)
library(peperr)
library(nlme)

### example data
load("example_data.Rdata")
####################################################################################
### Purpose: Dynamic prediction of LM and LM2, example_data is used for illustration
### 'data.raw.sim' in example data for training and 'data.raw.prediction' for prediction
###
### LM: landmark model, LM2: landmark model with slope
###
####################################################################################


####################################################################################
### 
###  Dynamic prediction of LM
###
####################################################################################

#difference between two landmark times
repeated.interval = 1
#number of landmark times
n.landmark.time = 6
#landmark times
landmark.time = repeated.interval * seq(0, n.landmark.time - 1) 

##landmark times from s = 1 to 6
s = 1
prediction.horizon = 1.5

data.raw.predict.LM = data.raw.predict[data.raw.predict$time == landmark.time[s], ]

fit <- coxph(Surv(Ts, status) ~ longitudinal1 + covariate1  + covariate2,
             data = data.raw.sim, 
             subset = time == (s - 1) * repeated.interval, ties = "breslow")
condition.prediction.LM = predictProb(fit, response = Surv(data.raw.predict.LM$Ts, 
                                                           data.raw.predict.LM$status),  x = data.raw.predict.LM,
                                      time  = c(prediction.horizon))

### ROC curve
data.raw.predict.LM$survivalProb = c(condition.prediction.LM)
data.raw.predict.LM$True = (data.raw.predict.LM$Ts > prediction.horizon)
roc.LM = roc(data.raw.predict.LM$True, data.raw.predict.LM$survivalProb, 
             levels=c("TRUE", "FALSE"))
AUC.LM = pROC::auc(roc.LM)
###Brier Score
BS.LM = mean((as.numeric(data.raw.predict.LM$True) - data.raw.predict.LM$survivalProb)^2)

###############################################################################################
###
### Dynamic prediction of LM-2 with slope 
###
###############################################################################################
set.seed(8)

### landmark times greater than 1, s greater or equal than 2
s= 2
#number of subjects
n.subjects = 500
### get the slope using linear mixed effects model
data.raw.LM2.fit =  data.raw.sim[data.raw.sim$time <= landmark.time[s],]
data.raw.LM2.predict =  data.raw.predict[data.raw.predict$time <= landmark.time[s], ]
ctrl <- lmeControl(10000000, 1000000, opt='optim')
fitLME.fit.slope = lme(longitudinal1  ~ time + covariate1 + covariate2, 
                       random = ~ time | num, 
                       data = data.raw.LM2.fit, control=ctrl)
fitLME.predict.slope = lme(longitudinal1  ~ time + covariate1 + covariate2, 
                           random = ~ time | num, 
                           data = data.raw.LM2.predict, control = ctrl)

slope.fit = c()
slope.predict = c()
for(i in 1 : n.subjects){
  slope.fit = cbind(slope.fit, cbind(t(rep(summary(fitLME.fit.slope)[[4]][[2]][[1]][i,2], 
                                           sum(data.raw.LM2.fit$num == i) ) )) )
}
data.raw.LM2.fit$slope = c(slope.fit)
for(i in 1 : n.subjects){
  slope.predict = cbind(slope.predict, cbind(t(rep(summary(fitLME.predict.slope)[[4]][[2]][[1]][i,2], 
                                                   sum(data.raw.LM2.predict$num == i) ) )) )
}
data.raw.LM2.predict$slope = c(slope.predict)
data.raw.LM2.predict =  data.raw.LM2.predict[data.raw.LM2.predict$time == landmark.time[s], ]

###model fitting using coxph
LM2.coxph.fit = coxph(Surv(Ts, status) ~ covariate1 + covariate2 + longitudinal1 + slope, 
                      data = data.raw.LM2.fit, subset = time == (s - 1) * repeated.interval, ties = "breslow")

###prediction using predictProb
condition.prediction.LM2 = predictProb(LM2.coxph.fit, 
                                       response = Surv(data.raw.LM2.predict$Ts, 
                                                       data.raw.LM2.predict$status),  
                                       x = data.raw.LM2.predict, 
                                       time  = c(prediction.horizon))

data.raw.LM2.predict$survivalProb = c(condition.prediction.LM2)
data.raw.LM2.predict$True = data.raw.LM2.predict$Ts > prediction.horizon
### AUC
AUC.LM2 = pROC::auc(roc(data.raw.LM2.predict$True, 
                        data.raw.LM2.predict$survivalProb, 
                        levels = c("TRUE", "FALSE")) )

### Brier score
BS.LM2  = mean((as.numeric(data.raw.LM2.predict$True) - data.raw.LM2.predict$survivalProb)^2)


