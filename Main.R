# Simulation exercise for "Simultaneous Inference for Conditional Average Treatment Effect and Other Structural Functions"

## Authors: Victor Chernozhukov, Vira Semenova
## Corresponding author: Vira Semenova vsemen@mit.edu

## Abstract: Monte Carlo simulations for two stage estimation of the Best Linear Predictor 
##           in a misspecified linear model Y = p (X)'beta_0 + r(X) + epsilon, 
##           whose outcome Y is partilly missing.

##           control vector Z
##           outcome: Y = g(Z) + eps =  Z %*%gamma + eps, gamma_j = 1/j^2, eps ~ N(0,1)
##           presence indicator: D = B(L(Z%*%delta )), L(t) is logistic function
##                                   B(p) is a bernoulli draw with probability p
##           The technical regressors p(X) = Z[,1:d] is the first d components of control vector Z

## Description: simulate data (Y,D,Z) using simulate_data
#               estimate the doubly robust signal Y.hat = g.hat + D/s.hat (y  - g.hat)
#               run OLS on complete data (D=1), IPW, and Orthogonal Series using cross-fitting with (K=2) folds
#               compute bias, st.error, mse, rejection frequency.

rm(list=ls())
# Set a path to your folder
my_folder<-"/Users/virasemenora/Dropbox (MIT)/18.657/Draft/Best Linear Predictor/R_code"
# Latex Table with output
filename<-"Results.tex"
# Libraries
library(glmnet)
library(plyr)
library(xtable)
library(plotrix)
library(gamlr)

setwd(my_folder)
source("Functions.R")
# Number of observations per draw
N<-1000
# Dimension of the controls Z in main regression
p<-1000
# Dimension of the controls Z in propensity score regression
p1<-100
# Dimension of BLP
p.OLS<-6

# Number of Monte Carlo repetitions
N_rep<-300
## Objects to compute:

# estimator
b.hat<-list()
# standard error
st.error.hat<-list()
# studentized value: (b.hat - true_b)/st.error.hat
norm.b.hat<-list()
# Methods to compare:
### OLS: Ordinary Least Squares on Complete Data only
### IPW: Weighting by Inverse Propensity Score
### OS:  Orthogonal Series (see paper) y.hat = g.hat + D/s.hat ( y - g.hat)
methods<-c("OLS","IPW","OS")
# Method Names
Name<-as.list(c("OLS",
                "IPW",
                "Orthogonal Series"))
names(Name)<-methods



Results<-list()
for (method in methods) {
  b.hat[[method]]<-matrix(0,p.OLS,N_rep)
  st.error.hat[[method]]<-matrix(0,p.OLS,N_rep)
  norm.b.hat[[method]]<-matrix(0,p.OLS,N_rep)
  
}

## c.gamma govers misspecification
## c.gamma = 0: no misspecification, OLS is unbiased
## c.gamma ->0: model is misspecified, OLS is biased
for (c.gamma in c(0.01,20)) {
  for (num in 1:N_rep) {
    show(num)
    set.seed(num) 
    ### Generate data
    data.s<-simulate(N=N,p=p,seed=num,
                     c.gamma=c.gamma,rho=0.8)
   
    ### Robust Estimate of Signal for Orthogonal Series
    est<-est_signal(data=data.s)
    ### Estimate of the propensity score
    s.hat<-est$s.hat
    ### Technical regressors p(X)
    p.X<-data.s$X[,1:p.OLS]
      
    method<-"OLS"
    # weighting matrix (D = 1: Y is observed; Y = 0: is unobserved )
    D<-matrix(rep(data.s$D,p.OLS),nrow=N,ncol=p.OLS)
    # observed covariates
    p.X.D<-p.X*D
    # OLS estimate on complete data
    b.hat[[method]][,num]<-ols_est(p.X.D,data.s$D*data.s$y)
    # st.error
    st.error.hat[[method]][,num]<-white_sd(p.X.D,data.s$D*data.s$y,b.hat =   b.hat[[method]][,num])
    # z-stat = (b.hat - b.true)/ st.error
    norm.b.hat[[method]][,num]<-sqrt(N)*student(b.hat[[method]][,num],data.s$beta.0[1:p.OLS],st.error.hat[[method]][,num])
 
    
    method<-"IPW"
    #weighting matrix (D = 1/sqrt(s.hat): Y is observed; Y = 0: is unobserved )
    D<-matrix(rep(data.s$D/sqrt(s.hat),p.OLS),nrow=N,ncol=p.OLS)
    p.X.D<-p.X*D
    y.D<-data.s$D*data.s$y/sqrt(s.hat)
    b.hat[[method]][,num]<-ols_est(p.X.D,y.D)
    st.error.hat[[method]][,num]<-white_sd(p.X.D,y.D,b.hat =   b.hat[[method]][,num])
    norm.b.hat[[method]][,num]<-sqrt(N)*student(b.hat[[method]][,num],data.s$beta.0[1:p.OLS],st.error.hat[[method]][,num])
    
 
    method<-"OS"
    y.est<-est$y.est
    b.hat[[method]][,num]<-ols_est(p.X,y.est)
    st.error.hat[[method]][,num]<-white_sd(p.X,y.est,b.hat =   b.hat[[method]][,num])
    norm.b.hat[[method]][,num]<-sqrt(N)*student(b.hat[[method]][,num],data.s$beta.0[1:p.OLS],st.error.hat[[method]][,num])
    
      
    
       
  }
}
  #### Output tex table with bias, mse, error, and rejection frequency
  
  # Info to be added to the output table
  add_info<-c(N,p,p.OLS,N_rep,data.s$R.sq,c.gamma)
  names(add_info)<-c("Number of observations (N) ","Number of covariates (d) ",
                     "Dimension of BLP (p)",
                     "Number of MC repetitions",
                     "R.squared",
                     "C.gamma")
  # Output result to Latex file
  output_result(
              true.value=data.s$beta.0[1:p.OLS],
              filename=filename,
              directoryname=directoryname,
              methods=methods,
              add_info=add_info,
              alpha = 0.1,
              digs=3,N_real_rep = N_real_rep)

