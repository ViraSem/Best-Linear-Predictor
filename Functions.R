est_signal<-function(data) {
  # X - covariates
  X<-data$X
  # outcome
  y<-data$y
  # Number of covariates
  p<-dim(X)[2]
  # Presence indicator
  D<-data$D
  # Train/test split
  inds.train<-data$inds.train
  inds.test<-setdiff(1:N,inds.train)
  # propensity score estimate
  # initialize s.hat variable
  s.hat<-rep(0,N)
  # estimate on the first half
  p.score.fit.train<-gamlr(x=X[inds.train,],y=D[inds.train],family="binomial",standardize=TRUE)
  # estimate on the second half
  p.score.fit.test<-gamlr(x=X[inds.test,],y=D[inds.test],family="binomial",standardize=TRUE)
  
  # predict using the first half
  s.hat[inds.test]<-predict(p.score.fit.train,newdata=X[inds.test,],type="response")
  # predict using the second half
  s.hat[inds.train]<-predict(p.score.fit.test,newdata=X[inds.train,],type="response")
  # regression function estimate
  g.hat<-rep(0,N)
  
  ## inds = the complete observations (D=1) in the first half
  inds<-intersect(inds.train,(1:N)[D==1])
  ## fit Lasso on inds
  regression.fit.train<-gamlr(x=X[inds.train,],y=y[inds.train],family="gaussian",standardize=TRUE)
  ## use first half estimate to predict the second half
  g.hat[inds.test]<-predict(regression.fit.train,X[inds.test,],type="response")
  ## inds = the complete observations (D=1) in the second half
  inds<-intersect(inds.test,(1:N)[D==1])
  ##  fit Lasso on inds
  regression.fit.test<-gamlr(x=X[inds.test,],y=y[inds.test],family="gaussian",standardize=TRUE)
  ## use first half estimate to predict the first half
  g.hat[inds.train]<-predict(regression.fit.test,X[inds.train,],type="response")
  
  ## Doubly Robust (Orthogonal) Estimate of the signal
  y.est<-g.hat+D/s.hat*(y-g.hat)
  est<-list(y.est=y.est,s.hat=s.hat,g.hat=g.hat)
  return(est) 
}
# least squares
ols_est<-function(X,y,...) {
  return(solve(t(X)%*%X)%*%t(X)%*%y)
}
# white standard error
white_sd<-function(X,y,b.hat,...) {
  # take subset of data
  X<-X[,1:p.OLS]
  # compute the estimated residuals
  error.sq<-(y-X%*%b.hat)^2
  # 
  Q<-t(X)%*%X/N
  Sigma<-matrix(0,p.OLS,p.OLS)
  for (i in 1:N) {
    Sigma<-Sigma+X[i,]%*%t(X[i,])*error.sq[i]
  }
  Sigma<-Sigma/N
  Omega.hat<-solve(Q)%*%Sigma%*%solve(Q)
  return( sqrt(diag(Omega.hat)))
}
# studentized(normalized) z-stat
student<-function(beta.hat,beta.0,st.error) {
  return((beta.hat-beta.0)/st.error)
}

simulate<-function(N,p,seed=888,c.gamma=1, p.OLS=6,rho=0.01) {
  # set random number generator
  set.seed(seed)
  # initialize R.square
  R.sq<-NULL
  #### generate covariates X ~ N(0,Sigma)
      # Sigma = Toeplitz(rho)
      # first row of Toeplitz(rho)
       first.row<-rho^(0:(p-1))
       # Covariance matrix 
      Sigma<-toeplitz(first.row)
      # Uncorrelated covariates N(0,I)
      X.raw<-matrix(rnorm(N*p),N,p)
      # Covariates N(0, \Sigma)
      X<-t(chol(Sigma)%*%t(X.raw))
  
  #### true nuisance parameter: conditional expectation f.X = X gamma and propensity score s.x
      #f. X = X%*%gamma
      gamma<-1/(1:p)
      # down/up scale misspecification error r_g(X) = f.X - p(X)%*%beta0
      gamma[(p.OLS+1):p]<-gamma[(p.OLS+1):p]*c.gamma
      #gamma[30:p]<-0
 
      # mu(Z) = Z%*%gamma: regression function
       f.X<-X%*%gamma
       # Missingness pattern
       delta<-1/(1:p1)^2
       s.x<-plogis(X[,1:p1]%*%delta)
       
  ### Generate presence indicator D and outcome Y
      xi<-rnorm(N)
      # Outcome variable
      y<-f.X+xi
      # Best Linear Predictor
      beta.0<-solve(Sigma[1:p.OLS,1:p.OLS])%*%Sigma[1:p.OLS,]%*%gamma
      # R squared
      R.sq<-  (t(beta.0)%*%Sigma[1:p.OLS,1:p.OLS]%*%beta.0)/(t(gamma)%*%Sigma%*%gamma+1)
      beta.0<-c(beta.0,rep(0,p-p.OLS))

      R.sq<-(t(beta.0[1:p.OLS])%*%Sigma[1:p.OLS,1:p.OLS]%*%beta.0[1:p.OLS])/sd(y)^2
      eps<-y-X%*%beta.0
 
     D<-sapply(s.x,rbinom,size=1,n=1)
    # Training/ Test split
    inds.train<-sample(1:N,floor(N/2))
 
  
  
  return(list(X=X,
              y=y,
              beta.0=beta.0,
              eps=eps,
              D=D,
              inds.train=inds.train,
              f.X=f.X,
              s.x=s.x,
              R.sq=R.sq
              ))
  
}

# AUXILIARY ---------------------------------------------------------------

output_result<-function(
  true.value,filename,directoryname,
  methods,add_info,
  alpha = 0.05,
  digs=3, measures=c("Bias","St.Error","RMSE","Rej.Freq"),
  N_real_rep=N_rep, runtime, which_coefs=1:p.OLS) {
  
  # b.hat<-res$b.hat
  # st.error.hat<-res$st.error.hat
  # norm.b.hat<-res$norm.b.hat
  # which_coefs<-res$which_coefs
  ####show(b.hat)
  
  Results<-array(0,c(length(true.value),length(methods),length(measures))) 
  dimnames(Results)<-list(names(true.value),methods,measures)
  
  
  for (met in methods) {
    # estimate
    b<-matrix(b.hat[[met]][which_coefs,1:N_real_rep],nrow=length(which_coefs))
    # st.error
    st.error<-matrix(st.error.hat[[met]][which_coefs,1:N_real_rep],nrow=length(which_coefs))
    # rejection frequency
    rej.freq<-matrix(norm.b.hat[[met]][which_coefs,1:N_real_rep],nrow=length(which_coefs))
    # Results table fill in
    Results[,met,"Bias"]<-apply(b,1,mean)-true.value
    Results[,met,"St.Error"]<-apply(b,1,sd)
    Results[,met,"RMSE"]<-sqrt((Results[,met,"Bias"])^2+(Results[,met,"St.Error"])^2)
    Results[,met,"Rej.Freq"]<-as.numeric(apply(abs(rej.freq)>qnorm(1-alpha/2),1,mean))
  }
    add_info_caption<-c()
    for (j in 1:length(add_info)) {
      add_info_caption<-c(add_info_caption,paste (names(add_info)[j],add_info[j][[1]]))
    
   }
    add_info_caption<-paste0(add_info_caption,collapse=".")  
    Results<-matrix(Results,nrow=length(true.value),ncol = length(measures)*length(methods))

    #rown<-paste0("n=", N, ",dim.Z=",dim.Z,",\rho=", rho,",p.cat=",p.cat,",t=",decay)
    #rown<-paste0("$(",rown,")$")
    rown<-""
    rnames<-true.value
    for (j in 1:length(true.value)) {
      if (j==1 ) {
        rnames[j]<-paste0(rown,",true ",j, " coef ",round(true.value[j],1))
      } else {
        rnames[j]<-paste0("true ",j, " coef ",round(true.value[j],1))
      }
    
  }
  
  rownames(Results)<-rnames
  colnames(Results)<-rep(methods,length(measures))
  caption<-paste0(paste0(measures,collapse=", "),'.',
                  paste0(methods,collapse=", "), '.', 
                  add_info_caption)
  label<-paste0(paste0(measures,collapse=", "))
  to.Latex(digs=digs,Results,cap=caption,lab=label,filename=filename,k=length(methods),step=length(measures),measures)

  
 
}

#### Auxiliary functions
#### Create latex table "Bias,St.error,RMSE,Rej.freq. of different estimators
to.Latex<-function (digs=3,matrix,cap="Caption here",lab="Description here",filename,k=length(methods),
                    step=4,measures=c("Bias","St.Error","RMSE","Rej.Freq")) {
  # step<-length(measures)
  #options(digits=12)
  align<-rep("r",k)
  align<-paste0(align,collapse="")
  align<-paste0(align,"|")
  align<-rep(align,step)
  align<-paste0(align,collapse="")
  align<-paste0("r|",align)
  tab<-xtable(matrix,digits=digs,caption=cap,label=lab,align=align)
  
  addtorow <- list()
  addtorow$pos <- list(0)
  addtorow$command <- paste0(paste0('& \\multicolumn{',as.character(k),'}{c}{', measures, '}', collapse=''), '\\\\')
  
  print(tab,file=filename,append=TRUE, 
        include.rownames=TRUE,
        include.colnames=TRUE,
        sanitize.text.function=function(x){x},
        add.to.row=addtorow)
  
  
}
