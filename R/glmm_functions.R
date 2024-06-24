
### functions needed to run GLMM for MWAS
#' calculate_GRM
#' 
#' calculate the GRM from the gene matrix
#' 
#' @param genematrix must have first column as sample names and other columns as gene names
#' @return GRM
#' @export
calculate_GRM<-function(gene_matrix){
  freq_mat_dist_man<-parDist(as.matrix(gene_matrix[,-1]), method="manhattan")/(ncol(gene_matrix[,-1])-1)
  freq_mat_dist_man<-as.matrix(freq_mat_dist_man)
  
  freq_mat_GRM_man<-1-freq_mat_dist_man
  GRM<-as.matrix(freq_mat_GRM_man)
  colnames(GRM)<- as.matrix(gene_matrix[,1])
  rownames(GRM)<-colnames(GRM)
  return(GRM)
  }

getCoefficients<-function(Y, X, W, tau, GRM){
  # Y is working vector Y=alpha X + b
  # n number of samples 
  # N number of covarites
  # inputs are Y (nx1) diagnosis
  #X (nxN) covarities and diagnosis
  # W coefficient of variation 
  # tau is the variance of the residual errors
  # adapted from SAIGE package
  simga=gen_sp_Sigma(W,tau,GRM)# V
  Y=as.vector(Y)
  Sigma_iY=solve(simga,Y) # V^-1 Y
  sigma_X=solve(simga,X) # V^-1 X
  cov_var=Matrix::solve(forceSymmetric(t(X) %*% sigma_X),sparse=TRUE,tol = 1e-10) # (Xt V^-1 X)^-1
  Sigma_iXt=t(sigma_X) # t(V^-1 X)
  SigmaiXtY= Sigma_iXt%*%Y # XtV^-1Y
  alpha = cov_var %*% SigmaiXtY # (Xt V X)^-1 XtVY 
  #b1= tau[2]*GRM %*%solve(simga)%*%(Y-X %*% alpha)
  #epsilon=Y-(X%*%alpha+b)
  epsilon=tau[1] * (t(Sigma_iY) - t(sigma_X %*% alpha)) / as.vector(W)#tau[1] to act on W
 
  eta = as.vector(Y - epsilon) # Y-tau \sigma (Y-X\alpha) 
  
  
  b= eta-X %*% alpha
  #eta=Y-epsilon
  re=list("Sigma_iY"=Sigma_iY,"Sigma_iX"=sigma_X,"cov_var"=cov_var,"alpha"=alpha,"eta"=eta,"b"=b,"epsilon"=epsilon)
}
gen_sp_Sigma<-function(W,tau,kinship){
  ### update kinship with W and tau
  ## value vector is kin
  #kinship is an (nxn) symetric matrix
  # adapted from SAIGE package
  kinship<-as.matrix(kinship)
  dtkin=W^-1 * (tau[1]) # inverse W 
  new_kin = kinship * tau[2]
  diag_new_kin =diag(new_kin)+dtkin
  diag(new_kin)<-diag_new_kin
  new_kin[new_kin<1e-4]=1e-4
  return(as.matrix(new_kin))
}

get_AI_score<-function(Y,X,GRM,W,tau,Sigma_iY,Sigma_iX,cov_var){
  ## get score for finding tau function from supplment of Saige paper
  ## Inputs Y, X, GRM, W, Tau, Sigma_Y, Sigma_X, cov_var
  # adapted from SAIGE package
  Sigma=gen_sp_Sigma(W,tau,GRM)
  Sigma_iXt = t(Sigma_iX) #transpose X
  P=solve(Sigma) - Sigma_iX %*% cov_var %*% Sigma_iXt
  PY1 = P %*% Y# \hat{Y}-\hat(X) (Xt V X)^-1 PY
  APY = GRM %*% PY1 # GRM (\hat{Y}-\hat(X) (Xt V X)^-1) 
  YPAPY = t(PY1) %*% APY# dot product
  YPAPY=YPAPY[1]# dot product
  PA= P %*% GRM
  Trace_P_GRM = sum(diag(PA))
  score1=YPAPY-Trace_P_GRM
  PAPY=P%*% APY
  
  AI = (t(PAPY) %*% APY)# AI=t(Y)%*%P%*%GRM%*%P%*%GRM%*%P%*%Y
  return(list(YPAPY=YPAPY,PY=PY1,Trace_P_GRM=Trace_P_GRM,score1=score1,AI=AI[1]))
}

get_AI_score_quant<-function(Y,X,GRM,W,tau,Sigma_iY,Sigma_iX,cov_var){
  ## get score for finding tau function from supplment of Saige paper
  ## Inputs Y, X, GRM, W, Tau, Sigma_Y, Sigma_X, cov_var
  # adapted from SAIGE package
  n=length(W)
  Sigma=gen_sp_Sigma(W,tau,GRM)
  Sigma_iXt = t(Sigma_iX) #transpose X
  P=solve(Sigma) - Sigma_iX %*% cov_var %*% Sigma_iXt
  diag_P=diag(P)/W
  PY1 = P %*% Y# \hat{Y}-\hat(X) (Xt V X)^-1 PY
  wPY=PY1/W
  YPwPY = t(PY1) %*% wPY
  YPwPY=YPwPY[1]
  APY=GRM %*% PY1
  YPAPY = t(PY1) %*% APY# dot product
  YPAPY=YPAPY[1]# dot product
  PA= P %*% GRM
  Trace_P_GRM =  (sum(solve(Sigma)*GRM)-sum(Sigma_iX*crossprod(GRM,t(cov_var %*% Sigma_iXt))))
  Trace_PW=sum(diag_P)
  score1=YPAPY-Trace_P_GRM#score 1
  score0=YPwPY-Trace_PW #score 0 good
  score_vector=as.matrix(c(score0[1],score1[1]))
  PwPY = P %*% wPY
  PAPY=P%*% APY
  
  
  AI_11 = (t(PAPY) %*% APY)#good
  AI_00=(t(PwPY) %*% wPY) 
  AI_01= (t(PAPY) %*% wPY)
  AI_mat=matrix(c(AI_00[1],AI_01[1],AI_01[1],AI_11[1]),2,2)
  Dtau <- solve(AI_mat, score_vector)
  
  
  
  return(list(YPAPY=YPAPY,PY=PY1,YPwPY=YPwPY,Trace_P_GRM=Trace_P_GRM,Trace_PW=Trace_PW,AI=AI_mat,score_vector=score_vector))
}




ScoreTest_NULL_Model = function(mu, y, X){
  ## score test for null model uses fitted mu, real y, and X
  # adapted from SAIGE package
  mu2=mu*(1-mu)
  V = as.vector(mu2)
  res = as.vector(y - mu)
  XV = t(X * V)
  XVX = t(X) %*% (t(XV))
  XVX_inv = solve(XVX)
  XXVX_inv = X %*% XVX_inv
  XVX_inv_XV = XXVX_inv * V
  S_a =  colSums(X * res)
  re = list(XV = XV, XVX = XVX, XXVX_inv = XXVX_inv, XVX_inv = XVX_inv, S_a = S_a, XVX_inv_XV = XVX_inv_XV, V = V)
  class(re) = "SA_NULL"
  return(re) 
}	

ScoreTest_NULL_Model_quant = function(mu,tau, y, X){
  # adapted from SAIGE package
  V = rep(1/tau[1], length(y))
  res = as.vector(y - mu)
  XV = t(X * V)
  XVX = t(X) %*% (t(XV))
  XVX_inv = solve(XVX)
  XXVX_inv = X %*% XVX_inv
  XVX_inv_XV = XXVX_inv * V
  S_a =  colSums(X * res)
  re = list(XV = XV, XVX = XVX, XXVX_inv = XXVX_inv, XVX_inv = XVX_inv, S_a = S_a, XVX_inv_XV = XVX_inv_XV, V = V)
  class(re) = "SA_NULL"
  return(re) 
}	



fitglmmaiRPCG<-function(Yvec, Xmat,GRM,wVec,  tauVec, Sigma_iY, Sigma_iX, cov_var,tol,quant=FALSE,verbose,write_log){
  # adapted from SAIGE package
  if(!quant){ 
    re.AI = get_AI_score(Yvec, Xmat,GRM,wVec,  tauVec, Sigma_iY, Sigma_iX, cov_var)
    score1 = re.AI$score1# this is equation 8 from paper 
    AI1 = re.AI$AI
    Dtau = score1/AI1
    
    tau0 = tauVec
    
    tauVec[2] = tau0[2] + Dtau
    step = 1.0
    while(tauVec[2]<0){
      
      step = step*0.5
      tauVec[2] = tau0[2] + step*Dtau
      
    }
    
  }else{
    re.AI = get_AI_score_quant(Yvec, Xmat,GRM,wVec,  tauVec, Sigma_iY, Sigma_iX, cov_var)
    YPAPY = re.AI$YPAPY
    YPwPY = re.AI$YPwPY
    Trace_PW=re.AI$Trace_PW
    Trace_P_GRM=re.AI$Trace_P_GRM
    Dtau = solve(re.AI$AI, re.AI$score_vector)
    tau0 = tauVec
    tauVec = tau0 + Dtau
    step = 1.0
    while(any(tauVec<0)){
      
      step = step*0.5
      tauVec = tau0 + step * Dtau
      
    }
  }
  
  
  
  
  if(any(tauVec < tol)){
    tauVec[which(tauVec < tol)]=0
  }
  
  return(list("tau" = tauVec))
}

#' fit_tau_test
#' 
#' Fit the base model for population structure finding random effects
#' 
#' @param glm_fit0 glm model. Model output with no sample relatedness accounted for
#' @param GRM Genetic Relatedness Matrix (from scripts or user) NxN matrix of sample relatedness
#' @param species_id for tracking species
#' @param tau vector for initial values for the variance component parameter estimates usually c(1,1)
#' @param maxiter maximum iterations to fit the glmm model
#' @param verbose whether outputting messages in the process of model fitting
#' @param log_file log file to write to
#' @return model output for the baseline structure glmm 
#' @export
fit_tau_test = function(glm_fit0, GRM,species_id,tau=c(1,1),maxiter =100, verbose = TRUE,tol=.0001,log_file=NA) {
  #Fits the null generalized linear mixed model for a binary trait
  # adapted from SAIGE package
  #Args:
  #  glm_fit0: glm model. Logistic model output (with no sample relatedness accounted for) 
  #  GRM Genetic Relatedness Matrix 
  # Species ID of the species for record
  #  tau: vector for initial values for the variance component parameter estimates
  #  maxiter: maximum iterations to fit the glmm model
  #  verbose: whether outputting messages in the process of model fitting
  #Returns:
  #  model output for the null glmm
  write_log=!is.na(log_file) & log_file!=FALSE
  t_begin = proc.time()
  if(verbose){
    cat("begining time ")
    cat(t_begin)	
    
  }
  if(write_log){
    put("begining time ",console=FALSE)
    put(t_begin,console=FALSE)
  }
  y = glm_fit0$y
  n = length(y)
  X = model.matrix(glm_fit0)
  Xorig= model.matrix(glm_fit0)
  offset = glm_fit0$offset
  if(is.null(offset)){
    offset = rep(0, n)
  }
  
  family = glm_fit0$family
  eta = glm_fit0$linear.predictors
  mu = glm_fit0$fitted.values
  mu.eta = family$mu.eta(eta)
  Y = eta - offset + (y - mu)/mu.eta
  alpha0 = glm_fit0$coef
  eta0 = eta
  sample_ids<-colnames(GRM)
  if(family$family %in% c("poisson", "binomial")) {
    tau[1] = 1
    quant=FALSE
  }else{
    quant=TRUE
  }

  if(verbose) cat(" Fixed-effect coefficients: ", glm_fit0$coef,"\n")
  if(write_log) put(paste(" Fixed-effect coefficients: ", glm_fit0$coef),console=FALSE)
  if(verbose) cat(" inital tau is ", tau,"\n")
  if(write_log) put(paste(" inital tau is ", tau),console=FALSE)
  tau0=tau
  if(tau[1]<=0){
      stop("ERROR! The first variance component parameter estimate is 0\n")
    }
  sqrtW_0 = mu.eta/sqrt(family$variance(mu))
  W_0 = sqrtW_0^2
  re.coef = Get_Coef(y, X, tau, GRM,family, alpha0, eta0,  offset, maxiter=maxiter,verbose=verbose,tol.coef = tol,write_log=write_log)
  
  ######
  if(quant){
    re = get_AI_score_quant(re.coef$Y, X, GRM,re.coef$W, tau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov_var)
    tau[2] = max(0, as.numeric(tau0[2] + tau0[2]^2 * (re$YPAPY - re$Trace_P_GRM)/n))
    tau[1] = max(0, as.numeric(tau0[1] + tau0[1]^2 * (re$YPwPY - re$Trace_PW)/n))
  }else{
    re = get_AI_score(re.coef$Y, X, GRM,re.coef$W, tau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov_var)
    
    tau[2] = max(0, as.numeric(tau0[2] + tau0[2]^2 * ((re$YPAPY - re$Trace_P_GRM))/n)) #tau + Dtau dumb way
  }

  for(i in seq_len(maxiter)){
    if(verbose) {
      cat(paste("\ni",i))
    }
    if(verbose) cat("\nIteration ", i, "tau is ", tau, "\n")
    if(write_log) put(paste(" Iteration ", i, "tau is: ", tau),console=FALSE)
    alpha0 = re.coef$alpha
    tau0 = tau
    #cat("tau0_v1: ", tau0, "\n")
    eta0 = eta
    rss_0=sum((y-mu)^2)
    # use Get_Coef before getAIScore       
    t_begin_Get_Coef = proc.time()
    re.coef = Get_Coef(y, X, tau, GRM,family, alpha0, eta0,  offset,verbose=verbose,maxiter=maxiter,tol.coef = tol,write_log=write_log)
    t_end_Get_Coef =  proc.time()
    if(verbose) {
    cat("\nt_end_Get_Coef - t_begin_Get_Coef\n")
    cat(t_end_Get_Coef - t_begin_Get_Coef)
    }
    if(write_log) {
      put("t_end_Get_Coef - t_begin_Get_Coef",console=FALSE)
      put(t_end_Get_Coef - t_begin_Get_Coef,console=FALSE)
    }
    ##update tau
   
    fit = fitglmmaiRPCG(re.coef$Y, X, GRM, re.coef$W, tau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov_var,tol=tol,verbose=verbose,write_log=write_log,quant=quant)
    t_end_fitglmmaiRPCG= proc.time()
    if(verbose) {
    cat("\nt_end_fitglmmaiRPCG - t_begin_fitglmmaiRPCG\n")
    cat(t_end_fitglmmaiRPCG - t_end_Get_Coef)
    }
    if(write_log) {
      put("t_end_fitglmmaiRPCG - t_begin_fitglmmaiRPCG",console=FALSE)
      put(t_end_fitglmmaiRPCG - t_end_Get_Coef,console=FALSE)
    }
    tau = as.numeric(fit$tau)
    cov_var = re.coef$cov_var
    alpha = re.coef$alpha
    eta = re.coef$eta
    Y = re.coef$Y
    mu = re.coef$mu
    res = y - mu
     if(verbose) {
       cat(paste("\nchange in tau",abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)))
    cat("\ntau: ", tau, "\n")
    cat("\ntau0: ", tau0, "\n")
     }
    if(write_log) {
      put(paste("change in tau",abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)),console=FALSE)
      put(paste("tau: ", tau),console=FALSE)
      put(paste("tau0: ", tau0),console=FALSE)
    }
    if(tau[1]<=0){
      stop("\nERROR! The first variance component parameter estimate is 0\n")
    }
    
    
    if(tau[2] == 0) break
    # Use only tau for convergence evaluation, because alpha was evaluated already in Get_Coef
    tau_condition=max(abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol
    
 
    rss=sum(res^2)
    rss_condition=rss_0-rss
    if(verbose){
      cat(paste("\nres",rss))
      cat(paste("\nrss change",rss_condition))
    }
    if(write_log){
      put(paste("res",rss),console=FALSE)
      put(paste("rss change",rss_condition),console=FALSE)
    }
    
    abs_condition=sum(res^2)
    if(tau_condition) break
    
    if(max(tau) > tol^(-2)) {
      i = maxiter
      break
    }
  }
  if(verbose) cat("\niter break at ",i)
  if(verbose) cat("\nFinal " ,tau, ":\n")
  if(write_log) put(paste("iter break at ",i),console=FALSE)
  if(write_log) put(paste("Final " ,tau, ":"),console=FALSE)
  if(max(tau) > tol^(-2)){
    cat("Model not converged")
    if(write_log){
      put("Model not converged",console=FALSE)
    }
    return(glm_fit0)
  }
  
  re.coef = Get_Coef(y, X, tau, GRM,family, alpha, eta,  offset,verbose=verbose, maxiter=maxiter,tol.coef = tol,write_log=write_log)
  if(quant){
    re.final = get_AI_score_quant(re.coef$Y, X, GRM,re.coef$W, tau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov_var)
    tau[2] = max(0, tau0[2] + tau0[2]^2 * (re.final$YPAPY - re.final$Trace_P_GRM)/n)
    tau[1] = max(0, tau0[1] + tau0[1]^2 * (re.final$YPwPY - re.final$Trace_PW)/n)
  }else{
    re.final = get_AI_score(re.coef$Y, X, GRM,re.coef$W, tau, re.coef$Sigma_iY, re.coef$Sigma_iX, re.coef$cov_var)
    
    tau[2] = max(0, as.numeric(tau0[2] + tau0[2]^2 * ((re.final$YPAPY - re.final$Trace_P_GRM))/n)) #tau + Dtau dumb way
  }
  cov_var = re.coef$cov_var
  
  alpha = re.coef$alpha
  eta = re.coef$eta
  Y = re.coef$Y
  mu = re.coef$mu
  mu.eta = family$mu.eta(eta)
  sqrtW = mu.eta/sqrt(family$variance(mu))
  W=sqrtW^2
  Sigma=gen_sp_Sigma(W,tau,GRM)
  converged = ifelse(i < maxiter, TRUE, FALSE)
  res = y - mu
  rss=sum(res^2)
  if(quant){
    mu2 = rep(1/tau[1],length(y))
    obj.noK = ScoreTest_NULL_Model_quant(mu,tau,y,Xorig)
  }else{
    mu2 = mu * (1-mu)
    obj.noK = ScoreTest_NULL_Model(mu, y, Xorig)
  }
  
  ss=sum((y-mean(y))^2)
  var_fixed=var(Xorig%*%alpha)
  var_random=var(as.vector(re.coef$b))
  var_error=var(res)
  #https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210x.2012.00261.x
  model_metrics=list("S"=sd(res),"R-sq"=(1-rss/ss),"R-sq(marginal)"=var_fixed/(var_fixed+var_random+var_error),"r-sq(conditional)"=(var_fixed+var_random)/(var_fixed+var_random+var_error))
  call=paste("pop_structure_test formula=",glm_fit0$formula,"family=",glm_fit0$family)
  Coefficients=paste("Coefficents:",alpha)
  ###make another option for qaunt
  if(!quant){
    AUC=paste("AUC",suppressMessages(suppressWarnings(auc(glm_fit0$y,as.vector(mu)))))
    tau_script=paste("tau is:",tau[2],"t value is",sum(re.coef$b^2))
    summary_text=paste(call,Coefficients,AUC,tau_script,sep="\n")
  }else{
    tau_script=paste("tau is:",tau[2],"t value is",sum(re.coef$b^2))
    summary_text=paste(call,Coefficients,sep="\n")
  }
  
  glmmResult = list(summary_text=summary_text,tau=tau,
                    coefficients=alpha, b=re.coef$b,t=sum(re.coef$b^2),
                    linear.predictors=eta, linear_model=Xorig%*%alpha+re.coef$b, 
                    fitted.values=mu, var_mu=mu2,Y=Y, residuals=res, 
                    cov_var=cov_var, converged=converged,
                    sampleID = sample_ids, 
                    obj.noK=obj.noK, 
                    y = y, X = Xorig, 
                    traitType=glm_fit0$family,
                    iter_finised=i,
                    model_metrics=model_metrics,species_id=species_id,GRM=GRM)
  
  t_end_null = proc.time()
  if(verbose) {
  cat("\nt_end_null - t_begin,  fitting the structure model took\n")
  cat(t_end_null - t_begin)
  cat("\n")
  }
  
  if(write_log) {
    put("t_end_null - t_begin, fitting the structure model took",console=FALSE)
    put(t_end_null - t_begin,console=FALSE)
  }
  
  return(glmmResult)
}


Get_Coef = function(y, X, tau, GRM,family, alpha0, eta0,  offset, verbose=FALSE,maxiter,tol.coef=tol,write_log=FALSE){
  # adapted from SAIGE package
  mu = family$linkinv(eta0)
  mu.eta = family$mu.eta(eta0)
  Y = eta0 - offset + (y - mu)/mu.eta
  
  sqrtW = mu.eta/sqrt(family$variance(mu))
  
  W = sqrtW^2
  for(i in 1:maxiter){
    re.coef = getCoefficients(Y, X, W, tau, GRM)
    
    alpha = as.matrix(re.coef$alpha)
    eta = as.matrix(re.coef$eta + offset)
    
    if(verbose) {
      cat("\n Tau:\n")
      cat(tau)
      cat("\n Fixed-effect coefficients:\n")
      cat(alpha)
    }
    if(write_log) {
      put(paste(" Tau:",tau," Fixed-effect coefficients:",alpha),console = FALSE)
    }
    mu = family$linkinv(eta)
    mu.eta = family$mu.eta(eta)
    
    Y = eta - offset + (y - mu)/mu.eta
    sqrtW = mu.eta/sqrt(family$variance(mu))
    W = sqrtW^2
    
    if( max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol.coef))< tol.coef){
      break
    }
    alpha0 = alpha
  }
  
  re = list(Y=Y, alpha=alpha, eta=eta, W=W, cov_var=re.coef$cov_var, sqrtW=sqrtW, Sigma_iY = re.coef$Sigma_iY, Sigma_iX = re.coef$Sigma_iX, mu=mu,eta_2=re.coef$eta_2,b=re.coef$b)
}

re_fit_glmm<-function(obj.pop.strut,glm_fit0,GRM){
  tau=obj.pop.strut$tau
  y=glm_fit0$y
  X=obj.pop.strut$X
  family = glm_fit0$family
  eta = obj.pop.strut$linear.predictors
  mu = obj.pop.strut$fitted.values
  mu.eta = family$mu.eta(eta)
  offset=rep(0,length(y))
  Y = eta - offset + (y - mu)/mu.eta
  alpha = obj.pop.strut$coefficients 
  re.coef = Get_Coef(y, X, tau, GRM,family, alpha, eta,  offset=offset,verbose=FALSE, maxiter=30,tol.coef = .01,write_log=FALSE)
  alpha = re.coef$alpha
  eta = re.coef$eta
  Y = re.coef$Y
  mu = re.coef$mu
  mu.eta = family$mu.eta(eta)
  obj.noK=ScoreTest_NULL_Model(mu,y,X)
  mu2=mu*(1-mu)
  res = y - mu
  converged=TRUE
  cov_var = re.coef$cov_var
  sample_ids=obj.pop.strut$sampleID
  rss=sum(res^2)
  ss=sum((y-mean(y))^2)
  var_fixed=var(X%*%alpha)
  var_random=var(as.vector(re.coef$b))
  var_error=var(res)
  #https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210x.2012.00261.x
  model_metrics=list("S"=sd(res),"R-sq"=(1-rss/ss),"R-sq(marginal)"=var_fixed/(var_fixed+var_random+var_error),"r-sq(conditional)"=(var_fixed+var_random)/(var_fixed+var_random+var_error))
  tau_script=paste("tau is:",tau[2],"t value is",sum(re.coef$b^2))
  glmmSNPResult = list(summary_text=tau_script,tau=tau,
                       coefficients=alpha, b=as.vector(re.coef$b),t=sum(re.coef$b^2),
                       linear.predictors=eta, linear_model=X%*%alpha+re.coef$b, 
                       fitted.values=mu, var_mu=mu2,Y=Y, residuals=res, 
                       cov_var=cov_var, converged=converged,
                       sampleID = sample_ids, 
                       obj.noK=obj.noK, 
                       y = y, X = X, 
                       traitType=glm_fit0$family,
                       iter_finised=30,
                       model_metrics=model_metrics,species_id=obj.pop.strut$species_id,GRM=GRM)
  return(glmmSNPResult)

}

#' fit_beta
#' 
#' Fit a beta for each genes including the population structure model with the random effects
#' 
#' @param obj.pop.strut output of fit_tau_test; GLMM of species with GRM accounted for
#' @param glm_fit0 glm model. Model output with no sample relatedness accounted for
#' @param GRM Genetic Relatedness Matrix (from scripts or user) NxN matrix of sample relatedness
#' @param gene_df long data frame with, gene_id, sample_name, and gene_value 
#' @param SPA whether to run Saddle point approximation for pvalues (will slow down output) 
#' @return model output for the baseline structure glmm 
#' @export
fit_beta = function(obj.pop.strut,
                    glm_fit0,GRM,
                    gene_df,SPA=FALSE){
  ## inputs: full null fitted model
  ## null model no random effect
  # GRM is kinship
  # gene_df is the copy number by gene 
  ## gene_df must have column for sample; column for gene_id; column for gene_value
  ## returns list of values for each gene examined
  # adapted from SAIGE package
  list_vec<-NULL
  t_begin = proc.time()

  
  obj.noK = obj.pop.strut$obj.noK
  family = glm_fit0$family
  
  eta = obj.pop.strut$linear.predictors
  mu = obj.pop.strut$fitted.values
  mu.eta = family$mu.eta(eta)
  sqrtW = mu.eta/sqrt(glm_fit0$family$variance(mu))
  W1 = sqrtW^2   ##(mu*(1-mu) for binary) theses are the same
  #tauVecNew = obj.pop.strut$tau
  #X = obj.pop.strut$X
  #Sigma=gen_sp_Sigma(W1,tauVecNew,GRM)
  #obj.pop.strut$Sigma<-Sigma
  #Sigma_iX<-solve(Sigma,X)
  #obj.pop.strut$Sigma_iX<-Sigma_iX
  #Sigma_iY<-solve(Sigma,obj.pop.strut$Y)
  #obj.pop.strut$Sigma_iY<-Sigma_iY
  sample_lookup<-data.frame(sampleID=obj.pop.strut$sampleID,index=seq(1,length(obj.pop.strut$sampleID)),y=obj.pop.strut$y)
  sample_genes<-unique(gene_df$gene_id)

  for(k in sample_genes){
    iter=which(sample_genes==k)
    if(iter %% 1000 == 0){
      cat(paste("number of genes done ",iter,"\n"))
      cat("time past:")
      t_now = proc.time()
      cat(t_now-t_begin)	
      cat("\n")
    }
    one_gene<-gene_df %>% ungroup %>% filter(gene_id==k)
    #one_gene_indexs<-sample_lookup %>% inner_join(one_gene,by=c("sampleID"="sample_name")) %>% select(sampleID,index)
    
    one_gene<-one_gene %>% inner_join(sample_lookup,by=c("sample_name"="sampleID"))
    ## filter obj to samples present in gene copy number
    #filtered_obj.pop.strut<-filter_null_obj(obj.pop.strut,one_gene)
    empty_mat<-matrix(0,nrow(one_gene),nrow(one_gene))
   
    G0<-as.vector(one_gene$gene_value)
    G_tilde = G0  -  obj.pop.strut$obj.noK$XXVX_inv %*%  (obj.pop.strut$obj.noK$XV %*% G0) # G1 is X adjusted
    #res=filtered_obj.pop.strut$residuals
    #eta = filtered_obj.pop.strut$linear.predictors
    
    #mu = filtered_obj.pop.strut$fitted.values
    #mu.eta = family$mu.eta(eta)
    #sqrtW = mu.eta/sqrt(glm_fit0$family$variance(mu))
    #W=sqrtW^2#mu*(1-mu)
    #Sigma_iG = solve(obj.pop.strut$Sigma,G_tilde)
    #PG_tilde<-Sigma_iG-obj.pop.strut$Sigma_iX%*%(solve(t(obj.pop.strut$X)%*%obj.pop.strut$Sigma_iX))%*%t(obj.pop.strut$X)%*%Sigma_iG
    Y = eta  + (obj.pop.strut$y - mu)/mu.eta
    
    #print(t(G_tilde)%*%PG_tilde)
    #t_score=t(PG_tilde)%*%(Y)/tauVecNew[1]
    #print(t_score)#these are the same!
    t_score=t(G_tilde)%*%(obj.pop.strut$y - mu)
    m1 = sum(mu * G_tilde)
    
    #qtilde=t(G_tilde)%*%filtered_obj.pop.strut$y
    var1<-sum(W1 * G_tilde^2)
    t_adj_2=(t_score^2)/var1
    beta=t_score/var1
    pval=(pchisq(t_adj_2, lower.tail = FALSE, df=1,log.p=FALSE))
    z=(qnorm(pval/2, log.p=F, lower.tail = F))
    se_beta=abs(beta)/sqrt(abs(z))
    new_eta=beta[1,1]*G0+obj.pop.strut$b+obj.pop.strut$X %*% obj.pop.strut$coefficients 
    qtilde=t_score+m1
    if(SPA){
      if(var1<0){
        list_vec=rbind(list_vec,data.frame("species_id"=obj.pop.strut$species_id,tau=obj.pop.strut$tau[2],"gene_id"=k,"cor"=cor(G0,obj.pop.strut$y),"cor_to_b"=cor(as.numeric(obj.pop.strut$b),G0),"z"=z,"var1"=var1,"beta"=NA,"se beta"=NA,"pvalue"=NA,
                                           "t_adj"=NA,"num_control"=sum(obj.pop.strut$y==0),
                                           "num_total"=length(G0),
                                           SPA_pvalue=NA,spa_score=NA,SPA_zvalue=NA,pvalue_noadj=NA,converged=FALSE))
      }else{
        out1 = Saddle_Prob(q=qtilde, mu = mu, g = G_tilde, var1,Cutoff = 2,log.p=FALSE)
        list_vec=rbind(list_vec,data.frame("species_id"=obj.pop.strut$species_id,tau=obj.pop.strut$tau[2],"gene_id"=k,"cor"=cor(G0,obj.pop.strut$y),"cor_to_b"=cor(as.numeric(obj.pop.strut$b),G0),"z"=z,"var1"=var1,"beta"=beta,"se beta"=se_beta,
                                           "t_adj"=t_adj_2,"num_control"=sum(obj.pop.strut$y==0),
                                           "num_total"=length(G0),
                                           SPA_pvalue=out1$p.value,spa_score=out1$Score,SPA_zvalue=out1$z_value,pvalue_noadj=out1$p.value.NA,converged=out1$Is.converge))
      }
   
    }else{
      #no SPA
      list_vec=rbind(list_vec,data.frame("species_id"=obj.pop.strut$species_id,tau=obj.pop.strut$tau[2],"gene_id"=k,"cor"=cor(G0,obj.pop.strut$y),"cor_to_b"=cor(obj.pop.strut$b,obj.pop.strut$y),"z"=z,"var1"=var1,"beta"=beta,"se beta"=se_beta,"pvalue"=pval,
                                         "t_adj"=t_adj_2,"num_control"=sum(obj.pop.strut$y==0),
                                         "num_total"=length(G0)))
      
    }
  }
  cat("total time past:")
  t_end = proc.time()
  cat(t_end-t_begin)
  cat("\n")
  return(list_vec)
}

filter_null_obj<-function(obj.pop.strut,sample_indexs){
  obj.pop.strut$gene_value<-sample_indexs$gene_value
  obj.pop.strut$residuals<-obj.pop.strut$residuals[sample_indexs$index]
  obj.pop.strut$b<-obj.pop.strut$b[sample_indexs$index]
  obj.pop.strut$linear.predictors<-obj.pop.strut$linear.predictors[sample_indexs$index]
  obj.pop.strut$fitted.values<-obj.pop.strut$fitted.values[sample_indexs$index]
  obj.pop.strut$obj.noK$XXVX_inv<-obj.pop.strut$obj.noK$XXVX_inv[sample_indexs$index,]
  obj.pop.strut$obj.noK$XV<-obj.pop.strut$obj.noK$XV[,sample_indexs$index]
  obj.pop.strut$Sigma<- obj.pop.strut$Sigma[sample_indexs$index,sample_indexs$index]
  obj.pop.strut$X<-obj.pop.strut$X[sample_indexs$index,]
  obj.pop.strut$y<-obj.pop.strut$y[sample_indexs$index]
  obj.pop.strut$Sigma_iX<-obj.pop.strut$Sigma_iX[sample_indexs$index,]
  obj.pop.strut$Sigma_iY<-obj.pop.strut$Sigma_iY[sample_indexs$index]
  return(obj.pop.strut)
}


simulate_tau_inner<-function(glm_fit0,GRM,species_id="s_id",tau0,phi0){
  family_to_fit<-glm_fit0$family
  data.new<-glm_fit0$data
  formulate_to_fit<-glm_fit0$formula
  data.new_shuffled<-data.new[sample(1:nrow(data.new),nrow(data.new)),]

  fit_logistic = glm(formulate_to_fit, data = data.new_shuffled, family = family_to_fit)
  fit_glmm_snp<-tryCatch(fit_tau_test(fit_logistic,GRM,tau=c(phi0,tau0),verbose = FALSE,species_id=species_id, log_file =NA),error=function(e) e)
  if(!is.na(fit_glmm_snp$t)){
    t=sum(fit_glmm_snp$b^2,na.rm=TRUE)
    tau=fit_glmm_snp$tau[2]
  }else{
    print("error")
    tau=0
    t=0
  }
  return(data.frame("tau"=tau,t))
}


#' run_tau_test
#' 
#' take output from population structure test and test how exterme the tau is for that set of data
#' 
#' @param glm_fit0 glm model. Model output with no sample relatedness accounted for
#' @param GRM Genetic Relatedness Matrix (from scripts or user) NxN matrix of sample relatedness
#' @param n_tau number of tau to simulate
#' @param species_id species id for bactrial species
#' @param tau0 starting tau
#' @param phi0 starting phi
#' @return df of values of T for tau for different runs
#' @export
run_tau_test<-function(glm_fit0,GRM,n_tau,species_id="s_id",tau0,phi0){
  list_of_tau<-lapply(seq(1,n_tau),function(x) simulate_tau_inner(glm_fit0,GRM,species_id=species_id,tau0,phi0))
  df_of_tau<-do.call(rbind, list_of_tau)
  return(df_of_tau)
}






Saddle_Prob<-function(q, mu, g, var1,Cutoff=2,output="P",log.p=FALSE)
{
  #### taken from ‘SPAtest’ with a few changes for use case
  m1<-sum(mu * g)
  var1<-sum(mu * (1-mu) * g^2)
  p1=NULL
  p2=NULL
  
  Score<- q-m1
  #
  
  qinv = -sign(q-m1) * abs(q-m1) + m1
  t_adj_sq=((q - m1)^2)/var1
  # Noadj
  pval.noadj<-pchisq(t_adj_sq, lower.tail = FALSE, df=1,log.p=FALSE)
  Is.converge=TRUE
  
  if(is.na(abs(q - m1)) || is.na(var1) || var1<0){
    Is.converge=FALSE
    pval=pval.noadj	
  }
  if(isTRUE(abs(q - m1)/sqrt(var1) < Cutoff )){
    
    pval=pval.noadj	
  } else {
    out.uni1<-getroot_K1(0, mu=mu, g=g, q=q)
    out.uni2<-getroot_K1(0, mu=mu, g=g, q=qinv)
    if(out.uni1$Is.converge==TRUE && out.uni2$Is.converge==TRUE)
    {
      p1<-tryCatch(Get_Saddle_Prob(out.uni1$root, mu, g, q,log.p=log.p),error=function(e) {
        if(log.p) return(pval.noadj-log(2))
        else return(pval.noadj/2)})
      p2<-tryCatch(Get_Saddle_Prob(out.uni2$root, mu, g, qinv,log.p=log.p),error=function(e) {
        if(log.p) return(pval.noadj-log(2))
        else return(pval.noadj/2)})
      if(log.p)
      {
        pval = add_logp(p1,p2)
      } else {

        pval = abs(p1)+abs(p2)
      }
      Is.converge=TRUE
    } else {
      cat("Error_Converge")
      pval<-pval.noadj
      Is.converge=FALSE	
    }				
  }
  z_value=qnorm(pval/2,log.p=F, lower.tail = F)*sign(q-m1)
  if(pval!=0 && pval.noadj/pval>10^3)
  {
    return(Saddle_Prob(q, mu, g,var1, Cutoff=Cutoff*2,output,log.p=log.p))
   }else if(pval==0){
    return(list(p.value=pval, p.value.NA=pval.noadj,z_value=z_value, Is.converge=FALSE, Score=Score))
  } 
  else {
    return(list(p.value=pval, p.value.NA=pval.noadj, z_value=z_value,Is.converge=Is.converge, Score=Score))
  }
}
getroot_K1<-function(init,mu,g,q,m1,tol=.Machine$double.eps^0.25,maxiter=1000)
{
  #### taken from ‘SPAtest package’
  g.pos<-sum(g[which(g>0)])
  g.neg<- sum(g[which(g<0)])
  if(q>=g.pos || q<=g.neg)
  {
    return(list(root=Inf,n.iter=0,Is.converge=TRUE))
  } else {
    t<-init
    K1_eval<-K1_adj(t,mu,g,q)
    prevJump<- Inf
    rep<-1
    repeat
    {
      K2_eval<-K2(t,mu,g)
      tnew<-t-K1_eval/K2_eval
      if(is.na(tnew))
      {
        conv=FALSE
        break
      }
      if(abs(tnew-t)<tol)
      {
        conv<-TRUE
        break
      }
      if(rep==maxiter)
      {
        conv<-FALSE
        break
      }
      
      newK1<-K1_adj(tnew,mu,g,q)
      if(sign(K1_eval)!=sign(newK1))
      {
        if(abs(tnew-t)>prevJump-tol)
        {
          tnew<-t+sign(newK1-K1_eval)*prevJump/2
          newK1<-K1_adj(tnew,mu,g,q)
          prevJump<-prevJump/2
        } else {
          prevJump<-abs(tnew-t)
        }
      }
      
      rep<-rep+1
      t<-tnew
      K1_eval<-newK1
    }
    return(list(root=t,n.iter=rep,Is.converge=conv))
  }
}

Korg<-function(t, mu, g){
  #### taken from ‘SPAtest’ 
  n.t<-length(t)
  out<-rep(0,n.t)
  
  for(i in 1:n.t){
    t1<-t[i]
    temp<-log(1 - mu + mu * exp(g* t1))
    out[i]<-sum(temp)
  }
  return(out)
}

Get_Saddle_Prob<-function(zeta, mu, g, q,log.p=FALSE) 
{
  #### taken from ‘SPAtest package’
  k1<-Korg(zeta, mu, g)
  k2<-K2(zeta, mu, g)

  if(is.finite(k1) && is.finite(k2))
  {
    temp1<-zeta * q - k1
    
    
    w<-sign(zeta) * (2 *temp1)^{1/2}
    v<- zeta * (k2)^{1/2}
    
    Z.test<-w + 1/w * log(v/w)	
    
    
    if(Z.test > 0){
      pval<-pnorm(Z.test, lower.tail = FALSE,log.p=log.p)
    } else {
      pval= -pnorm(Z.test, lower.tail = TRUE,log.p=log.p)
    }	
  } else {
    if(log.p)
    {
      pval<- -Inf
    } else {
      pval<-0
    }
  }
  
  return(pval)
}
K1_adj<-function(t, mu, g, q)
{
  #### taken from ‘SPAtest package’
  n.t<-length(t)	
  out<-rep(0,n.t)
  
  for(i in 1:n.t){
    t1<-t[i]
    temp1<-(1 - mu)* exp(-g * t1) + mu
    temp2<-mu *g
    out[i]<-sum(temp2/temp1)-q
  }
  return(out)
}
K2<-function(t, mu, g)
{
  n.t<-length(t)
  out<-rep(0,n.t)
  
  for(i in 1:n.t){
    t1<-t[i]
    temp1<-((1 - mu)* exp(-g * t1) + mu)^2
    temp2<-(1-mu) * mu * g^2 * exp(-g*t1)
    out[i]<-sum(temp2/temp1,na.rm=TRUE)
  }
  return(out)
}

