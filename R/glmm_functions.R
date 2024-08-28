### Functions needed to run mircoSLAM, software for performing metagenome-wide association studies (MWAS)
### The approach is to fit generalized linear mixed effects models (GLMMs) and use them to test for trait associations
###
### Code for fitting GLMMs (REML AI) is adapted from the SAIGE package (GNU GPL):
###   https://saigegit.github.io/SAIGE-doc/
###  Methodology and equations are described in the supplement of this paper: 
###   Zhou W, Nielsen JB, Fritsche LG, Dey R, Gabrielsen ME, Wolford BN, LeFaive J,
###   VandeHaar P, Gagliano SA, Gifford A, Bastarache LA, Wei WQ, Denny JC, Lin M, 
###   Hveem K, Kang HM, Abecasis GR, Willer CJ, Lee S. 
###   Efficiently controlling for case-control imbalance and sample relatedness in large-scale genetic association studies. 
###   Nature Genetics. 2018
###   doi: 10.1038/s41588-018-0184-y
### Functions for saddle point approximation (SPA) are from the SPAtest package (GNU GPL): 
###   https://github.com/leeshawn/SPAtest

#' calculate_grm
#'
#' calculate the genetic relatedness matrix from the gene matrix
#'
#' @param gene_matrix must have first column as sample names and other columns as gene names
#' @return grm genetic relatedness matrix for tau and beta test
#' @export
calculate_grm <- function(gene_matrix) {
  check_gene_matrix(gene_matrix)
  freq_mat_dist_man = parDist(as.matrix(gene_matrix[, -1]), method = "manhattan") / (ncol(gene_matrix[, -1]) - 1)
  freq_mat_dist_man = as.matrix(freq_mat_dist_man)
  freq_mat_grm_man = 1 - freq_mat_dist_man
  grm = as.matrix(freq_mat_grm_man)
  colnames(grm) = as.matrix(gene_matrix[, 1])
  rownames(grm) = colnames(grm)
  return(grm)
}

check_gene_matrix <- function(gene_matrix) {
 if(!colnames(gene_matrix)[1] == "sample_name"){
   warning("first column in gene matrix is not sample_name this might cause unexpected behavior")
 }
  if((ncol(gene_matrix)) < nrow(gene_matrix)){
    stop("gene matrix is not the expected dimensions, gene matrix should be a gene X sample matrix with more genes than samples")
  }
  if(!is.data.frame(gene_matrix)){
    warning("gene_matrix is not a data frame, could cause unexpected behavior")
  }
}


gen_sigma <- function(W, var_vec, grm) {
  ### update grm with W and var_vec
  # grm is an (nxn) symmetric matrix
  # adapted from SAIGE package
  grm = as.matrix(grm)
  dtkin = W^-1 * (var_vec[1]) 
  new_grm = grm * var_vec[2]
  diag_new_grm = diag(new_grm) + dtkin # update diag
  diag(new_grm) = diag_new_grm
  new_grm[new_grm < 1e-4] = 1e-4 #check nothing is 0
  return(as.matrix(new_grm))
}

get_coef_inner <- function(Y, X, W, var_vec, grm) {
  # Y is working vector Y=alpha X + b
  # n number of samples
  # N number of covariates
  # X (nxN) covariates
  # W coefficient of variation
  # var_vec is the variances of phi of the y varaible and tau of the population structure
  # adapted from SAIGE package
  simga = gen_sigma(W, var_vec, grm) # V
  Y = as.vector(Y)
  sigmai = Rfast::spdinv(simga) 
  if(all(is.na(sigmai))){
    print("While fitting Sigma, a transformation of the GRM became singular or was not symmetrical, this can happen when many columns of the GRM are highly correlated to the y, or when there is not enough varaiblity in the GRM.")
    stop()
  }
  sigmai_Y = sigmai %*% Y # V^-1 Y
  sigmai_X = sigmai %*% X # V^-1 X
  cov_var = Matrix::solve(forceSymmetric(t(X) %*% sigmai_X), sparse = TRUE, tol = 1e-10) # (Xt V^-1 X)^-1
  sigmai_Xt = t(sigmai_X) # t(V^-1 X)
  sigmaiXtY = sigmai_Xt %*% Y # XtV^-1Y
  alpha = cov_var %*% sigmaiXtY # (Xt V X)^-1 XtVY
  epsilon = var_vec[1] *(t(sigmai_Y) - t(sigmai_X %*% alpha)) / as.vector(W) # phi to act on W
  eta = as.vector(Y - epsilon) # Y-var_vec \sigma (Y-X\alpha)
  b = eta - X %*% alpha
  coef_list = list("sigmai_Y" = sigmai_Y, "sigmai_X" = sigmai_X, "cov_var" = cov_var, "alpha" = alpha, "eta" = eta, "b" = b, "epsilon" = epsilon)
  return(coef_list)
}

get_AI_score <- function(Y, X, grm, W, var_vec, sigmai_Y, sigmai_X, cov_var) {
  ## get score for finding var_vec function from supplement of Saige paper
  ## Inputs Y, X, grm, W, Tau, sigmai_Y (sigma times inverse Y), sigmai_X (sigma times inverse X), cov_var
  # adapted from SAIGE package
  sigma = gen_sigma(W, var_vec, grm)
  sigmai =  Rfast::spdinv(sigma)
  sigmai_Xt = t(sigmai_X) # transpose sigma inverse times X
  P = sigmai - sigmai_X %*% cov_var %*% sigmai_Xt # solve for P
  PY1 = P %*% Y # \hat{Y}-\hat(X) (Xt V X)^-1 PY
  APY = grm %*% PY1 # grm (\hat{Y}-\hat(X) (Xt V X)^-1)
  YPAPY = t(PY1) %*% APY # dot product
  YPAPY = YPAPY[1] 
  PA = P %*% grm
  trace_P_grm = sum(diag(PA))
  score1 = YPAPY - trace_P_grm
  PAPY = P %*% APY
  AI = (t(PAPY) %*% APY) # AI=t(Y)%*%P%*%grm%*%P%*%grm%*%P%*%Y
  return(list(YPAPY = YPAPY, PY = PY1, trace_P_grm = trace_P_grm, score1 = score1, AI = AI[1]))
}

get_AI_score_quant <- function(Y, X, grm, W, var_vec, sigmai_Y, sigmai_X, cov_var) {
  ## quanitative y verion
  ## get score for finding var_vec function from supplement of Saige paper 
  ## https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0184-y/MediaObjects/41588_2018_184_MOESM1_ESM.pdf
  ## Inputs Y, X, grm, W, Tau, sigmai_Y, sigmai_X, cov_var
  # adapted from SAIGE package
  n = length(W)
  sigma = gen_sigma(W, var_vec, grm)
  sigmai = Rfast::spdinv(sigma)
  sigmai_Xt = t(sigmai_X) # transpose X
  P = sigmai - sigmai_X %*% cov_var %*% sigmai_Xt
  diag_P = diag(P) / (W)
  PY1 = P %*% Y # \hat{Y}-\hat(X) (Xt V X)^-1 PY
  wPY = PY1 / (W)
  YPwPY = t(PY1) %*% wPY
  YPwPY = YPwPY[1]
  APY = grm %*% PY1
  YPAPY = t(PY1) %*% APY # dot product
  YPAPY = YPAPY[1] 
  PA = P %*% grm
  trace_P_grm = (sum(sigmai * grm) - sum(sigmai_X * crossprod(grm, t(cov_var %*% sigmai_Xt))))
  trace_PW = sum(diag_P)
  score1 = YPAPY - trace_P_grm # score 1
  score0 = YPwPY - trace_PW # score 0 
  score_vector = as.matrix(c(score0[1], score1[1]))
  PwPY = P %*% wPY
  PAPY = P %*% APY
  AI_11 = (t(PAPY) %*% APY)
  AI_00 = (t(PwPY) %*% wPY)
  AI_01 = (t(PAPY) %*% wPY)
  AI_mat = matrix(c(AI_00[1], AI_01[1], AI_01[1], AI_11[1]), 2, 2)
  Dtau = Matrix::solve(AI_mat, score_vector)
  
  return(list(YPAPY = YPAPY, PY = PY1, YPwPY = YPwPY, trace_P_grm = trace_P_grm, trace_PW = trace_PW, AI = AI_mat, score_vector = score_vector))
}

for_beta_bin <- function(mu, y, X) {
  ## score test for var_vec test used to fit beta test
  ## inputs:
  # fitted mu
  # original y
  # covariates X
  # adapted from SAIGE package
  mu2 = mu * (1 - mu)
  V = as.vector(mu2)
  res = as.vector(y - mu)
  XV = t(X * V)
  XVX = t(X) %*% (t(XV))
  XVX_inv = Matrix::solve(XVX)
  XXVX_inv = X %*% XVX_inv
  XVX_inv_XV = XXVX_inv * V
  S_a = colSums(X * res)
  for_beta = list(XV = XV, XVX = XVX, XXVX_inv = XXVX_inv, XVX_inv = XVX_inv, S_a = S_a, XVX_inv_XV = XVX_inv_XV, V = V)
  return(for_beta)
}

for_beta_quant <- function(mu, var_vec, y, X) {
  ## score test for var_vec test used to fit beta test
  ## inputs:
  # fitted mu
  # original y
  # covariates X
  # adapted from SAIGE package
  V = rep(1 / var_vec[1], length(y))
  res = as.vector(y - mu)
  XV = t(X * V)
  XVX = t(X) %*% (t(XV))
  XVX_inv = Matrix::solve(XVX)
  XXVX_inv = X %*% XVX_inv
  XVX_inv_XV = XXVX_inv * V
  S_a = colSums(X * res)
  for_beta = list(XV = XV, XVX = XVX, XXVX_inv = XXVX_inv, XVX_inv = XVX_inv, S_a = S_a, XVX_inv_XV = XVX_inv_XV, V = V)
  return(for_beta)
}

fit_vars <- function(Y_vec, X_mat, grm, w_vec, var_vec, sigmai_Y, sigmai_X, cov_var, tol, quant = FALSE, verbose, write_log) {
  ## fitting process for tau and phi
  # adapted from SAIGE package
  if (!quant) {
    re.AI = get_AI_score(Y_vec, X_mat, grm, w_vec, var_vec, sigmai_Y, sigmai_X, cov_var) 
    score1 = re.AI$score1 
    AI1 = re.AI$AI
    Dtau = score1 / AI1
    var_vec0 = var_vec
    var_vec[2] = var_vec0[2] + Dtau # update tau
    step = 1.0
    while (var_vec[2] < 0) {
      step = step * 0.5
      var_vec[2] = var_vec0[2] + step * Dtau
    }
  } else {
    re.AI = get_AI_score_quant(Y_vec, X_mat, grm, w_vec, var_vec, sigmai_Y, sigmai_X, cov_var)
    YPAPY = re.AI$YPAPY
    YPwPY = re.AI$YPwPY
    trace_PW = re.AI$trace_PW
    trace_P_grm = re.AI$trace_P_grm
    Dtau = Matrix::solve(re.AI$AI, re.AI$score_vector)
    var_vec0 = var_vec
    var_vec = var_vec0 + Dtau
    step = 1.0
    while (any(var_vec < 0)) { #if var_vec is less than zero step through
      step = step * 0.5
      var_vec = var_vec0 + step * Dtau
    }
  }

  if (var_vec[1] < tol) {
    print("Warning! The first variance component parameter estimate is set at the tolarance")
    var_vec[1] = tol*10
  }
  if (var_vec[2] < tol) {
    var_vec[2] = 0
  }
  
  return(list("var_vec" = var_vec))
}

get_alpha <- function(y, X, var_vec, grm, family, alpha0, eta0, offset, verbose = FALSE, maxiter, tol.coef = tol, write_log = FALSE) {
  # adapted from SAIGE package
  mu <- family$linkinv(eta0)
  mu_eta <- family$mu.eta(eta0)
  Y <- eta0 - offset + (y - mu) / mu_eta

  sqrt_W <- mu_eta / sqrt(family$variance(mu))

  W <- sqrt_W^2
  for (i in 1:maxiter) {
    alpha.obj <- get_coef_inner(Y, X, W, var_vec, grm)

    alpha <- as.matrix(alpha.obj$alpha)
    eta <- as.matrix(alpha.obj$eta + offset)

    if (verbose) {
      cat("\n Phi and Tau:\n")
      cat(var_vec)
      cat("\n Fixed-effect coefficients:\n")
      cat(alpha)
    }
    if (write_log) {
      put(paste(" Tau:", var_vec, " Fixed-effect coefficients:", alpha), console = FALSE)
    }
    mu <- family$linkinv(eta)
    mu_eta <- family$mu.eta(eta)

    Y <- eta - offset + (y - mu) / mu_eta
    sqrt_W <- mu_eta / sqrt(family$variance(mu))
    W <- sqrt_W^2

    if (max(abs(alpha - alpha0) / (abs(alpha) + abs(alpha0) + tol.coef)) < tol.coef) {
      break
    }
    alpha0 <- alpha
  }

  alpha_result <- list(Y = Y, alpha = alpha, eta = eta, W = W, cov_var = alpha.obj$cov_var, sqrt_W = sqrt_W, sigmai_Y = alpha.obj$sigmai_Y, sigmai_X = alpha.obj$sigmai_X, mu = mu, eta_2 = alpha.obj$eta_2, b = alpha.obj$b)
  return(alpha_result)
}

#' fit_tau_test
#'
#' Fit the base model for population structure finding random effects
#'
#' @param glm.fit0 glm model. Model output with no sample relatedness accounted for
#' @param grm Genetic Relatedness Matrix (from scripts or user) NxN matrix of sample relatedness
#' @param species_id for tracking species
#' @param tau0 inital tau estimate (variance on the population strucutre of genetic relatedness)
#' @param phi0 inital phi estimate (variance on means of outputs will be 1 for logit or binonimal dispersion)
#' @param maxiter maximum iterations to fit the glmm model
#' @param verbose whether outputting messages in the process of model fitting
#' @param log_file log file to write to
#' @return model output for the tau test on population structure 
#' @export
fit_tau_test <- function(glm.fit0, grm, species_id, tau0 = 1, phi0= 1, maxiter = 100, verbose = TRUE, tol = .0001, log_file = NA) {
  # Fits the null generalized linear mixed model for a binary trait
  # adapted from SAIGE package
  # Args:
  #  glm.fit0: glm model. Logistic model output (with no sample relatedness accounted for)
  #  grm: Genetic Relatedness Matrix  in the same sample order as glm.fit0!
  #  species_id: Species ID of the species for record
  #  tau0: initial values for the variance component parameter tau
  #  phi0: initial values for the variance component parameter phi (only for quant)
  #  maxiter: maximum iterations to fit the glmm model
  #  verbose: whether outputting messages in the process of model fitting
  #  tol: tolance for varaince estimates
  #  log_file: whether to write to a log file and what that would be
  # Returns:
  #  model object pop.struct.glmm
  write_log = !is.na(log_file) & log_file != FALSE
  t_begin = proc.time()
  if (verbose) {
    cat("begining time ")
    cat(t_begin)
  }
  if (write_log) {
    put("begining time ", console = FALSE)
    put(t_begin, console = FALSE)
  }
  check_inputs_tau(glm.fit0, grm, species_id, tau0, phi0, maxiter, tol, verbose, write_log, log_file)
  y = glm.fit0$y
  n = length(y)
  X = model.matrix(glm.fit0)
  Xorig = model.matrix(glm.fit0)
  offset = glm.fit0$offset
  if (is.null(offset)) {
    offset = rep(0, n)
  }

  family = glm.fit0$family
  eta = glm.fit0$linear.predictors
  mu = glm.fit0$fitted.values
  mu_eta = family$mu.eta(eta)
  Y = eta - offset + (y - mu) / mu_eta
  alpha0 = glm.fit0$coef
  eta0 = eta
  sample_ids = colnames(grm)
  
  ### check for quantitiative or not from glm obj
  if (family$family %in% c("poisson", "binomial")) {
    phi0 = 1
    quant = FALSE
  } else {
    quant = TRUE
  }

  if (verbose) cat(" Fixed-effect coefficients: \n", names(glm.fit0$coefficients), "\n", glm.fit0$coefficients)
  if (write_log) put(paste(" Fixed-effect coefficients: ", glm.fit0$coef), console = FALSE)
  if (verbose) cat(" inital tau is ", tau0, "\n")
  if (write_log) put(paste(" inital tau is ", tau0), console = FALSE)
  
  var_vec0 = c(phi0, tau0) ## var_vec is a vector of phi and tau for the variance estimates for tau test
  if (var_vec0[1] <= 0) {
    stop("\nERROR! The first variance component parameter estimate is 0\n")
  }
  alpha.obj = get_alpha(y, X, var_vec0, grm, family, alpha0, eta0, offset, maxiter = maxiter, verbose = verbose, tol.coef = tol, write_log = write_log)
  var_vec=var_vec0
 
  ### update var differently if quantitative 
  if (quant) {
    re = get_AI_score_quant(alpha.obj$Y, X, grm, alpha.obj$W, var_vec0, alpha.obj$sigmai_Y, alpha.obj$sigmai_X, alpha.obj$cov_var)
    var_vec[2] = max(0, as.numeric(var_vec0[2] + var_vec0[2]^2 * (re$YPAPY - re$trace_P_grm) / n))
    var_vec[1] = max(tol*10, as.numeric(var_vec0[1] + var_vec0[1]^2 * (re$YPwPY - re$trace_PW) / n))
  } else {
    re = get_AI_score(alpha.obj$Y, X, grm, alpha.obj$W, var_vec0, alpha.obj$sigmai_Y, alpha.obj$sigmai_X, alpha.obj$cov_var)
    var_vec[2] = max(0, as.numeric(var_vec0[2] + var_vec0[2]^2 * ((re$YPAPY - re$trace_P_grm)) / n)) # first fit vars
  }

  for (i in seq_len(maxiter)) {
    if (verbose) cat(paste("\ni", i))
    if (verbose) cat("\nIteration ", i, "tau is ", var_vec[2], "\n")
    if (write_log) put(paste(" Iteration ", i, "tau is: ", var_vec[2]), console = FALSE)
    tol_limt = FALSE
    alpha0 = alpha.obj$alpha
    var_vec0 = var_vec
    eta0 = eta
    rss_0 = sum((y - mu)^2)
    t_begin_alpha = proc.time()
    alpha.obj = get_alpha(y, X, var_vec, grm, family, alpha0, eta0, offset, verbose = verbose, maxiter = maxiter, tol.coef = tol, write_log = write_log)
    t_end_get_alpha = proc.time()
    if (verbose) {
      cat("\ntime to get alpha\n")
      cat(t_end_get_alpha - t_begin_alpha)
    }
    if (write_log) {
      put("time to get alpha", console = FALSE)
      put(t_end_get_alpha - t_begin_alpha, console = FALSE)
    }

    ## update tau and phi
    fit.obj = fit_vars(alpha.obj$Y, X, grm, alpha.obj$W, var_vec, alpha.obj$sigmai_Y, alpha.obj$sigmai_X, alpha.obj$cov_var, tol = tol, verbose = verbose, write_log = write_log, quant = quant)
    t_end_fit_tau = proc.time()
    if (verbose) {
      cat("\ntime to fit tau\n")
      cat(t_end_fit_tau - t_end_get_alpha)
    }
    if (write_log) {
      put("time to fit tau", console = FALSE)
      put(t_end_fit_tau - t_end_get_alpha, console = FALSE)
    }

    ## update all params
    tau = as.numeric(fit.obj$var_vec[2])
    var_vec = as.numeric(fit.obj$var_vec)
    cov_var = alpha.obj$cov_var
    alpha = alpha.obj$alpha
    eta = alpha.obj$eta
    Y = alpha.obj$Y
    mu = alpha.obj$mu
    res = y - mu
    if (verbose) {
      cat(paste("\nchange in tau", abs(tau - tau0) / (abs(tau) + abs(tau0) + tol)))
      cat("\ntau: ", tau, "\n")
    }
    if (write_log) {
      put(paste("change in tau", abs(tau - tau0) / (abs(tau) + abs(tau0) + tol)), console = FALSE)
      put(paste("tau: ", tau), console = FALSE)
    }
    #if(sum(res^2)/n < tol) break
    if (var_vec[1] <= 0) {
      stop("\nERROR! The first variance component parameter estimate is 0\n")
    }
    
    if (var_vec[2] == 0) break
    var_condition = max(abs(var_vec - var_vec0) / (abs(var_vec) + abs(var_vec0) + tol)) < tol 
    if (var_vec[1] <= tol*10 & var_vec0[1] <= tol*10) {
      tol_limt=TRUE
      print("Warning! The first variance component parameter estimate is set at the tolarance breaking iterations")
      break
    }
    rss = sum(res^2)
    rss_condition = rss_0 - rss
    if (verbose) {
      cat(paste("\nres", rss))
      cat(paste("\nrss change", rss_condition))
    }
    if (write_log) {
      put(paste("res", rss), console = FALSE)
      put(paste("rss change", rss_condition), console = FALSE)
    }

    if (var_condition) break
    if (max(var_vec) > tol^(-2)) {
      tol_limt = TRUE
      i = maxiter
      break
    }
  }
  if (verbose) cat("\niter break at ", i)
  if (verbose) cat("\nFinal ", var_vec, ":\n")
  if (write_log) put(paste("iter break at ", i), console = FALSE)
  if (write_log) put(paste("Final ", var_vec, ":"), console = FALSE)
  if (max(var_vec) > tol^(-2) | i == maxiter) {
    cat("Model not converged")
    if (write_log) {
      put("Model not converged", console = FALSE)
    }
  }
  alpha.obj = get_alpha(y, X, var_vec, grm, family, alpha, eta, offset, verbose = verbose, maxiter = maxiter, tol.coef = tol, write_log = write_log)
  if (quant) {
    if(tol_limt){
      fit.final = get_AI_score_quant(alpha.obj$Y, X, grm, alpha.obj$W, var_vec, alpha.obj$sigmai_Y, alpha.obj$sigmai_X, alpha.obj$cov_var)
      var_vec[2] = max(0, var_vec0[2] + var_vec0[2]^2 * (fit.final$YPAPY - fit.final$trace_P_grm) / n)
      var_vec[1] = max(tol*10, var_vec0[1] + var_vec0[1]^2 * (fit.final$YPwPY - fit.final$trace_PW) / n)
      i = maxiter
      
    }else{
      fit.final = get_AI_score_quant(alpha.obj$Y, X, grm, alpha.obj$W, var_vec, alpha.obj$sigmai_Y, alpha.obj$sigmai_X, alpha.obj$cov_var)
      var_vec[2] = max(0, var_vec0[2] + var_vec0[2]^2 * (fit.final$YPAPY - fit.final$trace_P_grm) / n)
      var_vec[1] = max(tol*10, var_vec0[1] + var_vec0[1]^2 * (fit.final$YPwPY - fit.final$trace_PW) / n)
    }
    
  } else {
    fit.final = get_AI_score(alpha.obj$Y, X, grm, alpha.obj$W, var_vec, alpha.obj$sigmai_Y, alpha.obj$sigmai_X, alpha.obj$cov_var)
    var_vec[2] = max(0, as.numeric(var_vec0[2] + var_vec0[2]^2 * ((fit.final$YPAPY - fit.final$trace_P_grm)) / n)) # tau + Dtau 
  }
  names(var_vec) = c("phi", "tau")
  cov_var = alpha.obj$cov_var
  alpha = alpha.obj$alpha
  names(alpha) = names(glm.fit0$coefficients)
  eta = alpha.obj$eta
  Y = alpha.obj$Y
  mu = alpha.obj$mu
  mu_eta = family$mu.eta(eta)
  sqrt_W = mu_eta / sqrt(family$variance(mu))
  W = sqrt_W^2
  sigma = gen_sigma(W, var_vec, grm)
  converged = ifelse(i < maxiter, TRUE, FALSE)
  res = y - mu
  rss = sum(res^2)
  if (quant) {
    mu2 = rep(1 / var_vec[1], length(y))
    forbeta.obj = for_beta_quant(mu, var_vec, y, Xorig)
  } else {
    mu2 = mu * (1 - mu)
    forbeta.obj = for_beta_bin(mu, y, Xorig)
  }

  ss = sum((y - mean(y))^2)
  var_fixed = var(Xorig %*% alpha)
  var_random = var(as.vector(alpha.obj$b))
  var_error = var(res)
  # https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210x.2012.00261.x

  Coefficients = paste("Coefficents:", alpha)
  if (!quant) {
    model_metrics = list(
      "S" = sd(res),
      "AUC" = suppressMessages(suppressWarnings(auc(glm.fit0$y, as.vector(mu)))),
      "R-sq" = (1 - rss / ss),
      "R-sq(marginal)" = var_fixed / (var_fixed + var_random + var_error),
      "r-sq(conditional)" = (var_fixed + var_random) / (var_fixed + var_random + var_error)
    )
  } else {
    model_metrics = list(
      "S" = sd(res),
      "R-sq" = (1 - rss / ss),
      "R-sq(marginal)" = var_fixed / (var_fixed + var_random + var_error),
      "r-sq(conditional)" = (var_fixed + var_random) / (var_fixed + var_random + var_error)
    )
  }

  glmm_result = list(
    tau = var_vec[2],
    var_vec = var_vec,
    coefficients = alpha, b = alpha.obj$b, 
    t = sum(alpha.obj$b^2)/length(sample_ids),
    linear_predictors = eta, linear_model = Xorig %*% alpha + alpha.obj$b,
    fitted_values = mu, var_mu = mu2, Y = Y, residuals = res,
    cov_var = cov_var, converged = converged,
    sample_names = sample_ids,
    forbeta.obj = forbeta.obj,
    y = y, X = Xorig,
    trait_type = glm.fit0$family,
    formula = paste0(glm.fit0$formula, " + b"),
    iter_finised = i,
    model_metrics = model_metrics, species_id = species_id, grm = grm
  )
  class(glmm_result) = "pop.struct.glmm"
  t_end = proc.time()
  if (verbose) {
    cat("\nfitting the structure model took\n")
    cat(t_end - t_begin)
    cat("\n")
  }
  if (write_log) {
    put("fitting the structure model took", console = FALSE)
    put(t_end - t_begin, console = FALSE)
  }
  return(glmm_result)
}

check_grm <- function(grm,glm.fit0,verbose){
  if(!is.matrix(grm)){
    stop("grm is expected to be a matrix")
  }
  if(nrow(grm) != ncol(grm)){
    stop("grm is expected to be a NxN matrix, with N as the number of samples")
  }
  if(!all(colnames(grm) == rownames(grm))){
    warning("column names do not match row names, grm should be a sample by sample matrix")
  }
  if(class(glm.fit0)[1] != "glm"){
    stop("Expected a glm object for glm.fit0, Example: glm_fit0=glm(`y ~ age  + 1`, data = exp_metadata, family = `binomial`)")
  }
  if(nrow(glm.fit0$data) != nrow(grm)){
    stop("number of samples from baseline glm does not match number of samples from grm, these must match")
  }
  if ("sample_name" %in% colnames(glm.fit0$data)) {
    if (verbose) {
      cat("\nchecking sample names of grm and glm fit match\n")
    }
    if (!all(glm.fit0$data$sample_name == rownames(grm))) {
      stop("\nERROR! the sample names for glm and grm do not match")
    } else {
      if (verbose) {
        cat("check complete")
      }
    }
  } else {
    warning("not running check on sample names ensure grm sample order and glm fit sample order match! Will only check if data from glm has column named sample_name")
  }
  
}

check_inputs_tau <- function(glm.fit0,grm,species_id,tau0,phi0,maxiter,tol,verbose,write_log,log_file){
  if(!is.logical(verbose)){
    stop("verbose should be a logical")
  }
  if(write_log){
    if(!file_test("-x",log_file)){
      warning("log file is not writing to a file, check path, running without log file")
      write_log = FALSE
    }
  }
  if(!is.numeric(tol)){
    stop("tolarance should be numeric")
  }
  if(tol>.01){
    warning("tolarance is usually smaller than .01")
  }
  if(!is.numeric(maxiter)){
    stop("maxiter should be numeric")
  }
  if(maxiter<3){
    warning("maxiter is usually larger than 3")
  }
  if(maxiter>100){
    warning("large maxiter, might take a long time to converge")
  }
  if(!is.numeric(tau0)){
    stop("tau0 must be numeric")
  }
  if(!is.numeric(phi0)){
    stop("phi0 must be numeric")
  }
  if(tau0 <= 0 | phi0 <= 0){
    stop("tau0 and phi0 must be a postive number")
  }
  check_grm(grm,glm.fit0,verbose)
  
}

check_beta <- function(pop.struct.glmm,
                       glm.fit0, grm,
                       gene_df, SPA = FALSE){
  if(!is.logical(SPA)){
    stop("SPA should be a logical")
  }
  if(!(class(pop.struct.glmm)=="pop.struct.glmm")){
    stop("pop.struct.glmm should be an output of fit_tau_test")
  }
  if(all(pop.struct.glmm$y != glm.fit0$y)){
    warning("y for pop.struct.glmm does not match y for glm.fit0, if this is not intentional double check data.")
  }
  if(!is.data.frame(gene_df)){
    stop("gene_df is expected to be a data frame")
  }
  if(!any("gene_id" == colnames(gene_df))){
    stop("Expected column named gene_id, please use a gene_df with a column named gene_id")
  }
  
  if(!any("sample_name" == colnames(gene_df))){
    stop("Expected column named sample_name, please use a gene_df with a column named sample_name with sample names")
  }
  if(!any("gene_value" == colnames(gene_df))){
    stop("Expected column named gene_value, please use a gene_df with a column named gene_value with gene value")
  }
  sample_genes = unique(gene_df$gene_id)
  one_gene =  gene_df[which(gene_df$gene_id == sample_genes[1]),]
  if(!all(one_gene$sample_name == pop.struct.glmm$sample_names)){
    warning("sample names for gene_df do not match sample names for  pop.struct.glmm, if this is not intentional double check data for order of samples and number of samples.")
  }
  
  if(ncol(gene_df)>3){
    warning("expected stacked data frame of 3 columns, gene_id, sample_name, gene_value, data frame has more than 3 columns, check that data frame is stacked.")
  }
  
}

#' fit_beta
#'
#' Fit a beta for each genes including the population structure model with the random effects
#'
#' @param pop.struct.glmm output of fit_tau_test; GLMM of species with grm accounted for
#' @param glm.fit0 glm model. Model output with no sample relatedness accounted for
#' @param grm Genetic Relatedness Matrix (from calculate_grm or user) NxN matrix of sample relatedness
#' @param gene_df long data frame with, gene_id, sample_name, and gene_value
#' @param SPA whether to run Saddle point approximation for pvalues (will slow down output)
#' @return dataframe of beta estimates for all genes tested
#' @export
fit_beta <- function(pop.struct.glmm,
                     glm.fit0, grm,
                     gene_df, SPA = FALSE) {
  # Args:
  # pop.struct.glmm: output from fit_tau_test
  # glm.fit0: output from glm model for data in same order
  # grm: genetic relatedness matrix
  # gene_df: the presance or absence for each gene in a long format 
  # gene_df must have column for sample; column for gene_id; column for gene_value
  # SPA: whether to run saddle point approximation for pvalues (mostly for binomial case)
  # Returns:
  # list of values for each gene examined
  # adapted from SAIGE package
  t_begin = proc.time()
  list_vec=NULL
  check_grm(grm,glm.fit0,verbose = TRUE)
  check_beta(pop.struct.glmm, glm.fit0, grm, gene_df, SPA)
  forbeta.obj = pop.struct.glmm$forbeta.obj
  family = glm.fit0$family

  eta = pop.struct.glmm$linear_predictors
  mu = pop.struct.glmm$fitted_values
  mu_eta = family$mu.eta(eta)
  sqrt_W = mu_eta / sqrt(glm.fit0$family$variance(mu))
  W1 = sqrt_W^2 
  sample_lookup = data.frame(sample_names = pop.struct.glmm$sample_names, index = seq(1, length(pop.struct.glmm$sample_names)), y = pop.struct.glmm$y)
  sample_genes = unique(gene_df$gene_id)

  for (k in sample_genes) {
    iter = which(sample_genes == k)
    if (iter %% 1000 == 0) {
      cat(paste("number of genes done ", iter, "\n"))
      cat("time past:")
      t_now = proc.time()
      cat(t_now - t_begin)
      cat("\n")
    }
    one_gene = gene_df[which(gene_df$gene_id == k),] %>%
      ungroup() 
    rownames(one_gene) = one_gene$sample_name
    one_gene = one_gene[sample_lookup$sample_names, ]
    G0 = as.vector(one_gene$gene_value)
    G_tilde = G0 - pop.struct.glmm$forbeta.obj$XXVX_inv %*% (pop.struct.glmm$forbeta.obj$XV %*% G0) # G1 is X adjusted
    Y = eta + (pop.struct.glmm$y - mu) / mu_eta
    t_score = t(G_tilde) %*% (pop.struct.glmm$y - mu)
    m1 = sum(mu * G_tilde)
    var1 = sum(W1 * G_tilde^2)
    t_adj_2 = (t_score^2) / var1
    beta = t_score / var1
    pval = (pchisq(t_adj_2, lower.tail = FALSE, df = 1, log.p = FALSE))
    z = (qnorm(pval / 2, log.p = F, lower.tail = F))
    se_beta = abs(beta) / sqrt(abs(z))
    qtilde = t_score + m1
    if (SPA) {
      if (var1 < 0) {
        list_vec = rbind(list_vec, data.frame(
          "species_id" = pop.struct.glmm$species_id, tau = as.numeric(pop.struct.glmm$var_vec[2]), "gene_id" = k, "cor" = cor(G0, pop.struct.glmm$y), "cor_to_b" = cor(as.numeric(pop.struct.glmm$b), G0), "z" = z, "var1" = var1, "beta" = NA, "se_beta" = NA, "pvalue" = NA,
          "t_adj" = NA,
          SPA_pvalue = NA, spa_score = NA, SPA_zvalue = NA, pvalue_noadj = NA, converged = FALSE
        ))
      } else {
        out1 = saddle_prob(q = qtilde, mu = mu, g = G_tilde, var1, cutoff = 2, log.p = FALSE)
        list_vec = rbind(list_vec, data.frame(
          "species_id" = pop.struct.glmm$species_id, tau = as.numeric(pop.struct.glmm$var_vec[2]), "gene_id" = k, "cor" = cor(G0, pop.struct.glmm$y), "cor_to_b" = cor(as.numeric(pop.struct.glmm$b), G0), "z" = z, "var1" = var1, "beta" = beta, "se_beta" = se_beta,
          "t_adj" = t_adj_2,
          SPA_pvalue = out1$p.value, spa_score = out1$score, SPA_zvalue = out1$z_value, pvalue_noadj = out1$pvalue_noadj, converged = out1$converged
        ))
      }
    } else {
      # no SPA
      list_vec = rbind(list_vec, data.frame(
        "species_id" = pop.struct.glmm$species_id, tau = as.numeric(pop.struct.glmm$var_vec[2]), "gene_id" = k, "cor" = cor(G0, pop.struct.glmm$y), "cor_to_b" = cor(as.numeric(pop.struct.glmm$b), G0), "z" = z, "var1" = var1, "beta" = beta, "se_beta" = se_beta, "pvalue" = pval,
        "t_adj" = t_adj_2
      ))
    }
  }
  cat("total time past:")
  t_end = proc.time()
  cat(t_end - t_begin)
  cat("\n")
  return(list_vec)
}

filter_pop_obj <- function(pop.struct.glmm, sample_indexs) {
  pop.struct.glmm$gene_value = sample_indexs$gene_value
  pop.struct.glmm$residuals = pop.struct.glmm$residuals[sample_indexs$index]
  pop.struct.glmm$b = pop.struct.glmm$b[sample_indexs$index]
  pop.struct.glmm$linear.predictors = pop.struct.glmm$linear.predictors[sample_indexs$index]
  pop.struct.glmm$fitted.values = pop.struct.glmm$fitted.values[sample_indexs$index]
  pop.struct.glmm$forbeta.obj$XXVX_inv = pop.struct.glmm$forbeta.obj$XXVX_inv[sample_indexs$index, ]
  pop.struct.glmm$forbeta.obj$XV = pop.struct.glmm$forbeta.obj$XV[, sample_indexs$index]
  pop.struct.glmm$sigma = pop.struct.glmm$sigma[sample_indexs$index, sample_indexs$index]
  pop.struct.glmm$X = pop.struct.glmm$X[sample_indexs$index, ]
  pop.struct.glmm$y = pop.struct.glmm$y[sample_indexs$index]
  pop.struct.glmm$sigmai_X = pop.struct.glmm$sigmai_Xt[sample_indexs$index, ]
  pop.struct.glmm$sigmai_Y = pop.struct.glmm$sigmai_Y[sample_indexs$index]
  return(pop.struct.glmm)
}

simulate_tau_inner <- function(glm.fit0, grm, species_id = "s_id", tau0, phi0) {
  family_to_fit = glm.fit0$family
  data_new = glm.fit0$data
  formulate_to_fit = glm.fit0$formula
  data_new_shuffled = data_new[sample(1:nrow(data_new), nrow(data_new)), ]
  data_new_shuffled$sample_name = rownames(grm)
  refit0 = glm(formulate_to_fit, data = data_new_shuffled, family = family_to_fit)
  fit.glmm = tryCatch(fit_tau_test(refit0, grm, tau0 = tau0, phi0 = phi0, verbose = FALSE, species_id = species_id, log_file = NA), error = function(e) e)
  if (length(fit.glmm$t)>0) {
    t = sum(fit.glmm$b^2, na.rm = TRUE)/length(fit.glmm$sample_names)
    tau = fit.glmm$var_vec[2]
    phi  = fit.glmm$var_vec[1]
  } else {
    print("error")
    tau = 0
    t = 0
    phi = 0
  }
  return(data.frame("tau" = tau, t,phi))
}

#' run_tau_test
#'
#' take output from population structure test and test how significant the tau is for that set of data
#'
#' @param glm.fit0 glm model. Model output with no sample relatedness accounted for
#' @param grm Genetic Relatedness Matrix (from scripts or user) NxN matrix of sample relatedness
#' @param n_tau number of tau to simulate
#' @param species_id species id for bacterial species
#' @param tau0 starting tau
#' @param phi0 starting phi
#' @return df of values of T for tau for different permutations of the covarites matrix
#' @export
run_tau_test <- function(glm.fit0, grm, n_tau, species_id = "s_id", tau0, phi0, seed=1) {
  set.seed(seed)
  list_of_tau = lapply(seq(1, n_tau), function(x) simulate_tau_inner(glm.fit0, grm, species_id = species_id, tau0, phi0))
  df_of_tau = do.call(rbind, list_of_tau)
  return(df_of_tau)
}




#' summary.pop.struct.glmm 
#'
#' Summarize the output of fit tau test
#'
#' @param x a pop.struct.glmm objecct the output of fit_tau_test; GLMM of species with grm accounted for
#' @export summary.pop.struct.glmm
#' @export 
summary.pop.struct.glmm <- function(x, ...) {
  cat("Species ID: ", x$species_id, "\n")
  cat("Formula: ", x$formula, "\n")
  cat("family: ", paste(x$trait_type[1:2]), "\n")
  cat("Fixed-effect covariates estimates: \n", names(x$coefficients), "\n", round(x$coefficients,3), "\n")
  cat("Converged: ", x$converged, "\n")
  cat("Number of iterations:",x$iter_finised,"\n")
  cat("Tau: ", round(x$var_vec[2],3), "\n")
  cat("Phi: ", round(x$var_vec[1],3), "if logit or binomail should be 1", "\n")
  cat("T value of tau:", round(x$t,3),"\n")
  cat("Number of Samples:", length(x$sample_names),"\n")
}

#' ... all the usual documentation for summary() ...
#' @export
summary <- function(x, ...) {
  UseMethod("summary")
}


saddle_prob <- function(q, mu, g, var1, cutoff = 2, log.p = FALSE) {
  #### taken from ‘SPAtest’ with a few changes for use case
  m1 = sum(mu * g)
  var1 = sum(mu * (1 - mu) * g^2)
  p1 = NULL
  p2 = NULL

  score = q - m1

  qinv = -sign(q - m1) * abs(q - m1) + m1
  t_adj_sq = ((q - m1)^2) / var1
  pval_noadj = pchisq(t_adj_sq, lower.tail = FALSE, df = 1, log.p = FALSE)
  converged = TRUE

  if (is.na(abs(q - m1)) || is.na(var1) || var1 < 0) {
    converged = FALSE
    pval = pval_noadj
  }
  if (isTRUE(abs(q - m1) / sqrt(var1) < cutoff)) {
    pval = pval_noadj
  } else {
    out_uni1 = get_root_K1(0, mu = mu, g = g, q = q)
    out_uni2 = get_root_K1(0, mu = mu, g = g, q = qinv)
    if (out_uni1$converged == TRUE && out_uni2$converged == TRUE) {
      p1 = tryCatch(get_saddle_prob(out_uni1$root, mu, g, q, log.p = log.p), error = function(e) {
        if (log.p) {
          return(pval_noadj - log(2))
        } else {
          return(pval_noadj / 2)
        }
      })
      p2 = tryCatch(get_saddle_prob(out_uni2$root, mu, g, qinv, log.p = log.p), error = function(e) {
        if (log.p) {
          return(pval_noadj - log(2))
        } else {
          return(pval_noadj / 2)
        }
      })
      if (log.p) {
        pval = add_logp(p1, p2)
      } else {
        pval = abs(p1) + abs(p2)
      }
      converged = TRUE
    } else {
      cat("Error_Converge")
      pval = pval_noadj
      converged = FALSE
    }
  }
  z_value = qnorm(pval / 2, log.p = F, lower.tail = F) * sign(q - m1)
  if (pval != 0 && pval_noadj / pval > 10^3) {
    return(saddle_prob(q, mu, g, var1, cutoff = cutoff * 2, log.p = log.p))
  } else if (pval == 0) {
    return(list(p.value = pval, pvalue_noadj = pval_noadj, z_value = z_value, converged = FALSE, score = score))
  } else {
    return(list(p.value = pval, pvalue_noadj = pval_noadj, z_value = z_value, converged = converged, score = score))
  }
}

get_root_K1 <- function(init, mu, g, q, m1, tol = .0001, maxiter = 1000) {
  #### taken from ‘SPAtest package’
  g_pos = sum(g[which(g > 0)])
  g_neg = sum(g[which(g < 0)])
  if (q >= g_pos || q <= g_neg) {
    return(list(root = Inf, n_iter = 0, converged = TRUE))
  } else {
    t = init
    K1_eval = K1_adj(t, mu, g, q)
    prevJump = Inf
    rep = 1
    repeat    {
      K2_eval = K2(t, mu, g)
      tnew = t - K1_eval / K2_eval
      if (is.na(tnew)) {
        conv = FALSE
        break
      }
      if (abs(tnew - t) < tol) {
        conv = TRUE
        break
      }
      if (rep == maxiter) {
        conv = FALSE
        break
      }

      newK1 = K1_adj(tnew, mu, g, q)
      if (sign(K1_eval) != sign(newK1)) {
        if (abs(tnew - t) > prevJump - tol) {
          tnew = t + sign(newK1 - K1_eval) * prevJump / 2
          newK1 = K1_adj(tnew, mu, g, q)
          prevJump = prevJump / 2
        } else {
          prevJump = abs(tnew - t)
        }
      }

      rep = rep + 1
      t = tnew
      K1_eval = newK1
    }
    return(list(root = t, n_iter = rep, converged = conv))
  }
}

Korg <- function(t, mu, g) {
  #### From ‘SPAtest’ package
  n.t = length(t)
  out = rep(0, n.t)

  for (i in 1:n.t) {
    t1 = t[i]
    temp = log(1 - mu + mu * exp(g * t1))
    out[i] = sum(temp)
  }
  return(out)
}

get_saddle_prob <- function(zeta, mu, g, q, log.p = FALSE) {
  #### From ‘SPAtest' package
  k1 = Korg(zeta, mu, g)
  k2 = K2(zeta, mu, g)

  if (is.finite(k1) && is.finite(k2)) {
    temp1 = zeta * q - k1


    w = sign(zeta) * (2 * temp1)^{
      1 / 2
    }
    v = zeta * (k2)^{
      1 / 2
    }

    Z.test = w + 1 / w * log(v / w)


    if (Z.test > 0) {
      pval = pnorm(Z.test, lower.tail = FALSE, log.p = log.p)
    } else {
      pval = -pnorm(Z.test, lower.tail = TRUE, log.p = log.p)
    }
  } else {
    if (log.p) {
      pval = -Inf
    } else {
      pval = 0
    }
  }

  return(pval)
}

K1_adj <- function(t, mu, g, q) {
  #### From ‘SPAtest' package
  n.t = length(t)
  out = rep(0, n.t)

  for (i in 1:n.t) {
    t1 = t[i]
    temp1 = (1 - mu) * exp(-g * t1) + mu
    temp2 = mu * g
    out[i] = sum(temp2 / temp1) - q
  }
  return(out)
}

K2 <- function(t, mu, g) {
  n.t = length(t)
  out = rep(0, n.t)

  for (i in 1:n.t) {
    t1 = t[i]
    temp1 = ((1 - mu) * exp(-g * t1) + mu)^2
    temp2 = (1 - mu) * mu * g^2 * exp(-g * t1)
    out[i] = sum(temp2 / temp1, na.rm = TRUE)
  }
  return(out)
}
