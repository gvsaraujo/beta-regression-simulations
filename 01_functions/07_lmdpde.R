# HERE RETURNS THE VALUES OF THE SCORE FUNCTION
#********************Score-type function*******************#
Score_LMDPDE <- function(theta, y, X, Z, alpha_const, linkmu, linkphi) {
  kk1 <- ncol(X)
  kk2 <- ncol(Z)
  beta <- theta[1:kk1]
  gama <- theta[(kk1+1.0):(kk1+kk2)]
  mean_info <- mean_link(beta = beta, y=y, X=X, linkmu=linkmu)
  mu  <- mean_info$mu
  T_1 <- mean_info$T_1
  precision_info <- precision_link(gama = gama, Z=Z, linkphi=linkphi)
  phi <- precision_info$phi
  T_2 <- precision_info$T_2

  a <- mu*phi
  b <- (1.0 - mu)*phi
  a_alpha <- mu*phi*(1 + alpha_const)
  b_alpha <- (1 - mu)*phi*(1 + alpha_const)
  K <- exp(lgamma(a_alpha) + lgamma(b_alpha) - lgamma(a_alpha + b_alpha) -  (1 + alpha_const)* (lgamma(a) + lgamma(b) - lgamma(a + b)))

  m_phi <- phi
  ystar <- log(y/(1.0 - y))
  ydagger <- log(1.0 - y)
  F_q <- ((y*(1-y))*dbeta(y, mu*phi, (1-mu)*phi))^alpha_const
  mustar <- psigamma(a, 0) - psigamma(b, 0)
  mudagger <-  psigamma(b, 0) - psigamma(a + b, 0)
  mustar_alpha <- psigamma(a_alpha, 0) - psigamma(b_alpha, 0)
  mudagger_alpha <- psigamma(b_alpha, 0) - psigamma(a_alpha + b_alpha, 0)

  E1q <- (1 + alpha_const)*as.vector(crossprod(X, K*m_phi*T_1*(mustar_alpha - mustar)))
  E2q <- (1 + alpha_const)*as.vector(crossprod(Z, K*T_2*(mu*(mustar_alpha - mustar) + mudagger_alpha - mudagger)))

  #****vector type Score*****#
  avScore <- numeric(kk1 + kk2)
  avScore[1:kk1] <- E1q - (1 + alpha_const)*as.vector(crossprod(X, m_phi*F_q*T_1*(ystar - mustar)))
  avScore[(kk1+1.0):(kk1+kk2)] <- E2q - (1 + alpha_const)*as.vector(crossprod(Z, F_q*T_2*(mu*(ystar - mustar) + ydagger - mudagger)))
  return(avScore)
}#ends score-type function

#********************Covariance matrix****************#
V_matrix_LMDPDE <- function(theta, y, X, Z, alpha_const, linkmu, linkphi){  
  kk1 <- ncol(X)
  kk2 <- ncol(Z)
  beta <- theta[1:kk1]
  gama <- theta[(kk1+1.0):(kk1+kk2)]
  
  mean_info <- mean_link(beta = beta, y=y, X=X, linkmu=linkmu)
  mu  <- mean_info$mu
  t_1 <- mean_info$T_1
  precision_info <- precision_link(gama = gama, Z=Z, linkphi=linkphi)
  phi <- precision_info$phi
  t_2 <- precision_info$T_2

  a <- mu*phi
  b <- (1 - mu)*phi

  mustar <- psigamma(a, 0) - psigamma(b, 0)
  mudagger <- psigamma(b, 0) - psigamma(a + b, 0)

  psi1 <- psigamma(a, 1.0)
  psi2 <- psigamma(b, 1.0) 
  psi3 <- psigamma(a + b, 1.0) 
  
  a_alpha <- mu*phi*(1 + alpha_const)			 
  b_alpha <- (1 - mu)*phi*(1 + alpha_const)	
  
  psi1_alpha <- psigamma(a_alpha, 1.0)
  psi2_alpha <- psigamma(b_alpha, 1.0) 
  psi3_alpha <- psigamma(a_alpha + b_alpha, 1.0) 
  
  a2_alpha <- mu*phi*(1 + 2*alpha_const)		
  b2_alpha <- (1 - mu)*phi*(1 + 2*alpha_const)	 
  psi1_2alpha <- psigamma(a2_alpha, 1.0)
  psi2_2alpha <- psigamma(b2_alpha, 1.0) 
  psi3_2alpha <- psigamma(a2_alpha + b2_alpha, 1.0)  
  
  mustar_alpha <- psigamma(a_alpha, 0) - psigamma(b_alpha, 0)
  mustar_2alpha <- psigamma(a2_alpha, 0) - psigamma(b2_alpha, 0)
  
  mudagger_alpha <-  psigamma(b_alpha, 0) - psigamma(a_alpha + b_alpha, 0)
  mudagger_2alpha <-  psigamma(b2_alpha, 0) - psigamma(a2_alpha + b2_alpha, 0)
  
  K <-  exp(lgamma(a_alpha) + lgamma(b_alpha) - lgamma(a_alpha + b_alpha) - (1 + alpha_const)*(lgamma(a) + lgamma(b) - lgamma(a + b)))							  
  K2 <-  exp(lgamma(a2_alpha) + lgamma(b2_alpha) - lgamma(a2_alpha + b2_alpha) - (1 + 2*alpha_const)*(lgamma(a) + lgamma(b) - lgamma(a + b)))   
  
  gama11_alpha <- K*(phi^2.0)*(t_1^2.0)*(psi1_alpha + psi2_alpha + (mustar_alpha - mustar)^2.0)
  gama12_alpha <- K*phi*t_1*t_2*(mu*(psi1_alpha +  psi2_alpha + (mustar_alpha - mustar)^2.0) - psi2_alpha + (mustar_alpha - mustar)*(mudagger_alpha - mudagger))
  gama22_alpha <- K*(t_2^2.0)*((mu^2.0)*psi1_alpha + ((1.0 - mu)^2.0)*psi2_alpha - psi3_alpha + (mu*(mustar_alpha - mustar) + mudagger_alpha - mudagger)^2.0)

  gama11_2alpha <- K2*(phi^2.0)*(t_1^2.0)*(psi1_2alpha + psi2_2alpha + (mustar_2alpha - mustar)^2.0)
  gama12_2alpha <- K2*phi*t_1*t_2*(mu*(psi1_2alpha + psi2_2alpha + (mustar_2alpha - mustar)^2.0) - psi2_2alpha + (mustar_2alpha - mustar)*(mudagger_2alpha - mudagger))
  gama22_2alpha <- K2*(t_2^2.0)*((mu^2.0)*psi1_2alpha + ((1.0 - mu)^2.0)*psi2_2alpha - psi3_2alpha + (mu*(mustar_2alpha - mustar) + mudagger_2alpha - mudagger)^2.0)

  gama1_alpha <- K*phi*t_1*(mustar_alpha - mustar)
  gama2_alpha <- K*t_2*(mu*(mustar_alpha - mustar) + mudagger_alpha - mudagger)

  Psin_betabeta <- (1 + alpha_const)*crossprod(X, X*gama11_alpha)
  Psin_betagamma <- (1 + alpha_const)*crossprod(X, Z*gama12_alpha)
  Psin_gammagamma <- (1 + alpha_const)*crossprod(Z, Z*gama22_alpha)
  Psin <- matrix(numeric(0), kk1+kk2, kk1+kk2)
  Psin[1:kk1,1:kk1] <- Psin_betabeta
  Psin[1:kk1,(kk1+1):(kk1+kk2)] <- Psin_betagamma
  Psin[(kk1+1):(kk1+kk2),1:kk1] <- t(Psin_betagamma)
  Psin[(kk1+1):(kk1+kk2),(kk1+1):(kk1+kk2)] <- Psin_gammagamma

  Omegan_betabeta <- ((1 + alpha_const)^2.0)*crossprod(X, X*(gama11_2alpha - gama1_alpha^2.0))
  Omegan_betagamma <- ((1 + alpha_const)^2.0)*crossprod(X, Z*(gama12_2alpha - gama1_alpha*gama2_alpha))
  Omegan_gammagamma <- ((1 + alpha_const)^2.0)*crossprod(Z, Z*(gama22_2alpha - gama2_alpha^2.0))
  Omegan <- matrix(numeric(0), kk1+kk2, kk1+kk2)
  Omegan[1:kk1,1:kk1] <- Omegan_betabeta
  Omegan[1:kk1,(kk1+1):(kk1+kk2)] <- Omegan_betagamma 
  Omegan[(kk1+1):(kk1+kk2),1:kk1] <- t(Omegan_betagamma)
  Omegan[(kk1+1):(kk1+kk2),(kk1+1):(kk1+kk2)] <- Omegan_gammagamma		
  Vq <- tryCatch(solve(Psin)%*%(Omegan)%*%t(solve(Psin)), error=function(e) {e})  #asymptotic covariance matrix
  if(inherits(Vq, "error")){
    Vq <- MASS::ginv(Psin)%*%(Omegan)%*%t(MASS::ginv(Psin))
  }
  return(list(Vq = Vq, H = Psin, J = Omegan))
}#ends Covariance matrix function

#**********LSMLE function**********#
LMDPDE_BETA <- function(y, X, Z, qoptimal=TRUE, q0=1, m=3, L=0.02, qmin=0.5,
                        spac=0.02, method="BFGS", startV ="CP", linkmu="logit", linkphi="log",
                        beta_fix = NULL, idx_beta = NULL,  
                        gamma_fix = NULL, idx_gamma = NULL){
  
  #**********Objective function**********#
  log_liksurrogate <- function(theta) 
  {
    kk1 <- ncol(X_full) # Number of columns in X_full
    kk2 <- ncol(Z_full) # Number of columns in Z_full
    beta <- numeric(kk1) # Initialize beta vector
    gama <- numeric(kk2) # Initialize gamma vector
    
    idx_beta_free <- setdiff(1:kk1, idx_beta) # Indices of free beta
    idx_gamma_free <- setdiff(1:kk2, idx_gamma) # Indices of free gamma
    
    beta[idx_beta] <- beta_fix # Assign fixed values to beta
    beta[idx_beta_free] <- theta[1:length(idx_beta_free)] # Assign free values to beta
    gama[idx_gamma] <- gamma_fix # Assign fixed values to gamma
    gama[idx_gamma_free] <- theta[(length(idx_beta_free) + 1):length(theta)] # Assign free values to gamma
    
    mu <- mean_link(beta = beta, y=y, X=X_full, linkmu=linkmu)$mu
    phi <- precision_link(gama = gama, Z=Z_full, linkphi=linkphi)$phi
    a <- mu*phi						 
    b <- (1.0 - mu)*phi	
    a_alpha <- mu*phi*(1 + alpha_const)			 
    b_alpha <- (1 - mu)*phi*(1 + alpha_const)	
    K <-  exp(lgamma(a_alpha) + lgamma(b_alpha) - lgamma(a_alpha + b_alpha) -  (1 + alpha_const)* (lgamma(a) + lgamma(b) - lgamma(a + b)))
    log_likS <- (1/n)*sum(K - ((1 + alpha_const)/alpha_const)*((y*(1-y))*dbeta(y, a, b))^alpha_const)
    return(log_likS)
  }
  
  #************Score function**********#
  Score <- function(theta) {
    kk1 <- ncol(X_full) # Number of columns in X_full
    kk2 <- ncol(Z_full) # Number of columns in Z_full
    beta <- numeric(kk1) # Initialize beta vector
    gama <- numeric(kk2) # Initialize gamma vector

    idx_beta_free <- setdiff(1:kk1, idx_beta) # Indices of free beta
    idx_gamma_free <- setdiff(1:kk2, idx_gamma) # Indices of free gamma
    beta[idx_beta] <- beta_fix # Assign fixed values to beta
    beta[idx_beta_free] <- theta[1:length(idx_beta_free)] # Assign free values to beta
    gama[idx_gamma] <- gamma_fix # Assign fixed values to gamma
    gama[idx_gamma_free] <- theta[(length(idx_beta_free) + 1):length(theta)] # Assign free values to gamma

    mean_info <- mean_link(beta = beta, y=y, X=X_full, linkmu=linkmu)
    mu  <- mean_info$mu
    T_1 <- mean_info$T_1
    precision_info <- precision_link(gama = gama, Z=Z_full, linkphi=linkphi)
    phi <- precision_info$phi
    T_2 <- precision_info$T_2

    a <- mu*phi
    b <- (1.0 - mu)*phi
    a_alpha <- mu*phi*(1 + alpha_const)
    b_alpha <- (1 - mu)*phi*(1 + alpha_const)
    K <- exp(lgamma(a_alpha) + lgamma(b_alpha) - lgamma(a_alpha + b_alpha) -  (1 + alpha_const)* (lgamma(a) + lgamma(b) - lgamma(a + b)))

    m_phi <- phi
    ystar <- log(y/(1.0 - y))
    ydagger <- log(1.0 - y)
    F_q <- ((y*(1-y))*dbeta(y, mu*phi, (1-mu)*phi))^alpha_const
    mustar <- psigamma(a, 0) - psigamma(b, 0)
    mudagger <-  psigamma(b, 0) - psigamma(a + b, 0)
    mustar_alpha <- psigamma(a_alpha, 0) - psigamma(b_alpha, 0)
    mudagger_alpha <- psigamma(b_alpha, 0) - psigamma(a_alpha + b_alpha, 0)

    E1q <- (1 + alpha_const)*as.vector(crossprod(X_full, K*m_phi*T_1*(mustar_alpha - mustar)))
    E2q <- (1 + alpha_const)*as.vector(crossprod(Z_full, K*T_2*(mu*(mustar_alpha - mustar) + mudagger_alpha - mudagger)))

    #****vector type Score*****#
    avScore_full <- numeric(kk1 + kk2)
    avScore_full[1:kk1] <- E1q - (1 + alpha_const)*as.vector(crossprod(X_full, m_phi*F_q*T_1*(ystar - mustar)))
    avScore_full[(kk1+1.0):(kk1+kk2)] <- E2q - (1 + alpha_const)*as.vector(crossprod(Z_full, F_q*T_2*(mu*(ystar - mustar) + ydagger - mudagger)))

    score_out <- c(avScore_full[idx_beta_free], avScore_full[kk1 + idx_gamma_free])
    return(score_out)
  }
  
  #****Covariance function V***#
  V_matrix <- function(theta) {  
    kk1 <- ncol(X_full) # Number of columns in X_full
    kk2 <- ncol(Z_full) # Number of columns in Z_full
    beta <- numeric(kk1) # Initialize beta vector
    gama <- numeric(kk2) # Initialize gamma vector
    
    idx_beta_free <- setdiff(1:kk1, idx_beta) # Indices of free beta
    idx_gamma_free <- setdiff(1:kk2, idx_gamma) # Indices of free gamma
    beta[idx_beta] <- beta_fix # Assign fixed values to beta
    beta[idx_beta_free] <- theta[1:length(idx_beta_free)] # Assign free values to beta
    gama[idx_gamma] <- gamma_fix # Assign fixed values to gamma
    gama[idx_gamma_free] <- theta[(length(idx_beta_free) + 1):length(theta)] # Assign free values to gamma
    
    
    mean_info <- mean_link(beta = beta, y=y, X=X_full, linkmu=linkmu)
    mu  <- mean_info$mu
    t_1 <- mean_info$T_1
    precision_info <- precision_link(gama = gama, Z=Z_full, linkphi=linkphi)
    phi <- precision_info$phi
    t_2 <- precision_info$T_2

    a <- mu*phi
    b <- (1 - mu)*phi

    mustar <- psigamma(a, 0) - psigamma(b, 0)
    mudagger <- psigamma(b, 0) - psigamma(a + b, 0)

    psi1 <- psigamma(a, 1.0)
    psi2 <- psigamma(b, 1.0) 
    psi3 <- psigamma(a + b, 1.0) 
    
    a_alpha <- mu*phi*(1 + alpha_const)			 
    b_alpha <- (1 - mu)*phi*(1 + alpha_const)	
    
    psi1_alpha <- psigamma(a_alpha, 1.0)
    psi2_alpha <- psigamma(b_alpha, 1.0) 
    psi3_alpha <- psigamma(a_alpha + b_alpha, 1.0) 
    
    a2_alpha <- mu*phi*(1 + 2*alpha_const)		
    b2_alpha <- (1 - mu)*phi*(1 + 2*alpha_const)	 
    psi1_2alpha <- psigamma(a2_alpha, 1.0)
    psi2_2alpha <- psigamma(b2_alpha, 1.0) 
    psi3_2alpha <- psigamma(a2_alpha + b2_alpha, 1.0)  
    
    mustar_alpha <- psigamma(a_alpha, 0) - psigamma(b_alpha, 0)
    mustar_2alpha <- psigamma(a2_alpha, 0) - psigamma(b2_alpha, 0)
    
    mudagger_alpha <-  psigamma(b_alpha, 0) - psigamma(a_alpha + b_alpha, 0)
    mudagger_2alpha <-  psigamma(b2_alpha, 0) - psigamma(a2_alpha + b2_alpha, 0)
    
    K <-  exp(lgamma(a_alpha) + lgamma(b_alpha) - lgamma(a_alpha + b_alpha) - (1 + alpha_const)*(lgamma(a) + lgamma(b) - lgamma(a + b)))							  
    K2 <-  exp(lgamma(a2_alpha) + lgamma(b2_alpha) - lgamma(a2_alpha + b2_alpha) - (1 + 2*alpha_const)*(lgamma(a) + lgamma(b) - lgamma(a + b)))   
    
    gama11_alpha <- K*(phi^2.0)*(t_1^2.0)*(psi1_alpha + psi2_alpha + (mustar_alpha - mustar)^2.0)
    gama12_alpha <- K*phi*t_1*t_2*(mu*(psi1_alpha +  psi2_alpha + (mustar_alpha - mustar)^2.0) - psi2_alpha + (mustar_alpha - mustar)*(mudagger_alpha - mudagger))
    gama22_alpha <- K*(t_2^2.0)*((mu^2.0)*psi1_alpha + ((1.0 - mu)^2.0)*psi2_alpha - psi3_alpha + (mu*(mustar_alpha - mustar) + mudagger_alpha - mudagger)^2.0)

    gama11_2alpha <- K2*(phi^2.0)*(t_1^2.0)*(psi1_2alpha + psi2_2alpha + (mustar_2alpha - mustar)^2.0)
    gama12_2alpha <- K2*phi*t_1*t_2*(mu*(psi1_2alpha + psi2_2alpha + (mustar_2alpha - mustar)^2.0) - psi2_2alpha + (mustar_2alpha - mustar)*(mudagger_2alpha - mudagger))
    gama22_2alpha <- K2*(t_2^2.0)*((mu^2.0)*psi1_2alpha + ((1.0 - mu)^2.0)*psi2_2alpha - psi3_2alpha + (mu*(mustar_2alpha - mustar) + mudagger_2alpha - mudagger)^2.0)

    gama1_alpha <- K*phi*t_1*(mustar_alpha - mustar)
    gama2_alpha <- K*t_2*(mu*(mustar_alpha - mustar) + mudagger_alpha - mudagger)

    Psin_betabeta <- (1 + alpha_const)*crossprod(X_full, X_full*gama11_alpha)
    Psin_betagamma <- (1 + alpha_const)*crossprod(X_full, Z_full*gama12_alpha)
    Psin_gammagamma <- (1 + alpha_const)*crossprod(Z_full, Z_full*gama22_alpha)
    Psin <- matrix(numeric(0), kk1+kk2, kk1+kk2)
    Psin[1:kk1,1:kk1] <- Psin_betabeta
    Psin[1:kk1,(kk1+1):(kk1+kk2)] <- Psin_betagamma
    Psin[(kk1+1):(kk1+kk2),1:kk1] <- t(Psin_betagamma)
    Psin[(kk1+1):(kk1+kk2),(kk1+1):(kk1+kk2)] <- Psin_gammagamma

    Omegan_betabeta <- ((1 + alpha_const)^2.0)*crossprod(X_full, X_full*(gama11_2alpha - gama1_alpha^2.0))
    Omegan_betagamma <- ((1 + alpha_const)^2.0)*crossprod(X_full, Z_full*(gama12_2alpha - gama1_alpha*gama2_alpha))
    Omegan_gammagamma <- ((1 + alpha_const)^2.0)*crossprod(Z_full, Z*(gama22_2alpha - gama2_alpha^2.0))
    Omegan <- matrix(numeric(0), kk1+kk2, kk1+kk2)
    Omegan[1:kk1,1:kk1] <- Omegan_betabeta
    Omegan[1:kk1,(kk1+1):(kk1+kk2)] <- Omegan_betagamma 
    Omegan[(kk1+1):(kk1+kk2),1:kk1] <- t(Omegan_betagamma)
    Omegan[(kk1+1):(kk1+kk2),(kk1+1):(kk1+kk2)] <- Omegan_gammagamma		
    Vq <- tryCatch(solve(Psin)%*%(Omegan)%*%t(solve(Psin)), error=function(e) {e})  #asymptotic covariance matrix
    if(inherits(Vq, "error")){
      Vq <- MASS::ginv(Psin)%*%(Omegan)%*%t(MASS::ginv(Psin))
    }
    idx_all <- c(idx_beta_free, kk1 + idx_gamma_free)
    Vq_restricted <- Vq[idx_all, idx_all, drop = FALSE]
    
    return(Vq_restricted)
  }
  
  #----------------------------------------------------------------------------------------------------------------#
  #------------------------------------->Estimation process<-----------------------------------------#
  #----------------------------------------------------------------------------------------------------------------#
  
  #***Evaluating MLE and starting values for SMLE***#
  kk1 <- ncol(X); kk2 <- ncol(Z); n <- nrow(X) # Important quantities
  
  # Free parameters
  idx_beta_free  <- setdiff(1:kk1, idx_beta) # Indices of free beta
  idx_gamma_free <- setdiff(1:kk2, idx_gamma) # Indices of free gamma
  idx_all_free   <- c(idx_beta_free, kk1 + idx_gamma_free)
  
  theta_hat_full <- numeric(kk1 + kk2) # Initialize full theta vector
  X_free <- X[, idx_beta_free, drop = FALSE] # X matrix with free beta columns
  Z_free <- Z[, idx_gamma_free, drop = FALSE] # Z matrix with free gamma columns
  X_full <- X; Z_full <- Z # Full X and Z matrices
  X <- X_free; Z <- Z_free # Matrices used in estimation process
  
  # Condition when all parameters are fixed
  if(ncol(X_full) == length(idx_beta) && ncol(Z_full) == length(idx_gamma)) {
    theta_hat_full <- c(beta_fix, gamma_fix)
    alpha_const <- 1 - q0 # Calculate q_alpha
    score <- Score_LMDPDE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const = alpha_const, linkmu=linkmu, linkphi=linkphi)
    Vq <- V_matrix_LMDPDE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const=alpha_const, linkmu=linkmu, linkphi=linkphi)
    
    final_results <- list(theta_hat = theta_hat_full, 
                          score = score, Vq = Vq, alpha_const = alpha_const)
    return(final_results)
  }
  # Function to complete theta
  fill_theta <- function() {
    kk1 <- ncol(X_full); kk2 <- ncol(Z_full) # Recalculate kk1 and kk2
    idx_beta_free  <- setdiff(1:kk1, idx_beta)
    idx_gamma_free <- setdiff(1:kk2, idx_gamma)
    
    theta_hat_full <- numeric(kk1 + kk2)
    
    theta_hat_full[idx_beta] <- beta_fix
    theta_hat_full[idx_beta_free] <- theta_hat[1:length(idx_beta_free)]
    
    theta_hat_full[kk1 + idx_gamma] <- gamma_fix
    theta_hat_full[kk1 + idx_gamma_free] <- theta_hat[(length(idx_beta_free) + 1):length(theta_hat)]
    return(theta_hat_full)
  }
  
  # Estimates via MLE
  SV <- MLE_BETA(y, X_full, Z_full,
                 beta_fix = beta_fix, idx_beta = idx_beta,
                 gamma_fix = gamma_fix, idx_gamma = idx_gamma
  )
  
  thetaStart <- SV$theta_hat # Initial parameters for optimization
  theta_mle <- SV$theta_hat # Parameters estimated via MLE
  theta_mle_full <- SV$theta_hat_full # Completing the theta vector 
  
  #----------------------REMOVE----------------------------# 
  vcov <- solve(SV$fisher_info)
  vcov <- vcov[idx_all_free, idx_all_free, drop = FALSE]
  se_mle <- suppressWarnings(sqrt(diag(vcov)))
  z_mle <- as.matrix(theta_mle / se_mle)
  #print(se_mle)
  #----------------------REMOVE----------------------------#
  
  #------------------------------->For fixed q<----------------------------------#
  if(qoptimal==FALSE){# Starting with qoptimal = FALSE
    
    # If q0 is equal to 1, we have MLE
    if(q0==1){
      alpha_const <- 0 # Fixed value of alpha
      theta_hat <- thetaStart # Initial parameters 
      theta_hat_full <- fill_theta() # Completing the theta vector
      
      # Important quantities
      score <- Score_LMDPDE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const =alpha_const, linkmu=linkmu, linkphi=linkphi)
      Vq <- V_matrix_LMDPDE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const=alpha_const, linkmu=linkmu, linkphi=linkphi)
      
      # Final result
      final_results <- list(theta_hat = theta_hat_full, 
                            score = score, 
                            Vq = Vq, 
                            alpha_const = alpha_const)
      
      return(final_results)
    }
    else{ # If q0 is less than 1
      q_const <- q0 # Value of q0
      alpha_const <- 1 - q_const
      
      #-------->Maximizing the objective function<------------#
      res <- suppressWarnings(tryCatch(optim(thetaStart, log_liksurrogate, hessian = F, control=list(fnscale=1,maxit=200), 
                                             gr = Score, method=method), error=function(e) {e}))#maximization
      
      # If optimization fails, return MLE results
      if(inherits(res, "error")){
        alpha_const <- 0
        theta_hat <- SV$theta_hat
        theta_hat_full <- fill_theta() # Completing the theta vector
        score <- Score_LMDPDE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const =alpha_const, linkmu=linkmu, linkphi=linkphi)
        Vq <- V_matrix_LMDPDE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const=alpha_const, linkmu=linkmu, linkphi=linkphi)
        
        final_results <- list(theta_hat = theta_hat_full, 
                              score = score, Vq = Vq, alpha_const = alpha_const)
        return(final_results)
      }
      else{ # If optimization converged
        if(res$convergence == 0||res$convergence == 1){ 
          alpha_const <- 1 - q_const
          theta_hat <- res$par # Estimated parameters
          theta_hat_full <- fill_theta() # Completing the theta vector
          
          score <- Score_LMDPDE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const =alpha_const, linkmu=linkmu, linkphi=linkphi)
          Vq <- V_matrix_LMDPDE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const=alpha_const, linkmu=linkmu, linkphi=linkphi)
          
          final_results <- list(theta_hat = theta_hat_full, 
                                score = score, Vq = Vq, 
                                alpha_const = alpha_const)
          return(final_results)
        }
        else { # If optimization did not converge, return MLE results
          
          alpha_const <- 0
          theta_hat <- SV$theta_hat
          theta_hat_full <- fill_theta() # Completing the theta vector
          
          # Important quantities
          score_full <- Score_LMDPDE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const =alpha_const, linkmu=linkmu, linkphi=linkphi)
          Vq_full <- V_matrix_LMDPDE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const=alpha_const, linkmu=linkmu, linkphi=linkphi)
          
          # Final result
          final_results <- list(theta_hat = theta_hat_full, 
                                score = score_full, Vq = Vq_full, 
                                alpha_const = alpha_const)
          
          return(final_results)
        }
      }
    }
    
    #*** Finding the optimal value for q***
  }else {# Condition for qoptimal = TRUE
    
    # Important to modify the values of kk1 and kk2
    kk1 <- ncol(X); kk2 <- ncol(Z); n <- nrow(X) # Important quantities
    
    #***************Searching for an optimal value of q****************#
    q_gridI <- sort(seq(from = 0.8, to = 1, by = spac), decreasing = TRUE)
    Size_gridI <- length(q_gridI)
    thetahat_I <- matrix(numeric(0), nrow = length(theta_mle), ncol = Size_gridI)
    se_smleI <- matrix(numeric(0), nrow = length(theta_mle), ncol = Size_gridI)
    trace_smleI <- matrix(numeric(0), nrow = 1, ncol = Size_gridI)
    SQVI <- numeric(Size_gridI - 1)
    thetahat_I[, 1] <- theta_mle
    se_smleI[, 1] <- se_mle
    
    for(j in 1:(Size_gridI-1)){
      q_const <- q_gridI[j+1]
      alpha_const <- 1-q_const
      resT <- suppressWarnings(tryCatch(optim(thetaStart, log_liksurrogate, hessian = F, control=list(fnscale=1, maxit=200), 
                                              gr = Score, method=method), error=function(e) {e}))
      if(inherits(resT, "error")){ # Error case
        thetahat_I[,j + 1.0] <- NA   
        se_smleI[,j + 1.0] <- NA
        trace_smleI[,j+1.0] <- NA
        SQVI[j] <- NA
      }
      else{# Success case
        
        estimatesI <- resT$par
        thetahat_I[,j + 1.0] <- estimatesI
        VqI <- V_matrix(thetahat_I[,j + 1.0])	
        se_smleI[,j + 1.0] =  suppressWarnings(t(sqrt(diag(VqI))))
        trace_smleI[,j+1.0] <- sum(diag(VqI)) 	  
        
        #**** Calculating SQV ****#
        SQVI[j] <- round((1.0/(kk1+kk2))*sqrt(sum( (thetahat_I[,j]/(sqrt(n)*se_smleI[,j]) - 
                                                      thetahat_I[,j + 1.0]/(sqrt(n)*se_smleI[,j + 1.0]))^2)),5); 
      }                                                     
    }  
    
    qvaluesout <- q_gridI[-length(q_gridI)][SQVI>=L]
    qstart <- suppressWarnings(min(qvaluesout,na.rm=T))  
    seq_pQS <- 1:Size_gridI
    position_qstart <- seq_pQS[q_gridI==qstart]
    
    if(qstart==Inf){# Instability in estimates q in [0.8, 1]
      q_optimal <- 1.0 # Return maximum likelihood
      alpha_const <- 0
      
      theta_hat <- theta_mle_full # Parameters estimated via MLE
      
      score <- Score_LMDPDE(theta_hat, y=y, X=X_full, Z=Z_full, alpha_const =alpha_const, linkmu=linkmu, linkphi=linkphi)
      Vq <- V_matrix_LMDPDE(theta_hat, y=y, X=X_full, Z=Z_full, alpha_const=alpha_const, linkmu=linkmu, linkphi=linkphi)
      
      final_results <- list(theta_hat = theta_hat, 
                            score = score, Vq = Vq, 
                            alpha_const = alpha_const)
      
    }else{ # If there is SQV>=L in [0.8, 1]
      
      #***Starting with the new value of q***#
      q_values <- sort(seq(from = qmin, to = qstart, by = spac), decreasing = TRUE)
      nq <- length(q_values)
      SQV <- numeric(nq - 1)
      thetahat_n <- matrix(numeric(0),nrow=kk1+kk2,ncol= nq)
      se_smle <- matrix(numeric(0),nrow=kk1+kk2,ncol= nq)
      trace_smle <- matrix(numeric(0),nrow=1,ncol= nq)	
      thetahat_n[,1] =  thetahat_I[,position_qstart]
      se_smle[,1] =  se_smleI[,position_qstart]
      trace_smle[,1] =  sum(se_smle[,1])	 
      counter <- 1
      grid <- m;  q_const <- qstart
      f1 <- 0.0; fc <- 0; f1e <- 0.0; cfailure <- 0.0
      q_valuesF <- sort(seq(from=qmin,to=1,by=spac), decreasing = TRUE)
      SQVF <- rep(NA, length(q_valuesF) - 1)
      SQVF[1:length(SQVI)] <- SQVI
      signal1 <- 0; signal2 <- 0
      
      thetahat_all <- matrix(numeric(0),nrow=kk1+kk2,ncol= length(q_valuesF))
      se_smle_all <- matrix(numeric(0),nrow=kk1+kk2,ncol= length(q_valuesF))
      trace_smle_all <- matrix(numeric(0),nrow=1,ncol= length(q_valuesF))	
      
      thetahat_all[,1:Size_gridI] <-  thetahat_I 
      se_smle_all[,1:Size_gridI] <-  se_smleI  
      trace_smle_all[,1:Size_gridI] <- trace_smleI  
      
      #***************Searching for the optimal q from qstart****************#    
      while(grid > 0 && q_const > qmin){
        q_const <- q_values[1.0 + counter]
        alpha_const <- 1-q_const
        res <- suppressWarnings(tryCatch(optim(thetaStart, log_liksurrogate, hessian = F, control=list(fnscale=1, maxit=200), 
                                               gr = Score, method=method), error=function(e) {e}))
        
        
        # print(thetaStart)
        
        if(inherits(res, "error")&& q_const >= qmin){
          counter <- counter + 1.0
          grid = m # Did not converge, so take the next q from the grid
          f1= f1 + 1.0  # Number of failures
          if(f1 == nq - 1.0) {# If the number of failures equals the grid size, use MLE
            grid = 0;	  # grid = 0 means no more q to test
            q_optimal =	  1.0	 # q_optimal = 1.0 means we use MLE
          }
          
          if(f1>=3){
            signal1 <- 1
            grid <- 0
            qM <- as.numeric(na.omit(q_valuesF[-length(q_values)][SQVF<L]))
            if(length(qM)<=3){ # Use MLE
              q_optimal <- 1.0 
              alpha_const <- 0
              
              theta_hat <- theta_mle_full # Parameters estimated via MLE
              score <- Score_LMDPDE(theta_hat, y=y, X=X_full, Z=Z_full, alpha_const =alpha_const, linkmu=linkmu, linkphi=linkphi)
              Vq <- V_matrix_LMDPDE(theta_hat, y=y, X=X_full, Z=Z_full, alpha_const=alpha_const, linkmu=linkmu, linkphi=linkphi)
              
              final_results <- list(theta_hat = theta_hat, 
                                    score = score, Vq = Vq, 
                                    alpha_const = alpha_const)
              
            }else{
              gride <- m
              qtest <- 1.0
              counter_test <- 1
              while(gride>0 && qtest>min(qM)){
                if(is.na(SQVF[counter_test])){
                  gride = m
                }else{
                  if(SQVF[counter_test] < L){
                    gride = gride - 1.0;
                    if(gride==0) q_optimal <- q_valuesF[counter_test-m+1.0] 
                  }  else{  
                    gride = m
                  } 
                } 
                counter_test <- counter_test + 1.0 
                qtest <- q_valuesF[counter_test]                       
              }
              if(gride>0 && qtest==min(qM)){ # Return MLE
                q_optimal <- 1.0     
                alpha_const <- 0
                
                theta_hat <- theta_mle_full # Parameters estimated via MLE
                score <- Score_LMDPDE(theta_hat, y=y, X=X_full, Z=Z_full, alpha_const =alpha_const, linkmu=linkmu, linkphi=linkphi)
                Vq <- V_matrix_LMDPDE(theta_hat, y=y, X=X_full, Z=Z_full, alpha_const=alpha_const, linkmu=linkmu, linkphi=linkphi)
                
                final_results <- list(theta_hat = theta_hat, 
                                      score = score, Vq = Vq, 
                                      alpha_const = alpha_const)
              }   
            }
            
          } # End condition for failures >= 3
        }
        else{# If estimation procedure converges
          
          estimates <- res$par
          #print(estimates)
          log.lik=res$value
          thetahat_n[,counter + 1.0] <- estimates
          Vq <- V_matrix(thetahat_n[,counter + 1.0])	
          se_smle[,counter + 1.0] =  suppressWarnings(t(sqrt(diag(Vq))))	  
          trace_smle[,counter + 1.0] =  sum(diag(Vq))	 
          
          #****Checking the stability of the estimates****#
          SQV[counter] <- round((1.0/(kk1+kk2))*sqrt(sum( (thetahat_n[,counter]/(sqrt(n)*se_smle[,counter]) - 
                                                             thetahat_n[,counter+1.0]/(sqrt(n)*se_smle[,counter + 1.0]))^2)),5) 
          
          SQVF[length(q_valuesF) - length(q_values) + 1.0 + counter] <- SQV[counter]
          thetahat_all[,length(q_valuesF)-length(q_values)+ 1.0 + counter] <- estimates
          se_smle_all[,length(q_valuesF)-length(q_values)+ 1.0 + counter] <- suppressWarnings(t(sqrt(diag(Vq))))
          trace_smle_all[,length(q_valuesF)-length(q_values)+ 1.0 + counter] <- sum(diag(Vq))
          
          #***If SQV is NaN***#
          if(is.nan(SQV[counter])){
            grid = m # If SQV is NaN, take the next q from the grid
            f1= f1 + 1.0  # Number of failures
          }	else{  	
            if(SQV[counter] < L){# If the stability condition is satisfied
              grid = grid - 1.0 # Decrease grid	
              if(grid==0) q_optimal = q_values[counter-m+1.0] # if grid = 0, take the maximum m (i.e., take the maximum q from the grid)
            }  else{  # Stability condition not satisfied  
              grid = m
              f1e = f1e + 1.0       
            }
          }                                                             
          counter = counter + 1.0
          
        }# end 'else' for 'estimation procedure converged
        
        if((grid>0)&&q_const==qmin){# if qmin is reached and there is no stability in the estimate
          grid <- 0; signal2 <- 1
          qME<-   as.numeric(na.omit(q_valuesF[-length(q_values)][SQVF<L]))
          if(length(qME)<=3){ # Return MLE
            q_optimal <- 1.0;  
            alpha_const <- 0
            
            theta_hat <- theta_mle_full # Parameters estimated via MLE
            score <- Score_LMDPDE(theta_hat, y=y, X=X_full, Z=Z_full, alpha_const =alpha_const, linkmu=linkmu, linkphi=linkphi)
            Vq <- V_matrix_LMDPDE(theta_hat, y=y, X=X_full, Z=Z_full, alpha_const=alpha_const, linkmu=linkmu, linkphi=linkphi)
            
            final_results <- list(theta_hat = theta_hat, 
                                  score = score, Vq = Vq, 
                                  alpha_const = alpha_const)
          }else{
            
            grideE <- m
            qtestE <- 1.0
            counter_testE <- 1
            while(grideE>0&&qtestE>min(qME)){
              if(is.na(SQVF[counter_testE])){
                grideE = m
              }else{
                if(SQVF[counter_testE] < L){
                  grideE = grideE - 1.0;
                  if(grideE==0) q_optimal <- q_valuesF[counter_testE-m+1.0] 
                }  else{  
                  grideE = m
                } 
              } 
              counter_testE <- counter_testE + 1.0 
              qtestE <- q_valuesF[counter_testE]                      
            }
            if(grideE>0&&qtestE==min(qME)){ # Return MLE
              q_optimal <- 1.0     
              alpha_const <- 0
              
              theta_hat <- theta_mle_full # Parameters estimated via MLE
              score <- Score_LMDPDE(theta_hat, y=y, X=X_full, Z=Z_full, alpha_const =alpha_const, linkmu=linkmu, linkphi=linkphi)
              Vq <- V_matrix_LMDPDE(theta_hat, y=y, X=X_full, Z=Z_full, alpha_const=alpha_const, linkmu=linkmu, linkphi=linkphi)
              
              final_results <- list(theta_hat = theta_hat, 
                                    score = score, Vq = Vq, 
                                    alpha_const = alpha_const)
            }  
          }
        }
      }
      
      #print(L)
      #print(SQV)
      #***************The search for the optimal q ends here****************# 
      
      #****Selecting the estimates corresponding to the optimal q****#
      
      if(q_optimal==1.0){# Return MLE
        alpha_const <- 0 # If q_optimal is 1, then alpha_const = 0
        theta_hat_full <- theta_mle_full # Parameters estimated via MLE
        score <- Score_LMDPDE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const =alpha_const, linkmu=linkmu, linkphi=linkphi)
        Vq_full <- V_matrix_LMDPDE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const=alpha_const, linkmu=linkmu, linkphi=linkphi)
        
        final_results <- list(theta_hat = theta_hat_full, 
                              score = score, Vq = Vq_full, 
                              alpha_const = alpha_const)
      }
      else
      {
        if(signal1==1||signal2==1){
          seq <- 1:length(q_valuesF)
          index_op <- seq[q_valuesF==q_optimal]
          thehat_n_optimal <- thetahat_all[,index_op]
          
        }else{
          seq <- 1:nq
          index_op <- seq[q_values==q_optimal]
          thehat_n_optimal <- thetahat_n[,index_op]
        } 
        
        theta_hat <- thehat_n_optimal # Final parameters
        # Building the complete theta vector
        theta_hat_full <- fill_theta()
        q_optimal <- as.numeric(q_optimal) # Convert q_optimal to numeric
        alpha_const <- 1 - q_optimal # Calculating alpha_const
        
        # Calculating the final results
        score <- Score_LMDPDE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const =alpha_const, linkmu=linkmu, linkphi=linkphi)
        Vq <- V_matrix_LMDPDE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const=alpha_const, linkmu=linkmu, linkphi=linkphi)
        
        final_results <- list(theta_hat = theta_hat_full, 
                              score = score, Vq = Vq, 
                              alpha_const = alpha_const)
        
        
      }
    } #closes else "#if there is SQV>=L in [0.8, 1]"
    #return(results)
  }# ends option "qoptimal==TRUE"
  
  return(final_results)
  
}# ends function
