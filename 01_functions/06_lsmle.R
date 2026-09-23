# HERE RETURNS THE VALUES OF THE SCORE FUNCTION
# Link functions: mean_link()/precision_link() from 01_link_functions.R
# (previously duplicated here as mean_link_functions_v2/precision_link_functions_v2)

#********************Score-type function*******************#
Score_LSMLE <- function(theta, y, X, Z, alpha_const, linkmu, linkphi) {

  kk1 <- ncol(X)
  kk2 <- ncol(Z)
  beta <- theta[1:kk1]
  gama <- theta[(kk1+1.0):(kk1+kk2)]
  mean_info <- mean_link(beta = beta, y=y, X=X, linkmu=linkmu)
  mu <- mean_info$mu
  T_1_alpha <- mean_info$T_1
  precision_info <- precision_link(gama = gama, Z=Z, linkphi=linkphi)
  phi_alpha <- precision_info$phi
  T_2_alpha <- precision_info$T_2
  phi <- phi_alpha/(1 - alpha_const)
  a <- mu*phi
  b <- (1.0 - mu)*phi
  m_phi <- phi
  ystar <- log(y/(1.0 - y))
  ydagger <- log(1.0 - y)
  F_q <- ((y*(1-y))*dbeta(y, mu*phi, (1-mu)*phi))^alpha_const
  mustar <- psigamma(a, 0) - psigamma(b, 0)
  mudagger <-  psigamma(b, 0) - psigamma(a + b, 0)
  #****vector type Score*****#
  avScore <- numeric(kk1 + kk2)
  avScore[1:kk1] <- as.vector(crossprod(X, m_phi*F_q*T_1_alpha*(ystar - mustar)))
  avScore[(kk1+1.0):(kk1+kk2)] <- ((1-alpha_const)^(-1))*as.vector(crossprod(Z, F_q*T_2_alpha*(mu*(ystar - mustar) + ydagger - mudagger)))
  return(avScore)
}#ends score-type function

#********************Covariance matrix****************#
V_matrix_LSMLE <- function(theta, y, X, Z, alpha_const, linkmu, linkphi) {
  kk1 <- ncol(X)
  kk2 <- ncol(Z)
  beta_p <- theta[1:kk1]
  gama_p <- theta[(kk1+1):(kk1+kk2)]
  
  mean_info <- mean_link(beta = beta_p, y=y, X=X, linkmu=linkmu)
  muhat <- mean_info$mu
  T_1   <- mean_info$T_1
  precision_info <- precision_link(gama = gama_p, Z=Z, linkphi=linkphi)
  phihat_alpha <- precision_info$phi
  T_2_alpha <- precision_info$T_2

  ahat_alpha <- muhat*phihat_alpha
  bhat_alpha <- (1.0 - muhat)*phihat_alpha

  phihat_n <- phihat_alpha/(1-alpha_const)
  muhat_n <- 	muhat

  ahat_n <- muhat_n*phihat_n
  bhat_n <- (1.0 - muhat_n)*phihat_n
  #
  phihat_1_alpha <- phihat_n*(1 + alpha_const)	#expression of phi_(1 + alpha)
  #
  a1_alphahat <- muhat_n*phihat_1_alpha
  b1_alphahat <- (1.0 - muhat_n)*phihat_1_alpha

  m_phi <- phihat_n
  psi1_n <- psigamma(ahat_n, 1.0) # psi(mu*phi)
  psi2_n <- psigamma(bhat_n, 1.0) # psi((1-mu)*phi)
  psi3_n <- psigamma(phihat_n, 1.0) #psi(phi)

  V <- psi1_n + psi2_n	# ok
  B1 <- exp((1 - alpha_const)*(lgamma(ahat_n) + lgamma(bhat_n) - lgamma(phihat_n)) - (lgamma(ahat_alpha) + lgamma(bhat_alpha) - lgamma(phihat_alpha))) #ok
  B2 <- exp(lgamma(a1_alphahat) + lgamma(b1_alphahat) - lgamma(phihat_1_alpha) - (2.0*(alpha_const)*(lgamma(ahat_n) + lgamma(bhat_n) - lgamma(phihat_n))  +
                                                                                     lgamma(ahat_alpha) + lgamma(bhat_alpha) - lgamma(phihat_alpha)))	# ok

  C <- phihat_n*(muhat_n*psi1_n - (1.0 - muhat_n)*psi2_n)	# ok
  D <- (muhat_n^2.0)*psi1_n + ((1.0 - muhat_n)^2.0)*psi2_n - psi3_n # ok

  psi1_1_alpha <- psigamma(a1_alphahat, 1.0)
  psi2_1_alpha <- psigamma(b1_alphahat, 1.0)
  psi3_1_alpha <- psigamma(phihat_1_alpha, 1.0)

  V_1_alpha <- psi1_1_alpha + psi2_1_alpha # ok
  C_1_alpha <- phihat_n*(muhat_n*psi1_1_alpha - (1.0 - muhat_n)*psi2_1_alpha) # ok
  D_1_alpha <- (muhat_n^2.0)*psi1_1_alpha + ((1.0 - muhat_n)^2.0)*psi2_1_alpha - psi3_1_alpha # ok

  Jalpha_betabeta <- (1 - alpha_const)*crossprod(X, X*(B1*T_1^2.0*m_phi^2.0*V))
  Jalpha_betagamma <- crossprod(X, Z*(B1*T_1*T_2_alpha*C))
  Jalpha_gammagamma <- ((1 - alpha_const)^(-1))*crossprod(Z, Z*(B1*T_2_alpha^2.0*D))

  Jalpha <- matrix(numeric(0), kk1+kk2, kk1+kk2)
  Jalpha[1:kk1,1:kk1] <- Jalpha_betabeta
  Jalpha[1:kk1,(kk1+1):(kk1+kk2)] <- Jalpha_betagamma
  Jalpha[(kk1+1):(kk1+kk2),1:kk1] <- t(Jalpha_betagamma)
  Jalpha[(kk1+1):(kk1+kk2),(kk1+1):(kk1+kk2)] <- t(Jalpha_gammagamma)
  Jalpha <- -Jalpha

  Kalpha_betabeta <- crossprod(X, X*(B2*T_1^2.0*m_phi^2.0*V_1_alpha))
  Kalpha_betagamma <- ((1 - alpha_const)^(-1))*crossprod(X, Z*(B2*T_1*T_2_alpha*C_1_alpha))
  Kalpha_gammagamma <- ((1 - alpha_const)^(-2))*crossprod(Z, Z*(B2*T_2_alpha^2.0*D_1_alpha))
  Kalpha <- matrix(numeric(0),kk1+kk2,kk1+kk2)
  Kalpha[1:kk1,1:kk1] <- Kalpha_betabeta
  Kalpha[1:kk1,(kk1+1):(kk1+kk2)] <- Kalpha_betagamma
  Kalpha[(kk1+1):(kk1+kk2),1:kk1] <- t(Kalpha_betagamma)
  Kalpha[(kk1+1):(kk1+kk2),(kk1+1):(kk1+kk2)] <- t(Kalpha_gammagamma)
  Vq <- tryCatch( solve(Jalpha)%*%Kalpha%*%t(solve(Jalpha)), error=function(e) {e})  #asymptotic covariance matrix
  if(inherits(Vq, "error")){
    Vq <- MASS::ginv(Jalpha)%*%Kalpha%*%t(MASS::ginv(Jalpha)) 
  }
  return(list(Vq=Vq, H = Jalpha, J = Kalpha))
}#ends Covariance matrix function


#**********LSMLE function**********#
LSMLE_BETA <- function(y, X, Z, qoptimal=TRUE, q0=1, m=3, L=0.02, qmin=0.5,
                       spac=0.02, method="BFGS", startV ="CP", linkmu="logit", linkphi="log",
                       beta_fix = NULL, idx_beta = NULL,  
                       gamma_fix = NULL, idx_gamma = NULL){
  
  #**********Objective function**********#
  log_liksurrogate <- function(theta) {
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
    
    # Calculate predicted values of mu and phi
    mu <- mean_link(beta = beta, y = y, X = X_full, linkmu = linkmu)$mu
    phi_alpha <- precision_link(gama = gama, Z = Z_full, linkphi = linkphi)$phi
    phi <- phi_alpha / (1 - alpha_const)
    y_star <- log(y/(1 - y))
    log_likS <- sum(((y * (1 - y)) * dbeta(y, mu * phi, (1 - mu) * phi))^alpha_const)
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

    # Calculate predicted values of mu and phi
    mean_info <- mean_link(beta = beta, y = y, X = X_full, linkmu = linkmu)
    mu <- mean_info$mu
    T_1_alpha <- mean_info$T_1
    precision_info <- precision_link(gama = gama, Z = Z_full, linkphi = linkphi)
    phi_alpha <- precision_info$phi
    T_2_alpha <- precision_info$T_2
    phi <- phi_alpha / (1 - alpha_const)
    a <- mu * phi
    b <- (1 - mu) * phi
    m_phi <- phi
    ystar <- log(y / (1 - y))
    ydagger <- log(1 - y)
    F_q <- ((y * (1 - y)) * dbeta(y, mu * phi, (1 - mu) * phi))^alpha_const
    mustar <- psigamma(a, 0) - psigamma(b, 0)
    mudagger <- psigamma(b, 0) - psigamma(a + b, 0)
    avScore_full <- numeric(kk1 + kk2)
    avScore_full[1:kk1] <- as.vector(crossprod(X_full, m_phi * F_q * T_1_alpha * (ystar - mustar)))
    avScore_full[(kk1 + 1):(kk1 + kk2)] <- ((1 - alpha_const)^(-1)) *
      as.vector(crossprod(Z_full, F_q * T_2_alpha * (mu * (ystar - mustar) + ydagger - mudagger)))
    score_out <- c(avScore_full[idx_beta_free], avScore_full[kk1 + idx_gamma_free])
    return(score_out)
  }
  
  #****Covariance function V***#
  V_matrix <- function(theta) {
    kk1 <- ncol(X_full) # Number of columns in X_full
    kk2 <- ncol(Z_full) # Number of columns in Z_full
    beta_p <- numeric(kk1) # Initialize beta vector
    gama_p <- numeric(kk2) # Initialize gamma vector
    
    idx_beta_free <- setdiff(1:kk1, idx_beta) # Indices of free beta
    idx_gamma_free <- setdiff(1:kk2, idx_gamma) # Indices of free gamma
    beta_p[idx_beta] <- beta_fix # Assign fixed values to beta
    beta_p[idx_beta_free] <- theta[1:length(idx_beta_free)] # Assign free values to beta
    gama_p[idx_gamma] <- gamma_fix # Assign fixed values to gamma
    gama_p[idx_gamma_free] <- theta[(length(idx_beta_free) + 1):length(theta)] # Assign free values to gamma
    #print(beta_p)
    #print(gama_p)
    # Calculate predicted values of mu and phi
    mean_info <- mean_link(beta = beta_p, y = y, X = X_full, linkmu = linkmu)
    muhat <- mean_info$mu
    T_1 <- mean_info$T_1
    precision_info <- precision_link(gama = gama_p, Z = Z_full, linkphi = linkphi)
    phihat_alpha <- precision_info$phi
    T_2_alpha <- precision_info$T_2
    ahat_alpha <- muhat * phihat_alpha
    bhat_alpha <- (1.0 - muhat) * phihat_alpha
    phihat_n <- phihat_alpha / (1 - alpha_const)
    muhat_n <- muhat
    ahat_n <- muhat_n * phihat_n
    bhat_n <- (1.0 - muhat_n) * phihat_n
    phihat_1_alpha <- phihat_n * (1 + alpha_const)
    a1_alphahat <- muhat_n * phihat_1_alpha
    b1_alphahat <- (1.0 - muhat_n) * phihat_1_alpha
    m_phi <- phihat_n
    psi1_n <- psigamma(ahat_n, 1.0)
    psi2_n <- psigamma(bhat_n, 1.0)
    psi3_n <- psigamma(phihat_n, 1.0)
    V <- psi1_n + psi2_n
    B1 <- exp((1 - alpha_const) * (lgamma(ahat_n) + lgamma(bhat_n) - lgamma(phihat_n)) -
                     (lgamma(ahat_alpha) + lgamma(bhat_alpha) - lgamma(phihat_alpha)))
    B2 <- exp(lgamma(a1_alphahat) + lgamma(b1_alphahat) - lgamma(phihat_1_alpha) -
                     (2.0 * alpha_const * (lgamma(ahat_n) + lgamma(bhat_n) - lgamma(phihat_n)) +
                        lgamma(ahat_alpha) + lgamma(bhat_alpha) - lgamma(phihat_alpha)))
    C <- phihat_n * (muhat_n * psi1_n - (1.0 - muhat_n) * psi2_n)
    D <- (muhat_n^2.0) * psi1_n + ((1.0 - muhat_n)^2.0) * psi2_n - psi3_n
    psi1_1_alpha <- psigamma(a1_alphahat, 1.0)
    psi2_1_alpha <- psigamma(b1_alphahat, 1.0)
    psi3_1_alpha <- psigamma(phihat_1_alpha, 1.0)
    V_1_alpha <- psi1_1_alpha + psi2_1_alpha
    C_1_alpha <- phihat_n * (muhat_n * psi1_1_alpha - (1.0 - muhat_n) * psi2_1_alpha)
    D_1_alpha <- (muhat_n^2.0) * psi1_1_alpha + ((1.0 - muhat_n)^2.0) * psi2_1_alpha - psi3_1_alpha
    Jalpha_betabeta <- (1 - alpha_const) * crossprod(X_full, X_full * (B1 * T_1^2 * m_phi^2 * V))
    Jalpha_betagamma <- crossprod(X_full, Z_full * (B1 * T_1 * T_2_alpha * C))
    Jalpha_gammagamma <- (1 - alpha_const)^(-1) * crossprod(Z_full, Z_full * (B1 * T_2_alpha^2 * D))
    Jalpha <- matrix(0, kk1 + kk2, kk1 + kk2)
    Jalpha[1:kk1, 1:kk1] <- Jalpha_betabeta
    Jalpha[1:kk1, (kk1 + 1):(kk1 + kk2)] <- Jalpha_betagamma
    Jalpha[(kk1 + 1):(kk1 + kk2), 1:kk1] <- t(Jalpha_betagamma)
    Jalpha[(kk1 + 1):(kk1 + kk2), (kk1 + 1):(kk1 + kk2)] <- Jalpha_gammagamma
    Jalpha <- -Jalpha
    Kalpha_betabeta <- crossprod(X_full, X_full * (B2 * T_1^2 * m_phi^2 * V_1_alpha))
    Kalpha_betagamma <- (1 - alpha_const)^(-1) * crossprod(X_full, Z_full * (B2 * T_1 * T_2_alpha * C_1_alpha))
    Kalpha_gammagamma <- (1 - alpha_const)^(-2) * crossprod(Z_full, Z_full * (B2 * T_2_alpha^2 * D_1_alpha))
    Kalpha <- matrix(0, kk1 + kk2, kk1 + kk2)
    Kalpha[1:kk1, 1:kk1] <- Kalpha_betabeta
    Kalpha[1:kk1, (kk1 + 1):(kk1 + kk2)] <- Kalpha_betagamma
    Kalpha[(kk1 + 1):(kk1 + kk2), 1:kk1] <- t(Kalpha_betagamma)
    Kalpha[(kk1 + 1):(kk1 + kk2), (kk1 + 1):(kk1 + kk2)] <- Kalpha_gammagamma
    Vq_full <- tryCatch(
      solve(Jalpha) %*% Kalpha %*% solve(Jalpha),
      error = function(e) MASS::ginv(Jalpha) %*% Kalpha %*% MASS::ginv(Jalpha)
    )
    idx_all <- c(idx_beta_free, kk1 + idx_gamma_free)
    Vq_restricted <- Vq_full[idx_all, idx_all, drop = FALSE]
    return(Vq_restricted)
  }
  
  
  #----------------------------------------------------------------------------------------------------------------#
  #------------------------------------->Estimation Process<-----------------------------------------#
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
    score <- Score_LSMLE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const = alpha_const, linkmu=linkmu, linkphi=linkphi)
    Vq <- V_matrix_LSMLE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const=alpha_const, linkmu=linkmu, linkphi=linkphi)
    
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
      score <- Score_LSMLE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const =alpha_const, linkmu=linkmu, linkphi=linkphi)
      Vq <- V_matrix_LSMLE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const=alpha_const, linkmu=linkmu, linkphi=linkphi)
      
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
      res <- suppressWarnings(tryCatch(optim(thetaStart, log_liksurrogate, hessian = F, control=list(fnscale=-1,maxit=200), 
                                             gr = Score, method=method), error=function(e) {e}))#maximization
      
      # If optimization fails, return MLE results
      if(inherits(res, "error")){
        alpha_const <- 0
        theta_hat <- SV$theta_hat
        theta_hat_full <- fill_theta() # Completing the theta vector
        score <- Score_LSMLE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const =alpha_const, linkmu=linkmu, linkphi=linkphi)
        Vq <- V_matrix_LSMLE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const=alpha_const, linkmu=linkmu, linkphi=linkphi)
        
        final_results <- list(theta_hat = theta_hat_full, 
                              score = score, Vq = Vq, alpha_const = alpha_const)
        return(final_results)
      }
      else{ # If optimization converged
        if(res$convergence == 0||res$convergence == 1){ 
          alpha_const <- 1 - q_const
          theta_hat <- res$par # Estimated parameters
          theta_hat_full <- fill_theta() # Completing the theta vector
          
          score <- Score_LSMLE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const =alpha_const, linkmu=linkmu, linkphi=linkphi)
          Vq <- V_matrix_LSMLE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const=alpha_const, linkmu=linkmu, linkphi=linkphi)
          
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
          score_full <- Score_LSMLE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const =alpha_const, linkmu=linkmu, linkphi=linkphi)
          Vq_full <- V_matrix_LSMLE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const=alpha_const, linkmu=linkmu, linkphi=linkphi)
          
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
      resT <- suppressWarnings(tryCatch(optim(thetaStart, log_liksurrogate, hessian = F, control=list(fnscale=-1, maxit=200), 
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
      
      score <- Score_LSMLE(theta_hat, y=y, X=X_full, Z=Z_full, alpha_const =alpha_const, linkmu=linkmu, linkphi=linkphi)
      Vq <- V_matrix_LSMLE(theta_hat, y=y, X=X_full, Z=Z_full, alpha_const=alpha_const, linkmu=linkmu, linkphi=linkphi)
      
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
        res <- suppressWarnings(tryCatch(optim(thetaStart, log_liksurrogate, hessian = F, control=list(fnscale=-1, maxit=200), 
                                               gr = Score, method=method), error=function(e) {e}))
        
        
        # print(thetaStart)
        
        if(inherits(res, "error")&& q_const >= qmin){
          counter <- counter + 1.0
          grid = m # Did not converge, so take the next q from the grid
          f1= f1 + 1.0  # Number of failures
          if(f1 == nq - 1.0) {# If the number of failures equals the grid size, use MLE
            grid = 0;	  # grid = 0 means no more q to test
            q_optimal =	  1.0	 # q_optimal = 1.0 means use MLE
          }
          
          if(f1>=3){
            signal1 <- 1
            grid <- 0
            qM <- as.numeric(na.omit(q_valuesF[-length(q_values)][SQVF<L]))
            if(length(qM)<=3){ # Use MLE
              q_optimal <- 1.0 
              alpha_const <- 0
              
              theta_hat <- theta_mle_full # Parameters estimated via MLE
              score <- Score_LSMLE(theta_hat, y=y, X=X_full, Z=Z_full, alpha_const =alpha_const, linkmu=linkmu, linkphi=linkphi)
              Vq <- V_matrix_LSMLE(theta_hat, y=y, X=X_full, Z=Z_full, alpha_const=alpha_const, linkmu=linkmu, linkphi=linkphi)
              
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
                score <- Score_LSMLE(theta_hat, y=y, X=X_full, Z=Z_full, alpha_const =alpha_const, linkmu=linkmu, linkphi=linkphi)
                Vq <- V_matrix_LSMLE(theta_hat, y=y, X=X_full, Z=Z_full, alpha_const=alpha_const, linkmu=linkmu, linkphi=linkphi)
                
                final_results <- list(theta_hat = theta_hat, 
                                      score = score, Vq = Vq, 
                                      alpha_const = alpha_const)
              }   
            }
            
          } # End of condition failures >= 3
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
          
          #***And in case of NaN in SQV***#
          if(is.nan(SQV[counter])){
            grid = m # If SQV is NaN, take the next q from the grid
            f1= f1 + 1.0  # Number of failures
          }	else{  	
            if(SQV[counter] < L){# If the stability condition is satisfied
              grid = grid - 1.0 # Decrease the grid	
              if(grid==0) q_optimal = q_values[counter-m+1.0] # if grid = 0, take the maximum m (i.e., take the maximum q from the grid)
            }  else{  # Stability condition not satisfied  
              grid = m
              f1e = f1e + 1.0       
            }
          }                                                             
          counter = counter + 1.0
          
        }# end of 'else' for 'estimation procedure converged'
        
        if((grid>0)&&q_const==qmin){# if qmin is reached and there is no stability in the estimate
          grid <- 0; signal2 <- 1
          qME<-   as.numeric(na.omit(q_valuesF[-length(q_values)][SQVF<L]))
          if(length(qME)<=3){ # Return MLE
            q_optimal <- 1.0;  
            alpha_const <- 0
            
            theta_hat <- theta_mle_full # Parameters estimated via MLE
            score <- Score_LSMLE(theta_hat, y=y, X=X_full, Z=Z_full, alpha_const =alpha_const, linkmu=linkmu, linkphi=linkphi)
            Vq <- V_matrix_LSMLE(theta_hat, y=y, X=X_full, Z=Z_full, alpha_const=alpha_const, linkmu=linkmu, linkphi=linkphi)
            
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
              score <- Score_LSMLE(theta_hat, y=y, X=X_full, Z=Z_full, alpha_const =alpha_const, linkmu=linkmu, linkphi=linkphi)
              Vq <- V_matrix_LSMLE(theta_hat, y=y, X=X_full, Z=Z_full, alpha_const=alpha_const, linkmu=linkmu, linkphi=linkphi)
              
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
        score <- Score_LSMLE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const =alpha_const, linkmu=linkmu, linkphi=linkphi)
        Vq_full <- V_matrix_LSMLE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const=alpha_const, linkmu=linkmu, linkphi=linkphi)
        
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
        q_optimal <- as.numeric(q_optimal) # Converting q_optimal to numeric
        alpha_const <- 1 - q_optimal # Calculating alpha_const
        
        # Calculating the final results
        score <- Score_LSMLE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const =alpha_const, linkmu=linkmu, linkphi=linkphi)
        Vq <- V_matrix_LSMLE(theta_hat_full, y=y, X=X_full, Z=Z_full, alpha_const=alpha_const, linkmu=linkmu, linkphi=linkphi)
        
        final_results <- list(theta_hat = theta_hat_full, 
                              score = score, Vq = Vq, 
                              alpha_const = alpha_const)
        
        
      }
    } #closes else "#if there is SQV>=L in [0.8, 1]"
    #return(results)
  }# ends option "qoptimal==TRUE"
  
  return(final_results)
  
}# ends function

