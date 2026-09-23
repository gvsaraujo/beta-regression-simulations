#**********link functions**********#
mean_link <- function(beta, y, X, linkmu="logit"){

  eta <- as.vector(X%*%beta)
  if(linkmu == "logit"){
    mu <- exp(eta)/(1.0+exp(eta))
    T_1 <- mu*(1.0-mu)
    yg1 <- log(y/(1-y))
  }
  if(linkmu == "probit"){
    mu <- pnorm(eta)
    T_1 <- dnorm(qnorm(mu))
    yg1 <- qnorm(y)
  }
  if(linkmu == "cloglog"){
    mu <- 1.0 - exp(-exp(eta))
    T_1 <- (mu - 1)*log(1 - mu)
    yg1 <- log(-log(1-y))
  }
  if(linkmu == "log"){
    mu <- exp(eta)
    T_1 <- mu
    yg1 <- log(y)
  }
  if(linkmu == "loglog"){
    mu <- exp(-exp(-eta))
    T_1 <- -mu*log(mu)
    yg1 <- -log(-log(y))
  }
  if(linkmu == "cauchit"){
    mu <- (pi^(-1))*atan(eta) + 0.5
    T_1 <- (1/pi)*((cospi(mu-0.5))^2)
    yg1 <- tan(pi*(y-0.5))
  }
  results <- list(mu= mu, T_1 = T_1, yg1=yg1)
  return(results)
}#ends mean link function

precision_link <- function(gama, Z, linkphi="log"){

  n <- nrow(Z);  delta <- as.vector(Z%*%gama)
  if(linkphi == "log"){
    phi <- exp(delta)
    T_2 <- phi
  }
  if(linkphi == "identify"){
    phi <- delta
    T_2 <- rep(1,n)
  }
  if(linkphi == "sqrt"){
    phi <- delta^2
    T_2 <- 2*sqrt(phi)
  }
  results <- list(phi=phi, T_2=T_2)
  return(results)
}