fit_dina_bmm_simple <- function(Y, Q, prior_type = "Unif", maxitr = 2000, 
                                epsilon = 1e-04, bound = 1e-04, verbose = 0){
  Y <- as.matrix(Y)
  Q <- as.matrix(Q)
  
  K <- ncol(Q)
  N <- nrow(Y)
  J <- ncol(Y)
  patterns <- attributepattern(K)
  
  R <- compute_R(patterns, Q)
  sg <- cbind(rep(0.2, J), rep(0.2, J))
  p <- rep(1/2^K, 2^K)
  
  ConstrMatrix <- matrix(c(-1, 1), 1, 2)
  
  prior_parms <- create_prior_parms_simple(J, K, prior_type)
  
  deltam1 <- prior_parms$delta - 1
  slip_am1 <- prior_parms$slip_a - 1
  slip_bm1 <- prior_parms$slip_b - 1
  guess_am1 <- prior_parms$guess_a - 1
  guess_bm1 <- prior_parms$guess_b - 1
  
  itr <- 0
  
  parm <- c(p, sg[1,1], sg[1,2])
  
  while(itr < maxitr){
    Z <- compute_z_post(Y, R, p, sg)
    Eta <- Z%*%R
    
    S1 <- sum((1-Y)*Eta)
    S2 <- sum(Y*Eta)
    G1 <- sum(Y*(1-Eta))
    G2 <- sum((1-Y)*(1-Eta))
    
    if((S1 + S2 + slip_am1 + slip_bm1) == 0 |
       (G1 + G2 + guess_am1 + guess_bm1) == 0){
      
      warning(paste("Nj contains 0"), call. = FALSE)
      
      return(list(success=FALSE))
      
    }
    
    #EM and BM estimates
    PJ <- c((G1 + guess_am1)/(G1 + G2 + guess_am1 + guess_bm1),
            1-(S1 + slip_am1)/(S1 + S2 + slip_am1 + slip_bm1))
    PJ[PJ <= bound] <- bound
    PJ[PJ >= 1 - bound] <- 1 - bound
    
    if(any(c(ConstrMatrix%*%PJ) < 0)){
      
      
      obj <- function(x0){
        -1*((G1+guess_am1)*log(x0[1]) + (G2+guess_bm1)*log(1-x0[1]) +
              (S1+slip_am1)*log(1-x0[2]) + (S2+slip_bm1)*log(x0[2]))
      }
      
      dev <- function(x0){
        -1*(c((G1+guess_am1)/x0[1] - (G2+guess_bm1)/(1-x0[1]),
              -(S1+slip_am1)/(1-x0[2]) + (S2+slip_bm1)/x0[2]))
      }
      
      ineq <- function(x0){
        c(ConstrMatrix %*% x0)
      }
      
      ineq.jac <- function(x0){
        ConstrMatrix
      }
      
      x00 <- c(sg[1,2], 1-sg[1,1])
      
      x00[x00 < bound] <- bound
      x00[x00 > 1-bound] <- 1-bound
      
      optims <- suppressMessages({try(slsqp(x0 = x00,
                                            fn = obj, 
                                            gr = dev,
                                            hin = ineq, 
                                            hinjac = ineq.jac,
                                            lower = rep(bound, length(x00)),
                                            upper = rep(1-bound,length(x00))), silent = TRUE)})
      
      if(inherits(optims,"try-error")){
        warning(paste("Optimization failed"),call. = FALSE)
        if(verbose==1)
          print(optims)
        return(list(success=FALSE))
      }
      
      sg <- cbind(rep(1-optims$par[2], J), rep(optims$par[1], J))
    }else{
      sg <- cbind(rep(1-PJ[2], J), rep(PJ[1], J))
    }
    
    p <- (colSums(Z) + deltam1)/(N + sum(deltam1))
    
    parm_new <- c(p, sg[1,1], sg[1,2])
    
    max_chg = max(abs(parm_new-parm),na.rm = TRUE)
    
    parm <- parm_new
    
    itr <- itr + 1
    
    if (max_chg < epsilon) break
  }
  
  Z <- compute_z_post(Y, R, p, sg)
  
  alpha_hat <- patterns[apply(Z, 1, which.max),]
  slip_hat <- sg[1,1]
  guess_hat <- sg[1,2]
  p_hat <- p
  
  results <- list(slip_hat = slip_hat, guess_hat = guess_hat, 
                  p_hat = p_hat, alpha_hat = alpha_hat)
  
  results
}

fit_dina_bmm <- function(Y, Q, prior_type = "Info", maxitr = 2000, 
                         epsilon = 1e-04, bound = 1e-04, verbose = 0){
  Y <- as.matrix(Y)
  Q <- as.matrix(Q)
  
  K <- ncol(Q)
  N <- nrow(Y)
  J <- ncol(Y)
  patterns <- attributepattern(K)
  
  R <- compute_R(patterns, Q)
  sg <- cbind(rep(0.2, J), rep(0.2, J))
  p <- rep(1/2^K, 2^K)
  
  ConstrMatrix <- matrix(c(-1, 1), 1, 2)
  
  prior_parms <- create_prior_parms_full(J, K, prior_type)
  
  deltam1 <- prior_parms$delta - 1
  slip_am1 <- prior_parms$slip_a - 1
  slip_bm1 <- prior_parms$slip_b - 1
  guess_am1 <- prior_parms$guess_a - 1
  guess_bm1 <- prior_parms$guess_b - 1

  itr <- 0
  
  parm <- c(p, c(sg))
  
  while(itr < maxitr){
    Z <- compute_z_post(Y, R, p, sg)
    Eta <- Z%*%R
    
    S1 <- colSums((1-Y)*Eta)
    S2 <- colSums(Y*Eta)
    G1 <- colSums(Y*(1-Eta))
    G2 <- colSums((1-Y)*(1-Eta))
    
    for(j in 1:J){
      
      S1_j <- S1[j]
      S2_j <- S2[j]
      G1_j <- G1[j]
      G2_j <- G2[j]
      
      slip_am1_j <- slip_am1[j]
      slip_bm1_j <- slip_bm1[j]
      guess_am1_j <- guess_am1[j]
      guess_bm1_j <- guess_bm1[j]
      
      if((S1_j + S2_j + slip_am1_j + slip_bm1_j) == 0 |
         (G1_j + G2_j + guess_am1_j + guess_bm1_j) == 0){
        
        warning(paste("Nj contains 0 for item", j), call. = FALSE)
        
        # cat("\nFor item ",j,"\n")
        # print(data.frame(Rj = Rj,Nj = Nj,Pj = Rj/Nj))
        
        return(list(success=FALSE))
        
      }
      
      #EM and BM estimates
      Pj <- c((G1_j + guess_am1_j)/(G1_j + G2_j + guess_am1_j + guess_bm1_j),
              1-(S1_j + slip_am1_j)/(S1_j + S2_j + slip_am1_j + slip_bm1_j))
      Pj[Pj <= bound] <- bound
      Pj[Pj >= 1 - bound] <- 1 - bound
      
      if(any(c(ConstrMatrix%*%Pj) < 0)){
        
        
        obj <- function(x0){
          -1*((G1_j+guess_am1_j)*log(x0[1]) + (G2_j+guess_bm1_j)*log(1-x0[1]) +
                (S1_j+slip_am1_j)*log(1-x0[2]) + (S2_j+slip_bm1_j)*log(x0[2]))
        }
        
        dev <- function(x0){
          -1*(c((G1_j+guess_am1_j)/x0[1] - (G2_j+guess_bm1_j)/(1-x0[1]),
                -(S1_j+slip_am1_j)/(1-x0[2]) + (S2_j+slip_bm1_j)/x0[2]))
        }
        
        ineq <- function(x0){
          c(ConstrMatrix %*% x0)
        }
        
        ineq.jac <- function(x0){
          ConstrMatrix
        }
        
        x00 <- c(sg[j,2], 1-sg[j,1])
        
        x00[x00 < bound] <- bound
        x00[x00 > 1-bound] <- 1-bound
        
        optims <- suppressMessages({try(slsqp(x0 = x00,
                                              fn = obj, 
                                              gr = dev,
                                              hin = ineq, 
                                              hinjac = ineq.jac,
                                              lower = rep(bound, length(x00)),
                                              upper = rep(1-bound,length(x00))), silent = TRUE)})
        
        if(inherits(optims,"try-error")){
          warning(paste("Optimization failed for item", j),call. = FALSE)
          if(verbose==1)
            print(optims)
          return(list(success=FALSE))
        }
        
        sg[j,] <- c(1-optims$par[2], optims$par[1])
      }else{
        sg[j,] <- c(1-Pj[2], Pj[1])
      }
    }
    
    p <- (colSums(Z) + deltam1)/(N + sum(deltam1))
    
    parm_new <- c(p, c(sg))
    
    max_chg = max(abs(parm_new-parm),na.rm = TRUE)
    
    parm <- parm_new
    
    itr <- itr + 1
    
    if (max_chg < epsilon) break
  }
  
  Z <- compute_z_post(Y, R, p, sg)
  
  alpha_hat <- patterns[apply(Z, 1, which.max),]
  slip_hat <- sg[,1]
  guess_hat <- sg[,2]
  p_hat <- p
  
  results <- list(slip_hat = slip_hat, guess_hat = guess_hat, 
                  p_hat = p_hat, alpha_hat = alpha_hat)
  
  results
}

fit_dina_cat_bmm <- function(Y, Q, tau_a = 2, tau_b = 1, tau_init = 4,
                             maxitr = 2000, epsilon = 1e-04, 
                             bound = 1e-04, verbose = 0){
  Y <- as.matrix(Y)
  Q <- as.matrix(Q)
  
  K <- ncol(Q)
  N <- nrow(Y)
  J <- ncol(Y)
  patterns <- attributepattern(K)
  
  R <- compute_R(patterns, Q)
  sg <- cbind(rep(0.2, J), rep(0.2, J))
  p <- rep(1/2^K, 2^K)

  tau <- tau_init

  prior_gen_mod <- fit_dina_bmm_simple(Y, Q)
  slip_hat0 <- prior_gen_mod$slip
  guess_hat0 <- prior_gen_mod$guess
  p_hat0 <- prior_gen_mod$p_hat
  p_hat0[p_hat0 == 0] <- .Machine$double.xmin
  pj_hat0 <- p_hat0%*%R
  
  deltam1 <- tau*p_hat0
  slip_am1 <- tau*slip_hat0*pj_hat0
  slip_bm1 <- tau*(1-slip_hat0)*pj_hat0
  guess_am1 <- tau*guess_hat0*(1-pj_hat0)
  guess_bm1 <- tau*(1-guess_hat0)*(1-pj_hat0)
  
  ConstrMatrix <- matrix(c(-1, 1), 1, 2)
  
  itr <- 0
  
  parm <- c(p, c(sg), tau)
  
  while(itr < maxitr){
    Z <- compute_z_post(Y, R, p, sg)
    Eta <- Z%*%R
    
    S1 <- colSums((1-Y)*Eta)
    S2 <- colSums(Y*Eta)
    G1 <- colSums(Y*(1-Eta))
    G2 <- colSums((1-Y)*(1-Eta))
    
    for(j in 1:J){
      
      S1_j <- S1[j]
      S2_j <- S2[j]
      G1_j <- G1[j]
      G2_j <- G2[j]
      
      slip_am1_j <- slip_am1[j]
      slip_bm1_j <- slip_bm1[j]
      guess_am1_j <- guess_am1[j]
      guess_bm1_j <- guess_bm1[j]
      
      if((S1_j + S2_j + slip_am1_j + slip_bm1_j) == 0 |
         (G1_j + G2_j + guess_am1_j + guess_bm1_j) == 0){
        
        warning(paste("Nj contains 0 for item", j), call. = FALSE)
        
        # cat("\nFor item ",j,"\n")
        # print(data.frame(Rj = Rj,Nj = Nj,Pj = Rj/Nj))
        
        return(list(success=FALSE))
        
      }
      
      #EM and BM estimates
      Pj <- c((G1_j + guess_am1_j)/(G1_j + G2_j + guess_am1_j + guess_bm1_j),
              1-(S1_j + slip_am1_j)/(S1_j + S2_j + slip_am1_j + slip_bm1_j))
      Pj[Pj <= bound] <- bound
      Pj[Pj >= 1 - bound] <- 1 - bound
      
      if(any(c(ConstrMatrix%*%Pj) < 0)){
        
        
        obj <- function(x0){
          -1*((G1_j+guess_am1_j)*log(x0[1]) + (G2_j+guess_bm1_j)*log(1-x0[1]) +
                (S1_j+slip_am1_j)*log(1-x0[2]) + (S2_j+slip_bm1_j)*log(x0[2]))
        }
        
        dev <- function(x0){
          -1*(c((G1_j+guess_am1_j)/x0[1] - (G2_j+guess_bm1_j)/(1-x0[1]),
                -(S1_j+slip_am1_j)/(1-x0[2]) + (S2_j+slip_bm1_j)/x0[2]))
        }
        
        ineq <- function(x0){
          c(ConstrMatrix %*% x0)
        }
        
        ineq.jac <- function(x0){
          ConstrMatrix
        }
        
        x00 <- c(sg[j,2], 1-sg[j,1])
        
        x00[x00 < bound] <- bound
        x00[x00 > 1-bound] <- 1-bound
        
        optims <- suppressMessages({try(slsqp(x0 = x00,
                                              fn = obj, 
                                              gr = dev,
                                              hin = ineq, 
                                              hinjac = ineq.jac,
                                              lower = rep(bound, length(x00)),
                                              upper = rep(1-bound,length(x00))), silent = TRUE)})
        
        if(inherits(optims,"try-error")){
          warning(paste("Optimization failed for item", j),call. = FALSE)
          if(verbose==1)
            print(optims)
          return(list(success=FALSE))
        }
        
        sg[j,] <- c(1-optims$par[2], optims$par[1])
      }else{
        sg[j,] <- c(1-Pj[2], Pj[1])
      }
    }
    
    p <- (colSums(Z) + deltam1)/(N + sum(deltam1))
    p[p == 0] <- .Machine$double.xmin
    
    tau_optims <- optimise(compute_tau_lp, lower = bound, upper = 100,
                           p = p, p_hat = p_hat0, 
                           slip = sg[,1], s_hat = slip_hat0,
                           guess = sg[,2], g_hat = guess_hat0,
                           pj_hat = pj_hat0, 
                           tau_a = tau_a, tau_b = tau_b, maximum = TRUE)
    
    tau <- tau_optims$maximum
    
    deltam1 <- tau*p_hat0
    slip_am1 <- tau*slip_hat0*pj_hat0
    slip_bm1 <- tau*(1-slip_hat0)*pj_hat0
    guess_am1 <- tau*guess_hat0*(1-pj_hat0)
    guess_bm1 <- tau*(1-guess_hat0)*(1-pj_hat0)
    
    parm_new <- c(p, c(sg), tau)
    
    max_chg = max(abs(parm_new-parm), na.rm = TRUE)
    
    parm <- parm_new
    
    itr <- itr + 1
    
    if (max_chg < epsilon) break
  }
  
  Z <- compute_z_post(Y, R, p, sg)
  
  alpha_hat <- patterns[apply(Z, 1, which.max),]
  slip_hat <- sg[,1]
  guess_hat <- sg[,2]
  p_hat <- p
  tau_hat <- tau
  
  results <- list(slip_hat = slip_hat, guess_hat = guess_hat, 
                  p_hat = p_hat, alpha_hat = alpha_hat, tau_hat = tau_hat,
                  gen_mod_p = p_hat0, gen_mod_slip = slip_hat0, 
                  gen_mod_guess = guess_hat0)
  
  results
}

fit_dina_ncpd <- function(Y, Q){
  npcd_fit <- AlphaNP(Y, Q, gate = "AND", method = "Hamming")
  
  npcd_fit
}

fit_dina_simple_mcmc <- function(Y, Q, init = NULL, prior_type = "Unif", 
                                 niter = 2000, nburn = 200){
  N <- nrow(Y)
  J <- ncol(Y)
  K <- ncol(Q)
  index_vector <- 1:2^K
  
  draws <- create_draws_container(N, J, K, niter, nburn)
  Z_draws <- draws$Z
  p_draws <- draws$p
  slip_draws <- draws$slip
  guess_draws <- draws$guess
  
  patterns <- unname(attributepattern(K))
  R <- compute_R(patterns, Q)
  
  init <- initialize_parms(N, J, K, init)
  
  Z <- init$Z
  p <- init$p
  sg <- init$sg
  
  prior <- create_prior_parms_simple(J, K, prior_type)
  
  delta <- prior$delta
  slip_a <- prior$slip_a
  slip_b <- prior$slip_b
  guess_a <- prior$guess_a
  guess_b <- prior$guess_b
  
  for (i in 1:niter){
    Z <- sample_Z(Y, R, p, sg)
    p <- sample_p(Z, delta)
    sg <- sample_sg_simple(Y, Z, R, sg, p, slip_a, slip_b, guess_a, guess_b)
    
    if (i > nburn){
      Z_draws[i-nburn,] <- Z%*%index_vector
      p_draws[i-nburn,] <- p
      slip_draws[i-nburn,] <- sg[,1]
      guess_draws[i-nburn,] <- sg[,2]
    }
  }
  
  list(Z = Z_draws, p = p_draws, slip = slip_draws, guess = guess_draws)
}

fit_dina_mcmc <- function(Y, Q, init = NULL, prior_type = "Unif", 
                          niter = 2000, nburn = 200){
  N <- nrow(Y)
  J <- ncol(Y)
  K <- ncol(Q)
  index_vector <- 1:2^K
  
  draws <- create_draws_container(N, J, K, niter, nburn)
  Z_draws <- draws$Z
  p_draws <- draws$p
  slip_draws <- draws$slip
  guess_draws <- draws$guess
  
  patterns <- unname(attributepattern(K))
  R <- compute_R(patterns, Q)
  
  init <- initialize_parms(N, J, K, init)
  
  Z <- init$Z
  p <- init$p
  sg <- init$sg
  
  prior <- create_prior_parms_full(J, K, prior_type)
  
  delta <- prior$delta
  slip_a <- prior$slip_a
  slip_b <- prior$slip_b
  guess_a <- prior$guess_a
  guess_b <- prior$guess_b
  
  for (i in 1:niter){
    Z <- sample_Z(Y, R, p, sg)
    p <- sample_p(Z, delta)
    sg <- sample_sg_full(Y, Z, R, sg, p, slip_a, slip_b, guess_a, guess_b)
    
    if (i > nburn){
      Z_draws[i-nburn,] <- Z%*%index_vector
      p_draws[i-nburn,] <- p
      slip_draws[i-nburn,] <- sg[,1]
      guess_draws[i-nburn,] <- sg[,2]
    }
  }
  
  list(Z = Z_draws, p = p_draws, slip = slip_draws, guess = guess_draws)
}

fit_dina_cat_mcmc <- function(Y, Q, init = NULL, niter = 2000, 
                              nburn = 200, tau_a = 2, tau_b = 1, 
                              tau_sd = 10, tau_init = 15){
  N <- nrow(Y)
  J <- ncol(Y)
  K <- ncol(Q)
  index_vector <- 1:2^K
  
  draws <- create_draws_container(N, J, K, niter, nburn)
  Z_draws <- draws$Z
  p_draws <- draws$p
  slip_draws <- draws$slip
  guess_draws <- draws$guess
  tau_draws <- rep(0, niter-nburn)
  
  patterns <- unname(attributepattern(K))
  R <- compute_R(patterns, Q)
  
  prior_gen_mod <- fit_dina_bmm_simple(Y, Q)
  slip_hat <- prior_gen_mod$slip
  guess_hat <- prior_gen_mod$guess
  p_hat <- prior_gen_mod$p_hat
  
  init <- initialize_parms(N, J, K, init)
  
  Z <- init$Z
  p <- init$p
  sg <- init$sg
  
  tau <- tau_init
  
  for (i in 1:niter){
    Z <- sample_Z(Y, R, p, sg)
    p <- sample_p_cat(Z, tau, p_hat)
    sg <- sample_sg_cat(Y, Z, R, sg, p, tau, slip_hat, guess_hat, p_hat)
    tau <- sample_tau_mh(tau, sg, p, slip_hat, guess_hat, p_hat, R, 
                         tau_a, tau_b, tau_sd = tau_sd)
    
    if (i > nburn){
      Z_draws[i-nburn,] <- Z%*%index_vector
      p_draws[i-nburn,] <- p
      slip_draws[i-nburn,] <- sg[,1]
      guess_draws[i-nburn,] <- sg[,2]
      tau_draws[i-nburn] <- tau
    }
  }
  
  list(Z = Z_draws, p = p_draws, slip = slip_draws, 
       guess = guess_draws, tau = tau_draws, 
       gen_mod_p = p_hat, gen_mod_slip = slip_hat, gen_mod_guess = guess_hat)
}
