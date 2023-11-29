compute_p_mvnorm <- function(R, K, patterns){
  Sigma <- diag(1-R, K)+R
  Mu <- rep(0, K)
  cutoffs <- 1:K/(K+1)
  
  lower <- t(apply(patterns, 1, function(x) ifelse(x == 0, -Inf, cutoffs)))
  upper <- t(apply(patterns, 1, function(x) ifelse(x == 0, cutoffs, Inf)))
  
  p_mvnorm <- sapply(1:2^K, function(x){
    pmvnorm(lower = lower[x,], upper = upper[x,], corr = Sigma)
  })
  
  p <- p_mvnorm/sum(p_mvnorm)
  
  p
}

generate_sim_data <- function(Q, N, R, slip, guess, p_list){
  K <- ncol(Q)
  patterns <- unname(attributepattern(K))
  
  p <- p_list[[paste(R, K, sep = "_")]]
  
  dina_data <- sim_dina(N, Q, slip, guess, p, patterns)
  
  Y <- dina_data$Y
  alpha <- dina_data$alpha
  
  sim_data <- list(Y = Y, Q = Q, alpha = alpha, p = p, slip = slip, guess = guess)
  
  sim_data
}

create_draws_container = function(N, J, K, niter, nburn){
  Z_draws <- matrix(0, niter-nburn, N)
  p_draws <- matrix(0, niter-nburn, 2^K)
  slip_draws <- matrix(0, niter-nburn, J)
  guess_draws <- matrix(0, niter-nburn, J)
  
  draws <- list(Z = Z_draws, p = p_draws, slip = slip_draws, guess = guess_draws)
  
  draws
}

initialize_parms <- function(N, J, K, init = NULL){
  if (is.null(init)){
    Z_diag <- diag(2^K)
    Z_init <- Z_diag[sample(1:2^K, N, TRUE), ]
    p_init <- rep(1/2^K, 2^K)
    sg_init <- matrix(0.1, J, 2)
    
    init <- list(Z = Z_init, p = p_init, sg = sg_init)
  } else{
    init
  }
  
  init
}

create_prior_parms_simple <- function(J, K, prior_type = "Unif"){
  
  if (prior_type == "Unif"){
    delta <- rep(1, 2^K)
    slip_a <- 1
    slip_b <- 1
    guess_a <- 1
    guess_b <- 1 
  } else if (prior_type == "Info"){
    delta <- rep(1, 2^K)
    slip_a <- 1.5
    slip_b <- 2.5
    guess_a <- 1.5
    guess_b <- 2.5
  }
  
  prior_parms <- list(delta = delta, slip_a = slip_a, slip_b = slip_b,
                      guess_a = guess_a, guess_b = guess_b)
  prior_parms
}

create_prior_parms_full <- function(J, K, prior_type = "Unif"){
  
  if (prior_type == "Unif"){
    delta <- rep(1, 2^K)
    slip_a <- rep(1, J)
    slip_b <- rep(1, J)
    guess_a <- rep(1, J)
    guess_b <- rep(1, J)
  } else if (prior_type == "Info"){
    delta <- rep(1, 2^K)
    slip_a <- rep(1.5, J)
    slip_b <- rep(2.5, J)
    guess_a <- rep(1.5, J)
    guess_b <- rep(2.5, J)
  }
  
  prior_parms <- list(delta = delta, slip_a = slip_a, slip_b = slip_b,
                      guess_a = guess_a, guess_b = guess_b)
  prior_parms
}

compute_z_post <- function(Y, R, p, sg){
  alpha_post_unnorm <- compute_z_post_unnorm(Y, R, p, sg)
  alpha_post <- alpha_post_unnorm/row_sums(alpha_post_unnorm)
  
  alpha_post
}

sample_Z <- function(Y, R, p, sg){
  alpha_post <- compute_z_post(Y, R, p, sg)
  Z <- sample_all_multinomial(alpha_post)
  
  Z
}

sample_p <- function(Z, delta){
  delta_post <- colSums(Z) + delta
  
  p <- rdirichlet(1, delta_post)
  
  p
}

sample_p_cat <- function(Z, tau, p_hat){
  delta_post <- colSums(Z) + tau*p_hat + 1
  
  p <- rdirichlet(1, delta_post)
  p[p == 0] <- .Machine$double.xmin
  
  p
}

compute_eta_Z <- function(Z, R){
  eta <- Z%*%R
  
  eta
}

rtbeta <- function(N, shape1, shape2, lower, upper){
  p_lower <- pbeta(lower, shape1, shape2)
  p_upper <- pbeta(upper, shape1, shape2)
  u <- runif(N, p_lower, p_upper)
  rv <- qbeta(u, shape1, shape2)
  
  rv
}

sample_sg_simple <- function(Y, Z, R, sg, p, slip_a, slip_b, guess_a, guess_b){
  N <- nrow(Y)
  J <- ncol(Y)
  slip <- sg[1,1]
  
  Eta <- compute_eta_Z(Z, R)
  
  n_par <- 1
  slip_post_a <- sum((1-Y)*Eta) + slip_a
  slip_post_b <- sum(Y*Eta) + slip_b
  guess_post_a <- sum(Y*(1-Eta)) + guess_a
  guess_post_b <- sum((1-Y)*(1-Eta)) + guess_b
  
  guess <- rtbeta(n_par, guess_post_a, guess_post_b, 0, 1-slip)
  slip <- rtbeta(n_par, slip_post_a, slip_post_b, 0, 1-guess)
  
  sg <- cbind(rep(slip, J), rep(guess, J))
  
  sg
}

sample_sg_full <- function(Y, Z, R, sg, p, slip_a, slip_b, guess_a, guess_b){
  N <- nrow(Y)
  J <- ncol(Y)
  slip <- sg[,1]
  
  Eta <- compute_eta_Z(Z, R)
  
  n_par <- length(slip)
  slip_post_a <- colSums((1-Y)*Eta) + slip_a
  slip_post_b <- colSums(Y*Eta) + slip_b
  guess_post_a <- colSums(Y*(1-Eta)) + guess_a
  guess_post_b <- colSums((1-Y)*(1-Eta)) + guess_b
  
  guess <- rtbeta(n_par, guess_post_a, guess_post_b, 0, 1-slip)
  slip <- rtbeta(n_par, slip_post_a, slip_post_b, 0, 1-guess)
  
  sg <- cbind(slip, guess)
  
  sg
}

sample_sg_cat <- function(Y, Z, R, sg, p, tau, slip_hat, guess_hat, p_hat){
  N <- nrow(Y)
  J <- ncol(Y)
  slip <- sg[,1]
  
  Eta <- compute_eta_Z(Z, R)
  
  pj_hat <- p_hat%*%R
  
  slip_a <- tau*slip_hat*pj_hat + 1
  slip_b <- tau*(1-slip_hat)*pj_hat + 1
  guess_a <- tau*guess_hat*(1-pj_hat) + 1
  guess_b <- tau*(1-guess_hat)*(1-pj_hat) + 1
  
  n_par <- length(slip)
  slip_post_a <- colSums((1-Y)*Eta) + slip_a
  slip_post_b <- colSums(Y*Eta) + slip_b
  guess_post_a <- colSums(Y*(1-Eta)) + guess_a
  guess_post_b <- colSums((1-Y)*(1-Eta)) + guess_b
  
  guess <- rtbeta(n_par, guess_post_a, guess_post_b, 0, 1-slip)
  slip <- rtbeta(n_par, slip_post_a, slip_post_b, 0, 1-guess)
  
  sg <- cbind(slip, guess)
  
  sg
}

compute_R = function(patterns, Q){
  num_patt <- nrow(patterns)
  PtQ <- patterns%*%t(Q)
  H <- rbind(row_sums(Q))[rep(1, num_patt), ]
  R_matrix <- (PtQ >= H)*1
  
  R_matrix
}

sim_dina <- function(N, Q, slip, guess, p, patterns){
  P <- rbind(p)[rep(1, N),]
  Z <- sample_all_multinomial(P)
  alpha <- patterns[apply(Z, 1, function(x) which(x == 1)),]
  R <- compute_R(patterns, Q)
  eta <- Z%*%R
  dina_probs <- compute_dina_probs(eta, slip, guess)
  Y <- sample_dina_data(dina_probs)
  
  dina_data <- list(Y = Y, alpha = alpha)
  
  dina_data
}

sample_tau_mh <- function(tau, sg, p, slip_hat, guess_hat, p_hat, R, 
                          tau_a, tau_b, tau_sd = 0.5){
  
  tau_star <- rnorm(1, tau, tau_sd)
  
  if (tau_star > 0){
    pj_hat <- p_hat%*%R
    
    lp_tau_cur <- compute_tau_lp(tau, p, p_hat, sg[,1], slip_hat, sg[,2], 
                                 guess_hat, pj_hat, tau_a, tau_b)
    
    lp_tau_star <- compute_tau_lp(tau_star, p, p_hat, sg[,1], slip_hat, sg[,2], 
                                  guess_hat, pj_hat, tau_a, tau_b)
    
    log_acc_prob <- lp_tau_star - lp_tau_cur
    
    if (log(runif(1)) < log_acc_prob){
      tau <- tau_star
    }
  }
  
  tau
}

get_p_from_alpha <- function(patterns, alpha_hat){
  pattern_table <- rep(0, nrow(patterns))
  names(pattern_table) <- apply(patterns, 1, paste, collapse = "")
  pattern_table_unnord <- table(apply(alpha_hat, 1, paste, collapse = ""))
  pattern_table[names(pattern_table_unnord)] <- pattern_table_unnord
  p_hat <- pattern_table/sum(pattern_table)
  
  p_hat
}

compute_eta_alpha <- function(A, Q, H){
  AQt <- A%*%t(Q)
  AQtH <- AQt >= H
  
  AQtH
}

create_H <- function(Y, Q){
  N <- nrow(Y)
  H <- t(replicate(N, rowSums(Q)))
  
  H
}

get_sg_from_alpha <- function(Y, Q, alpha){
  H <- create_H(Y, Q)
  eta <- compute_eta_alpha(alpha, Q, H)*1
  slip <- (colSums((1-Y)*eta)+1)/(colSums(eta)+1)
  guess <- (colSums(Y*(1-eta))+1)/(colSums(1-eta)+1)
  sg <- cbind(slip, guess)
  
  sg
}

get_Z_from_alpha <- function(alpha){
  K <- ncol(alpha)
  Z_base <- diag(2^K)
  patterns <- unname(attributepattern(K))
  patt_string <- apply(patterns, 1, paste, collapse = "")
  alpha_string <- apply(alpha, 1, paste, collapse = "")
  Z_ind <- unname(sapply(alpha_string, function(x)
    which(x == patt_string)))
  Z <- Z_base[Z_ind,]
  
  Z
}

extract_estimates <- function(Y, Q, model_object, model){
  if (model == "dina_simple_mcmc" |
      model == "dina_full_mcmc"){
    patterns <- unname(attributepattern(ncol(Q)))
    s_hat <- colMeans(model_object$slip)
    g_hat <- colMeans(model_object$guess)
    p_hat <- colMeans(model_object$p)
    alpha_hat <- patterns[apply(model_object$Z, 2, function(x){
      as.numeric(names(which.max(table(x))))
    }),]
  } else if (model == "dina_cat_mcmc"){
    patterns <- unname(attributepattern(ncol(Q)))
    s_hat <- colMeans(model_object$slip)
    g_hat <- colMeans(model_object$guess)
    p_hat <- colMeans(model_object$p)
    alpha_hat <- patterns[apply(model_object$Z, 2, function(x){
      as.numeric(names(which.max(table(x))))
    }),]
    tau_hat <- mean(model_object$tau)
  } else if (model == "dina_bmm"){
    s_hat <- model_object$slip_hat
    g_hat <- model_object$guess_hat
    p_hat <- model_object$p_hat
    alpha_hat <- model_object$alpha_hat
  }
  else if (model == "dina_cat_bmm"){
    s_hat <- model_object$slip_hat
    g_hat <- model_object$guess_hat
    p_hat <- model_object$p_hat
    alpha_hat <- model_object$alpha_hat
    tau_hat <- model_object$tau_hat
  } else if (model == "dina_npcd"){
    patterns <- unname(attributepattern(ncol(Q)))
    alpha_hat <- unname(model_object$alpha.est)
    p_hat <- unname(get_p_from_alpha(patterns, alpha_hat))
    sg <- get_sg_from_alpha(Y, Q, alpha_hat)
    s_hat <- sg[,1]
    g_hat <- sg[,2]
  }
  
  if (model == "dina_cat_mcmc" | model == "dina_cat_bmm"){
    estimates <- list(slip_hat = s_hat, guess_hat = g_hat, 
                      p_hat = p_hat, alpha_hat = alpha_hat,
                      tau_hat = tau_hat)
  } else {
    estimates <- list(slip_hat = s_hat, guess_hat = g_hat, 
                      p_hat = p_hat, alpha_hat = alpha_hat) 
  }
  
  estimates
}

get_metrics <- function(sim_datasets, estimates){
  num_metrics <- 4
  num_combs <- length(sim_datasets)
  num_models <- length(estimates[[1]])
  
  metrics <- vector("list", length(num_models))
  
  for (i in 1:num_models){
    metrics_i <- matrix(0, num_combs, num_metrics)
    
    for (j in 1:num_combs){
      K_j <- ncol(estimates[[j]][[i]]$alpha_hat)
      
      if (!is.null(K_j)){
        slip_mse <- mean((estimates[[j]][[i]]$slip_hat - sim_datasets[[j]]$slip)^2)
        guess_mse <- mean((estimates[[j]][[i]]$guess_hat - sim_datasets[[j]]$guess)^2)
        p_mse <- mean((estimates[[j]][[i]]$p_hat - sim_datasets[[j]]$p)^2)
        prof_accuracy <- mean(rowSums(estimates[[j]][[i]]$alpha_hat == sim_datasets[[j]]$alpha) == K_j) 
        metrics_i[j,] <- c(slip_mse, guess_mse, p_mse, prof_accuracy)
      } else {
        metrics_i[j,] <- rep(NA, 4)
      }
    }
    colnames(metrics_i) <- c("Slip_MSE", "Guess_MSE", "P_MSE", "Accuracy")
    metrics[[i]] <- metrics_i
  }
  
  metrics
}