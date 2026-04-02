#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rstan)
  library(tidyverse)
  library(posterior)   # for as_draws_matrix
  library(MASS)
  library(Matrix)
})

source("simulate_and_helpers.R")

simulate_cure_data <- function(
    N, M, K, P,
    seed = 2025,
    # pi_true = c(0.02, 0.01, 0.01, 0.96),   # [both, U-only, C-only, none]
    deterministic = TRUE,                  # divide P exactly, rounding
    tauU_both_true = 0.6,
    tauC_both_true = 0.6,
    rho_both_true  = 0.6,
    sigma_true = 0.5,                      # Weibull sigma (Stan: shape=1/sigma)
    target_median_uncured = 10,            # uncured median at x=0
    target_cure = 0.2,                     # target cure
    target_censor_overall = 0.4,           
    censor_dist = c("exp","unif"),         
    Cmax = 50,
    blocks = 5,
    rho = 0.5
) {
  stopifnot(P <= min(M, K))
  set.seed(seed)
  
  
  # counts <- .alloc_counts(P, pi_true)    
  counts <- c(5, 5, 5, P - 15)
  names(counts) <- c("both","Uonly","Conly","none")
  
  # block diagonal matrix for the covariates ####################################
  blocks <- blocks
  rho <- rho
  sizes <- rep(P %/% blocks, blocks)
  sizes[seq_len(P %% blocks)] <- sizes[seq_len(P %% blocks)] + 1
  blks <- lapply(sizes, function(b) outer(1:b, 1:b, function(i,j) rho^abs(i-j)))
  Sig  <- as.matrix(Matrix::bdiag(blks))
  L <- chol(Sig)
  sc <- matrix(rnorm(N * P), N , P) %*% L
  ################################################################################
  
  
  xonly <- if (M > P) as.matrix(scale(matrix(rnorm(N * (M-P)), N, (M-P))))
  else matrix(numeric(N*0), N, 0)
  zonly <- if (K > P) as.matrix(scale(matrix(rnorm(N * (K-P)), N, (K-P))))
  else matrix(numeric(N*0), N, 0)
  
  pat <- unlist(mapply(function(lbl, cnt) rep(lbl, cnt), 1:4, counts))
  if (length(pat) != P) stop("pattern length mismatch")
  
  
  betaU.sc <- numeric(P)
  betaC.sc <- numeric(P)
  
  wb <- which(pat == 1)
  if (length(wb) > 0) {
    Sigma_rho <- matrix(c(1, rho_both_true, rho_both_true, 1), 2, 2)
    Z <- matrix(MASS::mvrnorm(n = length(wb), mu = c(0, 0), Sigma = Sigma_rho), ncol = 2)
    
    U <- pnorm(Z)
    mag <- 0.5 + 0.3 * U
    sgn <- ifelse(Z >= 0, 1, -1)
    
    betaU.sc[wb] <- sgn[, 1] * mag[, 1]
    betaC.sc[wb] <- sgn[, 2] * mag[, 2]
  }
  
  
  
  # U-only: 
  wu <- which(pat == 2)
  if (length(wu) > 0) {
    # betaU.sc[wu] <- sample(c(-1.2, -1, 0.8, -0.8, 1, 1.2), length(wu), replace = TRUE)
    betaU.sc[wu] <- sample(c(0.8, -0.8, 0.5, -0.5), length(wu), replace = TRUE)
    betaC.sc[wu] <- 0
  }
  
  # C-only: 
  wc <- which(pat == 3)
  if (length(wc) > 0) {
    betaU.sc[wc] <- 0
    # betaC.sc[wc] <- sample(c(-1.2, -1, 0.8, -0.8, 1, 1.2), length(wc), replace = TRUE)
    betaC.sc[wc] <- sample(c(0.8, -0.8, 0.5, -0.5), length(wc), replace = TRUE)
  }
  
  ## -----------------------------
  ## 5. intercept 
  ## -----------------------------
  betaU.o <- if (M > P) numeric(M - P) else numeric(0)
  betaC.o <- if (K > P) numeric(K - P) else numeric(0)
  
  shape_true <- 1 / sigma_true
  lambda0 <- target_median_uncured / (log(2))^(1 / shape_true)
  alphaU  <- log(lambda0)
  
  linC   <- as.numeric(sc %*% betaC.sc + if (length(betaC.o)) zonly %*% betaC.o else 0)
  alphaC <- uniroot(
    f = function(a) mean(plogis(a + linC)) - target_cure,
    interval = c(-15, 15)
  )$root
  prop <- plogis(alphaC + linC)  # 개인별 cure prob (정확히 target에 맞춰짐)
  
  linU      <- as.numeric(sc %*% betaU.sc + if (length(betaU.o)) xonly %*% betaU.o else 0)
  scale_vec <- as.numeric(exp(alphaU + linU))
  
  is_cured <- rbinom(N, 1, prop)  # 1=cured
  Tevent   <- rep(Inf, N)
  idx_unc  <- which(is_cured == 0)
  if (length(idx_unc) > 0) {
    Tevent[idx_unc] <- rweibull(length(idx_unc), shape = shape_true, scale = scale_vec[idx_unc])
  }
  
  censor_dist <- match.arg(censor_dist)

  U_cens <- runif(N)
  make_censor <- function(par) {
    if (censor_dist == "exp") {
      C <- -log(U_cens) / pmax(par, 1e-12)
    } else {
      C <- 1 + (par - 1) * U_cens
    }
    pmin(C, Cmax)  
  }
  
  obj_overall <- function(par) {
    C <- make_censor(par)
    t_obs    <- pmin(Tevent, C)
    is_event <- as.integer(Tevent <= C)  # 1=event observed, 0=censored
    mean(1 - is_event) - target_censor_overall
  }
  
  if (censor_dist == "exp") {
    lambda_opt <- bracket_uniroot(obj_overall, lower = 1e-6, upper = 10)
    if (is.list(lambda_opt)) lambda_opt <- lambda_opt$root
    lambda_opt <- as.numeric(lambda_opt)
    Cens <- make_censor(lambda_opt)
  } else {
    Cmax_opt <- bracket_uniroot(obj_overall, lower = 1 + 1e-6, upper = max(Cmax, 5))
    if (is.list(Cmax_opt)) Cmax_opt <- Cmax_opt$root
    Cmax_opt <- as.numeric(Cmax_opt)
    Cens <- make_censor(Cmax_opt)
  }
  t_obs       <- numeric(N)
  is_event    <- integer(N)  # 1=event, 0=censored  (== is_censored)
  cured       <- integer(N)  # 1=cured
  
  
  t_obs    <- pmin(Tevent, Cens)
  is_event <- as.integer(Tevent <= Cens)
  cured    <- is_cured
  
  
  stan_data <- list(
    N = N,
    t = as.vector(pmax(t_obs, 1e-6)),
    is_censored = is_event,
    M = M, sc = matrix(sc, nrow = N), xonly = matrix(xonly, nrow = N), zonly = matrix(zonly, nrow = N),
    K = K, P = P
  )
  
  truth <- list(
    alphaU = alphaU, alphaC = alphaC,
    betaU_sc = betaU.sc, betaC_sc = betaC.sc,
    betaU_o = betaU.o,  betaC_o = betaC.o,
    sigma = sigma_true,
    pattern = pat,                    # 1=both, 2=U-only, 3=C-only, 4=none
    counts = counts,
    # pi_true = pi_true,
    realized_pi = counts / P
  )
  
  list(stan_data = stan_data, truth = truth,
       generated = list(t_obs = t_obs, is_event = is_event, cured = cured))
}


N = 200 ; M = 300; K = 300; P = 200; seed = 12
sim <- simulate_cure_data(N=N, M=M, K=K, P=P, seed=seed,
                          target_censor_overall = 0.6,
                          target_cure = 0.4,
                          censor_dist = "exp")
