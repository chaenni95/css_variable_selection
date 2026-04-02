# the simulation codes is run within HPC format - make sbatch file to run the code. 

suppressPackageStartupMessages({
  library(rstan)
  library(tidyverse)
  library(posterior)   
  library(MASS)
  library(Matrix)
})

# setwd("your directory")
source(dgp.R)
source(simulate_and_helpers.R)

STAN_FILE <- "your/path/css_dir4.stan"  # Stan 파일 경로
if (!file.exists(STAN_FILE)) stop("STAN_FILE not found: ", STAN_FILE)


# ---- set values -------------------------------------------------
N_vals <- 300
P_vals <- c(200, 500)   
REPS   <- 10L           

ITER_WARMUP   <- 5000L
ITER_SAMPLING <- 5000L
CHAINS        <- 2L
ADAPT_DELTA   <- 0.99
MAX_TREEDEPTH <- 14L
THIN          <- 1L

# ---- Stan code compile --------------------------------
message("Compiling Stan model ..."); flush.console()
sm <- rstan::stan_model(file = STAN_FILE)


make_cmdstan_like <- function(fit_rstan) {
  list(draws = function(variables = NULL, format = "draws_matrix") {
    posterior::as_draws_matrix(as.matrix(fit_rstan, pars = variables))
  })
}

has_par <- function(fit, name) {
  cn <- colnames(as.matrix(fit))
  any(cn == name | startsWith(cn, paste0(name, "[")))
}


for (P in P_vals) {
  
  N      <- N_vals
  M      <- P + 1L
  K      <- P + 1L
  
  # pair_id: P=200 → 1, P=500 → 2 
  pair_id    <- which(P_vals == P)
  block_id   <- 1L         
  base_seed  <- 2025L + 1000L * pair_id + 10000L * block_id
  
  message(sprintf("\n===== N=%d, P=%d | REPS=%d | base_seed=%d =====",
                  N, P, REPS, base_seed)); flush.console()
  
  # ---- set file path for output ------------------------------------------
  out_dir   <- file.path("simulation_out_g4_c03_ef05", sprintf("N%d_P%d", N, P))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  fileNZ    <- file.path(out_dir, "results_nonzero_fromPIP.csv")
  fileE     <- file.path(out_dir, "errors.csv")
  fileParam <- file.path(out_dir, "param_summary.csv")
  
  # 파일이 없을 때만 빈 파일 생성
  if (!file.exists(fileNZ))    readr::write_csv(tibble(), fileNZ)
  if (!file.exists(fileE))     readr::write_csv(
    tibble(N=integer(), P=integer(), replicate=integer(), error=character()), fileE)
  if (!file.exists(fileParam)) readr::write_csv(tibble(), fileParam)
  
  # ---- Replicate -----------------------------------------
  for (replicate in seq_len(REPS)) {
    seed <- base_seed + replicate
    
    message(sprintf("  [P=%d] replicate %d/%d (seed=%d)",
                    P, replicate, REPS, seed)); flush.console()
    
    tryCatch({
      sim <- simulate_cure_data(
        N = N, M = M, K = K, P = P, seed = seed,
        target_censor_overall = 0.4,
        target_cure          = 0.2,
        censor_dist          = "exp"
      )
      
      fit <- rstan::sampling(
        object  = sm,
        data    = sim$stan_data,
        chains  = CHAINS,
        iter    = ITER_WARMUP + ITER_SAMPLING,
        warmup  = ITER_WARMUP,
        control = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH),
        seed    = seed,
        refresh = 0,
        thin    = THIN,
        init    = "0"
      )
      
      if (!has_par(fit, "betaU_sc")) stop("Parameter 'betaU_sc' not found in fit.")
      if (!has_par(fit, "betaC_sc")) stop("Parameter 'betaC_sc' not found in fit.")
      
      # nonzero summary 저장
      nz_summary <- compare_nonzero_from_prnz(
        fit,
        betaU_true = sim$betaU_sc,
        betaC_true = sim$betaC_sc,
        thr = 0.5
      ) |>
        dplyr::mutate(
          N = N, P = P, replicate = replicate,
          family = "PIP_nonzero",
          .before = 1
        )
      
      readr::write_csv(
        nz_summary, fileNZ,
        append    = file.exists(fileNZ),
        col_names = !file.exists(fileNZ)
      )
      
      # param summary 저장
      draws_all  <- posterior::as_draws(fit)
      param_summ <- posterior::summarise_draws(draws_all)
      readr::write_csv(
        dplyr::mutate(param_summ, N=N, P=P, replicate=replicate, .before=1),
        fileParam, append = TRUE
      )
      
    }, error = function(e) {
      readr::write_csv(
        tibble(N=N, P=P, replicate=replicate, error=conditionMessage(e)),
        fileE, append = TRUE
      )
      message(sprintf("  ERROR at rep=%d: %s", replicate, conditionMessage(e)))
      flush.console()
    })
  } # end replicate loop
  
  message(sprintf("Done N=%d, P=%d. Outputs -> %s", N, P, out_dir))
  
} # end P loop

message("\n===== All done =====")














