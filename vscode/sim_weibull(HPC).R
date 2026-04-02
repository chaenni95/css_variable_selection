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


N_SCENARIOS <- 2L
# block number by scenario 
BLOCKS_PER_SCENARIO <- 10L

TASK_ID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
if (is.na(TASK_ID) || TASK_ID < 1L || TASK_ID > N_SCENARIOS * BLOCKS_PER_SCENARIO) {
  stop("SLURM_ARRAY_TASK_ID must be in 1..", N_SCENARIOS * BLOCKS_PER_SCENARIO,
       " (got ", TASK_ID, ")")
}

# scenario: 1..2  (1=P=200, 2=P=500)
pair_id  <- ((TASK_ID - 1L) %/% BLOCKS_PER_SCENARIO) + 1L
# block: 1..10 (ten replicates per block)
block_id <- ((TASK_ID - 1L) %%  BLOCKS_PER_SCENARIO) + 1L

N_vals <- 300
P_vals <- c(200, 500)
grid <- expand.grid(N = N_vals, P = P_vals, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
grid <- grid[order(grid$N, grid$P), ]  # row1: (300,200), row2: (300,500)

N <- grid$N[pair_id]
P <- grid$P[pair_id]
M <- P + 1L
K <- P + 1L

## ---------------- 반복 횟수 ----------------
# 10 array * 10 replicates
REPS <- 10L

## ---------------- seed 설계 ----------------
base_seed <- 2025L + 1000L * pair_id + 10000L * block_id

message(sprintf("TASK %d => pair_id=%d (N=%d,P=%d), block_id=%d, REPS=%d, base_seed=%d",
                TASK_ID, pair_id, N, P, block_id, REPS, base_seed)); flush.console()

## ---------------- Sampler controls ----------------
ITER_WARMUP   <- as.integer(Sys.getenv("ITER_WARMUP",   "5000"))
ITER_SAMPLING <- as.integer(Sys.getenv("ITER_SAMPLING", "5000"))
CHAINS        <- as.integer(Sys.getenv("CHAINS",        "2"))
ADAPT_DELTA   <- as.numeric(Sys.getenv("ADAPT_DELTA",   "0.99"))
MAX_TREEDEPTH <- as.integer(Sys.getenv("MAX_TREEDEPTH", "14"))
THIN          <- as.integer(Sys.getenv("THIN",          "1"))

rstan_options(auto_write = TRUE)
options(mc.cores = 1L)
Sys.setenv(OMP_NUM_THREADS = "1")


STAN_FILE <- Sys.getenv("STAN_FILE", "/vscode/css_dir4.stan") # set your directory
if (!file.exists(STAN_FILE)) stop("STAN_FILE not found: ", STAN_FILE)

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

## ---------------- Output paths ----------------
out_dir <- file.path("simulation_out_g4_c03_ef05", sprintf("N%d_P%d", N, P))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fileNZ <- file.path(out_dir, "results_nonzero_fromPIP.csv")
fileE <- file.path(out_dir, "errors.csv")
fileParam <- file.path(out_dir, "param_summary.csv")

if (!file.exists(fileNZ)) readr::write_csv(tibble(), fileNZ)
if (!file.exists(fileE)) readr::write_csv(tibble(N=integer(),P=integer(),
                                                 replicate=integer(), error=character()), fileE)
if (!file.exists(fileParam)) readr::write_csv(tibble(), fileParam)  
## ---------------- Replicate loop ONLY ----------------
for (replicate in seq_len(REPS)) {
  seed <- base_seed + replicate
  
  tryCatch({
    sim <- simulate_cure_data(N=N, P=P, seed=seed, 
                              blocks = 5, rho = 0.5, rho_both_true = 0.6, 
                              censoring_rate = 0.3, 
                              alpha = 1.5)
    
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
      init = "0"
    )
    if (!has_par(fit, "betaU_sc")) stop("Parameter 'betaU_sc' not found in fit.")
    if (!has_par(fit, "betaC_sc")) stop("Parameter 'betaC_sc' not found in fit.")
    
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
      append = file.exists(fileNZ),
      col_names = !file.exists(fileNZ)
    )
    
    draws_all <- posterior::as_draws(fit)
    param_summ <- posterior::summarise_draws(draws_all)
    param_summ_out <- dplyr::mutate(param_summ, N=N, P=P, replicate=replicate, .before=1)
    readr::write_csv(param_summ_out, fileParam, append = TRUE)
    
  }, error = function(e) {
    readr::write_csv(tibble(N=N, P=P, replicate=replicate, error=conditionMessage(e)), fileE, append=TRUE)
    message(sprintf("Error at N=%d, P=%d, rep=%d: %s", N, P, replicate, conditionMessage(e)))
    flush.console()
  })
}

message(sprintf("Done pair N=%d, P=%d. Outputs in %s\nErrors: %s", N, P, out_dir, fileE))
