# simulate_and_helpers.R
suppressPackageStartupMessages(library(posterior))  # for summarise_draws()

compare_nonzero <- function(fit, true_beta, param_name) {
  vars <- paste0(param_name, "[", seq_along(true_beta), "]")
  sm   <- summarise_draws(fit$draws(variables = vars))  # has q5, q95 by default
  
  sel      <- (sm$q5 > 0) | (sm$q95 < 0)   # estimated nonzero
  true_nz  <- (true_beta != 0)
  
  TP <- sum(sel & true_nz); FP <- sum(sel & !true_nz)
  FN <- sum(!sel & true_nz); TN <- sum(!sel & !true_nz)
  precision <- if ((TP+FP)>0) TP/(TP+FP) else NA_real_
  recall    <- if ((TP+FN)>0) TP/(TP+FN) else NA_real_
  f1        <- if (is.na(precision)||is.na(recall)||(precision+recall)==0) NA_real_ else
    2*precision*recall/(precision+recall)
  
  list(
    summary = data.frame(
      param = param_name,
      p_true_nonzero = mean(true_nz),
      p_est_nonzero  = mean(sel),
      TP, FP, FN, TN, precision, recall, F1 = f1
    ),
    selected_idx      = which(sel),
    missed_true_idx   = which(true_nz & !sel),
    false_pos_idx     = which(sel & !true_nz)
  )
}


bracket_uniroot <- function(f, lower, upper, ..., max_expand = 20, expand = 2) {
  fl <- f(lower); fu <- f(upper)
  if (is.nan(fl) || is.nan(fu)) stop("f(lower) or f(upper) is NaN")
  
  k <- 0
  while (sign(fl) == sign(fu) && k < max_expand) {
    # keep bounds valid for each parameter type
    lower <- lower / expand
    upper <- upper * expand
    # protect positivity
    if (lower <= 0) lower <- .Machine$double.eps
    fl <- f(lower); fu <- f(upper)
    k <- k + 1
  }
  
  if (sign(fl) != sign(fu)) {
    return(uniroot(f, c(lower, upper), ...)$root)
  }
  
  # Fallback: pick the point in [lower, upper] that minimizes |f|
  grid <- exp(seq(log(lower), log(upper), length.out = 201))
  vals <- vapply(grid, f, numeric(1))
  grid[which.min(abs(vals))]
}


# --- 작은 유틸 ---
.is_zero <- function(x, tol = 1e-12) abs(x) <= tol

.bin_metrics <- function(truth, pred) {
  stopifnot(length(truth) == length(pred))
  TP <- sum(truth == 1 & pred == 1)
  TN <- sum(truth == 0 & pred == 0)
  FP <- sum(truth == 0 & pred == 1)
  FN <- sum(truth == 1 & pred == 0)
  
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_   # = TPR
  f1        <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0)
    2 * precision * recall / (precision + recall) else NA_real_
  fpr       <- if ((FP + TN) > 0) FP / (FP + TN) else NA_real_
  fdr       <- if ((TP + FP) > 0) FP / (TP + FP) else NA_real_
  
  tibble::tibble(
    TP = TP, TN = TN, FP = FP, FN = FN,
    TPR = recall, FPR = fpr,
    precision = precision, recall = recall, F1 = f1, FDR = fdr
  )
}

# pr_*의 posterior mean을 뽑기
extract_pr_means <- function(fit) {
  dd <- posterior::as_draws_df(fit)
  keep_cols <- function(prefix) {
    sel <- grep(paste0("^", prefix, "\\["), names(dd), value = TRUE)
    if (length(sel) == 0) sel <- grep(paste0("^", prefix, "$"), names(dd), value = TRUE)
    sel
  }
  cols_both <- keep_cols("pr_both")
  cols_U    <- keep_cols("pr_Uonly")
  cols_C    <- keep_cols("pr_Conly")
  cols_none <- keep_cols("pr_none")
  if (any(lengths(list(cols_both, cols_U, cols_C, cols_none)) == 0))
    stop("pr_both / pr_Uonly / pr_Conly / pr_none not found in draws.")
  
  tibble::tibble(
    p        = seq_along(cols_both),
    pr_both  = sapply(dd[cols_both], mean, na.rm = TRUE),
    pr_Uonly = sapply(dd[cols_U],    mean, na.rm = TRUE),
    pr_Conly = sapply(dd[cols_C],    mean, na.rm = TRUE),
    pr_none  = sapply(dd[cols_none], mean, na.rm = TRUE)
  )
}

# 진실 라벨(비공유/공유 상관없이 U/C가 nonzero인지)
true_nonzero_UC <- function(betaU_true, betaC_true, tol = 1e-12) {
  stopifnot(length(betaU_true) == length(betaC_true))
  tibble::tibble(
    p = seq_along(betaU_true),
    y_nonzero_U = as.integer(!.is_zero(betaU_true, tol)),
    y_nonzero_C = as.integer(!.is_zero(betaC_true, tol))
  )
}


compare_nonzero_from_pr <- function(fit, betaU_true, betaC_true, thr = 0.5) {
  pr  <- extract_pr_means(fit)
  tru <- true_nonzero_UC(betaU_true, betaC_true)
  
  dat <- dplyr::left_join(pr, tru, by = "p") |>
    dplyr::mutate(
      pr_nonzero_U = pr_both + pr_Uonly,
      pr_nonzero_C = pr_both + pr_Conly,
      predU = as.integer(pr_nonzero_U >= thr),
      predC = as.integer(pr_nonzero_C >= thr)
    )
  
  mU <- .bin_metrics(dat$y_nonzero_U, dat$predU) |> dplyr::mutate(target = "nonzero_U", thr = thr, .before = 1)
  mC <- .bin_metrics(dat$y_nonzero_C, dat$predC) |> dplyr::mutate(target = "nonzero_C", thr = thr, .before = 1)
  
  dplyr::bind_rows(mU, mC)
}


## ---- helpers: posterior summaries (q025, mean, median, q975) ----
summarise_draws_any <- function(fit, var_prefix) {
  m <- as.matrix(fit, pars = var_prefix)
  if (is.null(m) || ncol(m) == 0) return(tibble())  
  tibble(
    name   = colnames(m),
    q025   = apply(m, 2, quantile, 0.025),
    mean   = colMeans(m),
    median = apply(m, 2, median),
    q975   = apply(m, 2, quantile, 0.975)
  )
}

summarise_scalars <- function(fit, vars) {
  purrr::map_dfr(vars, \(v) {
    s <- summarise_draws_any(fit, v)
    if (nrow(s) == 0) return(tibble())
    s |>
      mutate(param = v, p = NA_integer_, .before = 1) |>
      select(param, p, q025, mean, median, q975)
  })
}

summarise_vectors <- function(fit, vars) {
  purrr::map_dfr(vars, \(v) {
    s <- summarise_draws_any(fit, v)
    if (nrow(s) == 0) return(tibble())
    s |>
      tidyr::extract(name, into = c("param","p"),
                     regex = "^([A-Za-z_]+)\\[(\\d+)\\]$", convert = TRUE) |>
      select(param, p, q025, mean, median, q975)
  })
}

summarise_post_by_p <- function(fit) {
  # pr_*[p] 각각에 대해 q025/mean/median/q975
  pars <- c("pr_both","pr_Uonly","pr_Conly","pr_none")
  m <- as.matrix(fit, pars = pars)
  if (is.null(m) || ncol(m) == 0) return(tibble())
  dn <- colnames(m)
  long <- tibble(
    var    = dn,
    q025   = apply(m, 2, quantile, 0.025),
    mean   = colMeans(m),
    median = apply(m, 2, median),
    q975   = apply(m, 2, quantile, 0.975)
  ) |>
    tidyr::extract(var, into = c("cat","p"),
                   regex = "^(pr_[A-Za-z]+)\\[(\\d+)\\]$", convert = TRUE)
  
  # wide: cat별 열로 (mean_pr_both, median_pr_both, q025_pr_both, q975_pr_both, ...)
  long |>
    tidyr::pivot_wider(
      names_from = cat,
      values_from = c(q025, mean, median, q975)
    ) |>
    arrange(p)
}

has_par <- function(fit, name) {
  cn <- colnames(as.matrix(fit))
  any(cn == name | startsWith(cn, paste0(name, "[")))
}

save_param_and_post_summaries_q <- function(fit, out_dir, N, P, replicate_id) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 1) scalars
  scal_out <- summarise_scalars(fit, c("alphaU","alphaC"))
  
  # 2) vectors (shared)
  sc_out <- summarise_vectors(fit, c("betaU_sc","betaC_sc"))
  
  # 3) vectors (non-shared) — 존재할 때만
  vec_o <- c()
  if (has_par(fit, "betaU_o")) vec_o <- c(vec_o, "betaU_o")
  if (has_par(fit, "betaC_o")) vec_o <- c(vec_o, "betaC_o")
  o_out <- if (length(vec_o)) summarise_vectors(fit, vec_o) else tibble()
  
  param_tbl <- bind_rows(scal_out, sc_out, o_out) |>
    mutate(N = !!N, P = !!P, replicate = !!replicate_id, .before = 1)
  
  # 4) post (pr_*)
  post_tbl <- summarise_post_by_p(fit) |>
    mutate(N = !!N, P = !!P, replicate = !!replicate_id, .before = 1)
  
  fn_param <- file.path(out_dir, "params_summary.csv")  # 스키마: q025/mean/median/q975
  fn_post  <- file.path(out_dir, "post_summary.csv")    # 스키마: q025_/mean_/median_/q975_ 접두사
  
  readr::write_csv(param_tbl, fn_param, append = file.exists(fn_param))
  if (nrow(post_tbl)) {
    readr::write_csv(post_tbl, fn_post, append = file.exists(fn_post))
  }
}
