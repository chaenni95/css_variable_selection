// Mixture cure rate model with Weibull AFT and logistic incidence
// NOTE: is_censored == 1 means event observed
functions {
  real cen_weibull_lpdf(real t, int is_censored, real shape, real scale, real prop) {
    // prop: cure probability
    real log_lik_cen   = (1 - is_censored) *
                         log_sum_exp(log(prop), log1m(prop) + weibull_lccdf(t | shape, scale));
    real log_lik_uncen = is_censored * log1m(prop) + is_censored * weibull_lpdf(t | shape, scale);
    return log_lik_cen + log_lik_uncen;
  }
}

data {
  int N;
  vector<lower=0>[N] t;
  array[N] int<lower=0,upper=1> is_censored;

  int P; int M; int K;
  matrix[N, P] sc;
  matrix[N, M-P] xonly;
  matrix[N, K-P] zonly;
}

parameters {
  // shared features (P)
  vector[P] betaU_sc;
  vector[P] betaC_sc;

  // non-shared features
  vector[M-P] betaU_o;
  vector[K-P] betaC_o;

  // Weibull: shape = 1/sigma, scale = exp(x*betaU)
  real log_sigma;

  // --- Global mixture weights (shared across features) ---
  // order: [both, U-only, C-only, none]
  simplex[4] pi;

  // scales for "both" case (U, C)
  vector<lower=0.01>[2] tau_both;
  cholesky_factor_corr[2] L_Omega_both;

  // one-sided slabs
  real<lower=0.01> tauU_only;
  real<lower=0.01> tauC_only;

  // intercepts
  real alphaU;
  real alphaC;
}

transformed parameters {
  vector[N] muU   = exp(alphaU + sc * betaU_sc + xonly * betaU_o);           // Weibull scale
  real      sigma = exp(log_sigma);
  vector[N] prop  = inv_logit(alphaC + sc * betaC_sc + zonly * betaC_o);     // cure prob
}

model {
  // --- priors (more conservative for C to reduce FP) ---
  log_sigma ~ normal(log(0.5), 0.5);
  alphaU ~ normal(0, 1);
  alphaC ~ normal(0, 1);

  if (M-P > 0) betaU_o ~ normal(0, 1);
  if (K-P > 0) betaC_o ~ normal(0, 1);

  // global sparsity prior: favor "none" to control FP
  pi ~ dirichlet([1, 1, 1, 4]');       // [both, U-only, C-only, none]

  // correlation & slab scales (C slab slightly tighter than U)
  L_Omega_both ~ lkj_corr_cholesky(2); // η=4: weaker corr → fewer false "both"
  tau_both[1]  ~ student_t(3, 0, 2);   // U
  tau_both[2]  ~ student_t(3, 0, 2.5); // C
  tauU_only    ~ student_t(3, 0, 2);   // U-only slab
  tauC_only    ~ student_t(3, 0, 2.5); // C-only slab

  // likelihood
  for (i in 1:N)
    target += cen_weibull_lpdf(t[i] | is_censored[i], 1/sigma, muU[i], prop[i]);

  // mixture over feature pairs (betaU_sc[p], betaC_sc[p])
  {
    matrix[2,2] L_both = diag_pre_multiply(tau_both, L_Omega_both);
    // distinct spikes: C spike tighter to curb FP, U spike modest
    real eps = 0.01;

    for (p in 1:P) {
      vector[2] bpair;
      vector[4] lp;

      bpair[1] = betaU_sc[p];
      bpair[2] = betaC_sc[p];

      // 1) both nonzero: MVN(0, D * Omega * D)
      lp[1] = multi_normal_cholesky_lpdf(bpair | rep_vector(0, 2), L_both);

      // 2) U-only: U slab, C spike
      lp[2] = normal_lpdf(betaU_sc[p] | 0, tauU_only)
            + normal_lpdf(betaC_sc[p] | 0, eps);

      // 3) C-only: U spike, C slab
      lp[3] = normal_lpdf(betaU_sc[p] | 0, eps)
            + normal_lpdf(betaC_sc[p] | 0, tauC_only);

      // 4) none: both spikes
      lp[4] = normal_lpdf(betaU_sc[p] | 0, eps)
            + normal_lpdf(betaC_sc[p] | 0, eps);

      // global mixture weight (shared across features)
      target += log_sum_exp(log(pi) + lp);
    }
  }
}

generated quantities {
  // posterior feature-wise responsibilities
  vector[P] pr_both;
  vector[P] pr_Uonly;
  vector[P] pr_Conly;
  vector[P] pr_none;
  vector[P] pr_nonzero_U;
  vector[P] pr_nonzero_C;

  {
    matrix[2,2] L_both = diag_pre_multiply(tau_both, L_Omega_both);
    real eps = 0.01;

    for (p in 1:P) {
      vector[2] bpair;
      vector[4] lp;
      vector[4] post;

      bpair[1] = betaU_sc[p];
      bpair[2] = betaC_sc[p];

      lp[1] = multi_normal_cholesky_lpdf(bpair | rep_vector(0, 2), L_both);
      lp[2] = normal_lpdf(betaU_sc[p] | 0, tauU_only)
            + normal_lpdf(betaC_sc[p] | 0, eps);
      lp[3] = normal_lpdf(betaU_sc[p] | 0, eps)
            + normal_lpdf(betaC_sc[p] | 0, tauC_only);
      lp[4] = normal_lpdf(betaU_sc[p] | 0, eps)
            + normal_lpdf(betaC_sc[p] | 0, eps);

      post = softmax(log(pi) + lp);

      pr_both[p]      = post[1];
      pr_Uonly[p]     = post[2];
      pr_Conly[p]     = post[3];
      pr_none[p]      = post[4];
      pr_nonzero_U[p] = post[1] + post[2];
      pr_nonzero_C[p] = post[1] + post[3];
    }
  }
}
