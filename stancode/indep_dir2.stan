// Mixture cure rate model with Weibull AFT and logistic regression (logit link) log-likelihood structure 
// is_censored == 1 event observed 
functions{
  real cen_weibull_lpdf(real t, int is_censored, real shape, real scale, real prop){
                        //row_vector zrow, vector betaC, real alphaC) {
    // real prop = inv_logit(alphaC + dot_product(zrow, betaC)); // cure prob (logistic)
    real log_lik_cen   = (1 - is_censored) *
  log_sum_exp(log(prop), log1m(prop) + weibull_lccdf(t | shape, scale));
    real log_lik_uncen = is_censored * log1m(prop) + is_censored * weibull_lpdf(t | shape, scale);
    return log_lik_cen + log_lik_uncen;
  }
}

data{
  int N;
  vector<lower=0>[N] t;
  array[N] int<lower=0,upper=1> is_censored;
  
  int P; int M; int K;
  matrix[N, P] sc;
  matrix[N, M-P] xonly;
  matrix[N, K-P] zonly; 
} 


parameters{

  vector[P] betaU_sc;
  vector[P] betaC_sc;
  
  // real<lower = 0> gsc;
  // vector[P] betaC_sc_raw; 
  
  vector[M-P] betaU_o;
  vector[K-P] betaC_o; 
  
  real log_sigma;    // Weibull: shape = 1/sigma, scale = exp(x*betaU)
  
  // Independent spike–slab mixing weights for each shared feature
  array[P] simplex[2] piU;           // [slab, spike] for betaU_sc[p]
  array[P] simplex[2] piC;           // [slab, spike] for betaC_sc[p]

  // Slab scales (global, per part)
  real<lower=0.01> tauU;             // slab sd for betaU_sc
  real<lower=0.01> tauC;             // slab sd for betaC_sc
  
  real alphaU;
  real alphaC;
  
  // simplex[4] base;
  // real<lower=0> conc;
  
  
}

transformed parameters{
  vector[N] muU = exp(alphaU + sc * betaU_sc + xonly * betaU_o);  // Weibull scale per obs
  real sigma = exp(log_sigma);
  
  vector[N] prop = inv_logit(alphaC + sc * betaC_sc + zonly * betaC_o);
 
}

model{
  // ---------- Likelihood ----------
  log_sigma ~ normal(log(0.5), 0.5);
  
  alphaU ~ normal(0, 1);
  alphaC ~ normal(0, 1);
  real eps = 0.01;
  
  // nonshared covariates
  if (M-P > 0) betaU_o ~ normal(0, 1);
  if (K-P > 0) betaC_o ~ normal(0, 1);
  
  for(i in 1:N)
    target += cen_weibull_lpdf(t[i] | is_censored[i], 1/sigma, muU[i], prop[i]);
  
  piU ~ dirichlet([1, 2]');
  piC ~ dirichlet([1, 2]');
  // for (p in 1:P) {
  //   piU[p] ~ dirichlet(rep_vector(1.0, 2));
  //   piC[p] ~ dirichlet(rep_vector(1.0, 2));
  // }
  tauU ~ student_t(3, 0, 2);   
  tauC ~ student_t(3, 0, 2.5);
          
  {
    vector[2] lpU;
    vector[2] lpC;
    for (p in 1:P) {
      // U part: slab vs spike
      lpU[1] = normal_lpdf(betaU_sc[p] | 0, tauU);  // slab
      lpU[2] = normal_lpdf(betaU_sc[p] | 0, eps);   // spike
      target += log_sum_exp(log(piU[p]) + lpU);

      // C part: slab vs spike
      lpC[1] = normal_lpdf(betaC_sc[p] | 0, tauC);  // slab
      lpC[2] = normal_lpdf(betaC_sc[p] | 0, eps);   // spike
      target += log_sum_exp(log(piC[p]) + lpC);
    }
  }
}

generated quantities {
  // Posterior responsibilities per feature (independent parts)
  vector[P] prU_slab;
  vector[P] prU_spike;
  vector[P] prC_slab;
  vector[P] prC_spike;

  // (Optional) derived 4-way combos under independence assumption
  vector[P] pr_both;
  vector[P] pr_Uonly;
  vector[P] pr_Conly;
  vector[P] pr_none;
  
  vector[P] pr_nonzero_U;
  vector[P] pr_nonzero_C;
  
  real eps = 0.01;

  {
    vector[2] lpU;
    vector[2] lpC;

    for (p in 1:P) {
      // U responsibilities
      lpU[1] = normal_lpdf(betaU_sc[p] | 0, tauU);
      lpU[2] = normal_lpdf(betaU_sc[p] | 0, eps);
      {
        vector[2] postU = softmax(log(piU[p]) + lpU);
        prU_slab[p]  = postU[1];
        prU_spike[p] = postU[2];
      }

      // C responsibilities
      lpC[1] = normal_lpdf(betaC_sc[p] | 0, tauC);
      lpC[2] = normal_lpdf(betaC_sc[p] | 0, eps);
      {
        vector[2] postC = softmax(log(piC[p]) + lpC);
        prC_slab[p]  = postC[1];
        prC_spike[p] = postC[2];

        // 4-way (independence product)
        pr_both[p]  = prU_slab[p]  * prC_slab[p];
        pr_Uonly[p] = prU_slab[p]  * prC_spike[p];
        pr_Conly[p] = prU_spike[p] * prC_slab[p];
        pr_none[p]  = prU_spike[p] * prC_spike[p];
        
        
        pr_nonzero_U[p] = pr_both[p] + pr_Uonly[p];
        pr_nonzero_C[p] = pr_both[p] + pr_Conly[p];
      }
    }
  }
}