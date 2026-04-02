# css_variable_selection

This repository contains simulation code for Bayesian variable selection in **mixture cure rate (MCR) models**, comparing two spike-and-slab prior specifications: correlated and independent.

---

## Repository Structure
```
css_variable_selection/
├── stancode/
│   ├── css_dir1.stan
│   ├── css_dir2.stan
│   ├── css_dir4.stan
│   ├── css_dir8.stan
│   ├── indep_dir1.stan
│   ├── indep_dir2.stan
│   ├── indep_dir4.stan
│   └── indep_dir8.stan
└── vscode/
    ├── dgp.R
    ├── simulate_and_helpers.R
    ├── sim_weibull(HPC).R
    └── sim_weibull(local).R
```

---

## Folder Descriptions

### `stancode/`

Contains Stan model files implementing Bayesian variable selection for the mixture cure rate model via spike-and-slab priors. Each prior specification is fit under four Dirichlet prior settings (`dir1`, `dir2`, `dir4`, `dir8`) for sensitivity analysis.

**Correlated spike-and-slab (`css_*.stan`)** — Stan models using a correlated spike-and-slab prior, applied jointly to the latency and incidence components of the MCR model.

| File | Dirichlet Prior |
|------|----------------|
| `css_dir1.stan` | Dirichlet(1) |
| `css_dir2.stan` | Dirichlet(2) |
| `css_dir4.stan` | Dirichlet(4) |
| `css_dir8.stan` | Dirichlet(8) |

**Independent spike-and-slab (`indep_*.stan`)** — Stan models using an independent spike-and-slab prior, where inclusion indicators for each variable are treated as independent across components.

| File | Dirichlet Prior |
|------|----------------|
| `indep_dir1.stan` | Dirichlet(1) |
| `indep_dir2.stan` | Dirichlet(2) |
| `indep_dir4.stan` | Dirichlet(4) |
| `indep_dir8.stan` | Dirichlet(8) |

The `dir1`, `dir2`, `dir4`, and `dir8` suffixes correspond to the concentration parameter of the Dirichlet prior placed on the mixture weights, with larger values representing a more diffuse (uniform) prior. These variants are used to assess the sensitivity of posterior inclusion probabilities (PIPs) to the choice of prior.

### `vscode/`

Contains R scripts for generating data, running simulations, and summarizing results.

- **`dgp.R`** — Data generating process. Simulates survival data under the mixture cure rate model for use in the simulation study.
- **`simulate_and_helpers.R`** — Helper functions for running simulations and printing posterior inclusion probabilities (PIPs) for both the latency component (survival among susceptibles) and the incidence component (cure probability).
- **`sim_weibull(HPC).R`** — Main simulation script designed for high-performance computing (HPC) environments. Use this for large-scale replication studies.
- **`sim_weibull(local).R`** — Equivalent simulation script configured for local execution. Note: this script is computationally intensive and time-consuming; it is recommended only for small-scale testing or verification.

---

## Getting Started

1. Install [R](https://www.r-project.org/) and [Stan](https://mc-stan.org/) (via `rstan` or `cmdstanr`).
2. Clone the repository:
```bash
   git clone https://github.com/chaenni95/css_variable_selection.git
   cd css_variable_selection
```
3. Generate data using `dgp.R`, then run a simulation script:
   - For local testing: `sim_weibull(local).R`
   - For full-scale runs: submit `sim_weibull(HPC).R` to your HPC scheduler.
4. Use `simulate_and_helpers.R` to extract and print PIPs for both model components.

---

## Notes

- All Stan models assume a Weibull baseline hazard for the latency component.
- PIPs close to 1 indicate strong evidence for variable inclusion; PIPs close to 0 suggest the variable is not influential.
- The four Dirichlet prior settings (`dir1`, `dir2`, `dir4`, `dir8`) are intended for sensitivity analysis and can be compared across both prior specifications.
- For HPC users, adjust the number of replicates and parallelization settings within `sim_weibull(HPC).R` to match your cluster configuration.

---

## Author

**Chaeyeon Yoo**  
Ph.D. Candidate, Department of Statistics, University of Connecticut  
chaeyeon.yoo@uconn.edu
