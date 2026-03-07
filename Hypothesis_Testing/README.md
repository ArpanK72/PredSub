# Hypothesis testing

This folder contains code to reproduce the hypothesis testing experiments (level/power/p-values) in the paper.

## Repository Structure
```
Hypothesis_Testing/
├── Testing_Run_All.R       # Entry point: set parameters and run pipeline
├── Testing_Functions.R     # Core functions
├── Testing_Generation.R    # Generate B and X
├── Testing_Iterations.R    # Monte Carlo simulation loop
└── Testing_Tables.R        # Summarise results into table
```

## How to Run
Simply open and run `Testing_Run_All.R`. It sets all parameters and sources the four files in order:
```r
source("Testing_Functions.R")
source("Testing_Generation.R")
source("Testing_Iterations.R")
source("Testing_Tables.R")
View(tab)   # Recreates the table presented in the paper
```
The paper presents results for three network sizes. To replicate each, set the following parameters in `Testing_Run_All.R` before running:

**`n = 20000`**
```r
n        <- 20000
d        <- 10
rho_n    <- 0.01
r        <- 200
n_boot   <- 100
alist    <- seq(2.750, 3.000, by = 0.125)
epsilons <- c(0, 0.02, 0.04, 0.08, 0.16)
```
**`n = 30000`**
```r
n        <- 30000
d        <- 10
rho_n    <- 0.01
r        <- 200
n_boot   <- 100
alist    <- seq(2.875, 3.125, by = 0.125)
epsilons <- c(0, 0.02, 0.04, 0.08, 0.16)
```
**`n = 50000`**
```r
n        <- 50000
d        <- 10
rho_n    <- 0.01
r        <- 200
n_boot   <- 100
alist    <- seq(3.000, 3.250, by = 0.125)
epsilons <- c(0, 0.02, 0.04, 0.08, 0.16)
```
PureSub uses `b = a + 0.125` as the subsampling exponent, derived automatically from `alist` in `Testing_Iterations.R`.

All parameters can be freely modified in `Testing_Run_All.R` to suit different experimental settings.

---

## File Descriptions

### `Testing_Run_All.R`
The entry point for the full pipeline. Sets all parameters and sources the four files in order.

### `Testing_Functions.R`
Loads all required packages and defines all core functions used in the simulation, including graph generation, spectral embedding, predictive subsampling, distance computation, and bootstrap-based hypothesis tests.

### `Testing_Generation.R`
Generates the latent position matrix `X` via Dirichlet sampling and the block probability matrix `B` once, held in memory for all `r` iterations.

### `Testing_Iterations.R`
Runs `r` independent Monte Carlo iterations. In each iteration:

1. Generates a sparse adjacency matrix `A` under $H_0$ from $P = \rho_n X B \Pi^\top$.
2. For each `epsilon` in `epsilons`, generates `A_e` under $H_1$ from $P_\epsilon = \rho_n X (B + \epsilon J) X^\top$.
3. Runs **ASE test** on `(A, A_e)` and records p-value and runtime.
4. Runs **PredSub test** for each `a` in `alist` with `m = ceiling(log(n)^{1+a})` and records p-value and runtime.
5. Runs **PureSub test** for each `b = a + 0.125` and records p-value and runtime.

Results from all iterations are pooled into `combined_result`.

### `Testing_Tables.R`
Summarises `combined_result` into the final table `tab`:

- Computes level/power as the proportion of rejections at $\alpha = 0.05$ for each `(method, a, epsilon)` combination.
- Averages runtime in minutes across all `epsilon` values for each `(method, a)` combination.
- Computes `speedup = ASE_time / method_time` for each PredSub and PureSub configuration.

---

## Requirements

R (>= 4.0) with the following packages:
```r
install.packages(c("RSpectra", "Matrix", "doParallel", "foreach",
                   "doFuture", "doRNG", "matrixStats", "dplyr",
                   "tidyr", "igraph", "Rcpp", "RcppEigen", "gtools"))
```

---

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `n` | Network size | `30000` |
| `d` | Embedding dimension | `10` |
| `rho_n` | Sparsity scaling factor | `0.01` |
| `r` | Number of Monte Carlo replications | `200` |
| `n_boot` | Number of bootstrap resamples per test | `100` |
| `alist` | Subsampling exponents for PredSub; PureSub uses `a + 0.125` | `seq(2.875, 3.125, by = 0.125)` |
| `epsilons` | Signal strengths; `epsilon = 0` gives level, `epsilon > 0` gives power | `c(0, 0.02, 0.04, 0.08, 0.16)` |

All parameters can be freely modified in `Testing_Run_All.R` to suit different experimental settings.
