# Estimation

This repository contains the R code to reproduce the simulation results for the estimation experiments in the paper. The pipeline consists of four source files and a short driver script.

---

## Repository Structure
```
.
├── Run_All.R        # Entry point — set parameters and run full pipeline
├── Functions.R      # All core functions
├── Generation.R     # Generation of true probability matrices
├── Iterations.R     # 1000 Monte Carlo iterations
├── Tables.R         # Summarises results into final table
└── README.md
```

---

## How to Run

Simply open and run `Run_All.R`. It sets all parameters and sources the four files in order:
```r
source("Functions.R")
source("Generation.R")
source("Iterations.R")
source("Tables.R")
View(tab)   # Recreates the table presented in the paper
```

The paper presents results for three network sizes. To replicate each, set the following parameters in `Run_All.R` before running:

**`n = 60000`**
```r
model   <- "GRDPG"
n       <- 60000
alist   <- seq(2.625, 3.125, by = 0.25)
dlist   <- c(5, 10, 20)
rholist <- c(0.01, 0.04)
```

**`n = 100000`**
```r
model   <- "GRDPG"
n       <- 100000
alist   <- seq(2.75, 3.25, by = 0.25)
dlist   <- c(5, 10, 20)
rholist <- c(0.01, 0.04)
```

**`n = 150000`**
```r
model   <- "GRDPG"
n       <- 150000
alist   <- seq(2.875, 3.375, by = 0.25)
dlist   <- c(5, 10, 20)
rholist <- c(0.005, 0.01)
```

All parameters can be freely modified in `Run_All.R` to suit different experimental settings.

---

## File Descriptions

### `Run_All.R`
The entry point for the full pipeline. Sets all parameters (`model`, `n`, `alist`, `dlist`, `rholist`) and sources the four files in order.

### `Functions.R`
Loads all required packages and defines all core functions used in the simulation, including graph generation, spectral embedding, predictive subsampling, and distance computation.

### `Generation.R`
Generates the true latent positions `X` and block probability matrix `B` for each value of `d` in `dlist`, and stores them in `BX_list`. The underlying probability matrix is given by `P = rho_n * X B X^T`.

`B` and `X` are generated once and held in memory for all 1000 iterations.

### `Iterations.R`
Runs 1000 independent simulation iterations. In each iteration:

1. For each `d` in `dlist`, reads `B` and `X` from `BX_list`.
2. For each `rho` in `rholist`, computes `p = #{positive eigenvalues of B}` and generates a sparse adjacency matrix `A` via `sparsegraph_SBM`.
3. Runs **ASE** on the full graph and records runtime and error.
4. Runs **PredSub** for each `a` in `alist` with `m = ceiling(log(n)^{1+a})` and records runtime and error.

Results from all iterations are pooled into `combined_result`.

### `Tables.R`
Summarises `combined_result` into the final table `tab`:

- Averages `time` and `error` across 1000 iterations for each `(d, rho, method, a)` combination.
- Computes `speedup = ASE_time / PredSub_time` and `acc_ratio = PredSub_error / ASE_error` for each PredSub configuration.
- Computes standard deviations for time and error.

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
| `model`   | Underlying probability model | `"GRDPG"` |
| `n`       | Network size | `100000` |
| `alist`   | Subsampling exponents `a` | `seq(2.75, 3.25, by = 0.25)` |
| `dlist`   | Embedding dimensions | `c(5, 10, 20)` |
| `rholist` | Sparsity factors | `c(0.01, 0.04)` |
