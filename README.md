# PredSub

Supplementary codes for the paper **"Predictive Subsampling for Scalable Inference in Networks"** by Kumar, Tang, and Sengupta.

---

## Repository Structure
```
.
├── Estimation/           # Simulation study for latent position estimation
├── Hypothesis_Testing/   # Simulation study for hypothesis testing
├── Real_Data_Analysis/   # Real data analysis
└── README.md
```

---

## Contents

### `Estimation/`
Contains code to reproduce the estimation simulation results in the paper. Compares **ASE** (Adjacency Spectral Embedding) and **PredSub** (Predictive Subsampling) in terms of runtime and accuracy across various network sizes, embedding dimensions, and sparsity factors. See [`Estimation/README.md`](Estimation/README.md) for details.

### `Hypothesis_Testing/`
Contains code to reproduce the hypothesis testing experiments (level, power, and p-values) in the paper. See [`Hypothesis_Testing/README.md`](Hypothesis_Testing/README.md) for details.

### `Real_Data_Analysis/`
Contains code for the real data analysis in the paper. See [`Real_Data_Analysis/README.md`](Real_Data_Analysis/README.md) for details.

---

## Requirements

R (>= 4.0) with the following packages:
```r
install.packages(c("RSpectra", "Matrix", "doParallel", "foreach",
                   "doFuture", "doRNG", "matrixStats", "dplyr",
                   "tidyr", "igraph", "Rcpp", "RcppEigen", "gtools"))
```

---

## License

This repository is licensed under the MIT License. See [`LICENSE`](LICENSE) for details.
