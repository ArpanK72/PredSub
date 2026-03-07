# Real data analysis

This folder contains code to reproduce the real-data analysis in the paper.

## Data
 Two datasets are analysed: the **DBLP Co-Authorship** network and the **Cannes 2013 Twitter Multiplex** network. See `Real_Data_Analysis/Data/README.md`.
 
 In each case, the goal is to test whether two networks on the same node set have the same underlying latent structure.

---

## Repository Structure
```
.
├── RDA_DBLP_Run_All.R        # Entry point for DBLP analysis
├── RDA_Cannes_Run_All.R      # Entry point for Cannes analysis
├── RDA_Functions.R           # All core functions
├── RDA_Data_Gen_DBLP.R       # DBLP data processing
├── RDA_Data_Gen_Cannes.R     # Cannes data processing
├── RDA_H0.R                  # Observed test statistics (shared)
├── RDA_H1.R                  # Bootstrap iterations (shared)
├── RDA_Final_Calc.R          # p-values and final table (shared)
└── README.md
```

---

## How to Run

### DBLP Co-Authorship Network
Simply open and run `RDA_DBLP_Run_All.R`. It sets all parameters and sources the pipeline files in order:
```r
source("RDA_Functions.R")
source("RDA_Data_Gen_DBLP.R")
source("RDA_H0.R")
source("RDA_H1.R")
source("RDA_Final_Calc.R")
View(tab_rda)   # Recreates the table presented in the paper
```

Recommended parameters for the DBLP dataset:
```r
d     <- 10
B     <- 200
alist <- seq(2.75, 3.25, by = 0.25)
olist <- tail(alist, 2)
```

### Cannes 2013 Twitter Multiplex Network
Simply open and run `RDA_Cannes_Run_All.R`. It sets all parameters and sources the pipeline files in order:
```r
source("RDA_Functions.R")
source("RDA_Data_Gen_Cannes.R")
source("RDA_H0.R")
source("RDA_H1.R")
source("RDA_Final_Calc.R")
View(tab_rda)   # Recreates the table presented in the paper
```

Recommended parameters for the Cannes dataset:
```r
d     <- 10
B     <- 200
alist <- seq(2.75, 3.25, by = 0.25)
olist <- tail(alist, 2)
```

The paper presents results for `d = 5, 10, 20`. To replicate each, update `d` in the respective `Run_All` file and re-run the full pipeline. All other parameters remain the same.

All parameters can be freely modified in the respective `Run_All` files to suit different experimental settings.

---

## File Descriptions

### `RDA_DBLP_Run_All.R` and `RDA_Cannes_Run_All.R`
Entry points for the respective analyses. Set all parameters here and source the pipeline files in order.

### `RDA_Functions.R`
Loads all required packages and defines all core functions used in the pipeline.

### DBLP Co-Authorship Network

**`RDA_Data_Gen_DBLP.R`** — Processes the raw DBLP co-authorship hypergraph data. Converts publication records into a co-authorship edge list using a C++ function, filters edges into two time periods (Era 1: $2011--2014$, Era 2: $2015--2018$), restricts to authors active in both eras, removes low-degree nodes (degree < 5), and produces two sparse adjacency matrices `A1` and `A2` on the same node set.

### Cannes 2013 Twitter Multiplex Network

**`RDA_Data_Gen_Cannes.R`** — Processes the raw Cannes 2013 Twitter multiplex edge list. Separates edges into the Retweet (layer 1) and Mention (layer 2) layers, converts to undirected graphs on a common node set, removes isolated nodes, and produces two sparse adjacency matrices, `A1` (Retweet) and `A2` (Mention), on the same node set.

### Shared Pipeline Files

**`RDA_H0.R`** — Computes the observed test statistic $T_0$ and embeddings under $H_0$ for ASE, PredSub, PureSub, and estSS. Results are stored in `ASE_H0`, `PredSub_H0`, `PureSub_H0`, and `estSS_H0`.

**`RDA_H1.R`** — Runs `B` bootstrap iterations. In each iteration generates bootstrap graphs under $H_0$ using the estimated embeddings and computes bootstrap test statistics for each method. Results are accumulated in memory.

**`RDA_Final_Calc.R`** — Pools bootstrap results, computes p-values as the proportion of bootstrap statistics exceeding $T_0$, computes total runtime in minutes, and assembles the final table `tab_rda` with columns `method`, `a`, `m`, `pvalue`, `time`, and `speedup`.

---

## Requirements

R (>= 4.0) with the following packages:
```r
install.packages(c("RSpectra", "Matrix", "doParallel", "foreach",
                   "doFuture", "doRNG", "matrixStats", "dplyr",
                   "tidyr", "igraph", "Rcpp", "RcppEigen", "gtools",
                   "data.table"))
```

---

## Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `d` | Embedding dimension | `10` |
| `B` | Number of bootstrap resamples | `200` |
| `alist` | Subsampling exponents for PredSub; PureSub uses `a + 0.125` | `seq(2.75, 3.25, by = 0.25)` |
| `olist` | Subsampling exponents for estSS; derived from `tail(alist, 2)` | `c(3.00, 3.25)` |

All parameters can be freely modified in the respective `Run_All` files to suit different experimental settings.
