library(doParallel)
library(matrixStats)
library(dplyr)
library(tidyr)
library(RSpectra)
library(Matrix)
library(doFuture)
library(foreach)
library(doRNG)
library(igraph)
library(Rcpp)
library(RcppEigen)
library(gtools)

sourceCpp("UnderH0.cpp")
cat("Cpp Loaded\n")

# Sparsity of a matrix.
# Input  : A -> matrix.
# Output : proportion of non-zero entries.
sparsity <- function(A) {
  sum(A == 1) / (as.numeric(nrow(A)) * as.numeric(ncol(A)))
}

# Procrustes Transformation.
# Input  : X, Y -> matrices to align.
# Output : W -> orthogonal alignment matrix, d -> Frobenius distance after alignment.
procr <- function(X, Y) {
  tmp     <- t(X) %*% Y
  tmp.svd <- svd(tmp)
  W       <- tmp.svd$u %*% t(tmp.svd$v)
  return(list(d = norm(X %*% W - Y, type = "F"), W = W))
}

# Sparse Graph Generation under SBM.
# Input  : X -> membership matrix (n x d), B -> block probability matrix (d x d),
#          rho_n -> sparsity factor (default: 0.01).
# Output : A -> sparse symmetric adjacency matrix (n x n).
sparsegraph_SBM <- function(X, B, blocks = 100, rho_n = 0.01) {
  n    <- nrow(X)
  stor <- vector("list", n)
  
  block_size <- floor(n / blocks)
  remainder  <- n %% blocks
  
  # Compute block start/end indices, distributing remainder across first blocks
  starts <- c(1, cumsum(rep(block_size, blocks) + c(rep(1, remainder), rep(0, blocks - remainder))) + 1)
  ends   <- starts[-1] - 1
  starts <- starts[-(blocks + 1)]
  
  for (j in 1:blocks) {
    a <- starts[j]
    b <- ends[j]
    
    M <- (X[a:b, ] %*% B %*% t(X)) * rho_n
    
    # Last node in the entire matrix has no upper triangle entries
    upper_b <- ifelse(j == blocks, b - 1, b)
    
    for (i in a:upper_b) {
      num_trials <- n - i
      if (num_trials > 0) {
        indices_to_check <- (i + 1):n
        p_vec            <- M[i - a + 1, indices_to_check]
        u_vec            <- runif(num_trials)
        successes        <- u_vec < p_vec
        stor[[i]]        <- i + which(successes)
      }
    }
  }
  
  valid_indices <- which(!sapply(stor, is.null))
  size          <- sapply(stor[valid_indices], length)
  
  vec <- list(
    i = rep(valid_indices, size),
    j = unlist(stor[valid_indices])
  )
  
  A <- sparseMatrix(
    i    = c(vec$i, vec$j),
    j    = c(vec$j, vec$i),
    x    = 1,
    dims = c(n, n)
  )
  
  return(A)
}

# Matrix Distance.
# Input  : X -> (n x d) matrix, B -> (d x d) matrix (default: identity),
#          Y -> (n x d) matrix, C -> (d x d) matrix (default: identity),
#          rho_n -> sparsity scaling factor (default: 1),
#          type -> "relative" (default) or "absolute".
# Output : Frobenius norm distance between rho_n * X B X^T and Y C Y^T.
matrix_distance <- function(X, B = NULL, Y, C = NULL, rho_n = 1, type = "relative") {
  
  if (!is.matrix(X) || !is.matrix(Y))
    stop("X and Y must be matrices.")
  if (!all(dim(X) == dim(Y)))
    stop("X and Y must have the same dimensions (n x d).")
  if (!is.numeric(rho_n) || length(rho_n) != 1)
    stop("rho_n must be a single numeric value.")
  
  d <- ncol(X)
  
  if (is.null(B)) B <- diag(d)
  if (is.null(C)) C <- diag(d)
  
  if (!is.matrix(B) || !is.matrix(C))
    stop("B and C must be matrices.")
  if (!all(dim(B) == c(d, d)) || !all(dim(C) == c(d, d)))
    stop("B and C must be square d x d matrices with d = ncol(X).")
  
  
  Gx  <- crossprod(X)       # X^T X
  Gy  <- crossprod(Y)       # Y^T Y
  Cxy <- crossprod(X, Y)    # X^T Y
  
  Sx        <- B %*% Gx %*% B
  norm_A_sq <- (rho_n^2) * sum(Sx * Gx)
  
  Sy        <- C %*% Gy %*% C
  norm_D_sq <- sum(Sy * Gy)
  
  cross_term <- rho_n * sum(B * (Cxy %*% C %*% t(Cxy)))
  
  tolerance <- .Machine$double.eps * max(1, abs(norm_A_sq))
  
  if (abs(norm_A_sq) < tolerance) {
    warning("||rho * X B X^T||_F is close to zero.")
    return(if (abs(norm_D_sq) < tolerance) 0.0 else Inf)
  }
  
  diff_sq <- norm_A_sq + norm_D_sq - 2 * cross_term
  
  if (diff_sq < 0) {
    if (abs(diff_sq) < tolerance * max(1, abs(norm_A_sq), abs(norm_D_sq))) {
      diff_sq <- 0
      warning("Squared norm difference was slightly negative due to floating point; clamped to 0.")
    } else {
      warning(paste("Squared norm difference is significantly negative (", diff_sq, "). Returning NaN."))
      return(NaN)
    }
  }
  
  if (type == "absolute") sqrt(diff_sq) else sqrt(diff_sq / norm_A_sq)
}


# Adjacency Spectral Embedding (ASE).
# Input  : A -> symmetric adjacency matrix, d -> rank,
#          p_known -> whether p is known (default: TRUE).
#          p -> number of positive eigenvalues if known(default: d),
# Output : estimated latent position X, p -> number of positive eigenvalues, t -> runtime.
ASE <- function(A, d, p_known = TRUE, p = d) {
  st <- Sys.time()
  
  n <- nrow(A)
  if (nrow(A) != ncol(A)) stop("A must be a square matrix.")
  if (d >= n) stop("d must be less than the number of nodes n.")
  
  if (p_known) {
    if (p > d) stop("p must be less than or equal to d.")
    
    A.eig.top <- eigs_sym(A, p, which = "LA", opts = list(tol = 1e-6, maxitr = 5000L))
    
    if (p == d) {
      A.val <- abs(A.eig.top$values)
      A.vec <- A.eig.top$vectors
    } else {
      A.eig.low <- eigs_sym(A, d - p, which = "SA", opts = list(tol = 1e-6, maxitr = 5000L))
      A.val <- abs(c(A.eig.top$values, A.eig.low$values))
      A.vec <- cbind(A.eig.top$vectors, A.eig.low$vectors)
    }
  } else {
    A.eig <- eigs_sym(A, d, which = "LM", opts = list(tol = 1e-6, maxitr = 5000L))
    p     <- sum(A.eig$values > 0)
    
    ord   <- order(A.eig$values > 0, decreasing = TRUE)
    A.val <- abs(A.eig$values[ord])
    A.vec <- A.eig$vectors[, ord]
  }
  
  if(d == 1){
    A.coords <- A.vec[, 1, drop = FALSE] * sqrt(A.val[1])
  } else {
    A.coords <- A.vec %*% diag(sqrt(A.val))
  }
  
  time <- difftime(Sys.time(), st, units = "secs")
  return(list(Xhat = A.coords, p = p, t = time))
}


# Predictive Subsampling (PredSub).
# Input  : A -> symmetric adjacency matrix, m <- subsample size, d -> rank,
#          p_known -> whether p is known (default: TRUE).
#          p -> number of positive eigenvalues if known(default: d),
#          sampling -> sampling strategy (default: "random").
# Output : estimated latent position X, p -> number of positive eigenvalues, t -> runtime.
PredSub <- function(A, m, d, p_known = TRUE, p = d, sampling = "random") {
  st <- Sys.time()
  
  n <- nrow(A)
  if (nrow(A) != ncol(A)) stop("A must be a square matrix.")
  if (d >= n) stop("d must be less than the number of nodes n.")
  if (m >= n) stop("m must be less than the number of nodes n.")
  if (p_known && p > d) stop("p must be less than or equal to d.")
  
  X.fin <- matrix(0, nrow = n, ncol = d)
  
  if (sampling == "sparsity based") {
    nnzAi        <- rowSums(A)
    nnzA         <- sum(nnzAi)
    wt           <- (m * nnzAi) / nnzA
    prob         <- sapply(1:n, function(i) min(1, wt[i]))
    sample_index <- which(runif(n) < prob)
    if (length(sample_index) == 0) stop("Sparsity based sampling produced an empty sample. Consider increasing m.")
    
  } else if (sampling == "degree") {
    sample_index <- order(rowSums(A), decreasing = TRUE)[1:m]
    
  } else {
    sample_index <- sort(sample(1:n, size = m, replace = FALSE))
  }
  
  A.samp <- A[sample_index, sample_index]
  
  if (p_known) {
    A.hat <- ASE(A.samp, d, p_known = TRUE, p)$Xhat
  } else {
    est   <- ASE(A.samp, d, p_known = FALSE)
    A.hat <- est$Xhat
    p     <- est$p
  }
  
  X.fin[sample_index, ] <- A.hat
  
  if (p == d) {
    I_pq <- diag(d)
  } else {
    I_pq <- diag(c(rep(1, times = p), rep(-1, times = d - p)))
  }
  
  scaling <- A.hat %*% solve(t(A.hat) %*% A.hat) %*% I_pq
  oos     <- setdiff(1:n, sample_index)
  
  if (length(oos) > 0) {
    X.fin[oos, ] <- as.matrix(A[oos, sample_index] %*% scaling)
  }
  
  time <- difftime(Sys.time(), st, units = "secs")
  return(list(Xhat = X.fin, p = p, t = time))
}


# Sequential estSS Embedding (Chakraborty et al., 2025).
# Input  : A -> adjacency matrix, m -> overlap size, k -> number of subsamples,
#          d -> embedding dimension, overlap -> "dense" or "random" (default: "dense").
# Output : Xhat -> (n x d) estimated latent positions, t -> total runtime, loop -> loop runtime.
estSS_sequential <- function(A, m, k, d, overlap = "dense") {
  st <- Sys.time()
  n  <- nrow(A)
  b  <- (n - m) / k
  
  S <- rowSums(A)
  
  if (overlap == "dense") {
    common <- tail(order(S), m)
    x      <- head(order(S), n - m)
  } else if (overlap == "random") {
    common <- sample(1:n, size = m, replace = FALSE)
    x      <- setdiff(1:n, common)
  } else {
    stop("overlap must be 'dense' or 'random'.")
  }
  
  samp <- matrix(NA, nrow = k, ncol = b)
  for (i in 1:k) {
    samp[i, ] <- sample(x, size = b, replace = FALSE)
    x         <- setdiff(x, samp[i, ])
  }
  cat("Sampling done\n")
  
  ind1      <- c(common, samp[1, ])
  X.hat.ref <- ASE(A[ind1, ind1], d, p_known = TRUE, d)$Xhat
  X.ref.0   <- X.hat.ref[1:m, ]
  X.ref.1   <- X.hat.ref[(m + 1):(m + b), ]
  
  X.part_list <- vector("list", k - 1)
  timings     <- numeric(k - 1)
  
  for (i in 2:k) {
    t0        <- Sys.time()
    A.sub     <- A[c(common, samp[i, ]), c(common, samp[i, ])]
    X.hat.sub <- ASE(A.sub, d, p_known = TRUE, d)$Xhat
    X.sub.0   <- X.hat.sub[1:m, ]
    X.sub.i   <- X.hat.sub[(m + 1):(m + b), ]
    
    H                  <- procr(X.sub.0, X.ref.0)$W
    X.part_list[[i-1]] <- X.sub.i %*% H
    timings[i-1]       <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }
  
  X.part <- do.call(rbind, X.part_list)
  
  X.fin             <- matrix(0, n, d)
  X.fin[common, ]   <- X.ref.0
  X.fin[samp[1, ],] <- X.ref.1
  for (i in 2:k) {
    X.fin[samp[i, ], ] <- X.part[((i - 2) * b + 1):((i - 1) * b), ]
  }
  
  fulltime <- difftime(Sys.time(), st, units = "secs")
  return(list(Xhat = X.fin, t = fulltime, loop = sum(timings)))
}


# ASE-based Two-Sample Graph Test (Tang et al., 2017).
# Input  : A, A_e -> adjacency matrices, d -> embedding dimension,
#          p -> number of positive eigenvalues, B -> number of bootstrap resamples (default: 200).
# Output : pvalue -> bootstrap p-value, t1 -> runtime.
ASEtest_rcpp_sequential <- function(A, A_e, d, p, B = 200) {
  starttime <- Sys.time()
  
  Xhat1 <- ASE(A, d, p_known = TRUE, p)$Xhat
  Xhat2 <- ASE(A_e, d, p_known = TRUE, p)$Xhat
  
  I_pq <- if (p < d) diag(c(rep(1, p), rep(-1, d - p))) else diag(p)
  
  T0   <- matrix_distance(Xhat1, I_pq, Xhat2, I_pq, rho_n = 1, type = "absolute")
  Tvec <- numeric(B)
  
  for (b in 1:B) {
    s  <- Sys.time()
    A1 <- sparsegraph_under_H0_rcpp(Xhat1, Xhat2, 100, p)
    A2 <- sparsegraph_under_H0_rcpp(Xhat1, Xhat2, 100, p)
    
    Xhat1_new <- ASE(A1, d, p_known = TRUE, p)$Xhat; rm(A1)
    Xhat2_new <- ASE(A2, d, p_known = TRUE, p)$Xhat; rm(A2)
    
    Tvec[b] <- matrix_distance(Xhat1_new, I_pq, Xhat2_new, I_pq, rho_n = 1, type = "absolute")
    rm(Xhat1_new, Xhat2_new); gc()
    cat("ASE Bootstrap:", b, "done in", difftime(Sys.time(), s, units = "secs"), "\n")
  }
  
  pval     <- sum(Tvec > T0) / B
  testtime <- difftime(Sys.time(), starttime, units = "secs")
  return(list(pvalue = pval, t1 = testtime))
}

# PredSub-based Two-Sample Graph Test.
# Input  : A, A_e -> adjacency matrices, m -> subsample size, d -> embedding dimension,
#          p -> number of positive eigenvalues, B -> number of bootstrap resamples (default: 200).
# Output : pvalue -> bootstrap p-value, t1 -> runtime.
PredSubtest_rcpp_sequential <- function(A, A_e, m, d, p, B = 200) {
  n         <- nrow(A)
  starttime <- Sys.time()
  
  Xhat1 <- PredSub(A, m, d, p_known = TRUE, p, sampling = "random")$Xhat
  Xhat2 <- PredSub(A_e, m, d, p_known = TRUE, p, sampling = "random")$Xhat
  
  I_pq <- if (p < d) diag(c(rep(1, p), rep(-1, d - p))) else diag(p)
  
  T0   <- matrix_distance(Xhat1, I_pq, Xhat2, I_pq, rho_n = 1, type = "absolute")
  Tvec <- numeric(B)
  
  for (b in 1:B) {
    s  <- Sys.time()
    S1 <- sort(sample(1:n, size = m, replace = FALSE))
    A1 <- sparsegraph_predsub_under_H0_rcpp(Xhat1, Xhat2, 100, S = S1, p)
    
    A.hat           <- ASE(A1$B, d, p_known = TRUE, p)$Xhat
    Xhat1_new       <- matrix(0, nrow = n, ncol = d)
    Xhat1_new[S1, ] <- A.hat
    scaling         <- A.hat %*% solve(t(A.hat) %*% A.hat) %*% I_pq
    oos             <- setdiff(1:n, S1)
    if (length(oos) > 0) Xhat1_new[oos, ] <- as.matrix(A1$A %*% scaling)
    rm(A1, A.hat, scaling)
    
    S2 <- sort(sample(1:n, size = m, replace = FALSE))
    A2 <- sparsegraph_predsub_under_H0_rcpp(Xhat1, Xhat2, 100, S = S2, p)
    
    A.hat           <- ASE(A2$B, d, p_known = TRUE, p)$Xhat
    Xhat2_new       <- matrix(0, nrow = n, ncol = d)
    Xhat2_new[S2, ] <- A.hat
    scaling         <- A.hat %*% solve(t(A.hat) %*% A.hat) %*% I_pq
    oos             <- setdiff(1:n, S2)
    if (length(oos) > 0) Xhat2_new[oos, ] <- as.matrix(A2$A %*% scaling)
    rm(A2, A.hat, scaling)
    
    Tvec[b] <- matrix_distance(Xhat1_new, I_pq, Xhat2_new, I_pq, rho_n = 1, type = "absolute")
    rm(Xhat1_new, Xhat2_new); gc()
    cat("PredSub Bootstrap:", b, "for m:", m, "done in", difftime(Sys.time(), s, units = "secs"), "\n")
  }
  
  pval     <- sum(Tvec > T0) / B
  testtime <- difftime(Sys.time(), starttime, units = "secs")
  return(list(pvalue = pval, t1 = testtime))
}

# PureSub-based Two-Sample Graph Test.
# Input  : A, A_e -> adjacency matrices, m -> subsample size, d -> embedding dimension,
#          p -> number of positive eigenvalues, B -> number of bootstrap resamples (default: 200).
# Output : pvalue -> bootstrap p-value, t1 -> runtime.
PureSubtest_rcpp_sequential <- function(A, A_e, m, d, p, B = 200) {
  n         <- nrow(A)
  starttime <- Sys.time()
  
  S     <- sort(sample(1:n, size = m, replace = FALSE))
  Xhat1 <- ASE(A[S, S], d, p_known = TRUE, p)$Xhat
  Xhat2 <- ASE(A_e[S, S], d, p_known = TRUE, p)$Xhat
  
  I_pq <- if (p < d) diag(c(rep(1, p), rep(-1, d - p))) else diag(p)
  
  T0   <- matrix_distance(Xhat1, I_pq, Xhat2, I_pq, rho_n = 1, type = "absolute")
  Tvec <- numeric(B)
  
  for (b in 1:B) {
    s  <- Sys.time()
    A1 <- sparsegraph_under_H0_rcpp(Xhat1, Xhat2, 100, p)
    A2 <- sparsegraph_under_H0_rcpp(Xhat1, Xhat2, 100, p)
    
    Xhat1_new <- ASE(A1, d, p_known = TRUE, p)$Xhat; rm(A1)
    Xhat2_new <- ASE(A2, d, p_known = TRUE, p)$Xhat; rm(A2)
    
    Tvec[b] <- matrix_distance(Xhat1_new, I_pq, Xhat2_new, I_pq, rho_n = 1, type = "absolute")
    rm(Xhat1_new, Xhat2_new); gc()
    cat("PureSub Bootstrap:", b, "for m:", m, "done in", difftime(Sys.time(), s, units = "secs"), "\n")
  }
  
  pval     <- sum(Tvec > T0) / B
  testtime <- difftime(Sys.time(), starttime, units = "secs")
  return(list(pvalue = pval, t1 = testtime))
}

# estSS-based Two-Sample Graph Test (Chakraborty et al., 2025).
# Input  : A, A_e -> adjacency matrices, m -> overlap size, k -> number of subsamples,
#          d -> embedding dimension, B -> number of bootstrap resamples (default: 200).
# Output : pvalue -> bootstrap p-value, t1 -> adjusted runtime.
estSStest_rcpp_sequential <- function(A, A_e, m, k, d, B = 200) {
  starttime <- Sys.time()
  
  C1    <- estSS_sequential(A, m, k, d, overlap = "random")
  C2    <- estSS_sequential(A_e, m, k, d, overlap = "random")
  Xhat1 <- C1$Xhat
  Xhat2 <- C2$Xhat
  I_pq  <- diag(d)
  
  T0   <- matrix_distance(Xhat1, I_pq, Xhat2, I_pq, rho_n = 1, type = "absolute")
  Tmat <- matrix(nrow = B, ncol = 2)
  
  for (b in 1:B) {
    s  <- Sys.time()
    A1 <- sparsegraph_under_H0_rcpp(Xhat1, Xhat2, 100, d)
    A2 <- sparsegraph_under_H0_rcpp(Xhat1, Xhat2, 100, d)
    
    B1 <- estSS_sequential(A1, m, k, d, overlap = "random"); rm(A1)
    B2 <- estSS_sequential(A2, m, k, d, overlap = "random"); rm(A2)
    
    Xhat1_new <- B1$Xhat
    Xhat2_new <- B2$Xhat
    
    dist      <- matrix_distance(Xhat1_new, I_pq, Xhat2_new, I_pq, rho_n = 1, type = "absolute")
    rm(Xhat1_new, Xhat2_new); gc()
    Tmat[b, ] <- c(dist, B1$loop + B2$loop)
    cat("estSS Bootstrap:", b, "done in", difftime(Sys.time(), s, units = "secs"), "\n")
  }
  
  pval      <- sum(Tmat[, 1] > T0) / B
  testtime  <- difftime(Sys.time(), starttime, units = "secs")
  finaltime <- testtime - (k - 1) * (C1$loop + C2$loop) / k - (k - 1) * sum(Tmat[, 2]) / k
  return(list(pvalue = pval, t1 = finaltime))
}

cat("Functions Loaded\n")
