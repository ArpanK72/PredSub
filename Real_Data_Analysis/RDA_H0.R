if (nrow(A1) != nrow(A2)) stop("Dimensions of the two networks do not match.")
n <- nrow(A1)

cat("--- Network Summary ---\n")
cat(sprintf("Network 1: %d nodes, %d edges, sparsity = %.6f\n",
            n, sum(A1 == 1) / 2, sparsity(A1)))
cat(sprintf("Network 2: %d nodes, %d edges, sparsity = %.6f\n",
            n, sum(A2 == 1) / 2, sparsity(A2)))
cat("----------------------\n")

I_pq <- diag(d)

# ---- ASE ----
cat(sprintf("Running ASE | d = %d\n", d))
starttime <- Sys.time()
Xhat1     <- ASE(A1, d, p_known = TRUE, d)$Xhat
Xhat2     <- ASE(A2, d, p_known = TRUE, d)$Xhat
T0_ASE    <- matrix_distance(Xhat1, I_pq, Xhat2, I_pq, rho_n = 1, type = "absolute")
time_ASE  <- difftime(Sys.time(), starttime, units = "secs")
ASE_H0    <- list(Xhat1 = Xhat1, Xhat2 = Xhat2, T0 = T0_ASE, time = time_ASE, n = n, d = d)
cat(sprintf("ASE done | T0 = %.4f | time = %.2f secs\n", T0_ASE, as.numeric(time_ASE)))

# ---- PredSub ----
PredSub_H0 <- vector("list", length(alist))
for (i in seq_along(alist)) {
  a          <- alist[i]
  m          <- ceiling(log(n)^(1 + a))
  starttime  <- Sys.time()
  Xhat1      <- PredSub(A1, m, d, p_known = TRUE, d, sampling = "random")$Xhat
  Xhat2      <- PredSub(A2, m, d, p_known = TRUE, d, sampling = "random")$Xhat
  T0         <- matrix_distance(Xhat1, I_pq, Xhat2, I_pq, rho_n = 1, type = "absolute")
  time_taken <- difftime(Sys.time(), starttime, units = "secs")
  PredSub_H0[[i]] <- list(Xhat1 = Xhat1, Xhat2 = Xhat2, T0 = T0, time = time_taken, samplesize = m, a = a)
  cat(sprintf("PredSub done | a = %.3f | m = %d | T0 = %.4f | time = %.2f secs\n",
              a, m, T0, as.numeric(time_taken)))
}

# ---- PureSub ----
PureSub_H0 <- vector("list", length(alist))
for (i in seq_along(alist)) {
  a          <- alist[i]
  b          <- a + 0.125
  m          <- ceiling(log(n)^(1 + b))
  starttime  <- Sys.time()
  S          <- sort(sample(1:n, size = m, replace = FALSE))
  Xhat1      <- ASE(A1[S, S], d, p_known = TRUE, d)$Xhat
  Xhat2      <- ASE(A2[S, S], d, p_known = TRUE, d)$Xhat
  T0         <- matrix_distance(Xhat1, I_pq, Xhat2, I_pq, rho_n = 1, type = "absolute")
  time_taken <- difftime(Sys.time(), starttime, units = "secs")
  PureSub_H0[[i]] <- list(Xhat1 = Xhat1, Xhat2 = Xhat2, T0 = T0, time = time_taken, samplesize = m, a = b)
  cat(sprintf("PureSub done | b = %.3f | m = %d | T0 = %.4f | time = %.2f secs\n",
              b, m, T0, as.numeric(time_taken)))
}

# ---- estSS ----
estSS_H0 <- vector("list", length(olist))
k        <- 10
for (i in seq_along(olist)) {
  a          <- olist[i]
  s          <- ceiling(log(n)^(1 + a))
  m          <- find_o(n, k, s)
  starttime  <- Sys.time()
  C1         <- estSS_sequential(A1, m, k, d, overlap = "random")
  C2         <- estSS_sequential(A2, m, k, d, overlap = "random")
  Xhat1      <- C1$Xhat
  Xhat2      <- C2$Xhat
  T0         <- matrix_distance(Xhat1, I_pq, Xhat2, I_pq, rho_n = 1, type = "absolute")
  time_taken <- difftime(Sys.time(), starttime, units = "secs")
  estSS_H0[[i]] <- list(Xhat1 = Xhat1, Xhat2 = Xhat2, T0 = T0, time = time_taken,
                        looptime = C1$loop + C2$loop, overlapsize = m, k = k, a = a)
  cat(sprintf("estSS done | a = %.3f | m = %d | T0 = %.4f | time = %.2f secs\n",
              a, m, T0, as.numeric(time_taken)))
}