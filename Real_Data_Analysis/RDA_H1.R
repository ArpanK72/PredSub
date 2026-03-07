H1_ASE     <- vector("list", B)
H1_PredSub <- vector("list", length(alist))
H1_PureSub <- vector("list", length(alist))
H1_estSS   <- vector("list", length(olist))

for (i in seq_along(H1_PredSub))  H1_PredSub[[i]]  <- vector("list", B)
for (i in seq_along(H1_PureSub))  H1_PureSub[[i]]  <- vector("list", B)
for (i in seq_along(H1_estSS))    H1_estSS[[i]]    <- vector("list", B)

for (iter in 1:B) {
  cat(sprintf("Bootstrap iteration %d / %d\n", iter, B))
  
  # ---- ASE ----
  st <- Sys.time()
  A1_boot <- sparsegraph_under_H0_rcpp(ASE_H0$Xhat1, ASE_H0$Xhat2, blocks = 100, d)
  A2_boot <- sparsegraph_under_H0_rcpp(ASE_H0$Xhat1, ASE_H0$Xhat2, blocks = 100, d)
  
  Xhat1_new <- ASE(A1_boot, d, p_known = TRUE, d)$Xhat; rm(A1_boot)
  Xhat2_new <- ASE(A2_boot, d, p_known = TRUE, d)$Xhat; rm(A2_boot)
  dist      <- matrix_distance(Xhat1_new, I_pq, Xhat2_new, I_pq, rho_n = 1, type = "absolute")
  t_ASE     <- as.numeric(difftime(Sys.time(), st, units = "secs"))
  rm(Xhat1_new, Xhat2_new); gc()
  
  H1_ASE[[iter]] <- c(dist, t_ASE)
  cat(sprintf("  ASE done | time = %.2f secs\n", t_ASE))
  
  # ---- PredSub ----
  for (i in seq_along(alist)) {
    m     <- PredSub_H0[[i]]$samplesize
    st    <- Sys.time()
    
    S1    <- sort(sample(1:n, size = m, replace = FALSE))
    A1_boot <- sparsegraph_predsub_under_H0_rcpp(PredSub_H0[[i]]$Xhat1, PredSub_H0[[i]]$Xhat2, 100, S = S1, d)
    A.hat <- ASE(A1_boot$B, d, p_known = TRUE, d)$Xhat
    Xhat1_new       <- matrix(0, nrow = n, ncol = d)
    Xhat1_new[S1, ] <- A.hat
    scaling         <- A.hat %*% solve(t(A.hat) %*% A.hat)
    oos             <- setdiff(1:n, S1)
    if (length(oos) > 0) Xhat1_new[oos, ] <- as.matrix(A1_boot$A %*% scaling)
    rm(A1_boot, A.hat, scaling)
    
    S2    <- sort(sample(1:n, size = m, replace = FALSE))
    A2_boot <- sparsegraph_predsub_under_H0_rcpp(PredSub_H0[[i]]$Xhat1, PredSub_H0[[i]]$Xhat2, 100, S = S2, d)
    A.hat <- ASE(A2_boot$B, d, p_known = TRUE, d)$Xhat
    Xhat2_new       <- matrix(0, nrow = n, ncol = d)
    Xhat2_new[S2, ] <- A.hat
    scaling         <- A.hat %*% solve(t(A.hat) %*% A.hat)
    oos             <- setdiff(1:n, S2)
    if (length(oos) > 0) Xhat2_new[oos, ] <- as.matrix(A2_boot$A %*% scaling)
    rm(A2_boot, A.hat, scaling)
    
    dist  <- matrix_distance(Xhat1_new, I_pq, Xhat2_new, I_pq, rho_n = 1, type = "absolute")
    t_PS  <- as.numeric(difftime(Sys.time(), st, units = "secs"))
    rm(Xhat1_new, Xhat2_new); gc()
    
    H1_PredSub[[i]][[iter]] <- c(dist, t_PS)
    cat(sprintf("  PredSub | a = %.3f | m = %d | time = %.2f secs\n", alist[i], m, t_PS))
  }
  
  # ---- PureSub ----
  for (i in seq_along(alist)) {
    m     <- PureSub_H0[[i]]$samplesize
    st    <- Sys.time()
    
    A1_boot <- sparsegraph_under_H0_rcpp(PureSub_H0[[i]]$Xhat1, PureSub_H0[[i]]$Xhat2, blocks = 100, d)
    A2_boot <- sparsegraph_under_H0_rcpp(PureSub_H0[[i]]$Xhat1, PureSub_H0[[i]]$Xhat2, blocks = 100, d)
    
    Xhat1_new <- ASE(A1_boot, d, p_known = TRUE, d)$Xhat; rm(A1_boot)
    Xhat2_new <- ASE(A2_boot, d, p_known = TRUE, d)$Xhat; rm(A2_boot)
    dist      <- matrix_distance(Xhat1_new, I_pq, Xhat2_new, I_pq, rho_n = 1, type = "absolute")
    t_PS      <- as.numeric(difftime(Sys.time(), st, units = "secs"))
    rm(Xhat1_new, Xhat2_new); gc()
    
    H1_PureSub[[i]][[iter]] <- c(dist, t_PS)
    cat(sprintf("  PureSub | b = %.3f | m = %d | time = %.2f secs\n", PureSub_H0[[i]]$a, m, t_PS))
  }
  
  # ---- estSS ----
  for (i in seq_along(olist)) {
    m <- estSS_H0[[i]]$overlapsize
    k <- estSS_H0[[i]]$k
    st <- Sys.time()
    
    A1_boot <- sparsegraph_under_H0_rcpp(estSS_H0[[i]]$Xhat1, estSS_H0[[i]]$Xhat2, blocks = 100, d)
    A2_boot <- sparsegraph_under_H0_rcpp(estSS_H0[[i]]$Xhat1, estSS_H0[[i]]$Xhat2, blocks = 100, d)
    
    B1 <- estSS_sequential(A1_boot, m, k, d, overlap = "random"); rm(A1_boot)
    B2 <- estSS_sequential(A2_boot, m, k, d, overlap = "random"); rm(A2_boot)
    
    Xhat1_new <- B1$Xhat
    Xhat2_new <- B2$Xhat
    dist      <- matrix_distance(Xhat1_new, I_pq, Xhat2_new, I_pq, rho_n = 1, type = "absolute")
    t_SS      <- as.numeric(difftime(Sys.time(), st, units = "secs"))
    
    H1_estSS[[i]][[iter]] <- c(dist, t_SS, B1$loop + B2$loop)
    rm(Xhat1_new, Xhat2_new, B1, B2); gc()
    cat(sprintf("  estSS | a = %.3f | m = %d | time = %.2f secs\n", olist[i], m, t_SS))
  }
  
  cat(sprintf("Bootstrap iteration %d complete\n\n", iter))
}
