all_results <- vector("list", r)

for (iter in 1:r) {
  result <- matrix(nrow = length(rholist) * (length(alist) + 1) * length(dlist), ncol = 7)
  j      <- 1
  
  for (d in dlist) {
    B <- BX_list[[as.character(d)]]$B
    X <- BX_list[[as.character(d)]]$X
    
    for (rho in rholist) {
      A <- sparsegraph_SBM(X, B, blocks = 100, rho_n = rho)
      
      p    <- sum(eigen(B)$values > 0)
      
      # ASE
      est1          <- ASE(A, d, p_known = TRUE, p)
      p = est1$p
      if(p==d){
        I_pq = diag(d)
      } else{
        I_pq <- diag(c(rep(1, times = p), rep(-1, times = d - p)))
      }
      mat_dist1     <- matrix_distance(X, B, est1$Xhat, I_pq, rho_n = rho)
      result[j, ]   <- c(2, est1$t, mat_dist1, d, 0, rho, n)
      j             <- j + 1
      cat(sprintf("Iter %d | ASE done | rho = %.2f | d = %d\n", iter, rho, d))
      
      # PredSub
      for (k in 1:length(alist)) {
        m             <- ceiling(log(n)^(1 + alist[k]))
        est2          <- PredSub(A, m, d, p_known = TRUE, p)
        p = est2$p
        if(p==d){
          I_pq = diag(d)
        } else{
          I_pq <- diag(c(rep(1, times = p), rep(-1, times = d - p)))
        }
        mat_dist2     <- matrix_distance(X, B, est2$Xhat, I_pq, rho_n = rho)
        result[j, ]   <- c(alist[k], est2$t, mat_dist2, d, 1, rho, m)
        j             <- j + 1
        cat(sprintf("Iter %d | PredSub done | a = %.2f | rho = %.2f | d = %d\n", iter, alist[k], rho, d))
      }
      
      rm(A)
    }
  }
  
  colnames(result)    <- c("a", "time", "error", "d", "method", "rho", "m")
  all_results[[iter]] <- as.data.frame(result)
}

combined_result <- bind_rows(all_results) %>% na.omit()
