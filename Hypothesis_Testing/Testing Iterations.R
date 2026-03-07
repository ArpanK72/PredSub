all_results <- vector("list", r)

p    <- sum(eigen(B)$values > 0)
I_pq <- if (p < d) diag(c(rep(1, p), rep(-1, d - p))) else diag(p)

for (iter in 1:r) {
  result <- matrix(nrow = (2 * length(alist) + 1) * length(epsilons), ncol = 7)
  j      <- 1
  
  A <- sparsegraph_SBM(X, B, blocks = 100, rho_n = rho_n)
  cat(sprintf("Iter %d | A generated\n", iter))
  
  for (epsilon in epsilons) {
    B_e <- B + epsilon
    A_e <- sparsegraph_SBM(X, B_e, blocks = 100, rho_n = rho_n)
    cat(sprintf("Iter %d | epsilon = %.3f | A_e generated\n", iter, epsilon))
    
    # ASE test
    res1        <- ASEtest_rcpp_sequential(A, A_e, d, p, B = n_boot)
    result[j, ] <- c(n, n, epsilon, 1, res1$pvalue, res1$t1, -1)
    j           <- j + 1
    cat(sprintf("Iter %d | ASE test done | epsilon = %.3f\n", iter, epsilon))
    
    # PredSub and PureSub tests
    for (a in alist) {
      m_pred      <- ceiling(log(n)^(1 + a))
      res3        <- PredSubtest_rcpp_sequential(A, A_e, m_pred, d, p, B = n_boot)
      result[j, ] <- c(n, m_pred, epsilon, 3, res3$pvalue, res3$t1, a)
      j           <- j + 1
      cat(sprintf("Iter %d | PredSub done | a = %.3f | epsilon = %.3f\n", iter, a, epsilon))
      
      b_val       <- a + 0.125
      m_pure      <- ceiling(log(n)^(1 + b_val))
      res4        <- PureSubtest_rcpp_sequential(A, A_e, m_pure, d, p, B = n_boot)
      result[j, ] <- c(n, m_pure, epsilon, 4, res4$pvalue, res4$t1, b_val)
      j           <- j + 1
      cat(sprintf("Iter %d | PureSub done | a = %.3f | epsilon = %.3f\n", iter, b_val, epsilon))
    }
    
    rm(A_e); gc()
  }
  
  rm(A); gc()
  
  colnames(result)    <- c("n", "m", "epsilon", "method", "pvalue", "time", "a")
  all_results[[iter]] <- as.data.frame(result)
  cat(sprintf("Iteration %d complete\n", iter))
}

combined_result <- bind_rows(all_results) %>% na.omit()