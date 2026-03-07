B       <- matrix(0.5, nrow = d, ncol = d)
diag(B) <- runif(d)
X       <- rdirichlet(n, rep(0.5, times = d))