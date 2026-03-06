BX_list <- lapply(dlist, function(d) {
  if(model == "RDPG"){
    B       <- matrix(0, nrow = d, ncol = d)
  } else{
    B       <- matrix(0.5, nrow = d, ncol = d)
  }
  diag(B) <- runif(d)
  X       <- rdirichlet(n, rep(0.5, times = d))
  list(B = B, X = X)
})
names(BX_list) <- dlist