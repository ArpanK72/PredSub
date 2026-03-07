# ---- Pool bootstrap results ----
H1_ASE_mat     <- do.call(rbind, H1_ASE)
H1_PredSub_mat <- lapply(H1_PredSub, function(x) do.call(rbind, x))
H1_PureSub_mat <- lapply(H1_PureSub, function(x) do.call(rbind, x))
H1_estSS_mat   <- lapply(H1_estSS,   function(x) do.call(rbind, x))

# ---- Compute results ----
rows <- list()

# ASE
rows[[1]] <- data.frame(
  method   = "ASE",
  a        = NA_real_,
  m        = n,
  pvalue   = mean(H1_ASE_mat[, 1] > ASE_H0$T0),
  time     = (sum(H1_ASE_mat[, 2]) + as.numeric(ASE_H0$time)) / 60
)

# PredSub
for (i in seq_along(alist)) {
  rows[[1 + i]] <- data.frame(
    method = "PredSub",
    a      = alist[i],
    m      = PredSub_H0[[i]]$samplesize,
    pvalue = mean(H1_PredSub_mat[[i]][, 1] > PredSub_H0[[i]]$T0),
    time   = (sum(H1_PredSub_mat[[i]][, 2]) + as.numeric(PredSub_H0[[i]]$time)) / 60
  )
}

# PureSub
offset <- 1 + length(alist)
for (i in seq_along(alist)) {
  rows[[offset + i]] <- data.frame(
    method = "PureSub",
    a      = PureSub_H0[[i]]$a,
    m      = PureSub_H0[[i]]$samplesize,
    pvalue = mean(H1_PureSub_mat[[i]][, 1] > PureSub_H0[[i]]$T0),
    time   = (sum(H1_PureSub_mat[[i]][, 2]) + as.numeric(PureSub_H0[[i]]$time)) / 60
  )
}

# estSS
offset <- 1 + 2 * length(alist)
for (i in seq_along(olist)) {
  loop_adj <- (k - 1) * (as.numeric(estSS_H0[[i]]$looptime) + sum(H1_estSS_mat[[i]][, 3])) / k
  rows[[offset + i]] <- data.frame(
    method = "estSS",
    a      = olist[i],
    m      = estSS_H0[[i]]$overlapsize,
    pvalue = mean(H1_estSS_mat[[i]][, 1] > estSS_H0[[i]]$T0),
    time   = (sum(H1_estSS_mat[[i]][, 2]) + as.numeric(estSS_H0[[i]]$time) - loop_adj) / 60
  )
}

# ---- Final table ----
tab_rda <- bind_rows(rows) %>%
  mutate(
    speedup = round(time[method == "ASE"] / time, 1),
    speedup = ifelse(method == "ASE", NA_real_, speedup),
    time    = round(time, 2),
    pvalue  = round(pvalue, 3)
  )