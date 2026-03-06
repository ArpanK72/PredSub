summary_result <- combined_result %>%
  mutate(
    method = ifelse(method == 0, "ASE", "PredSub"),
    a      = ifelse(method == "ASE", NA_real_, a)
  ) %>%
  group_by(d, rho, method, a) %>%
  summarise(
    time_mean = round(mean(time), 1),
    time_sd   = round(sd(time),   2),
    error     = round(mean(error), 2),
    error_sd  = round(sd(error),   2),
    .groups   = "drop"
  )

ase_ref <- summary_result %>%
  filter(method == "ASE") %>%
  select(d, rho, ase_time = time_mean, ase_error = error)

tab <- summary_result %>%
  left_join(ase_ref, by = c("d", "rho")) %>%
  mutate(
    speedup   = round(ase_time  / time_mean, 2),
    acc_ratio = round(error / ase_error, 2),
    speedup   = ifelse(method == "ASE", NA_real_, speedup),
    acc_ratio = ifelse(method == "ASE", NA_real_, acc_ratio)
  ) %>%
  select(d, rho, method, a, time_mean, speedup, error, error_sd, acc_ratio) %>%
  arrange(d, rho, method == "PredSub", a)

View(tab)