summary_result <- combined_result %>%
  mutate(
    method = case_when(
      method == 1 ~ "ASE",
      method == 3 ~ "PredSub",
      method == 4 ~ "PureSub"
    ),
    a    = ifelse(a == -1, NA_real_, a),
    time = time / 60
  ) %>%
  group_by(method, a, epsilon) %>%
  summarise(
    level_power = round(mean(pvalue < 0.05), 2),
    time_mean   = round(mean(time), 2),
    .groups     = "drop"
  )

time_summary <- summary_result %>%
  group_by(method, a) %>%
  summarise(
    time_mean = round(mean(time_mean), 2),
    .groups   = "drop"
  )


ase_time <- time_summary %>%
  filter(method == "ASE") %>%
  pull(time_mean)


tab <- summary_result %>%
  select(method, a, epsilon, level_power) %>%
  pivot_wider(
    names_from   = epsilon,
    values_from  = level_power,
    names_prefix = "eps_"
  ) %>%
  left_join(time_summary, by = c("method", "a"), na_matches = "na") %>%
  mutate(
    speedup = round(ase_time / time_mean, 1),
    speedup = ifelse(method == "ASE", NA_real_, speedup)
  ) %>%
  arrange(method == "ASE", method == "PredSub", a) %>%
  rename_with(
    ~ gsub("eps_", "$\\epsilon$=", .x),
    starts_with("eps_")
  )