library(tidyverse)
source("simulation.R")

summary_bayes_tmle_hierarchical <- function(fit) {
  freq <- fit$group_index %>% select(g, fit) %>%
    mutate(freq_psi  = map_dbl(fit, \(fit) fit$estimates$ATE$psi),
           freq_low  = map_dbl(fit, \(fit) fit$estimates$ATE$CI[1]), 
           freq_high = map_dbl(fit, \(fit) fit$estimates$ATE$CI[2])) %>%
    select(-fit)
  
  bayes_hierarchical <- fit$bayes_hierarchical$draws("psi") %>%
    tidybayes::spread_draws(psi[g]) %>%
    tidybayes::median_qi() %>%
    left_join(fit$group_index %>% select(g, n, value))
  
  bayes_nonhierarchical <- fit$bayes_nonhierarchical$draws("psi") %>%
    tidybayes::spread_draws(psi[g]) %>%
    tidybayes::median_qi() %>%
    left_join(fit$group_index %>% select(g, n, value))
  
  left_join(bayes_hierarchical, freq, by = "g") %>%
    left_join(bayes_nonhierarchical, by = "g", suffix = c(".hierarchical", ".nonhierarchical")) %>%
    mutate(nuisance_method = fit$nuisance_method)
}

metrics <- function(fit, naive_ates, true_ates) {
  summary_bayes_tmle_hierarchical(fit) %>%
    left_join(true_ates, by = c(g = "group")) %>%
    left_join(naive_ates, by = c(g = "group")) %>%
    mutate(bayes_hierarchical_error = ate - psi.hierarchical, freq_error = ate - freq_psi, naive_error = ate - naive_ate,
           bayes_nonhierarchical_error = ate - psi.nonhierarchical,
           bayes_hierarchical_covered = .lower.hierarchical <= ate & .upper.hierarchical >= ate, freq_covered = freq_low <= ate & freq_high >= ate,
           bayes_nonhierarchical_covered = .lower.nonhierarchical <= ate & .upper.nonhierarchical >= ate) %>%
    summarize(bayes_hierarchical_mae = mean(abs(bayes_hierarchical_error)), freq_mae = mean(abs(freq_error)), naive_mae = mean(abs(naive_error)),
              bayes_nonhierarchical_mae = mean(abs(bayes_nonhierarchical_error)),
              bayes_nonhierarchical_coverage = mean(bayes_nonhierarchical_covered),
              bayes_hierarchical_coverage = mean(bayes_hierarchical_covered), freq_coverage = mean(freq_covered))
}

simulation_results <- read_rds("results/simulation_results.rds")

true_ates <- simulate(5e6, G = 20, seed = 5) %>%
  arrange(group) %>%
  group_by(group, group_effect) %>%
  nest() %>%
  mutate(ate = map_dbl(data, \(data) mean(data$Y1 - data$Y0))) %>%
  select(-data)

results <- simulation_results %>%
  mutate(metrics = map2(res, naive_ate, metrics)) %>%
  select(-res, -naive_ate) %>%
  unnest(c(metrics)) %>%
  group_by(N, G) %>%
  summarize_at(vars(ends_with("mae"), ends_with("coverage")), mean)
