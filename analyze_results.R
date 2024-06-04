library(tidyverse)
source("simulation.R")
source("bayes_tmle_hierarchical.R")

simulation_results <- read_rds("results/simulation_results.rds")

results <- simulation_results %>%
  select(-res, -naive_ate) %>%
  unnest(c(metrics)) %>%
  group_by(N, G) %>%
  summarize_at(vars(ends_with("mae"), ends_with("coverage")), mean)

tab_mae <- results %>%
  select(ends_with("mae")) %>%
  mutate_at(vars(ends_with("mae")), scales::number_format(accuracy = 0.001))

tab_coverage <- results %>%
  select(ends_with("coverage")) %>%
  mutate_at(vars(ends_with("coverage")), scales::percent_format(accuracy = 0.01))
