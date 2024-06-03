source("bayes_tmle_hierarchical.R")

simulate <- function(N, G = 5, seed = 5) {
  set.seed(seed)
  group_effects <- rnorm(G, mean = 0, sd = 1)
  
  tibble(
    W = rnorm(N),
    A = rbinom(W, 1, 0.5),
    group = sample(G, size = N, replace = TRUE, prob = rep(c(1, 2), length.out = G)),
    group_effect = group_effects[group],
    Y0 = rbinom(N, 1, plogis(W)),
    Y1 = rbinom(N, 1, plogis(W + 1 + group_effect)),
    Y = ifelse(A == 0, Y0, Y1)
  )
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

num_simulations <- 5
Q.SL.library <- c("SL.glm", "SL.mean")
g.SL.library <- c("SL.glm", "SL.mean")
parallel_chains <- 4
adapt_delta <- 0.999

true_ates <- simulate(5e6, G = 20, seed = 5) %>%
  arrange(group) %>%
  group_by(group, group_effect) %>%
  nest() %>%
  mutate(ate = map_dbl(data, \(data) mean(data$Y1 - data$Y0))) %>%
  select(-data)

simulations <- expand_grid(
  simulation_index = 1:num_simulations,
  N = c(5e2),
  G = c(10)
) %>%
  mutate(
    data = pmap(list(N, G, simulation_index), simulate),
    fit = map(data, \(data) {
      bayes_tmle_hierarchical(data$Y, data$A, data[, c("W"), drop = FALSE], data$group, nuisance_method = "combined", Q.SL.library = Q.SL.library, g.SL.library = g.SL.library, parallel_chains = parallel_chains, adapt_delta = adapt_delta)
    }),
    naive_ate = map(data, \(data) {
      data %>%
        group_by(group) %>%
        nest() %>%
        mutate(naive_ate = map_dbl(data, \(data) mean(data$Y[data$A == 1]) - mean(data$Y[data$A == 0])))
    }))

results <- simulations %>%
  mutate(metrics = map2(fit, naive_ate, metrics, true_ates = true_ates)) %>%
  select(-data, -fit, -naive_ate) %>%
  unnest(metrics)

results %>%
  group_by(N, G) %>%
  summarize_at(vars(ends_with("mae"), ends_with("coverage")), mean)
