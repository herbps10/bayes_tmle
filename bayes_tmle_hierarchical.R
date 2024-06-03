
bayes_tmle_hierarchical <- function(Y, A, W, group, nuisance_method = "combined", Q.SL.library = c("SL.glm", "SL.mean"), g.SL.library = c("SL.glm", "SL.mean"), ...) {
  o <- order(group)
  group <- group[o]
  Y <- Y[o]
  A <- A[o]
  W <- W[o,]
  
  model <- cmdstanr::cmdstan_model("./binomial_tmle_hierarchical.stan")
  
  group_index <- tibble(
    value = group
  )  %>%
    count(value) %>%
    arrange(value) %>%
    mutate(
      g = 1:n(),
      start = map_int(value, \(g) first(which(group == g))),
      end   = map_int(value, \(g) last(which(group == g)))
    )
  
  overall_fit <- tmle(Y, A, cbind(W, factor(group)), Q.SL.library = Q.SL.library, g.SL.library = g.SL.library, family = "binomial")
  
  group_index <- group_index %>% mutate(
    fit   = map2(start, end, function(start, end) {
      if(nuisance_method == "combined") {
        fit <- tmle(Y[start:end], A[start:end], W[start:end,], Q = overall_fit$Qinit$Q[start:end, ], g1W = overall_fit$g$g1W[start:end], family = "binomial")
      }
      else {
        fit <- tmle(Y[start:end], A[start:end], W[start:end,], Q.SL.library = Q.SL.library, g.SL.library = g.SL.library, family = "binomial")
      }
    }),
    Qbar  = pmap(list(fit, start, end), \(fit, start, end) ifelse(A[start:end] == 0, fit$Qinit$Q[, 1], fit$Qinit$Q[, 2])),
    Qbar0 = map(fit, \(fit) fit$Qinit$Q[, 1]),
    Qbar1 = map(fit, \(fit) fit$Qinit$Q[, 2]),
    H1    = pmap(list(fit, start, end), \(fit, start, end) ifelse(A[start:end] == 0, -1/(1 - fit$g$g1W), 1 / fit$g$g1W)),
    H1_0  = map(fit, \(fit) -1 / (1 - fit$g$g1W)),
    H1_1  = map(fit, \(fit) 1 / fit$g$g1W),
    H2    = map(fit, \(fit) fit$Qinit$Q[, 2] - fit$Qinit$Q[, 1] - mean(fit$Qinit$Q[, 2] - fit$Qinit$Q[, 1]))
  )
  
  Qbar  <- unlist(group_index$Qbar)
  Qbar0 <- unlist(group_index$Qbar0)
  Qbar1 <- unlist(group_index$Qbar1)
  H1    <- unlist(group_index$H1)
  H1_0  <- unlist(group_index$H1_0)
  H1_1  <- unlist(group_index$H1_1)
  H2    <- unlist(group_index$H2)
  
  stan_data <- list(
    N     = length(Y),
    G     = nrow(group_index),
    N_group = group_index$n,
    group = group_index$g[base::match(group, group_index$value)],
    group_start = group_index$start,
    group_end = group_index$end,
    y     = Y,
    Qbar  = Qbar,
    Qbar0 = Qbar0,
    Qbar1 = Qbar1,
    H1    = H1,
    H1_0  = H1_0,
    H1_1  = H1_1,
    H2    = H2,
    hierarchical = TRUE
  )
  
  stan_fit_hierarchical <- model$sample(
    data = stan_data,
    ...
  )
  
  stan_data$hierarchical = FALSE
  stan_fit_nonhierarchical <- model$sample(
    data = stan_data,
    ...
  )
  
  list(
    nuisance_method = nuisance_method,
    group_index = group_index,
    bayes_hierarchical = stan_fit_hierarchical,
    bayes_nonhierarchical = stan_fit_nonhierarchical,
    overall_fit = overall_fit
  )
}