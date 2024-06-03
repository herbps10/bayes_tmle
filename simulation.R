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
