simulate <- function(N, G = 5, d = 5, seed = 5) {
  set.seed(seed)
  group_effects <- runif(G, -0.5, 0.5)

  W <- matrix(rnorm(N * d), ncol = d, nrow = N)

  res <- as_tibble(W)
  names(res) <- paste0("W", 1:d)

  index <- 1
  repeat {
    res <- res %>% mutate(
      A = rbinom(N, 1, plogis(0.25 + 0.25 * W1 + 0.25 * W2)),
      group = sample(G, size = N, replace = TRUE, prob = rep(c(1, 2), length.out = G)),
      group_effect = group_effects[group],
      Y0 = rbinom(N, 1, plogis(-0.25 + 0.5 * W1 + 0.25 * W2)),
      Y1 = rbinom(N, 1, plogis(-0.25 + 0.5 * W1 + 0.25 * W2 + 1 + group_effect)),
      Y = ifelse(A == 0, Y0, Y1)
    )

    counts <- res %>%
      group_by(group) %>%
      summarize(cases = sum(A), controls = sum(1 - A))

    if(all(counts$cases > 10) && all(counts$controls > 10)) break
    index <- index + 1
    if(index > 100) stop("Could not generate dataset")
  }

  return(res)
}
