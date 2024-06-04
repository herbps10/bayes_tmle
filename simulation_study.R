library(tidyverse)
library(cmdstanr)
library(tmle)
source("simulation.R")
source("bayes_tmle_hierarchical.R")

num_simulations <- 1
Q.SL.library <- c("SL.knn", "SL.glm", "SL.mean", "SL.glm.interaction", "SL.glmnet")
g.SL.library <- c("SL.knn", "SL.glm", "SL.mean", "SL.glm.interaction", "SL.glmnet")
parallel_chains <- 4
adapt_delta <- 0.99

simulations <- expand_grid(
  simulation_index = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")),
  #d = c(2, 5, 10),
  d = c(5),
  #N = c(2.5e3, 5e3),
  N = c(1e3),
  G = c(20)
) %>%
  mutate(seed = N * 1e5 + 1e4 * G + 1e3 * d + simulation_index)

for(N in unique(simulations$N)) {
  dir <- glue::glue("/gpfs/scratch/susmah01/bayes_tmle/cache/{N}")
  if(!dir.exists(dir)) dir.create(dir)
}

simulations <- simulations %>%
  mutate(
    data = pmap(list(N, G, d, seed), simulate),
    naive_ate = map(data, \(data) {
      data %>%
        group_by(group) %>%
        nest() %>%
        mutate(naive_ate = map_dbl(data, \(data) mean(data$Y[data$A == 1]) - mean(data$Y[data$A == 0])))
    }),
    truth = pmap(list(simulation_index, G, seed), function(simulation_index, G, seed) {
      true_ates <- simulate(5e6, G = G, seed = seed) %>%
        arrange(group) %>%
	group_by(group, group_effect) %>%
	nest() %>%
	mutate(ate = map_dbl(data, \(data) mean(data$Y1 - data$Y0))) %>%
	select(-data)
    }),
    fit = pmap(list(simulation_index, N, G, d, data, naive_ate, truth), \(simulation_index, N, G, d, data, naive_ate, truth) {
      cache <- glue::glue("/gpfs/scratch/susmah01/bayes_tmle/cache/{N}/{d}_{G}_{simulation_index}.rds")
      if(file.exists(cache)) {
        return(read_rds(cache))
      }

      covars <- names(data)[str_detect(names(data), "W")]

      res <- bayes_tmle_hierarchical(
	data$Y, 
	data$A, 
	data[, covars, drop = FALSE], 
	data$group, 
	nuisance_method = "separate", 
	Q.SL.library = Q.SL.library, 
	g.SL.library = g.SL.library, 
	parallel_chains = parallel_chains, 
	adapt_delta = adapt_delta, 
	output_dir = "/gpfs/scratch/susmah01/bayes_tmle/fits/"
      )

      metrics <- bayes_tmle_hierarchical_metrics(res, naive_ate, truth)

      write_rds(tibble(
	simulation_index = simulation_index, 
        N                = N, 
	G                = G, 
	d                = d,
	res              = list(res), 
	metrics          = list(metrics), 
	naive_ate        = list(naive_ate)
     ), cache)

      return(res)
    })
  )
