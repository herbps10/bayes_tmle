library(tidyverse)
library(cmdstanr)
library(tmle)
source("simulation.R")
source("bayes_tmle_hierarchical.R")

num_simulations <- 1
Q.SL.library <- c("SL.glm", "SL.mean")
g.SL.library <- c("SL.glm", "SL.mean")
parallel_chains <- 4
adapt_delta <- 0.9

simulations <- expand_grid(
  simulation_index = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")),
  #N = c(1e3, 2.5e3, 5e3),
  N = 100,
  G = c(5)
)

for(N in unique(simulations$N)) {
  dir <- glue::glue("/gpfs/scratch/susmah01/bayes_tmle/cache/{N}")
  if(!dir.exists(dir)) dir.create(dir)
}

simulations <- simulations %>%
  mutate(
    data = pmap(list(N, G, simulation_index), simulate),
    naive_ate = map(data, \(data) {
      data %>%
        group_by(group) %>%
        nest() %>%
        mutate(naive_ate = map_dbl(data, \(data) mean(data$Y[data$A == 1]) - mean(data$Y[data$A == 0])))
    }),
    fit = pmap(list(simulation_index, N, G, data, naive_ate), \(simulation_index, N, G, data, naive_ate) {
      cache <- glue::glue("/gpfs/scratch/susmah01/bayes_tmle/cache/{N}/{simulation_index}_{G}.rds")
      if(file.exists(cache)) {
        return(read_rds(cache))
      }

      res <- bayes_tmle_hierarchical(
	data$Y, 
	data$A, 
	data[, c("W"), drop = FALSE], 
	data$group, 
	nuisance_method = "combined", 
	Q.SL.library = Q.SL.library, 
	g.SL.library = g.SL.library, 
	parallel_chains = parallel_chains, 
	adapt_delta = adapt_delta, 
	output_dir = "/gpfs/scratch/susmah01/bayes_tmle/fits/"
      )

      write_rds(tibble(N = N, G = G, simulation_index = simulation_index, res = list(res), naive_ate = list(naive_ate)), cache)

      return(res)
    })
  )
