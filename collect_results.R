library(tidyverse)
files <- Sys.glob("/gpfs/scratch/susmah01/bayes_tmle/cache/*/*.rds")

results <- map_df(files, read_rds) %>% dplyr::bind_rows()

write_rds(results, file = "/gpfs/home/susmah01/bayes_tmle/results/simulation_results.rds")
