# simulation for kociba data
library(ggplot2)
library(dplyr)
library(ToxicR)
library(simChef)
library(future)

source('simChef/dgps.R')
source('simChef/methods.R')
source('simChef/evaluators.R')
source('simChef/visualizers.R')
source('simChef/utility.R')

experiment = create_experiment(
  name = 'Kociba with quantal linear',
  save_dir = 'results/kociba_qlinear'
) %>%
  add_dgp(kociba1978_dgp) %>%
  add_method(uniform_design) %>%
  add_method(log_uniform_design) %>%
  add_method(lopt_D_design) %>%
  add_method(aug_D_design) %>%
  add_method(lopt_c_design) %>%
  add_method(aug_c_design) %>%
  add_evaluator(all_eval) %>%
  add_vary_across(
    .dgp = "Kociba et al. 1978",
    N1 = c(20, 30, 50, 75, 100, 117, 125, 150, 175, 200)
  )

plan(multisession, workers = 12)
set.seed(21225)
results = run_experiment(experiment, n_reps = 3000, save = T, use_cached = F)