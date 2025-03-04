# simulation for hill model with log-uniform design in the first stage
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
  name = 'JBRC data with Hill model and log-uniform first stage',
  save_dir = 'results/jbrc_hill_logunif'
) %>%
  add_dgp(JBRC_dgp) %>%
  add_method(uniform_design) %>%
  add_method(log_uniform_design) %>%
  add_method(lopt_D_design) %>%
  add_method(aug_D_design) %>%
  add_method(lopt_c_design) %>%
  add_method(aug_c_design) %>%
  add_evaluator(all_eval) %>%
  add_vary_across(
    .dgp = "JBRC",
    N1 = c(20, 30, 50, 75, 99, 100, 125, 150, 175)
  ) %>%
  add_vary_across(
    .method = c(lopt_D_design, aug_D_design, lopt_c_design, aug_c_design),
    first_stage = 'log-uniform',
    num_doses = 5
  ) %>%
  add_vary_across(
    .method = c(uniform_design, log_uniform_design),
    num_doses = 5
  )

plan(multisession, workers = 12)
set.seed(21125)
results = run_experiment(experiment, n_reps = 3000, save = T, use_cached = F)