# simulation study comparing naive designs
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
  name = 'Naive designs',
  save_dir = 'results/naive'
) %>%
  add_dgp(kociba1978_dgp) %>%
  add_dgp(JBRC_dgp) %>%
  add_dgp(deguelin_dgp) %>%
  add_dgp(kociba_loglogistic_dgp) %>%
  add_method(uniform_design) %>%
  add_method(log_uniform_design) %>%
  add_evaluator(all_eval) %>%
  add_vary_across(
    .method = "Single stage uniform design",
    num_doses = 4:10
  ) %>%
  add_vary_across(
    .method = "Single stage log-uniform design",
    num_doses = 4:10
  ) %>%
  add_vary_across(
    .dgp = "Kociba et al. 1978",
    N = c(20, 30, 50, 75, 100, 117, 125, 150, 175, 200, 234)
  ) %>%
  add_vary_across(
    .dgp = "JBRC",
    N = c(20, 30, 50, 75, 99, 100, 125, 150, 175, 198)
  ) %>%
  add_vary_across(
    .dgp = 'Deguelin',
    N = c(20, 30, 50, 75, 100, 125, 146, 150, 175, 200, 225, 250, 292)
  ) %>%
  add_vary_across(
    .dgp = 'Kociba log-logistic',
    N = c(20, 30, 50, 75, 100, 117, 125, 150, 175, 200)
  )



init_docs(experiment)

plan(multisession, workers = 12)
set.seed(162025)
results = run_experiment(experiment, n_reps = 3000, save = T, use_cached = F)

# render
#render_docs(experiment)
