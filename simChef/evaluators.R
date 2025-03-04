# evaluator functions

CI_length_eval = create_evaluator(
  .eval_fun = function(fit_results) {
    source('simChef/utility.R')
    fit_results %>%
      dplyr::mutate(
        `CI length` = sapply(mod, compute_BMD_CI_length)
      )
      
  },
  .name = 'BMD CI'
)

theta_mse_eval = create_evaluator(
  .eval_fun = function(fit_results) {
    source('simChef/utility.R')
    fit_results %>%
      dplyr::mutate(
        `Theta MSE` = as.numeric(purrr::map2(mod, theta, compute_theta_mse))
      )
  },
  .name = 'theta MSE'
)

# bias of theta
theta_bias_eval = create_evaluator(
  .eval_fun = function(fit_results) {
    source('simChef/utility.R')
    fit_results %>%
      mutate(
        theta_bias = purrr::map2(mod, theta, compute_theta_bias)
      )
  },
  .name = 'theta bias'
)


# collect all evaluators for two-stage designs
all_eval = create_evaluator(
  .eval_fun = function(fit_results) {
    source('simChef/utility.R')
    fit_results %>%
      mutate(
        CI_length = sapply(mod, compute_BMD_CI_length),
        bmd_diff = as.numeric(purrr::map2(mod, bmd, compute_bmd_diff)),
        coverage = as.numeric(purrr::map2(mod, bmd, bmd_coverage)),
        bmdl = sapply(mod, compute_BMDL)
      )
  },
  .name = 'all eval'
)
