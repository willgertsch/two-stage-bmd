# visualizer functions
CI_length_viz = create_visualizer(
  .viz_fun = function(eval_results) {
    eval_results[[1]] %>%
      mutate(method_name = reorder(as.factor(.method_name), `CI length`, FUN = median)) %>%
      ggplot(aes(x = method_name, y = log(`CI length`))) +
      geom_boxplot() +
      labs(y = 'Log confidence interval length', x = '') +
      theme_bw()
  }
)

CI_length_by_num_dose_viz = create_visualizer(
  .viz_fun = function(eval_results) {
    eval_results[[1]] %>%
      mutate(method_name = reorder(as.factor(.method_name), `CI length`, FUN = median)) %>%
      ggplot(aes(x = as.factor(num_doses), y = log(`CI length`))) +
      geom_boxplot() +
      labs(y = 'Log confidence interval length', x = '') +
      theme_bw() +
      facet_wrap(~.method_name)
  }
)

theta_mse_by_num_dose_viz = create_visualizer(
  .viz_fun = function(eval_results) {
    eval_results$`Compute MSE of theta` %>%
      ggplot(aes(x = as.factor(num_doses), y = `Theta MSE`)) +
      geom_point() + geom_smooth(se=F) +
      labs(y = 'Theta MSE', x = '') +
      theme_bw() +
      facet_wrap(~.method_name)
  }
)