# method functions
uniform_design = create_method(
  .method_fun = function(theta,
                         N1,
                         N,
                         dat,
                         max_dose,
                         model_type,
                         bmd,
                         num_doses = 4) {
    source('simChef/utility.R')
    x = seq(0, max_dose, length.out = num_doses)
    w = rep(1/num_doses, num_doses)
    dat = sim_data(N, theta, x, w, model_type)
    
    # fit model to the data
    mod = ToxicR::single_dichotomous_fit(
      D = dat$dose,
      Y = dat$events,
      N = dat$n,
      model_type = model_type,
      fit_type = 'mle',
      BMR = 0.1,
      alpha = 0.025
    )
    
    # create a smaller model object to save on memory
    small_mod = list(
      bmd = mod$bmd
    )
    
    return(list(mod = small_mod, theta = theta, theta_hat = mod$parameters, 
                bmd = bmd, stage1_x = dat$dose, stage1_n = dat$n/2,
                stage2_x = dat$dose, stage2_n = dat$n/2))
  },
  .name = 'Single stage uniform design'
)

# log-uniform designs have points that are equally far apart on the log scale
# always includes the max dose and the control 0 dose
log_uniform_design = create_method(
  .method_fun = function(theta,
                         N1,
                         N,
                         dat,
                         max_dose,
                         model_type,
                         bmd,
                         num_doses = 4) {
    source('simChef/utility.R')
    x = c(0, exp(seq(1, log(max_dose), by = (log(max_dose)-1)/(num_doses - 2))))
    w = rep(1/num_doses, num_doses)
    dat = sim_data(N, theta, x, w, model_type)
    
    # fit model to the data
    mod = ToxicR::single_dichotomous_fit(
      D = dat$dose,
      Y = dat$events,
      N = dat$n,
      model_type = model_type,
      fit_type = 'mle',
      BMR = 0.1,
      alpha = 0.025
    )
    
    # create a smaller model object to save on memory
    small_mod = list(
      bmd = mod$bmd
    )
    
    return(list(mod = small_mod, theta = theta, theta_hat = mod$parameters, 
                bmd = bmd, stage1_x = dat$dose, stage1_n = dat$n/2,
                stage2_x = dat$dose, stage2_n = dat$n/2))
  },
  .name = 'Single stage log-uniform design'
)

# locally optimal designs
lopt_D_design = create_method(
  .method_fun = function(theta,
                         N1,
                         N,
                         dat,
                         max_dose,
                         model_type,
                         bmd,
                         num_doses = 4,
                         first_stage = 'uniform',
                         swarm = 20,
                         iter = 500,
                         use_true_theta = F) {
    
    source('simChef/utility.R')
    source('design/aug_optimal.R')
    source('design/models.R')
    source('design/utility.R')
    
    # generate stage 1 data
    if (first_stage == 'uniform') {
      x1 = seq(0.00001, max_dose, length.out = num_doses)
      w1 = rep(1/num_doses, num_doses)
    }
    else if (first_stage == 'log-uniform') {
      x1 = c(0.00001, exp(seq(1, log(max_dose), by = (log(max_dose)-1)/(num_doses - 2))))
      w1 = rep(1/num_doses, num_doses)
    }
    else
      stop("invalid first stage design")
    
    dat1 = sim_data(N1, theta, x1, w1, model_type)
    
    # fit model to data
    mod1 = ToxicR::single_dichotomous_fit(
      D = dat1$dose,
      Y = dat1$events,
      N = dat1$n,
      model_type = model_type,
      fit_type = 'mle',
      BMR = 0.1,
      alpha = 0.025
    )
    if (use_true_theta) {
      theta1 = theta
    }
    else {
      theta1 = mod1$parameters
    }
    
    
    # look up model information
    mod_info = mod_info_lookup(model_type)
    
    # find the optimal design
    design = find_aug_opt(
      grad_fun = mod_info$grad_fun,
      dr_fun = mod_info$dr_fun,
      bmd_grad = mod_info$bmd_grad,
      obj = 'D',
      theta = theta1,
      num_doses = mod_info$pts,
      max_dose = max_dose,
      swarm = swarm,
      iter = iter,
      alg = 'DE',
      d0 = dat1$dose,
      n0 = dat1$n,
      N1 = N - N1,
      ignore_stage1 = T,
      show_progress = F
    )
    x2 = design$x
    w2 = design$w
    dat2 = sim_data(N - N1, theta, x2, w2, model_type)
    
    # combine data and fit final model
    dat12 = rbind(
      dat1, dat2
    )
    
    # fit model
    mod12 = ToxicR::single_dichotomous_fit(
      D = dat12$dose,
      Y = dat12$events,
      N = dat12$n,
      model_type = model_type,
      fit_type = 'mle',
      BMR = 0.1,
      alpha = 0.025
    )
    
    # create a smaller model object to save on memory
    small_mod = list(
      bmd = mod12$bmd
    )
    
    return(list(mod = small_mod, theta = theta, theta_hat = mod12$parameters, 
                bmd = bmd, stage1_x = dat1$dose, stage1_n = dat1$n,
                stage2_x = dat2$dose, stage2_n = dat2$n))
    
  },
  .name = 'Locally D-optimal design'
)

lopt_c_design = create_method(
  .method_fun = function(theta,
                         N1,
                         N,
                         dat,
                         max_dose,
                         model_type,
                         bmd,
                         num_doses = 4,
                         first_stage = 'uniform',
                         swarm = 20,
                         iter = 500,
                         use_true_theta = F) {
    source('simChef/utility.R')
    source('design/aug_optimal.R')
    source('design/models.R')
    source('design/utility.R')
    
    # generate stage 1 data
    if (first_stage == 'uniform') {
      x1 = seq(0.00001, max_dose, length.out = num_doses)
      w1 = rep(1/num_doses, num_doses)
    }
    else if (first_stage == 'log-uniform') {
      x1 = c(0.00001, exp(seq(1, log(max_dose), by = (log(max_dose)-1)/(num_doses - 2))))
      w1 = rep(1/num_doses, num_doses)
    }
    else
      stop("invalid first stage design")
    
    dat1 = sim_data(N1, theta, x1, w1, model_type)
    
    # fit model to data
    mod1 = ToxicR::single_dichotomous_fit(
      D = dat1$dose,
      Y = dat1$events,
      N = dat1$n,
      model_type = model_type,
      fit_type = 'mle',
      BMR = 0.1,
      alpha = 0.025
    )
    if (use_true_theta) {
      theta1 = theta
    }
    else {
      theta1 = mod1$parameters
    }
    
    # look up model information
    mod_info = mod_info_lookup(model_type)
    
    # find the optimal design
    design = find_aug_opt(
      grad_fun = mod_info$grad_fun,
      dr_fun = mod_info$dr_fun,
      bmd_grad = mod_info$bmd_grad,
      obj = 'c',
      theta = theta1,
      num_doses = mod_info$pts,
      max_dose = max_dose,
      swarm = swarm,
      iter = iter,
      alg = 'DE',
      d0 = dat1$dose,
      n0 = dat1$n,
      N1 = N - N1,
      ignore_stage1 = T,
      show_progress = F
    )
    x2 = design$x
    w2 = design$w
    dat2 = sim_data(N - N1, theta, x2, w2, model_type)
    
    # combine data and fit final model
    dat12 = rbind(
      dat1, dat2
    )
    
    # fit model
    mod12 = ToxicR::single_dichotomous_fit(
      D = dat12$dose,
      Y = dat12$events,
      N = dat12$n,
      model_type = model_type,
      fit_type = 'mle',
      BMR = 0.1,
      alpha = 0.025
    )
    
    # create a smaller model object to save on memory
    small_mod = list(
      bmd = mod12$bmd
    )
    
    return(list(mod = small_mod, theta = theta, theta_hat = mod12$parameters, 
                bmd = bmd, stage1_x = dat1$dose, stage1_n = dat1$n,
                stage2_x = dat2$dose, stage2_n = dat2$n))
    
  },
  .name = 'Locally c-optimal design'
)

aug_D_design = create_method(
  .method_fun = function(theta,
                         N1,
                         N,
                         dat,
                         max_dose,
                         model_type,
                         bmd,
                         num_doses = 4,
                         first_stage = 'uniform',
                         swarm = 20,
                         iter = 500,
                         use_true_theta = F) {
    source('simChef/utility.R')
    source('design/aug_optimal.R')
    source('design/models.R')
    source('design/utility.R')
    
    # generate stage 1 data
    if (first_stage == 'uniform') {
      x1 = seq(0.00001, max_dose, length.out = num_doses)
      w1 = rep(1/num_doses, num_doses)
    }
    else if (first_stage == 'log-uniform') {
      x1 = c(0.00001, exp(seq(1, log(max_dose), by = (log(max_dose)-1)/(num_doses - 2))))
      w1 = rep(1/num_doses, num_doses)
    }
    else
      stop("invalid first stage design")
    
    dat1 = sim_data(N1, theta, x1, w1, model_type)
    
    # fit model to data
    mod1 = ToxicR::single_dichotomous_fit(
      D = dat1$dose,
      Y = dat1$events,
      N = dat1$n,
      model_type = model_type,
      fit_type = 'mle',
      BMR = 0.1,
      alpha = 0.025
    )
    if (use_true_theta) {
      theta1 = theta
    }
    else {
      theta1 = mod1$parameters
    }
    
    # look up model information
    mod_info = mod_info_lookup(model_type)
    
    # find the optimal design
    design = find_aug_opt(
      grad_fun = mod_info$grad_fun,
      dr_fun = mod_info$dr_fun,
      bmd_grad = mod_info$bmd_grad,
      obj = 'D',
      theta = theta1,
      num_doses = mod_info$pts,
      max_dose = max_dose,
      swarm = swarm,
      iter = iter,
      alg = 'DE',
      d0 = dat1$dose,
      n0 = dat1$n,
      N1 = N - N1,
      ignore_stage1 = F,
      show_progress = F
    )
    x2 = design$x
    w2 = design$w
    dat2 = sim_data(N - N1, theta, x2, w2, model_type)
    
    # combine data and fit final model
    dat12 = rbind(
      dat1, dat2
    )
    
    # fit model
    mod12 = ToxicR::single_dichotomous_fit(
      D = dat12$dose,
      Y = dat12$events,
      N = dat12$n,
      model_type = model_type,
      fit_type = 'mle',
      BMR = 0.1,
      alpha = 0.025
    )
    
    # create a smaller model object to save on memory
    small_mod = list(
      bmd = mod12$bmd
    )
    
    return(list(mod = small_mod, theta = theta, theta_hat = mod12$parameters, 
                bmd = bmd, stage1_x = dat1$dose, stage1_n = dat1$n,
                stage2_x = dat2$dose, stage2_n = dat2$n))
    
  },
  .name = 'Augmented D-optimal design'
)

aug_c_design = create_method(
  .method_fun = function(theta,
                         N1,
                         N,
                         dat,
                         max_dose,
                         model_type,
                         bmd,
                         num_doses = 4,
                         first_stage = 'uniform',
                         swarm = 20,
                         iter = 500,
                         use_true_theta = F) {
    
    source('simChef/utility.R')
    source('design/aug_optimal.R')
    source('design/models.R')
    source('design/utility.R')
    
    # generate stage 1 data
    if (first_stage == 'uniform') {
      x1 = seq(0.00001, max_dose, length.out = num_doses)
      w1 = rep(1/num_doses, num_doses)
    }
    else if (first_stage == 'log-uniform') {
      x1 = c(0.00001, exp(seq(1, log(max_dose), by = (log(max_dose)-1)/(num_doses - 2))))
      w1 = rep(1/num_doses, num_doses)
    }
    else
      stop("invalid first stage design")
    
    dat1 = sim_data(N1, theta, x1, w1, model_type)
    
    # fit model to data
    mod1 = ToxicR::single_dichotomous_fit(
      D = dat1$dose,
      Y = dat1$events,
      N = dat1$n,
      model_type = model_type,
      fit_type = 'mle',
      BMR = 0.1,
      alpha = 0.025
    )
    if (use_true_theta) {
      theta1 = theta
    }
    else {
      theta1 = mod1$parameters
    }
    
    # look up model information
    mod_info = mod_info_lookup(model_type)
    
    # find the optimal design
    design = find_aug_opt(
      grad_fun = mod_info$grad_fun,
      dr_fun = mod_info$dr_fun,
      bmd_grad = mod_info$bmd_grad,
      obj = 'c',
      theta = theta1,
      num_doses = mod_info$pts,
      max_dose = max_dose,
      swarm = swarm,
      iter = iter,
      alg = 'DE',
      d0 = dat1$dose,
      n0 = dat1$n,
      N1 = N - N1,
      ignore_stage1 = F,
      show_progress = F
    )
    x2 = design$x
    w2 = design$w
    dat2 = sim_data(N - N1, theta, x2, w2, model_type)
    
    # combine data and fit final model
    dat12 = rbind(
      dat1, dat2
    )
    
    # fit model
    mod12 = ToxicR::single_dichotomous_fit(
      D = dat12$dose,
      Y = dat12$events,
      N = dat12$n,
      model_type = model_type,
      fit_type = 'mle',
      BMR = 0.1,
      alpha = 0.025
    )
    
    # create a smaller model object to save on memory
    small_mod = list(
      bmd = mod12$bmd
    )
    
    return(list(mod = small_mod, theta = theta, theta_hat = mod12$parameters, 
                bmd = bmd, stage1_x = dat1$dose, stage1_n = dat1$n,
                stage2_x = dat2$dose, stage2_n = dat2$n))
  },
  .name = 'Augmented c-optimal design'
)
