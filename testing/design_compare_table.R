# table for comparing different optimal designs
source('design/aug_optimal.R')
source('design/models.R')
source('design/utility.R')
source('design/equiv_theorem.R')

# Kociba data
dat = data.frame(
  dose = c(0, 1, 10, 100),
  events = c(9, 3, 18, 34),
  n = c(86, 50, 50, 48)
)

# fit logistic model
mod = ToxicR::single_dichotomous_fit(
  D = dat$dose, Y = dat$events, N = dat$n, model_type = 'log-logistic', 
  fit_type = 'mle', alpha = 0.025)

plot(mod)

# find locally D-optimal design
D_lopt = find_aug_opt(
  grad_fun = loglogistic.grad,
  dr_fun = loglogistic.fun,
  bmd_grad = loglogistic.bmdgrad,
  obj = 'D',
  theta = mod$parameters,
  num_doses = 3,
  max_dose = max(dat$dose),
  swarm = 100,
  iter = 500, 
  alg = 'DE',
  d0 = c(0.0001, 1, 10, 100),
  n0 = c(86, 50, 50, 48),
  N1 = 100,
  ignore_stage1 = T
)
D_lopt
check_eq(D_lopt)

# find augmented D-optimal design
D_aug = find_aug_opt(
  grad_fun = loglogistic.grad,
  dr_fun = loglogistic.fun,
  bmd_grad = loglogistic.bmdgrad,
  obj = 'D',
  theta = mod$parameters,
  num_doses = 3,
  max_dose = max(dat$dose),
  swarm = 100,
  iter = 500, 
  alg = 'DE',
  d0 = c(0.0001, 1, 10, 100),
  n0 = c(86, 50, 50, 48),
  N1 = 100,
  ignore_stage1 = F
)
D_aug
check_eq(D_aug)

# find locally c-optimal design
c_lopt = find_aug_opt(
  grad_fun = loglogistic.grad,
  dr_fun = loglogistic.fun,
  bmd_grad = loglogistic.bmdgrad,
  obj = 'c',
  theta = mod$parameters,
  num_doses = 3,
  max_dose = max(dat$dose),
  swarm = 100,
  iter = 500, 
  alg = 'DE',
  d0 = c(0.0001, 1, 10, 100),
  n0 = c(86, 50, 50, 48),
  N1 = 100,
  ignore_stage1 = T
)
c_lopt
check_eq(c_lopt)

# find augmented c-optimal design
c_aug = find_aug_opt(
  grad_fun = loglogistic.grad,
  dr_fun = loglogistic.fun,
  bmd_grad = loglogistic.bmdgrad,
  obj = 'c',
  theta = mod$parameters,
  num_doses = 3,
  max_dose = max(dat$dose),
  swarm = 100,
  iter = 500, 
  alg = 'DE',
  d0 = c(0.0001, 1, 10, 100),
  n0 = c(86, 50, 50, 48),
  N1 = 100,
  ignore_stage1 = F
)
c_aug
check_eq(c_aug)
